import copy
import os
import sys

from flask_restful import Resource, reqparse
from PIL import Image
import base64
import io
from app import app
import matplotlib as mpl
import networkx as nx
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import torch
from chebi_utils import LABEL_HIERARCHY, CHEBI_FRAGMENT, get_transitive_predictions
import hashlib

from chebifier.ensemble.weighted_majority_ensemble import WMVwithF1Ensemble

mpl.use("Agg")

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"

ENSEMBLE = WMVwithF1Ensemble(app.config["MODELS"])


def _build_node(ident, node, include_labels=True):
    d = dict(id=ident,
             color="#EEEEEC" if node.get("artificial") else "#729FCF")
    d["title"] = node["lbl"]
    if include_labels:
        d["label"] = node["lbl"]

    return d


def nx_to_graph(g: nx.Graph, colors=None):
    results = {
        "nodes": {n: g.nodes[n] for n in g.nodes},
        "edges": list(g.edges)
    }
    if colors is not None:
        results["node_colors"] = dict()
        for node, color in colors.items():
            if node in results["nodes"]:
                results["node_colors"][node] = color
    return results


class HierarchyAPI(Resource):
    def get(self):
        return {
            "available_models": [model.model_name for model in ENSEMBLE.models],
            "available_models_info_texts": [model.info_text for model in ENSEMBLE.models],
            "hierarchy": {r["ID"]: dict(label=r["LABEL"][0], children=r.get("SubClasses", [])) for r in LABEL_HIERARCHY}
        }


class ModelInfoAPI(Resource):

    def get(self):
        return {
            "available_models": [model.model_name for model in ENSEMBLE.models],
            "available_models_info_texts": [model.info_text for model in ENSEMBLE.models]
        }


def verify_disjointness(predicted_classes):
    disjoints = []
    with open(os.path.join("data", "disjoint_chebi.csv")) as f:
        for line in f:
            disjoints.append([f"CHEBI:{i}" for i in line.strip().split(",")])
    with open(os.path.join("data", "disjoint_additional.csv")) as f:
        for line in f:
            disjoints.append([f"CHEBI:{i}" for i in line.strip().split(",")])
    violations = []
    for sample in predicted_classes:
        violations_sample = []
        for disjoint in disjoints:
            if all(cls in sample for cls in disjoint):
                violations_sample.append(disjoint)
        violations.append(violations_sample)
    return violations


class BatchPrediction(Resource):
    def post(self):
        """
        Accepts a dictionary with the following structure
        {
            "smiles": [ ... list of smiles strings]
            "ontology": bool (Optional)
        }
        :return:
        A dictionary with the following structure
        {
            "predicted_parents": [ ... [... parent classes as predicted by the system] or None for each smiles ],
            "direct_parents": [ ... [... lowest possible predicted parents] or None for each smiles ] or None
            "ontology": Only returned if `ontology` is set. Returns a vis.js conform representation of the ontology containing all predicted classes.
        }

        If the system is unable to parse any smiles string, the respective entry in each list will be `None`.
        """
        parser = reqparse.RequestParser()
        parser.add_argument("smiles", type=str, action="append")
        parser.add_argument("ontology", type=bool, required=False, default=False)
        parser.add_argument("selectedModels", type=dict)
        args = parser.parse_args()
        smiles = args["smiles"]
        generate_ontology = args["ontology"]
        selected_models = args["selectedModels"]
        ensemble_models = ENSEMBLE.models
        ENSEMBLE.models = [model for model in ensemble_models if model.model_name in selected_models and selected_models[model.model_name]]

        if not smiles or len(smiles) == 0:
            result = {
                "predicted_parents": [],
                "direct_parents": [],
                "violations": []
            }
            if generate_ontology:
                # violation_colors = {violator: "#d92946" for violations_sample in violations for violation in violations_sample for violator in violation}
                result["ontology"] = nx_to_graph(CHEBI_FRAGMENT.subgraph([]))
            return result

        all_predicted = ENSEMBLE.predict_smiles_list(smiles, load_preds_if_possible=False)

        ENSEMBLE.models = ensemble_models  # restore original model list

        all_predicted = [[f"CHEBI:{cls}" for cls in sample if cls is not None] if sample is not None else None for sample in all_predicted]

        predicted_graph = CHEBI_FRAGMENT.subgraph([p for sample in all_predicted if sample is not None for p in sample])
        graphs_per_smiles = [CHEBI_FRAGMENT.subgraph(predicted) if predicted else None for predicted in all_predicted]

        direct_parents = [[cls for cls in graph.nodes if not any(predecessor in graph.nodes
                                                                 for predecessor in graph.predecessors(cls))] if graph is not None else None
                          for graph in graphs_per_smiles]

        # array with dimensions [n_samples, n_violations, 2]
        violations = verify_disjointness([p for p in all_predicted if p is not None])
        # replace the violation-causing classes with the direct predictions that are influenced by them
        violations_direct = [[[{v: [direct_pred for direct_pred in direct_preds
                                    if direct_pred == v or nx.has_path(graph_smiles, direct_pred, v)]
                                for v in violation}] for violation in violations_sample]
                             for violations_sample, direct_preds, graph_smiles in
                             zip(violations, direct_parents, graphs_per_smiles)]
        result = {
            "predicted_parents": all_predicted,
            "direct_parents": direct_parents,
            "violations": violations_direct
        }
        if generate_ontology:
            # violation_colors = {violator: "#d92946" for violations_sample in violations for violation in violations_sample for violator in violation}
            result["ontology"] = nx_to_graph(predicted_graph)
        return result


def to_color(model_name):
    # turn any string into a color
    h = hashlib.md5(model_name.encode()).hexdigest()
    return "#" + h[:6]


class PredictionDetailApiHandler(Resource):

    def load_image(self, path):
        im = Image.open(path)
        data = io.BytesIO()
        im.save(data, "PNG")
        encoded_img_data = base64.b64encode(data.getvalue())
        return encoded_img_data.decode("utf-8")

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument("type", type=str)
        parser.add_argument("smiles", type=str)
        parser.add_argument("selectedModels", type=dict)

        args = parser.parse_args()
        # note, the post req from frontend needs to match the strings here (e.g. 'type and 'message')
        request_type = args["type"]
        smiles = args["smiles"]
        selected_models = args["selectedModels"]

        predicted_classes = []
        predicted_by_model = {}
        explain_infos = {"models": {}}
        ensemble_models = ENSEMBLE.models
        ENSEMBLE.models = [model for model in ensemble_models if model.model_name in selected_models and selected_models[model.model_name]]
        for model in ENSEMBLE.models:
            pred = model.predict_smiles(smiles)
            if pred is not None:
                pred = [f"CHEBI:{cls}" for cls in pred if pred[cls] > ENSEMBLE.positive_prediction_threshold]
                predicted_classes += pred
                predicted_by_model[model.model_name] = pred
        for model in ENSEMBLE.models:
            explain_infos_model = model.explain_smiles(smiles)
            if explain_infos_model is not None:
                explain_infos_model["model_type"] = model.__class__.__name__
                explain_infos_model["model_info"] = model.info_text
                explain_infos["models"][model.model_name] = explain_infos_model

        # add ensemble predictions as own model
        pred = ENSEMBLE.predict_smiles_list([smiles], load_preds_if_possible=False)[0]
        ENSEMBLE.models = ensemble_models # restore original models
        if pred is not None:
            pred = [f"CHEBI:{cls}" for cls in pred if cls is not None]
            predicted_classes += pred
            predicted_by_model["Ensemble"] = pred

        all_predicted = get_transitive_predictions([predicted_classes])
        chebi = CHEBI_FRAGMENT.subgraph(all_predicted)

        predicted_by_per_node = {cls: [] for cls in all_predicted}
        for model_name, nodes in predicted_by_model.items():
            for node in nodes:
                predicted_by_per_node[node].append(model_name)
        node_colors = {node: to_color(str(predicted_by_per_node[node])) for node in all_predicted
                       }
        color_legend = {
            to_color(str(preds)): " and ".join(preds) if len(preds) > 0 else "Inferred based on predicted subclasses"
            for preds in predicted_by_per_node.values()}

        # todo replace with smiles (and paint with rdkit-js)
        with open("image.png", mode="wt") as svg1:
            rdmol = Chem.MolFromSmiles(smiles)
            d = rdMolDraw2D.MolDraw2DCairo(500, 500)
            rdMolDraw2D.PrepareAndDrawMolecule(d, rdmol)
            d.FinishDrawing()
            d.WriteDrawingText(svg1.name)
            mol_pic = self.load_image(svg1.name)

        explain_infos["figures"] = {"plain_molecule": mol_pic}
        explain_infos["classification"] = nx_to_graph(chebi, node_colors)
        explain_infos["color_legend"] = color_legend

        return explain_infos
