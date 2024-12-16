import queue
import sys

from flask_restful import Resource, reqparse
from chebai.result.molplot import AttentionMolPlot, AttentionNetwork
from prediction_models.electra_model import ElectraModel
from prediction_models.chemlog import ChemLog
from PIL import Image
import base64
import io
from app import app
import matplotlib as mpl
import json
import networkx as nx
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import torch
from chebi_utils import LABEL_HIERARCHY, CHEBI_FRAGMENT

mpl.use("Agg")

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"

AVAILABLE_MODELS = []
for model in app.config["MODELS"]:
    cls = getattr(sys.modules[__name__], model["class"])
    del model["class"]
    AVAILABLE_MODELS.append(cls(**model))

def _build_node(ident, node, include_labels=True):
    d = dict(id=ident,
             color="#EEEEEC" if node.get("artificial") else "#729FCF")
    d["title"] = node["lbl"]
    if include_labels:
        d["label"] = node["lbl"]

    return d


def nx_to_graph(g: nx.Graph):
    return dict(
        nodes={n: g.nodes[n] for n in g.nodes},
        edges=list(g.edges)
    )

class HierarchyAPI(Resource):
    def get(self):
        return {
            "available_models": [model.name for model in AVAILABLE_MODELS],
            "available_models_info_texts": [model.info_text for model in AVAILABLE_MODELS],
            "hierarchy": {r["ID"]: dict(label=r["LABEL"][0], children=r.get("SubClasses", [])) for r in LABEL_HIERARCHY}
        }


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

        If the system us unable to parse any smiles string, the respective entry in each list will be `None`.
        """
        parser = reqparse.RequestParser()
        parser.add_argument("smiles", type=str, action="append")
        parser.add_argument("ontology", type=bool, required=False, default=False)
        parser.add_argument("selectedModels", type=dict)
        args = parser.parse_args()
        smiles = args["smiles"]
        generate_ontology = args["ontology"]
        selected_models = args["selectedModels"]

        predicted_classes = [[] for _ in smiles]
        for model in AVAILABLE_MODELS:
            if model.name in selected_models.keys() and selected_models[model.name]:
                new_preds = model.predict(smiles)
                for i in range(len(smiles)):
                    if new_preds[i] is not None:
                        predicted_classes[i] += new_preds[i]

        # get all parents of predicted classes
        all_predicted = [cls for clss in predicted_classes for cls in clss]
        q = queue.Queue()
        for cls in all_predicted:
            q.put(cls)
        while not q.empty():
            cls = q.get()
            for parent in CHEBI_FRAGMENT.successors(cls):
                if parent not in all_predicted:
                    all_predicted.append(parent)
                    q.put(parent)

        # todo check for disjointness violations

        predicted_graph = CHEBI_FRAGMENT.subgraph(all_predicted)
        graphs_per_smiles = [CHEBI_FRAGMENT.subgraph(predicted) for predicted in predicted_classes]

        result = {
            "predicted_parents": predicted_classes,
            "direct_parents": [[cls for cls in graph.nodes if not any(predecessor in graph.nodes
                                                                      for predecessor in graph.predecessors(cls))]
                               for graph in graphs_per_smiles],
        }

        if generate_ontology:
            result["ontology"] = nx_to_graph(predicted_graph)

        return result


class PredictionDetailApiHandler(Resource):

    def load_image(self, path):
        im = Image.open(path)
        data = io.BytesIO()
        im.save(data, "PNG")
        encoded_img_data = base64.b64encode(data.getvalue())
        return encoded_img_data.decode("utf-8")

    def post(self):
        parser = reqparse.RequestParser()
        plotter = AttentionMolPlot()
        network_plotter = AttentionNetwork()
        parser.add_argument("type", type=str)
        parser.add_argument("smiles", type=str)

        args = parser.parse_args()

        # note, the post req from frontend needs to match the strings here (e.g. 'type and 'message')

        request_type = args["type"]
        smiles = args["smiles"]

        predicted_classes = []
        explain_infos = {}
        for model in AVAILABLE_MODELS:
            predicted_classes += model.predict([smiles])
            explain_infos.update(model.explain(smiles))

        chebi = CHEBI_FRAGMENT.subgraph(predicted_classes)

        with open("image.png", mode="wt") as svg1:
            rdmol = Chem.MolFromSmiles(smiles)
            d = rdMolDraw2D.MolDraw2DCairo(500, 500)
            rdMolDraw2D.PrepareAndDrawMolecule(d, rdmol)
            d.FinishDrawing()
            d.WriteDrawingText(svg1.name)
            mol_pic = self.load_image(svg1.name)

        explain_infos["figures"] = {"plain_molecule": mol_pic}
        explain_infos["classification"] = nx_to_graph(chebi)

        print("Prediction Detail Result: ", explain_infos)

        return explain_infos
