import torch.cuda
from flask_restful import Api, Resource, reqparse
from chebai.models.electra import Electra
from chebai.preprocessing.reader import ChemDataReader, EMBEDDING_OFFSET
from chebai.preprocessing.collate import RaggedCollater
from chebai.result.molplot import AttentionMolPlot, AttentionNetwork
from tempfile import NamedTemporaryFile
from PIL import Image
import base64
import io
import numpy as np
import sys
from app import app
import matplotlib as mpl
import json
import networkx as nx
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import torch
mpl.use("Agg")

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"

BATCH_SIZE = app.config.get("BATCH_SIZE", 100)
electra_model = None


def load_model():
    global electra_model
    electra_model = Electra.load_from_checkpoint(app.config["ELECTRA_CHECKPOINT"], map_location=torch.device(device), pretrained_checkpoint=None, criterion=None, strict=False, metrics=dict(train=dict(), test=dict(), validation=dict()))
    electra_model.eval()

try:
    import uwsgidecorators
except ModuleNotFoundError:
    load_model()
else:
    load_model = uwsgidecorators.postfork(load_model)


PREDICTION_HEADERS = ["CHEBI:" + r.strip() for r in open(app.config["CLASS_HEADERS"])]
LABEL_HIERARCHY = json.load(open(app.config["CHEBI_JSON"]))

def load_sub_ontology():
    g = nx.DiGraph()
    g.add_nodes_from([(c["ID"], dict(lbl=c["LABEL"][0])) for c in LABEL_HIERARCHY])
    g.add_edges_from([(t, c["ID"]) for c in LABEL_HIERARCHY for t in c.get("SubClasses",[])])
    return g

CHEBI_FRAGMENT = load_sub_ontology()


def get_relevant_chebi_fragment(predictions, smiles, labels=None):
    fragment_graph = CHEBI_FRAGMENT.copy()
    predicted_subsumtions = [(i, h) for i in range(predictions.shape[0]) for h, p in zip(PREDICTION_HEADERS, predictions[i].tolist()) if p >= 0]
    fragment_graph.add_edges_from(predicted_subsumtions)
    necessary_nodes = set()
    for i in range(predictions.shape[0]):
        fragment_graph.nodes[i]["lbl"] = labels[i] if labels else smiles[i]
        fragment_graph.nodes[i]["artificial"] = True
        necessary_nodes = necessary_nodes.union(set(nx.shortest_path(fragment_graph, i)))
    sub = nx.transitive_reduction(fragment_graph.subgraph(necessary_nodes))

    # Copy node data to subgraph
    sub.add_nodes_from(d for d in fragment_graph.nodes(data=True) if d[0] in sub.nodes)
    keys = set(i for i,_ in predicted_subsumtions)
    return sub, {k: [l for j,l in predicted_subsumtions if j == k] for k in keys}

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


def batchify(l):
    cache = []
    for r in l:
        cache.append(r)
        if len(cache) >= BATCH_SIZE:
            yield cache
            cache = []
    if cache:
        yield cache


def calculate_results(batch):
    collater = RaggedCollater()
    dat = electra_model._process_batch(collater(batch).to(electra_model.device), 0)
    dat["features"] = dat["features"].int()
    return electra_model(dat, **dat["model_kwargs"])


class HierarchyAPI(Resource):
    def get(self):
        return {r["ID"]: dict(label=r["LABEL"][0], children=r.get("SubClasses",[])) for r in LABEL_HIERARCHY}


class BatchPrediction(Resource):
    def post(self):
        """
        Accepts a dictionary with the following structure
        {
            "smiles": [ ... list of smiles strings]
            "ontology": bool (Optional)
        }
        :return:
        A dictionary wit hthe following structure
        {
            "predicted_parents": [ ... [... parent classes as predicted by the system] or None for each smiles ],
            "direct_parents": [ ... [... lowest possible predicted parents] or None for each smiles ] or None
            "ontology": Only returened if `ontology` is set. Returns a vis.js conform representation of the ontology containing all predicted classes.
        }

        If the system us unable to parse any smiles string, the respective entry in each list will be `None`.
        """

        parser = reqparse.RequestParser()
        parser.add_argument("smiles", type=str, action="append")
        parser.add_argument("ontology", type=bool, required=False, default=False)
        args = parser.parse_args()
        smiles = args["smiles"]
        generate_ontology = args["ontology"]

        reader = ChemDataReader()
        collater = RaggedCollater()
        token_dicts = []
        could_not_parse = []
        index_map = dict()
        for i, s in enumerate(smiles):
            try:
                # Try to parse the smiles string
                if not s:
                    raise ValueError()
                d = reader.to_data(dict(features=s, labels=None))
                # This is just for sanity checks
                rdmol = Chem.MolFromSmiles(s, sanitize=False)
            except Exception as e:
                # Note if it fails
                could_not_parse.append(i)
            else:
                if rdmol is None:
                    could_not_parse.append(i)
                else:
                    index_map[i] = len(token_dicts)
                    token_dicts.append(d)
        results = []
        if token_dicts:
            for batch in batchify(token_dicts):
                result = calculate_results(batch)
                results += result["logits"].cpu().detach().tolist()

            chebi, predicted_parents = get_relevant_chebi_fragment(np.stack(results, axis=0), smiles)
        else:
            chebi, predicted_parents = ([], [])

        result = {
            "predicted_parents": [(None if i in could_not_parse else predicted_parents[index_map[i]]) for i in range(len(smiles))],
            "direct_parents": [(None if i in could_not_parse else list(chebi.successors(index_map[i]))) for i in range(len(smiles))],

        }

        if generate_ontology:
            result["ontology"] = nx_to_graph(chebi)

        return result


class PredictionDetailApiHandler(Resource):
    def load_image(self, path):
        im = Image.open(path)
        data = io.BytesIO()
        im.save(data, "PNG")
        encoded_img_data = base64.b64encode(data.getvalue())
        return encoded_img_data.decode("utf-8")

    def build_graph_from_attention(self, att, node_labels, token_labels, threshold=0.0):
        n_nodes = len(node_labels)
        return dict(
            nodes=[
                dict(
                    label=token_labels[n],
                    id=f"{group}_{i}",
                    fixed=dict(x=True, y=True),
                    y=100 * int(group == "r"),
                    x=30 * i,
                    group=group,
                )
                for i, n in enumerate([0] + node_labels)
                for group in ("l", "r")
            ],
            edges=[
                {
                    "from": f"l_{i}",
                    "to": f"r_{j}",
                    "color": {"opacity": att[i, j].item()},
                    "smooth": False,
                    "physics": False,
                }
                for i in range(n_nodes)
                for j in range(n_nodes)
                if att[i, j] > threshold
            ],
        )





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

        reader = ChemDataReader()
        token_dict = reader.to_data(dict(features=smiles, labels=None))
        tokens = np.array(token_dict["features"]).astype(int).tolist()
        result = calculate_results([token_dict])

        token_labels = (
            ["[CLR]"] + [None for _ in range(EMBEDDING_OFFSET - 1)] + reader.cache
        )

        graphs = [
            [
                self.build_graph_from_attention(
                    a[0, i], tokens, token_labels, threshold=0.1
                )
                for i in range(a.shape[1])
            ]
            for a in result["attentions"]
        ]

        chebi, predicted_parents = get_relevant_chebi_fragment(result["logits"].detach().cpu().numpy(), [smiles])

        with NamedTemporaryFile(mode="wt", suffix=".png") as svg1:
            rdmol = Chem.MolFromSmiles(smiles)
            d = rdMolDraw2D.MolDraw2DCairo(500, 500)
            rdMolDraw2D.PrepareAndDrawMolecule(d, rdmol)
            d.FinishDrawing()
            d.WriteDrawingText(svg1.name)
            mol_pic = self.load_image(svg1.name)

        return {
            "figures": {"plain_molecule": mol_pic},
            "graphs": graphs,
            "classification": nx_to_graph(chebi)
        }
