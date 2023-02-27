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

mpl.use("TkAgg")

model_kwargs = dict(
    optimizer_kwargs=dict(lr=1e-4),
    pretrained_checkpoint=app.config["ELECTRA_CHECKPOINT"],
    out_dim=854,
    config=dict(
        vocab_size=1400,
        max_position_embeddings=1800,
        num_attention_heads=8,
        num_hidden_layers=6,
        type_vocab_size=1,
    ),
)

electra_model = Electra(**model_kwargs)
electra_model.eval()

PREDICTION_HEADERS = ["CHEBI:" + r.strip() for r in open("/home/glauer/dev/ChEBI_RvNN/data/ChEBI100/raw/classes.txt")]

def load_sub_ontology():
    d = json.load(open("/data/ontologies/chebi100.json"))
    g = nx.DiGraph()
    g.add_nodes_from([(c["ID"], dict(lbl=c["LABEL"][0])) for c in d])
    g.add_edges_from([(t, c["ID"]) for c in d for t in c.get("SubClasses",[])])
    return g

CHEBI_FRAGMENT = load_sub_ontology()


def get_relevant_chebi_fragment(predictions, smiles, labels=None):
    fragment_graph = CHEBI_FRAGMENT.copy()
    predicted_subsumtions = [(i, h) for i in range(predictions.shape[0]) for h, p in zip(PREDICTION_HEADERS, predictions[i].tolist()) if p >= 0]
    fragment_graph.add_edges_from(predicted_subsumtions)
    necessary_nodes = set()
    for i in range(predictions.shape[0]):
        fragment_graph.nodes[i]["lbl"] = labels[i] if labels else smiles[i]
        necessary_nodes = necessary_nodes.union(set(nx.shortest_path(fragment_graph, i)))
    sub = nx.transitive_reduction(fragment_graph.subgraph(necessary_nodes))

    # Copy node data to subgraph
    sub.add_nodes_from(d for d in fragment_graph.nodes(data=True) if d[0] in sub.nodes)
    return sub, dict(predicted_subsumtions)

def nx_to_graph(g: nx.Graph):
    return dict(
        nodes=[dict(id=n, label=g.nodes[n]["lbl"]) for n in g.nodes],
        edges=[{"from":a, "to":b, "arrows":dict(to=True)} for (a,b) in g.edges]
    )


class BatchPrediction(Resource):
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument("smiles", type=str, action="append")
        args = parser.parse_args()
        smiles = args["smiles"]

        reader = ChemDataReader()
        collater = RaggedCollater()
        token_dict = [reader.to_data(dict(features=s, labels=None)) for s in smiles]
        tokenised_input = electra_model._get_data_and_labels(collater(token_dict), 0)
        result = electra_model(tokenised_input)

        chebi, predicted_parents = get_relevant_chebi_fragment(result["logits"], smiles)

        return {"predicted_parents": [predicted_parents[i] for i in range(len(smiles))], "chebi": nx_to_graph(chebi), "direct_parents": [list(chebi.successors(i)) for i in range(len(smiles))]}


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

        print(args)
        # note, the post req from frontend needs to match the strings here (e.g. 'type and 'message')

        request_type = args["type"]
        smiles = args["smiles"]

        reader = ChemDataReader()
        collater = RaggedCollater()
        token_dict = reader.to_data(dict(features=smiles, labels=None))
        tokens = token_dict["features"]
        tokenised_input = electra_model._get_data_and_labels(collater([token_dict]), 0)
        result = electra_model(tokenised_input)

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

        _, chebi = get_relevant_chebi_fragment(result["logits"])

        atts = np.concatenate([a.detach().numpy() for a in result["attentions"]])

        with NamedTemporaryFile(mode="wt", suffix=".png") as svg1:
            d = plotter.draw_attention_molecule(
                smiles, np.max(np.max(atts, axis=2), axis=1)
            )
            d.WriteDrawingText(svg1.name)
            mol_pic = self.load_image(svg1.name)

        return {
            "figures": {"attention_mol": mol_pic},
            "graphs": graphs,
            "classification": nx_to_graph(chebi)
        }
