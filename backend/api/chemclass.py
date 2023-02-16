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


class PredictionApiHandler(Resource):
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

        atts = np.concatenate([a.detach().numpy() for a in result["attentions"]])

        with NamedTemporaryFile(mode="wt", suffix=".png") as svg1:
            d = plotter.draw_attention_molecule(
                smiles, np.max(np.max(atts, axis=2), axis=1)
            )
            d.WriteDrawingText(svg1.name)
            mol_pic = self.load_image(svg1.name)

        return {
            "logits": result["logits"].tolist(),
            "figures": {"attention_mol": mol_pic},
            "graphs": graphs,
        }
