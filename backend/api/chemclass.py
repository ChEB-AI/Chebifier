from flask_restful import Api, Resource, reqparse
from chebai.models.electra import Electra
from chebai.preprocessing.reader import ChemDataReader
from chebai.preprocessing.collate import RaggedCollater
from chebai.result.molplot import AttentionMolPlot
from tempfile import NamedTemporaryFile
from PIL import Image
import base64
import io
import numpy as np

model_kwargs = dict(
    optimizer_kwargs=dict(lr=1e-4),
    pretrained_checkpoint="Electra_ChEBI100.ckpt",
    out_dim=854,
    config=dict(
        vocab_size=1400,
        max_position_embeddings=1800,
        num_attention_heads=8,
        num_hidden_layers=6,
        type_vocab_size=1,
    )
)

electra_model = Electra(**model_kwargs)


class PredictionApiHandler(Resource):
    def post(self):
        parser = reqparse.RequestParser()
        plotter = AttentionMolPlot()
        parser.add_argument('type', type=str)
        parser.add_argument('smiles', type=str)

        args = parser.parse_args()

        print(args)
        # note, the post req from frontend needs to match the strings here (e.g. 'type and 'message')

        request_type = args['type']
        smiles = args['smiles']

        reader = ChemDataReader()
        collater = RaggedCollater()
        tokenised_input = electra_model._get_data_and_labels(collater([reader.to_data(dict(features=smiles, labels=None))]), 0)
        result = electra_model(tokenised_input)
        atts = np.concatenate([a.detach().numpy() for a in result["attentions"]])
        d = plotter.draw_attention_molecule(smiles, np.max(np.max(atts, axis=2), axis=1))
        with NamedTemporaryFile(mode="wt", suffix=".png") as svg1:
            d.WriteDrawingText(svg1.name)
            im = Image.open(svg1.name)
            data = io.BytesIO()
            im.save(data, "JPEG")
            encoded_img_data = base64.b64encode(data.getvalue())
        return {"logits": result["logits"].tolist(), "attention_fig": encoded_img_data.decode('utf-8')}

