from typing import Optional, Union

import numpy as np
import torch
from chebai.models.electra import Electra
from chebai.preprocessing.datasets.base import XYBaseDataModule
from chebai.preprocessing.datasets.chebi import ChEBIOver100, ChEBIOver50
from chebai.preprocessing.reader import EMBEDDING_OFFSET

from prediction_models.base import NNPredictionModel

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"


def build_graph_from_attention(att, node_labels, token_labels, threshold=0.0):
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


class ElectraModel(NNPredictionModel):

    def __init__(self, checkpoint_path: str, data_class: Union[XYBaseDataModule, str],
                 prediction_headers_path: str, batch_size: Optional[int] = 8, name: Optional[str] = None):
        super().__init__(prediction_headers_path, batch_size, name)
        self.model = Electra.load_from_checkpoint(
            checkpoint_path, map_location=torch.device(device), criterion=None, strict=False,
            metrics=dict(train=dict(), test=dict(), validation=dict()), pretrained_checkpoint=None)
        self.model.eval()
        if isinstance(data_class, str):
            if data_class == "ChEBI100":
                data_class = ChEBIOver100()
            elif data_class == "ChEBI50":
                data_class = ChEBIOver50()
            else:
                raise ValueError(f"Unknown data class {data_class}")
        self.data_class = data_class

    @property
    def default_name(self) -> str:
        return "ELECTRA"

    @property
    def info_text(self) -> str:
        return ("Transformer model for predicting arbitrary ChEBI classes that have a certain number of instances, "
                "see [1].")

    def explain(self, smiles) -> dict:
        reader = self.data_class.READER()
        token_dict = reader.to_data(dict(features=smiles, labels=None))
        tokens = np.array(token_dict["features"]).astype(int).tolist()
        result = self.calculate_results([token_dict])

        token_labels = (
                ["[CLR]"] + [None for _ in range(EMBEDDING_OFFSET - 1)] + reader.cache
        )

        graphs = [
            [
                build_graph_from_attention(
                    a[0, i], tokens, token_labels, threshold=0.1
                )
                for i in range(a.shape[1])
            ]
            for a in result["attentions"]
        ]
        return {"graphs": graphs}
