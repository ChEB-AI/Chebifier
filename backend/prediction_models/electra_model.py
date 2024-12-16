from typing import Optional, Union

import numpy as np
import torch
from chebai.models.electra import Electra
from chebai.preprocessing.datasets.base import XYBaseDataModule
from chebai.preprocessing.reader import EMBEDDING_OFFSET
from rdkit import Chem
from chebai.preprocessing.datasets.chebi import ChEBIOver100, ChEBIOver50


from prediction_models.base import PredictionModel

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"


class ElectraModel(PredictionModel):

    def __init__(self, checkpoint_path: str, data_class: Union[XYBaseDataModule, str],
                 prediction_headers_path: str, batch_size: Optional[int] = 8, name: Optional[str] = None):
        super().__init__(name)
        self.electra_model = Electra.load_from_checkpoint(
            checkpoint_path, map_location=torch.device(device), criterion=None, strict=False,
            metrics=dict(train=dict(), test=dict(), validation=dict()), pretrained_checkpoint=None)
        self.electra_model.eval()
        self.batch_size = batch_size
        if isinstance(data_class, str):
            if data_class == "ChEBI100":
                data_class = ChEBIOver100()
            elif data_class == "ChEBI50":
                data_class = ChEBIOver50()
            else:
                raise ValueError(f"Unknown data class {data_class}")
        self.data_class = data_class
        self.prediction_headers = ["CHEBI:" + r.strip() for r in open(prediction_headers_path, encoding="utf-8")]

    @property
    def default_name(self) -> str:
        return "ELECTRA"

    @property
    def info_text(self) -> str:
        return ("Transformer model for predicting arbitrary ChEBI classes that have a certain number of instances, "
                "see [1].")

    def calculate_results(self, batch):
        collator = self.data_class.READER.COLLATOR()
        dat = self.electra_model._process_batch(collator(batch).to(self.electra_model.device), 0)
        dat["features"] = dat["features"].int()
        return self.electra_model(dat, **dat["model_kwargs"])

    def batchify(self, batch):
        cache = []
        for r in batch:
            cache.append(r)
            if len(cache) >= self.batch_size:
                yield cache
                cache = []
        if cache:
            yield cache

    def predict(self, smiles_list) -> list:
        reader = self.data_class.READER()
        token_dicts = []
        could_not_parse = []
        index_map = dict()
        for i, smiles in enumerate(smiles_list):
            try:
                # Try to parse the smiles string
                if not smiles:
                    raise ValueError()
                d = reader.to_data(dict(features=smiles, labels=None))
                # This is just for sanity checks
                rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
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
            for batch in self.batchify(token_dicts):
                result = self.calculate_results(batch)
                results += result["logits"].cpu().detach().tolist()
            results = np.stack(results, axis=0)
            positive_preds = [[self.prediction_headers[j] for j, p in enumerate(results[index_map[i]]) if p > 0.5]
                              if i not in could_not_parse else None for i in range(len(smiles_list))]
            return positive_preds
        else:
            return [None for _ in smiles_list]

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
                self.build_graph_from_attention(
                    a[0, i], tokens, token_labels, threshold=0.1
                )
                for i in range(a.shape[1])
            ]
            for a in result["attentions"]
        ]
        return {"graphs": graphs}
