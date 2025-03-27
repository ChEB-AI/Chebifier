import abc
from typing import Optional

import numpy as np
from rdkit import Chem


class PredictionModel(abc.ABC):

    def __init__(self, name: Optional[str] = None, description: Optional[str] = None):
        self._name = name or self.default_name
        self._description = description

    @property
    def default_name(self) -> str:
        return "PredictionModel"

    @property
    def name(self) -> str:
        return self._name

    @property
    def info_text(self) -> str:
        if self._description is None:
            return "No description is available for this model."
        return self._description

    def predict(self, smiles_list) -> list:
        pass

    def explain(self, smiles) -> Optional[dict]:
        return None


class NNPredictionModel(PredictionModel):

    def __init__(self, prediction_headers_path: str, batch_size: Optional[int] = 32, name: Optional[str] = None,
                 description: Optional[str] = None):
        super().__init__(name, description)
        self.batch_size = batch_size
        self.prediction_headers = ["CHEBI:" + r.strip() for r in open(prediction_headers_path, encoding="utf-8")]
        self.model = None
        self.data_class = None

    def calculate_results(self, batch):
        collator = self.data_class.READER.COLLATOR()
        dat = self.model._process_batch(collator(batch).to(self.model.device), 0)
        return self.model(dat, **dat["model_kwargs"])

    def batchify(self, batch):
        cache = []
        for r in batch:
            cache.append(r)
            if len(cache) >= self.batch_size:
                yield cache
                cache = []
        if cache:
            yield cache

    def read_smiles(self, smiles):
        reader = self.data_class.READER()
        d = reader.to_data(dict(features=smiles, labels=None))
        return d

    def predict(self, smiles_list) -> list:
        token_dicts = []
        could_not_parse = []
        index_map = dict()
        for i, smiles in enumerate(smiles_list):
            try:
                # Try to parse the smiles string
                if not smiles:
                    raise ValueError()
                d = self.read_smiles(smiles)
                # This is just for sanity checks
                rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
            except Exception as e:
                # Note if it fails
                could_not_parse.append(i)
                print(f"Failing to parse {smiles} due to {e}")
            else:
                if rdmol is None:
                    could_not_parse.append(i)
                else:
                    index_map[i] = len(token_dicts)
                    token_dicts.append(d)
        results = []
        print(f"Predicting {len(token_dicts), token_dicts} out of {len(smiles_list)}")
        if token_dicts:
            for batch in self.batchify(token_dicts):
                result = self.calculate_results(batch)
                if isinstance(result, dict) and "logits" in result:
                    result = result["logits"]
                results += result.cpu().detach().tolist()
            results = np.stack(results, axis=0)
            positive_preds = [[self.prediction_headers[j] for j, p in enumerate(results[index_map[i]]) if p > 0]
                              if i not in could_not_parse else None for i in range(len(smiles_list))]
            return positive_preds
        else:
            return [None for _ in smiles_list]
