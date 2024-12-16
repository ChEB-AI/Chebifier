import abc
from typing import Optional


class PredictionModel(abc.ABC):


    def __init__(self, name: Optional[str] = None):
        self._name = name or self.default_name

    @property
    def default_name(self) -> str:
        return "PredictionModel"

    @property
    def name(self) -> str:
        return self._name

    @property
    def info_text(self) -> str:
        return "A prediction model."

    def predict(self, smiles_list) -> list:
        pass


    def explain(self, smiles) -> dict:
        pass


