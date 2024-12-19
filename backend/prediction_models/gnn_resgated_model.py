from typing import Union, Optional

import chebai_graph.preprocessing.properties as p
import torch
from chebai_graph.models.graph import ResGatedGraphConvNetGraphPred
from chebai_graph.preprocessing.datasets.chebi import GraphPropertiesMixIn, ChEBI100GraphProperties, \
    ChEBI50GraphProperties
from chebai_graph.preprocessing.property_encoder import IndexEncoder, OneHotEncoder
from torch_geometric.data.data import Data as GeomData

from prediction_models.base import NNPredictionModel

if torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"


class GNNResGated(NNPredictionModel):

    def __init__(self, checkpoint_path: str, data_class: Union[GraphPropertiesMixIn, str],
                 prediction_headers_path: str, batch_size: Optional[int] = 32, name: Optional[str] = None):
        super().__init__(prediction_headers_path, batch_size, name)
        self.model = ResGatedGraphConvNetGraphPred.load_from_checkpoint(
            checkpoint_path, map_location=torch.device(device), criterion=None, strict=False,
            metrics=dict(train=dict(), test=dict(), validation=dict()), pretrained_checkpoint=None,
            config={"in_length": 256, "hidden_length": 512, "dropout_rate": 0.1, "n_conv_layers": 3,
                    "n_linear_layers": 3, "n_atom_properties": 158, "n_bond_properties": 7,
                    "n_molecule_properties": 200})
        self.model.eval()
        if isinstance(data_class, str):
            properties = [p.AtomType(), p.NumAtomBonds(), p.AtomCharge(), p.AtomAromaticity(), p.AtomHybridization(),
                          p.AtomNumHs(), p.BondType(), p.BondInRing(), p.BondAromaticity(), p.RDKit2DNormalized()]
            if data_class == "ChEBI100":
                data_class = ChEBI100GraphProperties(properties=properties)
            elif data_class == "ChEBI50":
                data_class = ChEBI50GraphProperties(properties=properties)
            else:
                raise ValueError(f"Unknown data class {data_class}")
        self.data_class = data_class

    @property
    def default_name(self) -> str:
        return "ResGatedGraphConvNet"

    @property
    def info_text(self) -> str:
        return ("Residual gated graph convolutional network (see [3]) for predicting arbitrary ChEBI classes that "
                "have a certain number of instances.")

    def read_smiles(self, smiles):
        reader = self.data_class.READER()
        d = reader.to_data(dict(features=smiles, labels=None))
        geom_data = d["features"]
        edge_attr = geom_data.edge_attr
        x = geom_data.x
        molecule_attr = torch.empty((1, 0))
        for property in self.data_class.properties:
            property_values = reader.read_property(smiles, property)
            encoded_values = []
            for value in property_values:
                # cant use standard encode for index encoder because model has been trained on a certain range of values
                # use default value if we meet an unseen value
                if isinstance(property.encoder, IndexEncoder):
                    if str(value) in property.encoder.cache:
                        index = property.encoder.cache.index(str(value)) + property.encoder.offset
                    else:
                        index = 0
                        print(f"Unknown property value {value} for property {property} at smiles {smiles}")
                    if isinstance(property.encoder, OneHotEncoder):
                        encoded_values.append(torch.nn.functional.one_hot(
                            torch.tensor(index), num_classes=property.encoder.get_encoding_length()
                        ))
                    else:
                        encoded_values.append(torch.tensor([index]))

                else:
                    encoded_values.append(property.encoder.encode(value))
            if len(encoded_values) > 0:
                encoded_values = torch.stack(encoded_values)

            if isinstance(encoded_values, torch.Tensor):
                if len(encoded_values.size()) == 0:
                    encoded_values = encoded_values.unsqueeze(0)
                if len(encoded_values.size()) == 1:
                    encoded_values = encoded_values.unsqueeze(1)
            else:
                encoded_values = torch.zeros(
                    (0, property.encoder.get_encoding_length())
                )
            if isinstance(property, p.AtomProperty):
                x = torch.cat([x, encoded_values], dim=1)
            elif isinstance(property, p.BondProperty):
                edge_attr = torch.cat([edge_attr, encoded_values], dim=1)
            else:
                molecule_attr = torch.cat([molecule_attr, encoded_values[0]], dim=1)

        d["features"] = GeomData(
            x=x,
            edge_index=geom_data.edge_index,
            edge_attr=edge_attr,
            molecule_attr=molecule_attr,
        )
        return d
