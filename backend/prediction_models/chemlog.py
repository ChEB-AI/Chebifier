import queue
from typing import Optional

from chemlog.preprocessing.peptide_reader import PeptideCountMinAARReader
from chemlog.solving_strategies.peptide_chunk_strategy import PeptideChunkStrategy
from chemlog.solving_strategies.strategy import IsA
from rdkit import Chem

from chebi_utils import CHEBI_FRAGMENT
from prediction_models.base import PredictionModel

PEPTIDE_CLASSIFIER = IsA(limit_to_superclasses=None, available_strategies=[
    PeptideChunkStrategy(reader=PeptideCountMinAARReader(only_one_aar_per_chunk=True),
                         use_running_allocations=False)
])


class ChemLog(PredictionModel):

    def __init__(self, name: Optional[str] = None, strategy_collection: Optional[IsA] = None):
        super().__init__(name)
        if strategy_collection is None:
            strategy_collection = PEPTIDE_CLASSIFIER
        self.is_a = strategy_collection
        self.peptide_hierarchy = CHEBI_FRAGMENT.subgraph(
            [f"CHEBI:{ident}" for ident in self.is_a.strategies_by_class.keys()])

    @property
    def default_name(self):
        return "ChemLog"

    @property
    def info_text(self):
        return "A rule-based model for predicting peptides and peptide-like molecules, see [2]."

    def get_chemlog_results(self, smiles_list) -> list:
        all_preds = []
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None or not smiles:
                all_preds.append(None)
                continue
            mol.UpdatePropertyCache()
            preds = {}

            q = queue.Queue()
            for node in self.peptide_hierarchy.nodes:
                if self.peptide_hierarchy.out_degree(node) == 0:
                    q.put(node)
            while not q.empty():
                supercls = q.get()
                preds[supercls] = PEPTIDE_CLASSIFIER(mol, int(supercls.split(":")[1]))
                for child in self.peptide_hierarchy.predecessors(supercls):
                    if all(parent in preds and preds[parent] in [0, 4] for parent in
                           self.peptide_hierarchy.successors(child)):
                        q.put(child)
            all_preds.append(preds)
        return all_preds

    def predict(self, smiles_list):
        # only return positive predictions
        results = self.get_chemlog_results(smiles_list)
        return [[cls for cls, pred in result.items() if pred in [0, 4]] if result is not None else None for result in results]

    def build_explain_blocks(self, highlights, preds):
        blocks = []
        blocks.append(("text", f"The molecule has been analyzed for peptide properties. The full classification "
                               f"results are {preds}."))
        if not highlights["is_connected"]:
            blocks.append(("text", "The molecule is not connected, so it cannot be a peptide."))
            return blocks
        if len(highlights["excluded_atoms"]) > 0:
            blocks.append(("text", "The molecule contains atoms that are parts of 'functional units' "
                                   "(e.g. Coenzyme A, penicillin). "
                                   "The following atoms will be ignored:"))
            blocks.append(("single", highlights["excluded_atoms"]))
        if "carboxy_groups" not in highlights:
            blocks.append(("text", "The molecule does not contain any amide. Therefore, it cannot be a peptide."))
            return blocks
        blocks.append(("text", "The molecule contains the following functional groups:"))
        blocks.append(("tabs", {"Amide bonds": highlights["amide_bonds"],
                                "Carboxy groups / residues": highlights["carboxy_groups"],
                                "Amino groups / residues": highlights["amino_groups"]}))
        blocks.append(("text", "To divide up the molecule into potential amino acids, it has been split into the "
                               "following 'chunks' (based on amide bonds and disulfide bridges):"))
        blocks.append(("tabs", {"Chunks": highlights["chunks"]}))
        if len(highlights["aars"]) == 0:
            blocks.append(
                ("text", "No amino acids have been identified. Therefore, the molecule cannot be a peptide, "
                         "peptide anion, zwitterion or cation."))
            return blocks
        # todo check which amino acids actually form the peptide
        blocks.append(("text", "The following amino acids have been identified:"))
        blocks.append(("tabs", {"Amino acids": highlights["aars"]}))
        if highlights["net_charge"] > 0:
            blocks.append(("text", "The molecule has a net positive charge, therefore it is a cation."))
            return blocks
        if highlights["net_charge"] < 0:
            blocks.append(("text", "The molecule has a net negative charge, therefore it is an anion."))
            return blocks
        if len(highlights["zwitterion_atoms"]) > 0:
            blocks.append(("text", "The molecule is overall neutral, but a zwitterion. "
                                   "The following atoms are responsible:"))
            blocks.append(("single", highlights["zwitterion_atoms"]))
            return blocks
        blocks.append(("text", "The molecule is overall neutral and not a zwitterion. Therefore, it is a peptide."))
        return blocks
        # todo more classes

    def explain(self, smiles) -> dict:
        highlights_reader = PeptideCountMinAARReader(only_one_aar_per_chunk=True, return_highlights=True)
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache()
        domain_length, dense_extensions, highlights = highlights_reader.read_mol_as_subcls(mol)
        preds = self.get_chemlog_results([smiles])

        return {
            "smiles": smiles,
            "highlights": self.build_explain_blocks(highlights, preds),
        }
