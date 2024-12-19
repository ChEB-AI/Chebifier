import queue
from typing import Optional

from chemlog.preprocessing.peptide_reader import PeptideCountMinAARReader
from chemlog.solving_strategies.peptide_chunk_strategy import PeptideChunkStrategy
from chemlog.solving_strategies.strategy import IsA
from chemlog.solving_strategies.depsipeptide_strategy import DepsipeptideStrategy
from chemlog.preprocessing.build_rdkit import ChemDataExtensions
from chemlog.solving_strategies.fixed_fol_strategy import FixedFOLStrategy

from rdkit import Chem

from chebi_utils import CHEBI_FRAGMENT, get_transitive_predictions
from prediction_models.base import PredictionModel

CHEBI_VERSION = 237
PEPTIDE_CLASSIFIER = IsA(limit_to_superclasses=None, available_strategies=[
    FixedFOLStrategy(data_module=ChemDataExtensions(chebi_version=CHEBI_VERSION)),
    PeptideChunkStrategy(reader=PeptideCountMinAARReader(only_one_aar_per_chunk=True)),
    DepsipeptideStrategy(),
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
                preds[supercls] = self.is_a(mol, int(supercls.split(":")[1]))
                for child in self.peptide_hierarchy.predecessors(supercls):
                    if all(parent in preds and preds[parent] in [0, 4] for parent in
                           self.peptide_hierarchy.successors(child)):
                        q.put(child)
            all_preds.append(preds)
        return all_preds

    def predict(self, smiles_list):
        # only return positive predictions
        results = self.get_chemlog_results(smiles_list)
        positive = [[cls for cls, pred in result.items() if pred in [0, 4]] if result is not None else None for result in
                results]
        return [get_transitive_predictions([positive_i]) for positive_i in positive]

    def build_explain_blocks_atom_allocations(self, allocations, cls_name):
        print(allocations)
        return [
            ("heading", cls_name),
            ("text", f"The peptide has been identified as an instance of '"
                     f"{cls_name}'. This was decided based on the presence of the following structure:"),
            ("single", [alloc for var, alloc in allocations])
        ]

    def build_explain_blocks_depsi(self, highlights, preds, allocations):
        blocks = []
        blocks.append(("heading", "Depsipeptide"))
        if not highlights["is_connected"]:
            blocks.append(("text", "The molecule is not connected. Therefore, it cannot be a depsipeptide."))
            return blocks

        blocks.append(("text", "In addition to amides, carboxy and amino groups / residues, the molecule contains the "
                               "following functional groups:"))
        blocks.append(("tabs", {"Ester": highlights["ester_bonds"],
                                "Hydroxy carboxylic group / residue": highlights["carboxy_groups"],
                                }))
        blocks.append(("text", "To divide up the molecule into potential amino acids, it has been split into the "
                               f"{len(highlights['chunks'])} 'chunks' (based on amide bonds, ester bonds and"
                               f" disulfide bridges)."))
        blocks.append(("text", "For each, we have checked if it constitutes an amino acid residue "
                               "or hydroxy carboxylic acid residue."))
        if len(highlights["aars"]) == 0:
            blocks.append(
                ("text", "No amino acid residues have been identified. "
                         "Therefore, the molecule cannot be a depsipeptide."))
            blocks.append(("tabs", {"Chunks": highlights["chunks"]}))
            return blocks
        else:
            blocks.append(("tabs", {"Chunk": highlights["chunks"],
                                    "Amino acid residue": highlights["aars"],
                                    "Hydroxy carboxylic acid residue": highlights["hars"]}))

        if allocations is not None:
            blocks.append(("text", f"Out of the identified residues, we have identified a depsipeptide structure. "
                                   f"Note that this might not be the only possible configuration."
                                   f" Other configurations with the same number of amino acids might be possible. The "
                                   f"search was conducted up to a maximum of 10 amino acid residues."))
            depsipeptide_allocation = [highlights["aars"][alloc - len(highlights["amide_bonds"]) -
                                                          len(highlights["chunks"]) - len(highlights["ester_bonds"])]
                                       for var, alloc in allocations if var.startswith("a")]
            depsipeptide_allocation += [highlights["hars"][alloc - len(highlights["aars"]) -
                                                           len(highlights["amide_bonds"]) - len(highlights["chunks"]) -
                                                           len(highlights["ester_bonds"])]
                                        for var, alloc in allocations if var.startswith("h")]
            blocks.append(("single", depsipeptide_allocation))
        else:
            blocks.append(("text", "The identified residues can not be combined into a depsipeptide via amide or "
                                   "ester bonds."))
        return blocks

    def build_explain_blocks_peptides(self, highlights, preds, allocations):
        blocks = []
        #blocks.append(("text", f"The molecule has been analyzed for peptide properties. The full classification "
        #                       f"results are {preds}."))
        #blocks.append(("text", f"The following allocations have been made: {allocations}"))
        if not highlights["is_connected"]:
            blocks.append(("text", "The molecule is not connected. Therefore, it cannot be a peptide, "
                                   "peptide anion, peptide zwitterion or peptide cation."))
            return blocks
        blocks.append(("heading", "Functional groups"))
        if len(highlights["excluded_atoms"]) > 0:
            blocks.append(("text", "The molecule contains atoms that are parts of 'functional units' "
                                   "(e.g. Coenzyme A, penicillin) which are usually not considered as peptides. "
                                   "The following atoms will be ignored:"))
            blocks.append(("single", highlights["excluded_atoms"]))
        if "carboxy_groups" not in highlights:
            blocks.append(("text", "The molecule does not contain any amide. Therefore, it cannot be a peptide, "
                                   "peptide anion, peptide zwitterion or peptide cation."))
            return blocks
        blocks.append(("text", "The molecule contains the following functional groups:"))
        blocks.append(("tabs", {"Amide": highlights["amide_bonds"],
                                "Carboxy group / residue": highlights["carboxy_groups"],
                                "Amino group / residue": highlights["amino_groups"]}))
        blocks.append(("heading", "Identifying the peptide structure"))
        blocks.append(("text", "To divide up the molecule into potential amino acids, it has been split into the "
                               f"{len(highlights['chunks'])} 'chunks' (based on amide bonds and disulfide bridges)."))
        blocks.append(("text", "For each, we have checked if it constitutes an amino acid residue."))
        if len(highlights["chunks"]) == len(highlights["aars"]):
            blocks.append(("text", "All chunks have been identified as amino acid residues:"))
            blocks.append(("tabs", {"Amino acid residue": highlights["aars"]}))
        elif len(highlights["aars"]) == 0:
            blocks.append(
                ("text", "In these chunks, no amino acids have been identified. "
                         "Therefore, the molecule cannot be a peptide, "
                         "peptide anion, peptide zwitterion or peptide cation."))
            blocks.append(("tabs", {"Chunks": highlights["chunks"]}))
            return blocks
        else:
            blocks.append(("text", f"{len(highlights['aars'])} chunks have been identified as amino acid residues:"))
            blocks.append(("tabs", {"Chunks": highlights["chunks"],
                                    "Amino acid residue": highlights["aars"]}))
        if len(highlights["aars"]) == 0:
            blocks.append(
                ("text", "No amino acids have been identified. Therefore, the molecule cannot be a peptide, "
                         "peptide anion, peptide zwitterion or peptide cation."))
            return blocks
        if len(highlights["aars"]) < 2:
            blocks.append(("text", "Only one amino acid has been identified. Therefore, the molecule cannot be a "
                                   "peptide, peptide anion, peptide zwitterion or peptide cation."))
            return blocks
        if allocations is None:
            blocks.append(("text", "The identified amino acid residues can not be combined into a peptide via amide "
                                   "bonds."))
            return blocks
        n_total_aars = (len(allocations) + 1) // 2
        peptide_allocation = [highlights["aars"][alloc - len(highlights["amide_bonds"]) - len(highlights["chunks"])]
                              for var, alloc in allocations if var.startswith("a")]
        peptide_allocation = [p for pa in peptide_allocation for p in pa]
        blocks.append(("text", f"Out of the identified amino acid residues, "
                               f"{'all ' if n_total_aars == len(highlights['aars']) else ''}{n_total_aars} "
                               f"form a peptide structure. Note that this might not be the only possible configuration."
                               f" Other configurations with the same number of amino acids might be possible. The "
                               f"search was conducted up to a maximum of 10 amino acid residues. The complete "
                               f"identified structure is the following:"))
        blocks.append(("single", peptide_allocation))

        blocks.append(("heading", "Charge-based classification"))
        if highlights["net_charge"] > 0:
            blocks.append(("text", "The molecule has a net positive charge, therefore it is a peptide cation."))
            return blocks
        if highlights["net_charge"] < 0:
            blocks.append(("text", "The molecule has a net negative charge, therefore it is a peptide anion."))
            return blocks
        if len(highlights["zwitterion_atoms"]) > 0:
            blocks.append(("text", "The molecule is overall neutral, but a zwitterion. "
                                   "The following atoms are responsible:"))
            blocks.append(("single", highlights["zwitterion_atoms"]))
            if n_total_aars == 2:
                blocks.append(("text", "Since we have identified 2 amino acid residues, the final classification is "
                                       "'dipeptide zwitterion'."))
            if n_total_aars == 3:
                blocks.append(("text", "Since we have identified 3 amino acid residues, the final classification is "
                                       "'tripeptide zwitterion'."))
            return blocks
        subclasses_dict = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "oligo", 7: "oligo", 8: "oligo", 9: "oligo",
                           10: "poly"}
        blocks.append(("text", "The molecule is overall neutral and not a zwitterion. Therefore, it is a peptide."))
        blocks.append(("text", f"More specifically, since we have identified "
                               f"{n_total_aars if n_total_aars < 10 else str(n_total_aars) + '+'} amino acid residues,"
                               f"the final classification is '{subclasses_dict[n_total_aars]}peptide'."))
        return blocks

    def explain(self, smiles) -> dict:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache()

        highlights_reader = PeptideCountMinAARReader(only_one_aar_per_chunk=True, return_highlights=True)
        domain_length, dense_extensions, highlights = highlights_reader.read_mol_as_subcls(mol)
        preds = self.get_chemlog_results([smiles])
        allocs_peptide = self.is_a.get_strategy(16670).last_allocations
        highlight_blocks = self.build_explain_blocks_peptides(highlights, preds, allocs_peptide)

        if self.is_a.get_strategy(23643) is not None:
            allocs_depsi = self.is_a.get_strategy(23643).last_allocations
            highlights_reader_depsi = PeptideCountMinAARReader(only_one_aar_per_chunk=True, return_highlights=True,
                                                           depsi_mode=True)
            _, _, highlights_depsi = highlights_reader_depsi.read_mol_as_subcls(mol)
            highlight_blocks += self.build_explain_blocks_depsi(highlights_depsi, preds, allocs_depsi)

        for chebi_id in [64372, 65061]:
            chebi_id_long = f"CHEBI:{chebi_id}"
            if chebi_id_long in preds[0] and preds[0][chebi_id_long] in [0, 4]:
                label = CHEBI_FRAGMENT.nodes[chebi_id_long]["lbl"]
                strategy = self.is_a.get_strategy(chebi_id)
                # run strategy again to make sure that last_allocations correspond to the chebi_id
                strategy.set_mol(mol)
                strategy(chebi_id)
                highlight_blocks += self.build_explain_blocks_atom_allocations(strategy.last_allocations, label)


        return {
            "smiles": smiles,
            "highlights": highlight_blocks,
        }
