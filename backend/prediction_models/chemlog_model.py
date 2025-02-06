from typing import Optional

from chemlog2.classification.charge_classifier import get_charge_category, ChargeCategories
from chemlog2.classification.peptide_size_classifier import get_n_amino_acid_residues
from chemlog2.classification.proteinogenics_classifier import get_proteinogenic_amino_acids
from chemlog2.classification.substructure_classifier import is_emericellamide, is_diketopiperazine
from chemlog2.cli import resolve_chebi_classes

from rdkit import Chem

from chebi_utils import CHEBI_FRAGMENT, get_transitive_predictions
from prediction_models.base import PredictionModel

AA_DICT = {
        "A": "L-alanine",
        "C": "L-cysteine",
        "D": "L-aspartic acid",
        "E": "L-glutamic acid",
        "F": "L-phenylalanine",
        "G": "glycine",
        "H": "L-histidine",
        "I": "L-isoleucine",
        "K": "L-lysine",
        "L": "L-leucine",
        "M": "L-methionine",
        "fMet": "N-formylmethionine",
        "N": "L-asparagine",
        "O": "L-pyrrolysine",
        "P": "L-proline",
        "Q": "L-glutamine",
        "R": "L-arginine",
        "S": "L-serine",
        "T": "L-threonine",
        "U": "L-selenocysteine",
        "V": "L-valine",
        "W": "L-tryptophan",
        "Y": "L-tyrosine",
    }

class ChemLog(PredictionModel):

    def __init__(self, name: Optional[str] = None):
        super().__init__(name)

    @property
    def default_name(self):
        return "ChemLog Peptides"

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
            charge_category = get_charge_category(mol)
            n_amino_acid_residues, add_output = get_n_amino_acid_residues(mol)
            r = {"charge_category": charge_category.name, "n_amino_acid_residues": n_amino_acid_residues}
            if n_amino_acid_residues == 5:
                r["emericellamide"] = is_emericellamide(mol)[0]
            if n_amino_acid_residues == 2:
                r["2,5-diketopiperazines"] = is_diketopiperazine(mol)[0]

            chebi_classes = [f"CHEBI:{c}" for c in resolve_chebi_classes(r)]

            all_preds.append(chebi_classes)
        return all_preds

    def get_chemlog_result_info(self, smiles):
        """Get classification for single molecule with additional information."""
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None or not smiles:
            return {"error": "Failed to parse SMILES"}
        mol.UpdatePropertyCache()
        try:
            Chem.Kekulize(mol)
        except Chem.KekulizeException as e:
            pass

        charge_category = get_charge_category(mol)
        n_amino_acid_residues, add_output = get_n_amino_acid_residues(mol)
        if n_amino_acid_residues > 1:
            proteinogenics, proteinogenics_locations, _ = get_proteinogenic_amino_acids(
                mol,
                add_output["amino_residue"],
                add_output["carboxy_residue"])
        else:
            proteinogenics, proteinogenics_locations, _ = [], [], []
        results = {
            'charge_category': charge_category.name,
            'n_amino_acid_residues': n_amino_acid_residues,
            'proteinogenics': proteinogenics,
            'proteinogenics_locations': proteinogenics_locations,
        }

        if n_amino_acid_residues == 5:
            emericellamide = is_emericellamide(mol)
            results["emericellamide"] = emericellamide[0]
            if emericellamide[0]:
                results["emericellamide_atoms"] = emericellamide[1]
        if n_amino_acid_residues == 2:
            diketopiperazine = is_diketopiperazine(mol)
            results["2,5-diketopiperazines"] = diketopiperazine[0]
            if diketopiperazine[0]:
                results["2,5-diketopiperazines_atoms"] = diketopiperazine[1]
        print(results)

        return {**results, **add_output}

    def predict(self, smiles_list):
        return [get_transitive_predictions([positive_i]) for positive_i in self.get_chemlog_results(smiles_list)]

    def build_explain_blocks_atom_allocations(self, atoms, cls_name):
        return [
            ("heading", cls_name),
            ("text", f"The peptide has been identified as an instance of '"
                     f"{cls_name}'. This was decided based on the presence of the following structure:"),
            ("single", atoms)
        ]

    def build_explain_blocks_peptides(self, info):
        blocks = []
        if "error" in info:
            blocks.append(("text", f"An error occurred while processing the molecule: {info['error']}"))
            return blocks
        blocks.append(("heading", "Functional groups"))
        if len(info["amide_bond"]) == 0:
            blocks.append(("text", "The molecule does not contain any amide. Therefore, it cannot be a peptide, "
                                   "peptide anion, peptide zwitterion or peptide cation."))
            return blocks
        blocks.append(("text", "The molecule contains the following functional groups:"))
        blocks.append(("tabs", {"Amide": info["amide_bond"],
                                "Carboxylic acid derivative": info["carboxy_residue"],
                                "Amino group": [[a] for a in info["amino_residue"]]}))
        blocks.append(("heading", "Identifying the peptide structure"))
        if len(info["chunks"]) == 0:
            blocks.append(("text", "All atoms in the molecule are connected via a chain of carbon atoms. "
                                   "Therefore, the molecule cannot be a peptide, peptide anion, peptide zwitterion "
                                   "or peptide cation."))
            return blocks
        blocks.append(("text", "To divide up the molecule into potential amino acids, it has been split into the "
                               f"{len(info['chunks'])} 'building blocks' (based on heteroatoms)."))
        blocks.append(("text", "For each, we have checked if it constitutes an amino acid residue."))
        if len(info["chunks"]) == len(info["longest_aa_chain"]):
            blocks.append(("text", "All chunks have been identified as amino acid residues that are connected "
                                   "via amide bonds:"))
            blocks.append(("tabs", {"Amino acid residue": info["longest_aa_chain"]}))
        elif len(info["longest_aa_chain"]) == 0:
            blocks.append(("tabs", {"Chunks": info["chunks"]}))
            blocks.append(
                ("text", "In these chunks, no amino acids have been identified. "
                         "Therefore, the molecule cannot be a peptide, "
                         "peptide anion, peptide zwitterion or peptide cation."))
            return blocks
        else:
            blocks.append(("text", f"{len(info['longest_aa_chain'])} of these chunks have been identified as amino acid "
                                   f"residues and are connected via amide bonds:"))
            blocks.append(("tabs", {"Chunks": info["chunks"],
                                    "Amino acid residue": info["longest_aa_chain"]}))
        if len(info["longest_aa_chain"]) < 2:
            blocks.append(("text", "Only one amino acid has been identified. Therefore, the molecule cannot be a "
                                   "peptide, peptide anion, peptide zwitterion or peptide cation."))
            return blocks

        blocks.append(("heading", "Charge-based classification"))
        if info["charge_category"] == "SALT":
            blocks.append(("text", "The molecule consists of disconnected anionic and cationic fragments. "
                                   "Therefore, we classify it as a peptide salt. Since there is no class 'peptide salt'"
                                   "in ChEBI, no prediction is made."))
            return blocks
        elif info["charge_category"] == "CATION":
            blocks.append(("text", "The molecule has a net positive charge, therefore it is a 'peptide cation'."))
            return blocks
        elif info["charge_category"] == "ANION":
            blocks.append(("text", "The molecule has a net negative charge, therefore it is a 'peptide anion'."))
            return blocks
        elif info["charge_category"] == "ZWITTERION":
            blocks.append(("text", "The molecule is overall neutral, but a zwitterion, i.e., it contains connected "
                                   "(but non-adjacent) atoms with opposite charges."))
            if info["n_amino_acid_residues"] == 2:
                blocks.append(("text", "Since we have identified 2 amino acid residues, the final classification is "
                                       "'dipeptide zwitterion'."))
            if info["n_amino_acid_residues"] == 3:
                blocks.append(("text", "Since we have identified 3 amino acid residues, the final classification is "
                                       "'tripeptide zwitterion'."))
            return blocks
        subclasses_dict = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "oligo", 7: "oligo", 8: "oligo", 9: "oligo",
                           10: "poly"}
        blocks.append(("text", "The molecule is overall neutral and not a zwitterion. Therefore, it is a peptide."))
        blocks.append(("text", f"More specifically, since we have identified "
                               f"{info["n_amino_acid_residues"]} amino acid residues,"
                               f"the final classification is '{subclasses_dict[min(10, info["n_amino_acid_residues"])]}peptide'."))
        return blocks

    def build_explain_blocks_proteinogenics(self, proteinogenics, atoms):
        blocks = [("heading", "Proteinogenic amino acids")]
        if len(proteinogenics) == 0:
            blocks.append(("text", "No proteinogenic amino acids have been identified."))
            return blocks
        blocks.append(("text", "In addition to the classification, we have searched for the residues of 23 "
                           "proteinogenic amino acids in the molecule."))
        blocks.append(("text", "The following proteinogenic amino acids have been identified:"))
        proteinogenics_dict = {AA_DICT[aa]: [] for aa in proteinogenics}
        for aa, atoms_aa in zip(proteinogenics, atoms):
            proteinogenics_dict[AA_DICT[aa]].append(atoms_aa)
        blocks.append(("tabs", proteinogenics_dict))
        return blocks

    def explain(self, smiles) -> dict:
        info = self.get_chemlog_result_info(smiles)
        highlight_blocks = self.build_explain_blocks_peptides(info)

        for chebi_id, internal_name in [(64372, "emericellamide"), (65061, "2,5-diketopiperazines")]:
            chebi_id_long = f"CHEBI:{chebi_id}"
            if f"{internal_name}_atoms" in info:
                label = CHEBI_FRAGMENT.nodes[chebi_id_long]["lbl"]
                highlight_blocks += self.build_explain_blocks_atom_allocations(info[f"{internal_name}_atoms"], label)
        highlight_blocks += self.build_explain_blocks_proteinogenics(info["proteinogenics"], info["proteinogenics_locations"])
        return {
            "smiles": smiles,
            "highlights": highlight_blocks,
        }


