import base64
import io

from flask_restful import Api, Resource, reqparse
from PIL import Image

from chemlog.solving_strategies.strategy import IsA
from chemlog.preprocessing.peptide_reader import PeptideCountMinAARReader
from chemlog.solving_strategies.peptide_chunk_strategy import PeptideChunkStrategy
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


PEPTIDE_CLASSIFIER = IsA(limit_to_superclasses=None, available_strategies=[
    PeptideChunkStrategy(reader=PeptideCountMinAARReader(only_one_aar_per_chunk=True), use_running_allocations=False)
])

def predict_single_peptide(mol, cls):
    outcome = PEPTIDE_CLASSIFIER(mol, cls)
    return outcome in [0, 4]
def predict_peptides(smiles_list):
    superclasses = [16670, 25676, 46761, 47923]
    all_preds, all_direct = [], []
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            all_preds.append(None)
            all_direct.append(None)
            continue
        mol.UpdatePropertyCache()
        preds = []
        preds_direct = []

        # peptide zwitterion
        if predict_single_peptide(mol, 60466):
            preds.append(60466)
            if predict_single_peptide(mol, 90799):
                preds.append(90799)
                preds_direct.append(90799)
            elif predict_single_peptide(mol, 155837):
                preds.append(155837)
                preds_direct.append(155837)
            else:
                preds_direct.append(60466)
        # peptide anion
        elif predict_single_peptide(mol, 60194):
            preds.append(60194)
            preds_direct.append(60194)
        # peptide cation
        elif predict_single_peptide(mol, 60334):
            preds.append(60334)
            preds_direct.append(60334)
        # peptide
        elif predict_single_peptide(mol, 16670):
            preds.append(16670)
            # oligo
            if predict_single_peptide(mol, 25676):
                preds.append(25676)
                # di
                if predict_single_peptide(mol, 46761):
                    preds.append(46761)
                    # 2,5-diketopiperazines
                    if predict_single_peptide(mol, 65061):
                        preds.append(65061)
                    else:
                        preds_direct.append(46761)
                # tri
                if predict_single_peptide(mol, 47923):
                    preds.append(47923)
                    preds_direct.append(47923)
                # tetra
                elif predict_single_peptide(mol, 48030):
                    preds.append(48030)
                    preds_direct.append(48030)
                # penta
                elif predict_single_peptide(mol, 48545):
                    preds.append(48545)
                else:
                    preds_direct.append(25676)
            else:
                # poly
                preds.append(15841)
                preds_direct.append(15841)
        # depsipeptide
        if predict_single_peptide(mol, 23643):
            preds.append(23643)

        # emericellamide depends on depsi and penta
        if 23643 in preds and 48545 in preds:
            # emericellamide
            if predict_single_peptide(mol, 64372):
                preds.append(64372)
                preds_direct.append(64372)
            else:
                preds_direct.append(23643)
                preds_direct.append(64372)
        elif 23643 in preds:
            preds_direct.append(23643)
        elif 48545 in preds:
            preds_direct.append(48545)

        all_direct.append(preds_direct)
        all_preds.append(preds)

    return all_direct, all_preds


class PredictionDetailChemlog(Resource):
    def load_image(self, path):
        im = Image.open(path)
        data = io.BytesIO()
        im.save(data, "PNG")
        encoded_img_data = base64.b64encode(data.getvalue())
        return encoded_img_data.decode("utf-8")

    def build_explain_blocks(self, highlights):
        blocks = []
        if not highlights["is_connected"]:
            return [("text", "The molecule is not connected, so it cannot be a peptide.")]
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
            blocks.append(("text", "No amino acids have been identified. Therefore, the molecule cannot be a peptide, "
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



    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument("type", type=str)
        parser.add_argument("smiles", type=str)

        args = parser.parse_args()

        # note, the post req from frontend needs to match the strings here (e.g. 'type and 'message')

        request_type = args["type"]
        smiles = args["smiles"]

        highlights_reader = PeptideCountMinAARReader(only_one_aar_per_chunk=True, return_highlights=True)
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        mol.UpdatePropertyCache()
        domain_length, dense_extensions, highlights = highlights_reader.read_mol_as_subcls(mol)

        return {
            "smiles": smiles,
            "highlights": self.build_explain_blocks(highlights),
        }


if __name__ == "__main__":
    highlights_reader = PeptideCountMinAARReader(only_one_aar_per_chunk=True, return_highlights=True)
    mol = Chem.MolFromSmiles("O=C(N[C@@H](CCCN=C(N)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)C)CCC(O)=O", sanitize=False)
    mol.UpdatePropertyCache()
    print(highlights)