import base64
import io
import queue

from flask_restful import Api, Resource, reqparse
from PIL import Image

from chemlog.solving_strategies.strategy import IsA
from chemlog.preprocessing.peptide_reader import PeptideCountMinAARReader
from chemlog.solving_strategies.peptide_chunk_strategy import PeptideChunkStrategy
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from api.chebi_utils import PREDICTION_HEADERS, LABEL_HIERARCHY, CHEBI_FRAGMENT


PEPTIDE_CLASSIFIER = IsA(limit_to_superclasses=None, available_strategies=[
    PeptideChunkStrategy(reader=PeptideCountMinAARReader(only_one_aar_per_chunk=True), use_running_allocations=False)
])

PEPTIDE_HIERARCHY = CHEBI_FRAGMENT.subgraph([f"CHEBI:{ident}" for ident in PEPTIDE_CLASSIFIER.strategies_by_class.keys()])
print(PEPTIDE_HIERARCHY.nodes)
print(PEPTIDE_HIERARCHY.edges)
def predict_peptides(smiles_list):
    all_preds = []
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            all_preds.append(None)
            continue
        mol.UpdatePropertyCache()
        preds = {}

        q = queue.Queue()
        for node in PEPTIDE_HIERARCHY.nodes:
            if PEPTIDE_HIERARCHY.out_degree(node) == 0:
                q.put(node)
                print(f"put init {node}")

        while not q.empty():
            supercls = q.get()
            preds[supercls] = PEPTIDE_CLASSIFIER(mol, int(supercls.split(":")[1]))
            for child in PEPTIDE_HIERARCHY.predecessors(supercls):
                if all(parent in preds and preds[parent] in [0, 4] for parent in PEPTIDE_HIERARCHY.successors(child)):
                    q.put(child)
        all_preds.append(preds)

    return all_preds


class PredictionDetailChemlog(Resource):
    def load_image(self, path):
        im = Image.open(path)
        data = io.BytesIO()
        im.save(data, "PNG")
        encoded_img_data = base64.b64encode(data.getvalue())
        return encoded_img_data.decode("utf-8")

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
        preds = predict_peptides([smiles])


        return {
            "smiles": smiles,
            "highlights": self.build_explain_blocks(highlights, preds),
        }


if __name__ == "__main__":
    print("CHEBI:16670" in CHEBI_FRAGMENT.nodes)

    print(len(list(CHEBI_FRAGMENT.nodes)))
    preds = predict_peptides(["C[C@H]([NH3+])C(=O)NCC([O-])=O"])
    print(preds)