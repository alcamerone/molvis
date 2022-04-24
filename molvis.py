from msilib.schema import Component
from typing import Union
from pathlib import Path
from nglview.component import ComponentViewer
from Bio.PDB import PDBList
from chembl_webresource_client.new_client import new_client as chembl

class Macromolecule:
    def __init__(self, pdbid: str) -> Macromolecule:
        pdbl = PBDList()
        self.pdb_file = pdbl.retrieve_pdb_file(pdbid)

    def show(self, existing_viewer: Union[ComponentViewer, None] = None):
        pass

    def save(self, file_path: Union[str, Path]) -> Path:
        pass

class ChemicalMolecule:
    def __init__(self, chemblid: str) -> ChemicalMolecule:
        molecule = chembl.molecule
        m1 = molecule.filter(chembl_id=chemblid).only(['molecule_structures'])
        self.structure = m1[0].molecule_structures.molfile

    def show(self, existing_viewer: Union[ComponentViewer, None] = None):
        pass

    def save(self, file_path: Union[str, Path]) -> Path:
        pass
