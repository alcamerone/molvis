import nglview
import os
from Bio.PDB import PDBList
from chembl_webresource_client.new_client import new_client as chembl
from nglview.component import ComponentViewer
from nglview.show import StringIO
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Union

class Macromolecule:
    def __init__(self, pdbid: str):
        pdbl = PDBList()
        # TODO: Clean up after ourselves
        file_path = pdbl.retrieve_pdb_file(pdbid, pdir='.', file_format='pdb')
        self.molecule = Chem.MolFromPDBFile(file_path)
        self.molecule = Chem.AddHs(self.molecule)
        AllChem.EmbedMolecule(self.molecule)

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> ComponentViewer:
        if existing_viewer:
            file_handle = StringIO(Chem.MolToPDBBlock(self.molecule))
            existing_viewer.add_component(file_handle, ext='pdb')
            return existing_viewer
        return nglview.show_rdkit(self.molecule)

    def save(self, file_path: Union[str, Path]) -> Path:
        # `Path`s are not accepted as an argument, so we need to ensure the arg is a string
        writer = Chem.PDBWriter(str(file_path))
        writer.write(self.molecule)
        writer.close()
        if type(file_path) == 'Path':
            return file_path
        return Path(file_path)

class ChemicalMolecule:
    def __init__(self, chemblid: str):
        mol_structs = chembl.molecule.filter(chembl_id=chemblid).only(['molecule_structures'])
        tmp_file_name = f'tmp_{chemblid}.rdkit'
        # This file write is necessary because the argument to `MolFromMolFile` *must* be a string
        # representing a file name, a file handle is not allowed.
        tmp = open(tmp_file_name, 'w')
        tmp.write(mol_structs[0]['molecule_structures']['molfile'])
        tmp.close()
        self.molecule = Chem.MolFromMolFile(tmp_file_name)
        os.remove(tmp_file_name)
        self.molecule = Chem.AddHs(self.molecule)
        AllChem.EmbedMolecule(self.molecule)

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> ComponentViewer:
        if existing_viewer:
            file_handle = StringIO(Chem.MolToPDBBlock(self.molecule))
            existing_viewer.add_component(file_handle, ext='pdb')
            return existing_viewer
        return nglview.show_rdkit(self.molecule)

    def save(self, file_path: Union[str, Path]) -> Path:
        # `Path`s are not accepted as an argument, so we need to ensure the arg is a string
        writer = Chem.SDWriter(str(file_path))
        writer.write(self.molecule)
        writer.close()
        if type(file_path) == 'Path':
            return file_path
        return Path(file_path)
