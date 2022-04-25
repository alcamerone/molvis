import nglview
from chembl_webresource_client.new_client import new_client as chembl
from nglview.component import ComponentViewer
from nglview.show import StringIO
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Union

class ChemicalMolecule:
    def __init__(self, chemblid: str):
        mol_structs = chembl.molecule.filter(chembl_id=chemblid).only(['molecule_structures'])
        self.molecule = Chem.MolFromMolBlock(mol_structs[0]['molecule_structures']['molfile'])
        self.molecule = Chem.AddHs(self.molecule)
        AllChem.EmbedMolecule(self.molecule)

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> Union[ComponentViewer, None]:
        if existing_viewer:
            file_handle = StringIO(Chem.MolToPDBBlock(self.molecule))
            existing_viewer.add_component(file_handle, ext='pdb')
            # Do not return the ComponentViewer a second time,
            # just add to the provided viewer
            return
        return nglview.show_rdkit(self.molecule)

    def save(self, file_path: Union[str, Path]) -> Path:
        # `Path`s are not accepted as an argument, so we need to ensure the arg is a string
        writer = Chem.SDWriter(str(file_path))
        writer.write(self.molecule)
        writer.close()
        if type(file_path) == 'Path':
            return file_path
        return Path(file_path)