import nglview
from chembl_webresource_client.new_client import new_client as chembl
from nglview.widget import NGLWidget
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Union

class ChemicalMolecule:
    def __init__(self, chemblid: str):
        # Retrieve the molecular structure
        mol_structs = chembl.molecule.filter(chembl_id=chemblid).only(['molecule_structures'])
        # Convert to a RDKit molecule
        self.molecule = Chem.MolFromMolBlock(mol_structs[0]['molecule_structures']['molfile'])
        # Add Hs, for a more accurate conformer
        self.molecule = Chem.AddHs(self.molecule)
        # Create conformer
        AllChem.EmbedMolecule(self.molecule)

    def show(self, existing_viewer: Union[NGLWidget, None] = None) -> Union[NGLWidget, None]:
        if existing_viewer:
            # If a viewer is provided, add the molecule as a component to that viewer
            existing_viewer.add_component(self.molecule)
            # Do not return the NGLWidget a second time,
            # just add to the provided viewer
            return
        # If no viewer is provided, return the viewer
        return nglview.show_rdkit(self.molecule)

    def save(self, file_path: Union[str, Path]) -> Path:
        # `Path`s are not accepted as an argument, so we need to ensure the arg is a string
        writer = Chem.SDWriter(str(file_path))
        writer.write(self.molecule)
        writer.close()

        # We must return an instance of `pahtlib.Path`
        if issubclass(type(file_path), Path):
            # If the input is such an instance, return it directly
            return file_path
        # If the input is a string, convert it to a path
        return Path(file_path)