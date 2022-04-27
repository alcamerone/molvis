import nglview
import os
from Bio.PDB import PDBList
from nglview.widget import NGLWidget
from nglview.show import StringIO
from pathlib import Path
from typing import Union

class Macromolecule:
    def __init__(self, pdbid: str):
        pdbl = PDBList()
        file_path = pdbl.retrieve_pdb_file(pdbid, pdir='.', file_format='pdb')
        with open(file_path, 'r') as pdb_file:
            # Design decision: copying the file into memory makes it easier to tidy up after ourselves,
            # as we can delete the downloaded file once we're done.
            # However, for larger molecules on machines with less memory, it may be advantageous
            # to leave the files on disk.
            # I've opted for the tidier version for now.
            self.pdb_data = pdb_file.read()

        try:
            # Tidy up downloaded file
            os.remove(file_path)
            # Remove `obsolete` dir, which is created when `PDBList` is instantiated
            os.rmdir('./obsolete')
        except Exception as e:
            # Not the end of the world if this doesn't work.
            # Don't crash the app, but tell the user about it.
            print("Exception during cleanup:", e)

    def show(self, existing_viewer: Union[NGLWidget, None] = None) -> Union[NGLWidget, None]:
        file_handle = StringIO(self.pdb_data)
        if existing_viewer:
            # If a viewer is provided, add the molecule as a component to that viewer
            existing_viewer.add_component(file_handle, ext='pdb')
            # Do not return the NGLWidget a second time,
            # just add to the provided viewer
            return
        # If no viewer is provided, return the viewer
        return nglview.show_file(file_handle, ext='pdb')

    def save(self, file_path: Union[str, Path]) -> Path:
        # Dump the PDB data we downloaded on instantiation into a file
        with open(file_path, 'w') as outfile:
            outfile.write(self.pdb_data)

        # We must return an instance of `pahtlib.Path`
        if issubclass(type(file_path), Path):
            # If the input is such an instance, return it directly
            return file_path
        # If the input is a string, convert it to a path
        return Path(file_path)