import nglview
import os
from Bio.PDB import PDBList
from nglview.component import ComponentViewer
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

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> Union[ComponentViewer, None]:
        if existing_viewer:
            existing_viewer.add_component(StringIO(self.pdb_data), ext='pdb')
            # Do not return the ComponentViewer a second time,
            # just add to the provided viewer
            return
        return nglview.show_file(StringIO(self.pdb_data), ext='pdb')

    def save(self, file_path: Union[str, Path]) -> Path:
        with open(file_path, 'w') as outfile:
            outfile.write(self.pdb_data)
        if type(file_path) == 'Path':
            return file_path
        return Path(file_path)