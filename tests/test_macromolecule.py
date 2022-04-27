# Ensure we can access the `molvis` package
import sys
sys.path.append("..")

import os
import shutil
import unittest
from molvis.Macromolecule import Macromolecule
from nglview.widget import NGLWidget
from pathlib import Path
from unittest.mock import patch

class MockPDBList:
    def __init__(self):
        pass

    def retrieve_pdb_file(self, pdbid, pdir=".", file_format="pdb"):
        # Make a copy of the test file, as Macromolecule automatically cleans up after itself
        shutil.copyfile("./fixtures/pdb4hhb.ent", "./pdb4hhb.ent")
        return Path("./pdb4hhb.ent")

class TestMacromolecule(unittest.TestCase):
    @patch('molvis.Macromolecule.PDBList', new=MockPDBList)
    def test_init(self):
        # Smoke test, just ensure no exceptions are thrown
        Macromolecule("4HHB")

    @patch("molvis.Macromolecule.PDBList", new=MockPDBList)
    def test_show(self):
        haemoglobin = Macromolecule("4HHB")
        viewer = haemoglobin.show()
        # Ensure we got a NGLWidget back
        self.assertTrue(issubclass(type(viewer), NGLWidget))
        
        # This time, pass the viewer we just got, as an argument.
        # We shouldn't return anything
        nothing = haemoglobin.show(viewer)
        self.assertIsNone(nothing)

    @patch("molvis.Macromolecule.PDBList", new=MockPDBList)
    def test_save(self):
        haemoglobin = Macromolecule("4HHB")
        
        # Provide filepath as a string
        path = haemoglobin.save("./haemoglobin_str.pdb")
        # Assert with get a subclass of Path back
        # (Will be PosixPath or WindowsPath, depending on the OS)
        self.assertTrue(issubclass(type(path), Path))
        self.assertEqual(str(path), "haemoglobin_str.pdb")
        self.assertTrue(os.path.exists(path))
        # The file should be identical to the test fixture
        with open(path, 'r') as retrieved_file:
            retrieved_data = retrieved_file.read()
        with open("./fixtures/pdb4hhb.ent", 'r') as fixture_file:
            fixture_data = fixture_file.read()
        self.assertEqual(retrieved_data, fixture_data)
        try:
            os.remove(path)
        except Exception as e:
            print("Error removing temporary file:", path)

        # Now try with a Path instance
        path = haemoglobin.save("./haemoglobin_path.pdb")
        # Assert with get a subclass of Path back
        # (Will be PosixPath or WindowsPath, depending on the OS)
        self.assertTrue(issubclass(type(path), Path))
        self.assertEqual(str(path), "haemoglobin_path.pdb")
        self.assertTrue(os.path.exists(path))
        # The file should be identical to the test fixture
        with open(path, 'r') as retrieved_file:
            retrieved_data = retrieved_file.read()
        self.assertEqual(retrieved_data, fixture_data)
        try:
            os.remove(path)
        except Exception as e:
            print("Error removing temporary file:", path)

if __name__ == "__main__":
    unittest.main()