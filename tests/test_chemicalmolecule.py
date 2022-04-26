import os
import sys
import unittest
from unittest.mock import patch
sys.path.append("..")
from molvis.ChemicalMolecule import ChemicalMolecule
from nglview.widget import NGLWidget
from pathlib import Path

sildenafil_rdkit = """
     RDKit          2D

 33 36  0  0  0  0  0  0  0  0999 V2000
    2.1000   -0.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1000    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5375   -0.0042    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.4917   -0.3667    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.8792   -0.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8042    0.9083    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4917    1.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8792    0.6833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.2042    0.3458    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.8042   -0.2417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2875   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1583   -0.3750    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9333   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3208   -0.0333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1875    0.6083    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8958    0.6083    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3958   -1.0917    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7833   -0.0042    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1583   -1.0917    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2875   -1.1125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4917    1.7708    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9333   -1.1125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3208   -1.4542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3958   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7833   -1.4417    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0750    1.5750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8042   -0.9500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8792   -1.4542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9958   -1.4292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4958   -1.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4167   -1.3125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1125   -1.4500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0375   -0.9542    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  2  0
  3 13  1  0
  4  1  1  0
  5  4  2  0
  6  2  1  0
  7  2  1  0
  8  5  1  0
  9 10  2  0
 10  1  1  0
 11  5  1  0
 12  3  1  0
 13 14  2  0
 14 11  1  0
 15  3  2  0
 16  3  2  0
 17 25  1  0
 18 12  1  0
 19 12  1  0
 20 11  2  0
 21  7  2  0
 22 23  2  0
 23 20  1  0
 24 18  1  0
 25 19  1  0
 26  6  1  0
 27 10  1  0
 28 20  1  0
 29 17  1  0
 30 28  1  0
 31 27  1  0
 32 30  1  0
 33 31  1  0
  9  6  1  0
  8  7  1  0
 22 13  1  0
 17 24  1  0
M  END

> <chembl_id>
CHEMBL192

> <chembl_pref_name>
SILDENAFIL


"""

class MockChemblReturn:
    def filter(chembl_id):
        return MockChemblReturn()
    
    def only(self, attrs):
        return [
            {
                'molecule_structures': {
                    'molfile': sildenafil_rdkit
                }
            }
        ]

class TestChemicalMolecule(unittest.TestCase):
    @patch('molvis.ChemicalMolecule.chembl.molecule', new=MockChemblReturn)
    def test_init(self):
        # Smoke test, just ensure no exceptions are thrown
        ChemicalMolecule("CHEMBL192")
    
    @patch('molvis.ChemicalMolecule.chembl.molecule', new=MockChemblReturn)
    def test_show(self):
        sildenafil = ChemicalMolecule("CHEMBL192")
        viewer = sildenafil.show()
        # Ensure we got a NGLWidget back
        self.assertTrue(issubclass(type(viewer), NGLWidget))
        
        # This time, pass the viewer we just got, as an argument.
        # We shouldn't return anything
        nothing = sildenafil.show(viewer)
        self.assertIsNone(nothing)
    
    @patch('molvis.ChemicalMolecule.chembl.molecule', new=MockChemblReturn)
    def test_save(self):
        sildenafil = ChemicalMolecule("CHEMBL192")
        
        # Provide filepath as a string
        path = sildenafil.save("./sildenafil_str.sdf")
        # Assert with get a subclass of Path back
        # (Will be PosixPath or WindowsPath, depending on the OS)
        self.assertTrue(issubclass(type(path), Path))
        self.assertEqual(str(path), "sildenafil_str.sdf")
        self.assertTrue(os.path.exists(path))
        # Ensure the correct file is produced
        with open(path, 'r') as saved_file:
            saved_data = saved_file.read()
        with open("./fixtures/chembl192.sdf", 'r') as fixture_file:
            fixture_data = fixture_file.read()
        self.assertEqual(saved_data, fixture_data)
        try:
            os.remove(path)
        except Exception as e:
            print("Error removing temporary file:", path)

        # Now try with a Path instance
        path = sildenafil.save("./sildenafil_path.sdf")
        # Assert with get a subclass of Path back
        # (Will be PosixPath or WindowsPath, depending on the OS)
        self.assertTrue(issubclass(type(path), Path))
        self.assertEqual(str(path), "sildenafil_path.sdf")
        self.assertTrue(os.path.exists(path))
        with open(path, 'r') as saved_file:
            saved_data = saved_file.read()
        with open("./fixtures/chembl192.sdf", 'r') as fixture_file:
            fixture_data = fixture_file.read()
        self.assertEqual(saved_data, fixture_data)
        try:
            os.remove(path)
        except Exception as e:
            print("Error removing temporary file:", path)

if __name__ == "__main__":
    unittest.main()