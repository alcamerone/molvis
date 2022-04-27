# molvis

molvis provides utility classes for downloading, displaying (via NGLView), and saving PDB and CHEMBL structures.

The `Macromolecule` class can be used to work with PDB structures as follows:
```
protease = Macromolecule(pdbid="6LU7") # Get the crystal structure of the main COVID-19 protease
viewer = protease.show()
viewer # Will display a 3D representation of `protease`
protease.save("path/to/protease.pdb") # Saves the structure of `protease` in PDB format to `path/to/protease.pdb`
```

The `ChemicalMolecule` class can be used to work with CHEMBL structures:
```
lopinavir = ChemicalMolecule(chemblid="CHEMBL729") # Get the 3D structure of lopinavir
viewer = lopinavir.show()
viewer # Will display a 3D representation of `lopinavir`
lopinavir.save("path/to/lopinavir.sdf") # Saves the structure of `lopinavir` in SDF format to `path/to/lopinavir.sdf`
```

See the included Jupyter Notebook for a working demonstration.

# API

## molvis.Macromolecule.Macromolecule(pdbid: str)
Class representing a macromolecule structure, downloaded from PDB.

### Macromolecule.show([viewer: NGLWidget]) -> [NGLWidget | None]
Either adds a 3D representation of the Macromolecule to the provided widget, or returns a new widget displaying the representation.
If the `viewer` argument is provided, no value will be returned.

### Macromolecule.save(filepath: [pathlib.Path | str]) -> pathlib.Path
Saves the Macromolecule structure in PDB format to the provided path.

## molvis.ChemicalMolecule.ChemicalMolecule(pdbid: str)
Class representing a chemical's molecular structure, downloaded from CHEMBL.

### ChemicalMolecule.show([viewer: NGLWidget]) -> [NGLWidget | None]
Either adds a 3D representation of the ChemicalMolecule to the provided widget, or returns a new widget displaying the representation.
If the `viewer` argument is provided, no value will be returned.

### ChemicalMolecule.save(filepath: [pathlib.Path | str]) -> pathlib.Path
Saves the ChemicalMolecule structure in SDF format to the provided path.
