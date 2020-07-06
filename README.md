# pyrmsdmat
Superimpose a set of protein structures and report a RSMD matrix, in CSV and Mega-compatible formats, using Pymol as a module

# Rationale
Given a set of protein structures the pairwise RMSD values between equivalent CA atoms are
calculated. For the RMSD calculation each pair is superimposed, as implemented in Pymol, following a
refinement procedure that exclude equivalent atoms that are considered to far away. This refinement
procedure will led to smaller RMSD values that the one obtained using all the equivalent CA atoms
(as defined by residue equivalence in a sequence alignment). Use
[rmsdmat](https://github.com/amaurypm/rmsdmat) if you prefer to avoid the refinement.

It is up to the user to select sets of structures that make sense to compare. Remember that for proteins with non-identical sequences the reported RMSD value is for equivalent residues only, any non-conserved sequence is excluded, as there is not equivalent residue to compare with.

# Usage
```
pyrmsdmat [-h] [-o OUTPUT] [-v] structure [structure ...]

positional arguments:
  structure             Any format supported by Pymol.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output files base name [default: rmsd_matrix].
  -v, --version         Show program's version number and exit.
```

# Installation
This is a Python script, so, you can just run the `pyrmsdmat.py` file or put a symbolic link in any
directory of your PATH.

# Dependencies
* Python3
* argparse
* pymol
* numpy

