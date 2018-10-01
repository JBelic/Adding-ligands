# Adding-ligands
Does monosubstitution and disubstitution of ligands to different core molecules. 

Script is using functions from Python Library for Automating Molecular Simulations (PLAMS).

It connects XYZ coordinate files of given molecules resulting in new XYZ files of all possible combinations for given list of molecules.

Substitution is done on the first atom of coordinate list, which should be hydrogen atom. Hydrogen that is supposed to be substituted should be on the top of the coordinate list of both, core and ligand. There is a script for automating this process.
Program defines vector between hydrogen atom and one other atom that is connected to it (X atom), for core and ligand. Rotates ligand with aim to make H-X vectors parallel. Translates ligand to the position of H(core) and deletes redundant H atoms on both molecules.

Stores monosubstituted and disubstituted in specific folder. Molecules with problematic geometry are places in error folder.
