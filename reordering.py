from scm.plams import *
import os


def swap(mol, i, j):
	"""
	Takes three arguments: 1. list of atoms (molecule coordinates), 2. number of the present position of an atom
	and 3. number of desired new position for an atom.
	"""
	mol.atoms[i-1], mol.atoms[j] = mol.atoms[j], mol.atoms[i-1]


### User input ###

input_molecules = read_molecules('/home/jbelic/PLAMS/reordering/cor_test')

working_folder_path = '/home/jbelic/PLAMS/reordering/'
new_directory = 'reordered_molecules'

######

# Making new folder with specified path and given name for new folder
working_folder = os.path.join(working_folder_path,new_directory)
os.makedirs(working_folder)

# Uses the comment from a xyz file where the number of atoms that will be substituted should be written
# Specified atoms will be swapped with atoms on the top. First on the comment list gets the fist place,
# second gets second place and so on.
for name, mol in input_molecules.items():
	sb_str = mol.properties.comment
	sb_list = sb_str.split()
	sb_int = []
	for i in sb_list:
		sb_int.append(int(i))
	sb_number = len(sb_int)
	for i in range(sb_number):
                if sb_int[i] != i:
		    swap(mol,sb_int[i],i)
	mol.write(os.path.join(working_folder, name + '.xyz'))

a = [22, 35, 11]

sb_list = [int(item) for item in a.split()]
for i, index in enumerate(sb_list):
    mol[index].properties.bob = i

