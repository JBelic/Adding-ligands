from scm.plams import *
import numpy as np
import math
import copy
import time
from addlig_library import *
import os

# Time measuring
start = time.time()

##### Input #####

# Path to the working folder where are prepared molecules and where folder with new coordinares
# will be made with the specific name
working_folder_path = '/Users/jelena/Desktop/Adding-ligands'

input_ligands = read_molecules(os.path.join(working_folder_path,'LIGANDS'))
input_ligands = [bob(input_ligands[ligand]) for ligand in input_ligands]
input_cores = read_molecules(os.path.join(working_folder_path,'CORES'))
input_cores = [bob(input_cores[core]) for core in input_cores]


new_directory = 'new_molecules'
distance_limit_squared = 2.89 # square of 1.7; if atoms are closer, geometry is rejected


##############################################          new folder               ############################################

# Makes new folder
if not os.path.exists(os.path.join(working_folder_path,new_directory)):
    os.makedirs(os.path.join(working_folder_path,new_directory))

############################################ Monosubstitution and disubstitution ############################################

new_molecules = mono_di_substitution(input_ligands, input_cores, distance_limit_squared)
# mono_subs contains sublists with plams molecule ('molecule') and binary ('check')
# depending on binary value, False or True, name of the coordinate file is with or without 'err'.
for item in new_molecules: 
    print (new_molecules.index(item))
    n= ['mono','di']
    for molecule, check in item:
        if check == True:
            molecule.write(os.path.join(working_folder_path,new_directory, n[new_molecules.index(item)] + '_' + molecule.properties.name +'.xyz'))
        else:
            molecule.write(os.path.join(working_folder_path,new_directory, 'err_' + n[new_molecules.index(item)]+ molecule.properties.name + '.xyz'))


# The End
end = time.time()
print("Elapsed wall-clock time:  ", end - start)

