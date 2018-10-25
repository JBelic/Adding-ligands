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

working_folder_path = r'/Users/basvanbeek/Documents/GitHub/Adding-ligands'

input_ligands = read_molecules(os.path.join(working_folder_path,'LIGANDS'))
input_ligands = [bob(input_ligands[ligand]) for ligand in input_ligands]
input_cores = read_molecules(os.path.join(working_folder_path,'CORES'))
input_cores = [bob(input_cores[core]) for core in input_cores]
# Define bob, the greatest of all properties


new_directory = 'new_molecules'
distance_limit_squared = 2.89 # square of 1.7; if atoms are closer, geometry is rejected


##############################################          new folders              ############################################

# Arranging folders
#working_folder = os.path.join(working_folder_path,new_directory)
#folders = ["monosubstituted", "disupstituted", "error_folder"]
#for i in range(len(folders)):
if not os.path.exists(os.path.join(working_folder_path,new_directory)):
    os.makedirs(os.path.join(working_folder_path,new_directory))

##############################################         Monosubstitution      ################################################

##### Preparing coordinates #####


##### Monosubstitution #####

# To every list of cores one type of ligand is added
# mono_subs contaions of key = name of molecule, value = (coordinates of new molecule, shortest distance between core and ligand after its connection)
mono_subs = []
for ligand in input_ligands:
    # In each list c goes through all cores. New copies of ligands are needed in every loop
    for core in input_cores:
       	mono_subs.append(connect_two_molecules(core,ligand,distance_limit_squared))

#For diatomic molecules, distance is specified by user and it's given as a string.
#For the other molecule it's given as a float and it serves like indicator for steric clash, if distance is too small, molecule is stored in separate folder (list of folders made on the beggining)

for name, molecule in mono_subs.items():
    if molecule[1] == True:
        molecule[0].write(os.path.join(working_folder_path,new_directory, 'mono_' + name +'.xyz'))
    else:
        molecule[0].write(os.path.join(working_folder_path,new_directory, 'err_' + name + '.xyz'))
        print ("Molecule %s has problematic geometry!" % name)





##############################################         Disubstitution      ################################################

##### Preparing coordinates #####

# Copies of coordinates of monosupstituted molecules
monosub_coords = [x[0] for x in list(mono_subs.values())]
# Making list of lists with monosubs coordinates. Core molecules are copied as many times as there are ligands.
monosub_copies = [copy.deepcopy(monosub_coords) for i in input_ligands]

ligand_coords2 = list(input_ligands.values())
monosub_names = list(mono_subs.keys())

####### Disustitution #####

# Takes copies of monosubstituted cores and again adds ligands
di_subs = {}
for d in range(len(input_ligands)):
    for c in range(d*len(input_cores),len(mono_subs)):
        ligando = copy.deepcopy(ligand_coords2)
        di_subs [monosub_names[c]+ "_" + lig_names[d]] = (connect_two_molecules(monosub_copies[d][c],ligando[d],distance_limit_squared))

#For diatomic molecules, distance is specified by user and it's given as a string.
#For the other molecule it's given as a float and it serves like indicator for steric clash, if distance (float number) is too small molecule is stored in separate folder


for name, molecule in di_subs.items():
    if molecule[1] == True:
        molecule[0].write(os.path.join(working_folder_path,new_directory, 'di_' + name +'.xyz'))
    else:
        molecule[0].write(os.path.join(working_folder_path,new_directory, 'err_' + name + '.xyz'))
        print ("Molecule %s has problematic geometry!" % name)


# The End
end = time.time()
print("Elapsed wall-clock time:  ", end - start)

