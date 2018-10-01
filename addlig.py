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

working_folder_path = "/Users/jelena/Desktop/scripts/adding_ligands/1.11.square_dist"

input_ligands = read_molecules(os.path.join(working_folder_path,'lig_test'))
input_cores = read_molecules(os.path.join(working_folder_path,'cor_test'))

new_directory = 'new_molecules'
distance_limit_squared = 2.89 # square of 1.7; if atoms are closer, geometry is rejected


##############################################          new folders              ############################################

# Arranging folders
working_folder = os.path.join(working_folder_path,new_directory)
folders = ["monosubstituted", "disupstituted", "error_folder"]
for i in range(len(folders)):
    if not os.path.exists(os.path.join(working_folder,folders[i])):
        os.makedirs(os.path.join(working_folder, folders[i]))

##############################################         Monosubstitution      ################################################

##### Preparing coordinates #####

#Extracting core coordinates from imported dictionary
coords_c = mol_coords(input_cores)
#Extracting ligand coordinates from imported dictionary
ligand_coords = mol_coords(input_ligands)
#Making list of lists with core coordinates. Core molecules are copied as many times as there are ligands.
cores_list = []
for i in range(len(input_ligands)):
    kopi = copy.deepcopy(coords_c)
    cores_list.append(kopi)
    
#Exctracting core and ligand names from imported dictionary
lig_names = mol_names(input_ligands)
cor_names = mol_names(input_cores)

##### Monosubstitution #####

# To every list of cores one type of ligand is added
# mono_subs contaions of key = name of molecule, value = (coordinates of new molecule, shortest distance between core and ligand after its connection)
mono_subs = {}
for d in range(len(input_ligands)):
    # In each list c goes through all cores. New copies of ligands are needed in every loop
    for c in range(len(input_cores)): 
       	ligando = copy.deepcopy(ligand_coords)
       	mono_subs[cor_names[c]+ "_" + lig_names[d]] = (connect_two_molecules(cores_list[d][c],ligando[d],distance_limit_squared))

#For diatomic molecules, distance is specified by user and it's given as a string. 
#For the other molecule it's given as a float and it serves like indicator for steric clash, if distance is too small, molecule is stored in separate folder (list of folders made on the beggining) 

for name, molecule in mono_subs.items():
    if (molecule[1]) == 'diatomic':
        molecule[0].write(os.path.join(working_folder,folders[0],name+'.xyz'))
    elif float(molecule[1]) >= distance_limit_squared:
        molecule[0].write(os.path.join(working_folder,folders[0],name+'.xyz'))

    elif float(molecule[1]) < distance_limit_squared:
        if molecule[2] == 'HH':
            molecule[0].write(os.path.join(working_folder,folders[0],name+'.xyz'))
        else:
            molecule[0].write(os.path.join(working_folder,folders[2]+'opt.xyz'))
            print ("Molecule %s has problematic geometry!" % name)
        



##############################################         Disubstitution      ################################################

##### Preparing coordinates #####

# Copies of coordinates of monosupstituted molecules   
monosub_coords = [x[0] for x in mol_coords(mono_subs)]
# Making list of lists with monosubs coordinates. Core molecules are copied as many times as there are ligands.
monosub_copies = []
for i in range(len(input_ligands)):
    kopi = copy.deepcopy(monosub_coords)
    monosub_copies.append(kopi)

ligand_coords2 = mol_coords(input_ligands)
monosub_names = mol_names(mono_subs)

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
    if (molecule[1]) =='diatomic':
        molecule[0].write(os.path.join(working_folder,folders[1],name+'.xyz'))
    elif float(molecule[1]) >= distance_limit_squared:
        molecule[0].write(os.path.join(working_folder,folders[1],name+'.xyz'))
    elif float(molecule[1]) < distance_limit_squared:
        if molecule[2] == 'HH':
            molecule[0].write(os.path.join(working_folder,folders[1],name+'.xyz'))
        else:
            molecule[0].write(os.path.join(working_folder,folders[2],name+'opt.xyz'))
            print ("Molecule %s has problematic geometry!" % name)
        


# The End 
end = time.time()
print("Elapsed wall-clock time:  ", end - start)

