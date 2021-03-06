from scm.plams import *
import numpy as np
import math
import copy
import time
import os


@add_to_class(Molecule)
def restr_distance_to_mol(self, other, excluded=(None, None), result_unit='angstrom', return_atoms=False):
        """Modified function:
        Calculate the distance between this molecule and some *other* molecule.

        The distance is measured as the smallest distance between a pair of atoms, one belonging to each of the molecules. Returned distance is expressed in *result_unit*.

        If *return_atoms* is ``False``, only a single number is returned.  If *return_atoms* is ``True``, this method returns a tuple ``(distance, atom1, atom2)`` where ``atom1`` and ``atom2`` are atoms fulfilling the minimal distance, with atom1 belonging to this molecule and atom2 to *other*.

        Has arguments excluded1 and excluded2 for atoms in this molecule which distance from *other* molecule is not calculated.
        """
        dist = float('inf')
        for at1 in self.atoms:
            if at1 != excluded[0] and at1 != excluded[1]:
                for at2 in other.atoms:
                    newdist = (at1.x-at2.x)**2 + (at1.y-at2.y)**2 + (at1.z-at2.z)**2
                    if newdist < dist:
                        dist = newdist
                        atom1 = at1
                        atom2 = at2
        res = Units.convert(math.sqrt(dist), 'angstrom', result_unit)
        if return_atoms:
            return res, atom1, atom2
        return res


def rotation_matrix(vec1, vec2):
    """
    Calculates the rotation matrix rotating *vec1* to *vec2*. Vectors can be any containers with 3 numerical values.
    They don't need to be normalized. Returns 3x3 numpy array.
    """
    a = np.array(vec1)/np.linalg.norm(vec1)
    b = np.array(vec2)/np.linalg.norm(vec2)
    v1,v2,v3 = np.cross(a,b)
    M = np.array([[0, -v3, v2], [v3, 0, -v1], [-v2, v1, 0]])
    return (np.identity(3) + M + np.dot(M,M)/(1+np.dot(a,b)))

def rotation_check(ligand,lig_h,lig_other,core,dist_limit):
    """
    Function takes four arguments, list of atoms and atoms. It takes 1. ligand coordinates, 2. H atom that will be substituted on ligand,
    3. X atom on the ligand, 4. core coordinates and 5. distance limit between atoms.
    It rotates ligand around the bond between the H atom and atom connected to H. Criteria for good geometry is distance between ligand and core atoms.
    The H atom and atom connected to it are excuded from the distance check.
    Returns coordinates of rotated ligand and and binaries True or False; True for good and False for bad geometry, wow!
    """

    # Checks if ligand is already in good position
    distance = ligand.restr_distance_to_mol(core, excluded=(lig_h, lig_other), return_atoms=True)
    if distance[0] >= dist_limit:
        # If molecule has good geometry final_grade is True and it will serve as indicator for good or bad geometry
        final_grade = True

    # If distance is shorter than required, molecule is rotated for 360 degrees in specified 'angle_steps' and for every new position
    # shortest distance between two atoms and symbols of two atoms are placed in list 'atoms_distances'
    else:
        angle_step = 0.1745329252  # equals 10 degrees
        atoms_distances = []
        new_angles=[]
        # List of shorthest distances for different angles
        for alfa in range(37):
            ligand.rotate_bond(lig_h.bonds[0],lig_other, angle_step)
            atoms_distances.append(ligand.restr_distance_to_mol(core, excluded=(lig_h, lig_other), return_atoms=True)[:])

        # atoms_dstances returns list = [distance, ligand atom, core atom]
        # From that list extracts distances
        all_distances = []
        for sublist in atoms_distances:
            all_distances.append(sublist[0])

        # From all distances only ones that meet the criteria are extracted (larger distance than distance limit)
        matches = [all_distances.index(x) for x in all_distances if x >= dist_limit]
        # If the list of distances, that meet the criteria, is not empty - takes first distance that met the criteria
        if len(matches) != 0:
            good_angle = matches[0]*angle_step
            ligand.rotate_bond(lig_h.bonds[0],lig_other, good_angle)

            final_grade = True

        # If there are no atoms that meet the distance criteria, script does another search trough the simbols of closest atoms. 
        # If two hydrogens are closest atoms, distance between them can be shorter than for other atoms.
        else:
            ligand_atoms = []
            core_atoms = []
            for sublist in atoms_distances:
                ligand_atoms.append(sublist[1].symbol)
                core_atoms.append(sublist[2].symbol)

            h_matches = [all_distances.index(x) for x in all_distances if x >= 2.25] # 2.25 = 1.5^2

            # If there are two H atoms that are on distance higher than 1.5 A, geomtetry is accepted
            if len(h_matches) != 0:
                h_match = [x for x in h_matches if ligand_atoms[x] == core_atoms[x] and ligand_atoms[x] == "H"]
                h_good_angle = h_match[0]*angle_step
                ligand.rotate_bond(lig_h.bonds[0],lig_other, h_good_angle)

                final_grade = True

            else:
                final_grade = False

    return ligand, final_grade


def connect_two_molecules(core,ligand,dist_limit,rot_check=True):
    """
    Takes molecule coordinates as arguments.
    Connects two molecules substituting atoms mentioned is coordinate file
    Returns list of three items: 1. coordinates of new molecule, 2. binary - True for good geometry, False for bad.
    """
    core = copy.deepcopy(core)
    ligand = copy.deepcopy(ligand)

    # Lenght of C-C and C-N bonds for adjustment of new bond between core and ligand
    cc_lenght = 1.54
    nc_lenght = 1.469

    # Guess
    ligand.guess_bonds()
    core.guess_bonds()
    
    # Good guy Bob.
    # Bob searches for atoms that have properties 'bob'. Those atoms will be substituted.
    bobed_core = {}
    for atom in core:
        if isinstance(atom.properties.bob, int):
            bobed_core[atom.properties.bob] = atom
    core_h = bobed_core[sorted(bobed_core)[0]]

    bobed_lig = {}
    for atom in ligand:
        if isinstance(atom.properties.bob, int):
            bobed_lig[atom.properties.bob] = atom
    lig_h = bobed_lig[sorted(bobed_lig)[0]]
    
    # Defining atom connected to marked atom and vector between them
    core_other = core_h.bonds[0].other_end(core_h)
    core_vector = core_h.vector_to(core_other)
    lig_other = lig_h.bonds[0].other_end(lig_h)
    lig_vector = lig_other.vector_to(lig_h)

    # Rotation of ligand - aligning diresctions of two vectors
    rotmat = rotation_matrix(lig_vector, core_vector)
    ligand.rotate(rotmat)
    ligand.translate(lig_other.vector_to(core_h))

    # Deletes atom in core molecule
    core.delete_atom(core_h)

    # Resizing the distance for new bond considering it's not always C-C bond
    if core_other.symbol == "C":
        bond_lenght = 1.54
    elif core_other.symbol == "N":
        bond_lenght = 1.469
    else:
        bond_lenght = 1.5

    hc_vec = np.array(core_other.vector_to(lig_other))
    cc_vec = hc_vec*(bond_lenght/np.linalg.norm(hc_vec))
    diff_vec = cc_vec - hc_vec
    ligand.translate(diff_vec)

    # If ligand is not diatomic, rotation check is done to confirm the good geometry or search for the new one.
    # Function rotation_check returns list [coordinates, True(or False)]
    # Coordinates from input_ligand are changed by using rotation_check fuction

    if rot_check == True:
        if len(ligand) > 2:
            final_check = rotation_check(ligand,lig_h,lig_other,core,dist_limit)[1]
        else:
            # If it is diatomic, len(ligand) < 2, distance is set on 1.5 and there is no rotation check
            final_check = True


    # Deletes atom in ligand
    ligand.delete_atom(lig_h)

    # Adds two coordinates together
    new_molecule = core + ligand
    new_molecule.properties.name = core.properties.name + "_" + ligand.properties.name
    
    return  new_molecule, final_check


def bob(plams_mol):
    """
    Marks plams atoms with a property called bob.
    The argument bob (<int>) represents the order in which atoms will be substituted.

    plams_mol <plams.Molecule>: An input molecule with the plams_mol.properties.comment property.

    return <plams.Molecule>: A plams molecule with specific atoms marked with the
        plams_atom.properties.bob property.
    """
    comment = plams_mol.properties.comment
    comment = comment.split()
    comment_list = []
    for item in comment:
        try:
            index = int(item)
            comment_list.append(index)
        except ValueError:
            pass
    for i, index in enumerate(comment_list):
        plams_mol[index].properties.bob = i
    return plams_mol

def monosubstitution(input_ligands, input_cores, dist_limit):
    # To every list of cores one type of ligand is added
    # mono_subs contaions of key = name of molecule, value = (coordinates of new molecule, shortest distance between core and ligand after its connection)
    mono_subs = []
    for ligand in input_ligands:
        # In each list c goes through all cores. New copies of ligands are needed in every loop
        for core in input_cores:
            mono_subs.append(connect_two_molecules(core,ligand, dist_limit))
    
    return mono_subs

def mono_di_substitution(input_ligands, input_cores, dist_limit):

    mono_subs = monosubstitution(input_ligands, input_cores, dist_limit)
    
    # Exctracts coordinates for list of mono_subs
    monosub_coords = [x[0] for x in mono_subs]
    
    # Takes monosubstituted core and adds ligands, without duplicas
    di_subs = []
    for d in range(len(input_ligands)):
        ligand = input_ligands[d]
        for c in range(d*len(input_cores),len(mono_subs)):
            core = mono_subs[c][0]
            di_subs.append(connect_two_molecules(core,ligand, dist_limit))
    
    return mono_subs, di_subs





