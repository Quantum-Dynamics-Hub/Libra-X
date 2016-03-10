#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file extract.py
# This module implements the functions that extract
# atomic forces , molecular energies, molecular orbitals, and atomic basis information
# written in gamess output file.
#
# Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file
#        : main.py/main/unpack_file

#************************************************************
# This program extracting parameters from gamess output file.
#************************************************************

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_qchem_path"])
from libmmath import *
from libqchem import *


def atomic_basis_set(l_gam,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which contains extracted data from l_gam file.
    # This function returns the atomic orbital basis as "expo_" and "coef_" of param.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file/extract
    #        : main.py/main/unpack_file/extract

    #  atomic species
    ab_start = params["ab_start"]
    ab_end = params["ab_end"]
    l_atom_spec = []
    atom_spec = []
    for i in range(ab_start,ab_end+1):
        spline = l_gam[i].split()
        if len(spline) == 1:
            l_atom_spec.append(i)
            atom_spec.append(spline[0])
    params["atom_spec"] = atom_spec

    # atomic basis sets
    basis_type = []
    basis_expo = []
    basis_coef = []
    for i in range(0,len(atom_spec)):
        stmp = atom_spec[i]
        type_tmp = []
        expo_tmp = []
        coef_tmp = []
        if i < len(atom_spec) - 1:
            tmp_start = l_atom_spec[i] + 2
            tmp_end = l_atom_spec[i+1] - 2
        else:
            tmp_start = l_atom_spec[i] + 2
            tmp_end = ab_end
        for j in range(tmp_start,tmp_end+1):
            spline = l_gam[j].split()
            coef_tmp1 = []
            if len(spline) > 1:
                type_tmp.append(spline[1])
                expo_tmp.append(float(spline[3]))
                if spline[1] == "L":
                    coef_tmp1.append(float(spline[4]))
                    coef_tmp1.append(float(spline[5]))
                else:
                    coef_tmp1.append(float(spline[4]))
                coef_tmp.append(coef_tmp1)
        basis_type.append(type_tmp)
        basis_expo.append(expo_tmp)
        basis_coef.append(coef_tmp)

    params["basis_type"] = basis_type

    # ******* put into gaussian basis parameter

    l_atoms = params["l_atoms"]
    Ngbf = params["Ngbf"]
    expo_s = []
    expo_p = []
    expo_d = []

    coef_s = []
    coef_p = []
    coef_d = []

    for la in l_atoms: # all atoms
        for j in range(0,len(atom_spec)): # specify the kind of the atom
            if la == atom_spec[j]:
                i = j
        expo_stmp = []
        expo_ptmp = []
        expo_dtmp = []

        coef_stmp = []
        coef_ptmp = []
        coef_dtmp = []
        for j in range(0,len(basis_type[i])): # basis number of atoms
            b_tmp = basis_type[i][j]
            if b_tmp == "S":
                expo_stmp.append(basis_expo[i][j])
                coef_stmp.append(basis_coef[i][j][0])
            elif b_tmp == "P":
                expo_ptmp.append(basis_expo[i][j])
                coef_ptmp.append(basis_coef[i][j][0])
            elif b_tmp == "D":
                expo_dtmp.append(basis_expo[i][j])
                coef_dtmp.append(basis_coef[i][j][0])
            elif b_tmp == "L":
                expo_stmp.append(basis_expo[i][j])
                expo_ptmp.append(basis_expo[i][j])
                coef_stmp.append(basis_coef[i][j][0])
                coef_ptmp.append(basis_coef[i][j][1])

            # f orbitals are not taken into account, so should add them.
        expo_s.append(expo_stmp)
        expo_p.append(expo_ptmp)
        expo_d.append(expo_dtmp)
        coef_s.append(coef_stmp)
        coef_p.append(coef_ptmp)
        coef_d.append(coef_dtmp)

    params["expo_s"] = expo_s
    params["expo_p"] = expo_p
    params["expo_d"] = expo_d
    params["coef_s"] = coef_s
    params["coef_p"] = coef_p
    params["coef_d"] = coef_d


def molecular_orbitals(l_gam,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which contains extracted data from l_gam file.
    # This function returns the molecular orbitals info. as "E" and "C" of params.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file/extract
    #        : main.py/main/unpack_file/extract

    mo_start = params["mo_start"]
    mo_end = params["mo_end"]
    Ngbf = params["Ngbf"]
    stat_span = 4 + Ngbf
    #num_orb_colm = math.ceil(stat_s)

    mol_ene = []
    l_tmp = []
    mol_coef = []

    for i in range(mo_start,mo_end+1):
        spline = l_gam[i].split()
        di = i - mo_start
        # molecular energy
        if di % stat_span == 1 :
            l_tmp.append(i+2)
            for j in range(0,len(spline)):
                mol_ene.append(float(spline[j]))

    # molecular coefficients
    for i in l_tmp:
        for j in range(i,i+Ngbf):
            coef_tmp = []
            spline = l_gam[j].split()
            if len(mol_coef) < Ngbf:
                for k in range(4,len(spline)):
                    coef_tmp.append(float(spline[k]))
                mol_coef.append(coef_tmp)
            else:
                for k in range(4,len(spline)):
                    mol_coef[j-i].append(float(spline[k]))

    # define MATRIX of E, C

    E = MATRIX(Ngbf,Ngbf)
    C = MATRIX(Ngbf,Ngbf)

    for i in range(0,Ngbf):
        E.set(i,i,mol_ene[i])
        for j in range(0,Ngbf):
            C.set(i,j,mol_coef[i][j])

    params["E"] = E
    params["C"] = C
    

def coordinates_of_atoms(l_gam,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which contains extracted data from l_gam file.
    # This function returns the coordinates of atoms info. as "l_atoms" and "coor_atoms" of params.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file/extract
    #        : main.py/main/unpack_file/extract

    coor_start = params["coor_start"]
    coor_end = params["coor_end"]

    l_atoms = []
    l_charges = []
    coor_atoms = []
    for i in range(coor_start,coor_end+1):
        spline = l_gam[i].split()
        # explicit atoms
        l_atoms.append(spline[0])
        # atomic charges
        l_charges.append(float(spline[1]))
        # coordinates of atoms
        coor_tmp = []
        for j in range(len(spline)-3,len(spline)):
            coor_tmp.append(float(spline[j]))
        coor_atoms.append(coor_tmp)

    params["l_atoms"] = l_atoms
    params["l_charges"] = l_charges
    params["coor_atoms"] = coor_atoms


def gradient(l_gam,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which contains extracted data from l_gam file.
    # This function returns the gradients as "gradient" of params.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file/extract
    #        : main.py/main/unpack_file/extract


    grad_start = params["grad_start"]
    grad_end = params["grad_end"]

    gradient = []
    for i in range(grad_start,grad_end+1):
        spline = l_gam[i].split()
        grad_tmp = []
        for j in range(len(spline)-3,len(spline)):
            grad_tmp.append(float(spline[j]))
        gradient.append(grad_tmp)

    params["gradient"] = gradient


def extract(l_gam,params,flag):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which contains extracted data from l_gam file.
    # \param[in] flag : a flag for debugging detect module
    # This function returns the coordinates of atoms, gradients, atomic orbital basis,
    # and molecular orbitals extracted from the file, in the form of dictionary
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file
    #        : main.py/main/unpack_file

    coordinates_of_atoms(l_gam,params)

    gradient(l_gam,params)

    molecular_orbitals(l_gam,params)

    atomic_basis_set(l_gam,params)
    
    if flag == 1:
        print "expo_s=",params["expo_s"]
        print "expo_p=",params["expo_p"]
        print "expo_d=",params["expo_d"]
        print "coef_s=",params["coef_s"]
        print "coef_p=",params["coef_p"]
        print "coef_d=",params["coef_d"]

        print "E Matrix is";    params["E"].show_matrix()
        print "C Matrix is";    params["C"].show_matrix()

        print "l_atoms=", params["l_atoms"]
        print "l_charges=", params["l_charges"]
        print "coor_atoms=", params["coor_atoms"]

        print "gradient=",params["gradient"]

    print "********************************************"
    print "extract program ends"
    print "********************************************"
    print 

