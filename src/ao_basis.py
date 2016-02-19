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

## \file ao_basis.py
# This module implements the functions that constructs atomic orbital basis
# by n Gaussian Type Orbitals (nGTO)
#
# Used in: main.py/main/nve/nve_MD/gamess_to_libra/unpack_file
#        : main.py/main/initial_gamess_exe/unpack_file

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_qchem_path"])

from libmmath import *
from libqchem import *

def input_AO_name(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : The list which contains extracted data from the file.
    # This function returns the list of atomic orbital type (s, px, py, pz, etc...) 
    # in params.
    #
    # Used in: main.py/main/nve/nve_MD/gamess_to_libra/unpack_file/ao_basis
    #        : main.py/main/initial_gamess_exe/unpack_file/ao_basis

    atom_spec = params["atom_spec"]
    basis_type = params["basis_type"]
    l_atoms = params["l_atoms"]
    orb_name = []

    for la in l_atoms: # all atoms
        for j in range(0,len(atom_spec)): # specify the kind of the atom
            if la == atom_spec[j]:
                i = j
        aoa_tmp = []
        orb_name1 = []
        scount = 0
        pcount = 0
        dcount = 0
        lcount = 0
        for j in range(0,len(basis_type[i])): # basis number of atoms
            b_tmp = basis_type[i][j]
            if b_tmp == "S":
                scount += 1
                if scount < 2:
                    orb_name1.append("s")
            if b_tmp == "P":
                pcount += 1
                if pcount < 2:
                    orb_name1.append("px")
                    orb_name1.append("py")
                    orb_name1.append("pz")
            if b_tmp == "D":
                dcount += 1
                if dcount < 2:
                    orb_name1.append("dxy")
                    orb_name1.append("dyz")
                    orb_name1.append("dzx")
                    orb_name1.append("dx^2-y^2")
                    orb_name1.append("dz^2")
            elif b_tmp == "L":
                lcount += 1
                if lcount < 2:
                    orb_name1.append("s")
                    orb_name1.append("px")
                    orb_name1.append("py")
                    orb_name1.append("pz")
        orb_name.append(orb_name1)

    print "nGTO=",scount
    params["nGTO"] = scount
    print "orb_name=",orb_name
    params["orb_name"] = orb_name

def construct_ao_basis(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : The list which contains extracted data from the file.
    # This function returns the list of atomic orbital basis as "ao_basis".
    #
    # Used in: main.py/main/nve/nve_MD/gamess_to_libra/unpack_file/ao_basis
    #        : main.py/main/initial_gamess_exe/unpack_file/ao_basis

    l_atoms = params["l_atoms"]
    coor_atoms = params["coor_atoms"]
    expo_s = params["expo_s"]
    expo_p = params["expo_p"]
    expo_d = params["expo_d"]
    coef_s =  params["coef_s"]
    coef_p =  params["coef_p"]
    coef_d =  params["coef_d"]
    nGTO = params["nGTO"]
    orb_name = params["orb_name"]

    ao_basis = []
    Bohr_to_Angs = 0.529177208
    # define atomic orbitals by using the Gaussian basis
    for i in range(0,len(l_atoms)): # all atoms

        k_s = 0 # add nGTO for s orbital
        k_p = 0 #              p orbital

        for j in range(0,len(orb_name[i])):

            ao = AO()  # this is the AO we create - the j-th AO on the i-th atom

            if orb_name[i][j][0] == "s":
                expo_tmp = expo_s[i][k_s:k_s+nGTO];   coef_tmp = coef_s[i][k_s:k_s+nGTO]
                nx,ny,nz = 0,0,0
                k_s += nGTO
            elif orb_name[i][j][0] == "p":
                expo_tmp = expo_p[i][k_p:k_p+nGTO];   coef_tmp = coef_p[i][k_p:k_p+nGTO]
                if orb_name[i][j][1] == "x":
                    nx,ny,nz = 1,0,0
                elif orb_name[i][j][1] == "y":
                    nx,ny,nz = 0,1,0
                elif orb_name[i][j][1] == "z":
                    nx,ny,nz = 0,0,1
                    k_p += nGTO

            elif orb_name[i][j][0] == "d":
                expo_tmp = expo_p[i][k_p:k_p+nGTO];   coef_tmp = coef_p[i][k_p:k_p+nGTO]
                if orb_name[i][j][1:3] == "xy":
                    nx,ny,nz = 1,1,0
                elif orb_name[i][j][1:3] == "yz":
                    nx,ny,nz = 0,1,1
                elif orb_name[i][j][1:3] == "zx":
                    nx,ny,nz = 1,0,1
                elif orb_name[i][j][1:3] == "z^":
                    nx,ny,nz = 0,0,2
                # in the case of dx^2-y^2, i should add AO later. 

            # Construct AO from the primitive Gaussians
            for k in range(0,len(expo_tmp)):         

                # Contraction coefficients correspond to the Gaussian primitives as they are
                R = VECTOR(coor_atoms[i][0], coor_atoms[i][1], coor_atoms[i][2])

                g = PrimitiveG(nx, ny, nz, expo_tmp[k], R) # single point
                # g = PrimitiveG(nx, ny, nz, expo_tmp[k], R/Bohr_to_Angs) : optimization

                # Contraction coefficients correspond to the Gaussian primitives as they are
                ao.add_primitive(coef_tmp[k], g )
                   
                # Contraction coefficients correspond to normalized Gaussian primitives
                # this looks like a more probable case
                #ao.add_primitive(g.normalization_factor() * coef_tmp[k], g )

                # One more option:
                #ao.add_primitive(g.norm1() * coef_tmp[k], g )


            # Normalize the overall contraction
            ao.normalize()

            ao_basis.append(ao)
    
    return ao_basis

def ao_basis(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : The list which contains extracted data from the file.
    # This function returns the list of atomic orbital basis as "ao".
    #
    # Used in: main.py/main/nve/nve_MD/gamess_to_libra/unpack_file
    #        : main.py/main/initial_gamess_exe/unpack_file

    input_AO_name(params)
    
    ao = construct_ao_basis(params)

    #params["ao"] = ao

    return ao
