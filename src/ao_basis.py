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

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
print "Using the Libra installation at", cwd
sys.path.insert(1,cwd+"/_build/src/mmath")
sys.path.insert(1,cwd+"/_build/src/qchem")
from libmmath import *
from libqchem import *



def file_names(params):
# The name of this function is not informative - choose a more appropriate one

# Need a description of what this function does
# what are the inputs and what are the outputs

    l_atoms = params["l_atoms"]
    atom_spec = params["atom_spec"]
    basis_type = params["basis_type"]
    basis_expo = params["basis_expo"]
    basis_coef = params["basis_coef"]
    nGTO = params["nGTO"]
    Ngbf = params["Ngbf"]

    expo_1 = []
    expo_2 = []
    expo_3 = []

    coef_1s = []
    coef_2s = []
    coef_2p = []
    coef_3s = []
    coef_3p = []
    coef_3d = []

    for la in l_atoms: # all atoms
        for j in range(0,len(atom_spec)): # specify the kind of the atom
            if la == atom_spec[j]:
                i = j
        expo_1tmp = []
        expo_2tmp = []
        #expo_3tmp = []
        coef_1stmp = []
        coef_2stmp = []
        coef_2ptmp = []
        #coef_3stmp = []
        #coef_3ptmp = []
        #coef_3dtmp = []
        for j in range(0,len(basis_type[i])): # basis number of atoms
            b_tmp = basis_type[i][j]
            if b_tmp == "S":
                expo_1tmp.append(basis_expo[i][j])
                coef_1stmp.append(basis_coef[i][j][0])
            elif b_tmp == "L":
                expo_2tmp.append(basis_expo[i][j])
                coef_2stmp.append(basis_coef[i][j][0])
                coef_2ptmp.append(basis_coef[i][j][1])
            #elif b_tmp == "D":
            #    expo_3tmp.append(basis_expo[i][j])
            #    coef_3dtmp.append(basis_coef[i][j][0])
            #    print "you inputed illegal character (or D), so exit"
            #    sys.exit
            # f orbitals are not taken into account, so should add them.
        expo_1.append(expo_1tmp)
        expo_2.append(expo_2tmp)
        coef_1s.append(coef_1stmp)
        coef_2s.append(coef_2stmp)
        coef_2p.append(coef_2ptmp)

    print "expo_1=",expo_1
    print "expo_2=",expo_2
    print "coef_1s=",coef_1s
    print "coef_2s=",coef_2s
    print "coef_2p=",coef_2p

    params["expo_1"] = expo_1
    params["expo_2"] = expo_2
    #params["expo_3"] = expo_3
    params["coef_1s"] = coef_1s
    params["coef_2s"] = coef_2s
    params["coef_2p"] = coef_2p

    return 1.0

def construct_ao_basis(params): # old add_PrimitiveG
# The function takes the parameters in the format adopted for
# storing info read from the Gamess output file
# The function returns the list of AO objects - the basis
    
    l_atoms = params["l_atoms"]
    coor_atoms = params["coor_atoms"]
    expo_1 = params["expo_1"]
    expo_2 = params["expo_2"]
    coef_1s =  params["coef_1s"]
    coef_2s =  params["coef_2s"]
    coef_2p =  params["coef_2p"]
    orb_name = params["orb_name"]
    aoa = params["aoa"]

    ao_basis = []
    # define atomic orbitals using the Gaussian basis
    for i in range(0,len(l_atoms)): # all atoms
        for j in range(0,len(aoa[i])):

            ao = AO()  # this is the AO we create - the j-th AO on the i-th atom

            if orb_name[i][j] == "1s":
                expo_tmp = expo_1[i];   coef_tmp = coef_1s[i]
                nx,ny,nz = 0,0,0

            elif orb_name[i][j] == "2s":
                expo_tmp = expo_2[i];   coef_tmp = coef_2s[i]
                nx,ny,nz = 0,0,0

            elif orb_name[i][j][0:2] == "2p":
                expo_tmp = expo_2[i];   coef_tmp = coef_2p[i]

                if orb_name[i][j][2] == "x":
                    nx,ny,nz = 1,0,0
                elif orb_name[i][j][2] == "y":
                    nx,ny,nz = 0,1,0
                elif orb_name[i][j][2] == "z":
                    nx,ny,nz = 0,0,1

            # Construct AO from the primitive Gaussians
            for k in range(0,len(expo_tmp)):         
                # Contraction coefficients correspond to the Gaussian primitives as they are
                R = VECTOR(coor_atoms[i][0], coor_atoms[i][1], coor_atoms[i][2])
                g = PrimitiveG(nx,ny,nz, expo_tmp[k], R)

                # Contraction coefficients correspond to the Gaussian primitives as they are
                #ao.add_primitive(coef_tmp[k], g )
                   
                # Contraction coefficients correspond to normalized Gaussian primitives
                # this looks like a more probable case
                ao.add_primitive(g.normalization_factor() * coef_tmp[k], g )


            # Normalize the overall contraction
            ao.normalize()

            ao_basis.append(ao)
    
    
    return ao_basis

