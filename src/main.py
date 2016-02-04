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


from detect import *
from extract import *
from ao_basis import *
#from test_AO import *

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
print "Using the Libra installation at", cwd
sys.path.insert(1,cwd+"/_build/src/mmath")
sys.path.insert(1,cwd+"/_build/src/qchem")


print "\nTest 1: Importing the library and its content"
from libmmath import *
from libqchem import *


# !!!!!!!!!!!!!!!!
# The code below must be a part of a bigger function - to run the code, just call the function
# in the end of the file


# *************** define vacant parameter list ******************
params = {}

# *************** open the GAMESS output file  ******************
f_gam = open("../input/exam01.out","r")
l_gam = f_gam.readlines()
f_gam.close()

# *************** read columns **********************************

read_parameters(l_gam,params)

show_parameters(l_gam,params)

# *************** extract parameters ****************************

atomic_basis_set(l_gam,params)

molecular_orbitals(l_gam,params)

coordinates_of_atoms(l_gam,params)

gradient(l_gam,params)

# *************** the constructor of the molecular coefficients.

mol_coef = params["mol_coef"]  
Ngbf = params["Ngbf"]
C_mat = MATRIX(Ngbf,Ngbf)
for a in range(0,Ngbf): # orbital number
    for i in range(0,Ngbf):
        C_mat.set(a,i,mol_coef[a][i])
print "C_matrix has below matrices"
C_mat.show_matrix()
# *************** the constructor of the atomic orbitals ***************

file_names(params)

atom_spec = params["atom_spec"]
basis_type = params["basis_type"]
l_atoms = params["l_atoms"]
aoa = []
orb_name = []

for la in l_atoms: # all atoms
    for j in range(0,len(atom_spec)): # specify the kind of the atom
        if la == atom_spec[j]:
            i = j
    aoa_tmp = []
    orb_name1 = []
    scount = 0
    lcount = 0
    for j in range(0,len(basis_type[i])): # basis number of atoms
        b_tmp = basis_type[i][j]
        if b_tmp == "S":
            scount += 1
            if scount < 2:
                aoa_1s = AO()
                aoa_tmp.append(aoa_1s)
                orb_name1.append("1s")
        elif b_tmp == "L":
            lcount += 1
            if lcount < 2:
                aoa_2s = AO()
                aoa_2px = AO()
                aoa_2py = AO()
                aoa_2pz = AO()
                aoa_tmp.append(aoa_2s)
                aoa_tmp.append(aoa_2px)
                aoa_tmp.append(aoa_2py)
                aoa_tmp.append(aoa_2pz)
                orb_name1.append("2s")
                orb_name1.append("2px")
                orb_name1.append("2py")
                orb_name1.append("2pz")
    aoa.append(aoa_tmp)
    orb_name.append(orb_name1)

# define other constructor
print "aoa=",aoa
print "orb_name=",orb_name
params["orb_name"] = orb_name
params["aoa"] = aoa


ao = construct_ao_basis(params)
Norb = len(ao)

S = MATRIX(Norb,Norb)

# overlap matrix 
for i in range(0,Norb): # all orbitals
    for j in range(0,Norb):
        S.set(i,j,gaussian_overlap(ao[i],ao[j]))              

print "The overlap matrix is"
S.show_matrix()


print "C_matrix is"
C_mat.show_matrix()

print "Orthogonality check!!!"
(C_mat.T()*S*C_mat).show_matrix()

# overlap matrix
#C_mat.Transpose()
#CS_mat = C_mat * S_mat
#CS_mat.show_matrix()
#CS_mat.Transpose()
#C_mat.Transpose()
#D_mat = CS_mat * C_mat
#print "D_matrix is"
#D_mat.show_matrix()
#gb_i.add_primitive(coeff_2s[i],PrimitiveG(0,0,0, ksi*ksi*alp_2s[i], VECTOR(0,0,0)))



#for o in orb_name1:
#    if o not in orb_name:
#        orb_name.append(o)

# the constructor of atomic basis
#for o in orb_name:
#    ltmp = "vec_" + o
#    ltmp = AO()
# lists for exponent and contraction coefficient
#for o in orb_name:
#    if "1s" in o:
#        print o
    


# ksi = 1.3
#for i in range(0,nGTO):
#    aoa_2s_v1.add_primitive(coeff_2s[i],PrimitiveG(0,0,0, ksi*ksi*alp_2s[i], VECTOR(0,0,0)))

#print "\nTest 5 : Confirm 2s orbital is orthogonalized to 2p orbital"
#print "<2s|2px>=",gaussian_overlap(aoa_2s_v1,aoa_2p_v1)
