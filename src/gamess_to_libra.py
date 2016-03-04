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

## \file gamess_to_libra.py 
# This module implements the functions that extract parameters from the gamess output file:
# atomic forces , molecular energies, molecular orbitals, and atomic basis information.
# The forces are used for simulating Classical MD on Libra 
# and the others for calculating time-averaged energies and Non-Adiabatic Couplings(NACs).



from detect import *
from extract import *
from ao_basis import *
from overlap import *
from Ene_NAC import *
from reduce_matrix import *
from moment import *

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_qchem_path"])

#print "\nTest 1: Importing the library and its content"
from libmmath import *
from libqchem import *

def unpack_file(filename):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] filename  GAMESS output file name
    # This function returns the data extracted from the file, in the form of dictionary :
    # atomic basis sets, molecular energies, molecular coefficients, gradients, respectively.
    #
    # Used in:  main.py/main
    #           main.py/main/nve_MD/gamess_to_libra

    f_gam = open(filename,"r")
    l_gam = f_gam.readlines()
    f_gam.close()

    data = {}

    # detect the columns showing parameters
    detect(l_gam,data)    

    # extract the parameters from the columns detected
    extract(l_gam,data)

    # Construct the AO basis
    data["ao_basis"] = ao_basis(data) 

    return data["ao_basis"], data["E"], data["C"], data["gradient"], data


def gamess_to_libra(params, ao, E, C, ite):
    ## 
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params :  contains input parameters , in the directory form
    # \param[in,out] ao :  atomic orbital basis at "t" old
    # \param[in,out] E  :  molecular energies at "t" old
    # \param[in,out] C  :  molecular coefficients at "t" old
    # \param[in] ite : The number of iteration
    # This function outputs the files for excited electron dynamics
    # in "res" directory.
    # It returns the forces which act on the atoms.
    # Also, it returns new atomic orbitals, molecular energies, and
    # molecular coefficients used for calculating time-averaged
    # molecular energies and Non-Adiabatic Couplings(NACs).
    #
    # Used in: main.py/nve_MD/

    # 2-nd file - time "t+dt"  new
    ao2, E2, C2, Grad, data = unpack_file(params["gms_out"])

    # calculate overlap matrix of atomic and molecular orbitals
    P11, P22, P12, P21 = overlap(ao,ao2,C,C2,params["basis_option"])

    # calculate transition dipole moment matrices in the MO basis:
    # mu_x = <i|x|j>, mu_y = <i|y|j>, mu_z = <i|z|j>
    # this is done for the "current" state only    
    mu_x, mu_y, mu_z = transition_dipole_moments(ao2,C2)
    data["mu_x"] = mu_x
    data["mu_y"] = mu_y
    data["mu_z"] = mu_z

    if params["debug_mu_output"]==1:
        print "mu_x:";    mu_x.show_matrix()
        print "mu_y:";    mu_y.show_matrix()
        print "mu_z:";    mu_z.show_matrix()
 
    if params["debug_densmat_output"]==1:
        print "P11 and P22 matrixes should show orthogonality"
        print "P11 is";    P11.show_matrix()
        print "P22 is";    P22.show_matrix()
        print "P12 and P21 matrixes show overlap of MOs for different molecular geometries "
        print "P12 is";    P12.show_matrix()
        print "P21 is";    P21.show_matrix()


    # calculate molecular energies and Non-Adiabatic Couplings(NACs) on MO basis
    E_mol = average_E(E,E2)
    D_mol = NAC(P12,P21,params["dt_nucl"])
    # reduce the matrix size
    E_mol_red = reduce_matrix(E_mol,params["excitations"],params["HOMO"])
    D_mol_red = reduce_matrix(D_mol,params["excitations"],params["HOMO"])

    if params["print_mo_ham"]==1:
        E_mol.show_matrix(params["mo_ham"] + "full_re_Ham_" + str(ite))
        D_mol.show_matrix(params["mo_ham"] + "full_im_Ham_" + str(ite))
        E_mol_red.show_matrix(params["mo_ham"] + "reduced_re_Ham_" + str(ite))
        D_mol_red.show_matrix(params["mo_ham"] + "reduced_im_Ham_" + str(ite))

    # store "t+dt"(new) parameters on "t"(old) ones
    for i in range(0,len(ao2)):
        ao[i] = AO(ao2[i])
    E = MATRIX(E2)
    C = MATRIX(C2)

    return Grad, data, E_mol, D_mol, E_mol_red, D_mol_red
