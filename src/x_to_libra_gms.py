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

## \file x_to_libra_gms.py 
# This module implements the functions that extract parameters from the gamess output file:
# atomic forces , molecular energies, molecular orbitals, and atomic basis information.
# The forces are used for simulating Classical MD on Libra 
# and the others for calculating time-averaged energies and Non-Adiabatic Couplings(NACs).

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from extract_gms import *
from overlap import *
from hamiltonian_el import *
from moment import *
from misc import *


def exe_gamess(params):
    ##
    # This is a function that call GAMESS execution on the compute node
    # \param[in] params Input data containing all manual settings and some extracted data.
    #
    # Used in main.py/main and md.py/run_MD

    inp = params["gms_inp"]
    out = params["gms_out"]
    nproc = params["nproc"]

    scr_dir = params["scr_dir"]
    rungms = params["rungms"]
    VERNO = params["VERNO"]

    # set environmental variables for GAMESS execution
    os.environ["SCR"] = scr_dir
    os.environ["USERSCR"] = scr_dir
    os.environ["GMSPATH"] = params["GMSPATH"]

    #os.system("/usr/bin/time rungms.slurm %s 01 %s > %s" % (inp,nproc,out))
    os.system("/usr/bin/time %s %s %s %s > %s" % (rungms,inp,VERNO,nproc,out))

    # delete the files except input and output ones to do another GAMESS calculation.
    os.system("rm *.dat")
    os.system("rm -r %s/*" %(scr_dir))


def gamess_to_libra(params, ao, E, sd_basis, active_space,suff):
    ## 
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params         contains input parameters , in the directory form
    # \param[in,out] ao         atomic orbital basis at "t" old
    # \param[in,out] E          total excitation energies at "t" old
    # \param[in,out] sd_basis   Basis of Slater determinants at "t" old (list of CMATRIX object).
    #                           In the present implementation, it contains a single determinant
    # \param[in] active_space   A list of indices (starting from 1) of the MOs to include in calculations (and to read from the QE output files)  
    # \param[in] suff           A suffix to add to the name of the output files; this suffix is now considered to be of a string type 
    #                           - so you can actually encode both the iteration number (MD timestep),
    #                           the nuclear cofiguration (e.g. trajectory), and any other related information
    # Returned datas are explained above the return line.
    #
    # Used in: md.py/run_MD

    # 2-nd file - time "t+dt"  new
    label, Q, R, Grad, E2, sd_basis2, ao2 = gms_extract(params["gms_out"],params["excitations"],params["min_shift"],active_space,params["debug_gms_unpack"])

    # Gradients
    # in this implementation (CPA), the gradients on all excited states are the same
    nstates = len(params["excitations"])
    all_grads = []
    for i in xrange(nstates):
        grd = []
        for g in Grad:
            grd.append(VECTOR(g))
        all_grads.append(grd)
    
    t = Timer()
    # calculate overlap matrix of atomic and molecular orbitals
    #P11, P22, P12, P21 = overlap(ao,ao2,sd_basis[0],sd_basis2,params["basis_option"])
    # In semi-empirical case, AO basis isn't used and last input(d2_max) can be set 0.0.
    sz = len(active_space)
    act = Py2Cpp_int(active_space)
    P11, P22, P12, P21 = CMATRIX(sz,sz), CMATRIX(sz,sz), CMATRIX(sz,sz), CMATRIX(sz,sz) 
    tt = Timer()
    tt.start()
    MO_overlap(P11,sd_basis,sd_basis,act,act,0.0)
    MO_overlap(P22,sd_basis2,sd_basis2,act,act,0.0)
    MO_overlap(P12,sd_basis,sd_basis2,act,act,0.0)
    tt.stop() # 3 components same as dipole moment
    MO_overlap(P21,sd_basis2,sd_basis,act,act,0.0)
    print "Time to compute in MO_overlap= ",tt.show(),"sec"

    # calculate transition dipole moment matrices in the MO basis:
    # mu_x = <i|x|j>, mu_y = <i|y|j>, mu_z = <i|z|j>
    # this is done for the "current" state only    
    t.start()
    mu_x, mu_y, mu_z = transition_dipole_moments(ao2,sd_basis2)
    mu = [mu_x, mu_y, mu_z] # now mu is defined as a CMATRIX list.
    t.stop()
    print "Time to compute in dipole moment= ",t.show(),"sec"

    sys.exit(0)

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


    ### TO DO: In the following section, we need to avoid computing NAC matrices in the full
    # basis. We will need the information on cropping, in order to avoid computations that 
    # we do not need (the results are discarded anyways)
    # calculate molecular energies and Non-Adiabatic Couplings(NACs) on MO basis
    E_ave = average_E(E,E2)
    #nac = compute_nac_sd(sd_basis[0], sd_basis2, params["dt_nucl"])
    nac = NAC(P12,P21,params["dt_nucl"])

    # reduce the matrix size
    #E_mol_red = reduce_matrix(E_mol,params["min_shift"], params["max_shift"],params["HOMO"])
    #D_mol_red = reduce_matrix(D_mol,params["min_shift"], params["max_shift"],params["HOMO"])

    ### END TO DO

    #if params["print_mo_ham"]==1:
    #E_mol.show_matrix(params["mo_ham"] + "full_re_Ham_" + suff)
    #D_mol.show_matrix(params["mo_ham"] + "full_im_Ham_" + suff)
    #E_mol_red.show_matrix(params["mo_ham"] + "reduced_re_Ham_" + suff)
    #D_mol.show_matrix(params["mo_ham"] + "reduced_im_Ham_" + suff)
    # ********** "CMATRIX.show_matrix(filename)" is not exported ****** 

    # store "t+dt"(new) parameters on "t"(old) ones
    for i in range(0,len(ao2)):
        ao[i] = AO(ao2[i])
    E = MATRIX(E2)  # at time t+dt


    #sd_basis = [sd_basis2] #******* modified ******
    sd_basis = []
    for ex_st in xrange(nstates):
            sd_basis.append(CMATRIX(sd_basis2))

    # useless lines: nac is already defined as CMATRIX.
    #nac = CMATRIX(D_mol.num_of_rows, D_mol.num_of_cols)
    #for i in xrange(D_mol.num_of_rows):
    #    for j in xrange(D_mol.num_of_cols):
    #        nac.set(i,j,D_mol.get(i,j),0.0)

    # Returned data:
    ### Grad: Grad[k] - the gradient w.r.t. to k-th nucleus
    ### data: a dictionary containing transition dipole moments
    ### E_mol: the matrix of the 1-el orbital energies in the full space of the orbitals
    ### D_mol: the matrix of the NACs computed with 1-el orbitals. Same dimension as E_mol
    ### E_mol_red (MATRIX): the matrix of the 1-el orbital energies in the reduced (active) space
    # E_ave : the matrix of the total excitation energy averaged over energies at "t" and "t+dt"
    # nac (CMATRIX): the matrix of the NACs computed with 1-el orbital. Same dimension as E_mol_red
    # sd_basis : (list of CMATRIX, only 1 element): the SD of the present calculation - in the full dimension
    # all_grads: all_grads[i][k] - the gradient w.r.t. to k-th nucleus of i-th excitation state
    # mu : mu[i] transition dipole moment of i-th DOF. (mu_x, mu_y, mu_z)

    return E_ave, nac, sd_basis, all_grads, mu

