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

## \file main.py
# This module sets initial parameters from GAMESS output, creates initial system, 
# and executes runMD script.
# 
# It returns the data from runMD for debugging the code.

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


from create_input_gms import *
from create_input_qe import *
from x_to_libra_gms import *
from x_to_libra_qe import *
from md import *
from spin_indx import *


def construct_active_space(params):
##
# This function constructs the active space using:
# params["excitations"] - a list of user-provided excitations to include in the TD-SE basis
# params["nel"] - the number of electrons in the system - used to determine the index of HOMO
#
    print "In construct_active_space..."

    active_space = []

    homo = params["nel"]/2 +  params["nel"] % 2  # the index of HOMO
    print "nel = ", params["nel"]
    print "homo = ", homo

    # Find the lowest orbital from which excitation occurs and find the highest orbital to where
    # electron is sent
    min_from = 0
    max_to = 0
    for ex in params["excitations"]:
        print "ex = ", ex.from_orbit[0], ex.to_orbit[0]

        if min_from>ex.from_orbit[0]:
            min_from = ex.from_orbit[0]

        if max_to<ex.to_orbit[0]:
            max_to = ex.to_orbit[0]

    print "min_from = ", min_from
    print "max_to = ", max_to

    # This definition may be a bit excessive, but it is general and correct (at least for single
    # excitations)
    active_space = range(min_from + homo, max_to + homo + 1)

    print "active_space = ", active_space

    return active_space




def main(params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params  A list of  input parameters from {gms,qe}_run.py 
    #### Returned data:
    #### test_data - the output data for debugging, in the form of dictionary
    #### data - the data extracted from gamess output file, in the form of dictionary
    #
    # This function prepares lists of initial parameters from GAMESS output file
    # and executes run_MD function where NA-MD calculation is done.
    #
    # Used in:  run.py

    dt_nucl = params["dt_nucl"]
    nstates = len(params["excitations"])
    ninit = params["nconfig"]  
    SH_type = params["tsh_method"]

    num_SH_traj = 1
    if SH_type >= 1: # calculate no SH probs.  
        num_SH_traj = params["num_SH_traj"]

    ntraj = nstates*ninit*num_SH_traj

    #######
    active_space = []
    mo_pool_alp, mo_pool_bet = None, None

    if params["interface"]=="QE":
        pass
#        active_space = [5,6,7]  # For C2H4 
#    #print "Implement the algorithm to define the active space"
#    #********** active space is defined here *****************
    elif params["interface"]=="GAMESS":
        for i in range(params["min_shift"],params["max_shift"]+1):
            active_space.append(i+params["HOMO"]+1) # Here MO order start from 1, not 0.
    #*********************************************************
    #sys.exit(0)
    ######

    ################# Step 0: Use the initial file to create a working input file ###############

    if params["interface"]=="GAMESS": 
        os.system("cp %s %s" %(params["gms_inp0"], params["gms_inp"]))

    elif params["interface"]=="QE":
        for ex_st in xrange(nstates):
            os.system("cp x%i.scf.in x%i.scf_wrk.in" % (ex_st, ex_st))

    #### Step 1: Read initial input, run first calculation, and initialize the "global" variables ####

    t = Timer()
    t.start()
    # Initialize variables for a single trajectory first!
    label, Q, R, ao, tot_ene = None, [], None, [], None
    sd_basis = SDList()
    all_grads = []
    e = MATRIX(nstates,nstates)  # contains total energy of excited states!
    
    if params["interface"]=="GAMESS":

        params["gms_inp_templ"] = read_gms_inp_templ(params["gms_inp"])
        exe_gamess(params)
        label, Q, R, grads, E, c, ao, params["nel"] = gms_extract(params["gms_out"],params["excitations"],params["min_shift"],\
                                                            active_space,params["debug_gms_unpack"])
        e = MATRIX(E)
        homo = params["nel"]/2 +  params["nel"] % 2

        for ex_st in xrange(nstates): 
            mo_pool_alp = CMATRIX(c)
            mo_pool_bet = CMATRIX(c)
            alp,bet = index_spin(params["excitations"][ex_st],active_space, homo)

            # use excitation object to create proper SD object for different excited state
            sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet))
            sd_basis.append(sd)
            all_grads.append(copy.deepcopy(grads)) # newly defined
        #print "sd_basis=",sd_basis
        #print "all_grads=",all_grads
        #sys.exit(0)

    elif params["interface"]=="QE":
        params["qe_inp_templ"] = []
        for ex_st in xrange(nstates):

            params["qe_inp_templ"].append( read_qe_inp_templ("x%i.scf_wrk.in" % ex_st) )
            exe_espresso(ex_st)
            flag = 0

            # Basically, here we automatically determine the position of HOMO and will construct
            # the active space before actually using it
            if ex_st==0:
                tot_ene, params["norb"], params["nel"], params["nat"], params["alat"], icoord, iforce = qe_extract_info("x%i.scf.out" % ex_st, flag)
                active_space = construct_active_space(params)

            tot_ene, label, R, grads, mo_pool, params["norb"], params["nel"], params["nat"], params["alat"] = qe_extract("x%i.scf.out" % ex_st, flag, active_space, ex_st)

            mo_pool_alp = CMATRIX(mo_pool)
            mo_pool_bet = CMATRIX(mo_pool)

            homo = params["nel"]/2 +  params["nel"] % 2

            alp,bet = index_spin(params["excitations"][ex_st], active_space, homo)
            # use excitation object to create proper SD object for different excited state
            sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet) )

            #sd_basis2 = SDList()
            sd_basis.append(sd)

            all_grads.append(grads)
            e.set(ex_st, ex_st, tot_ene)


    # Now, clone the single-trajectory variables, to initialize the bunch of such parameters
    # for all trajectories
    sd_basis_list = []
    ham_el_list = []    
    ao_list = []
    e_list = []
    grad_list = []
    label_list = []
    Q_list = []
    R_list = []

    for i in xrange(ntraj):
        # AO
        ao_tmp = []
        for x in ao:
            ao_tmp.append(AO(x))
        ao_list.append(ao_tmp)

        # E and C
        e_list.append(MATRIX(e))

        # Slater determinants
        # eventually, the ordering is this: sd_basis_list[traj][ex_st] - type CMATRIX
        #sd_basis_tr = []
        #for sd in sd_basis:        
            #sd_basis_tr.append(CMATRIX(sd))
            #sd_basis_tr.append(sd)
        #sd_basis_list.append(sd_basis_tr)
        sd_basis_list.append(sd_basis)

        # Gradients
        # eventually, the ordering is this: grad_list[traj][ex_st][n_atom] - type VECTOR
        grd2 = []
        for grad in all_grads: # over all excited states
            grd1 = []
            for g in grad:     # over all atoms
                grd1.append(VECTOR(g))
            grd2.append(grd1)
        grad_list.append(grd2)
    

        # Coords
        rr = []
        for r in R:
            rr.append(VECTOR(r))
        R_list.append(rr)

        # Labels
        lab = []
        for i in xrange(len(label)):
            lab.append(label[i])
        label_list.append(lab)

        # Q
        qq  = []
        for q in Q:
            qq.append(q)
        Q_list.append(qq)

    t.stop()
    print "Step 1 in main takes",t.show(),"sec"

    #sys.exit(0)

    ################## Step 2: Initialize molecular system and run MD part with TD-SE and SH####


    t.start()
    print "Initializing nuclear configuration and electronic variables..."
    rnd = Random() # random number generator object
    syst = []
    el = []

    # all excitations for each nuclear configuration
    for i in xrange(ninit):
        print "init_system..." 
        for i_ex in xrange(nstates):
            for itraj in xrange(num_SH_traj):
                print "Create a copy of a system"  
                df = 0 # debug flag
                # Here we use libra_py module!
                # Utilize the gradients on the ground (0) excited state
                x = init_system.init_system(label_list[i], R_list[i], grad_list[i][0], rnd, params["Temperature"], params["sigma_pos"], df, "elements.txt")
                syst.append(x)    

                print "Create an electronic object"
                el.append(Electronic(nstates,i_ex))
    
    # set list of SH state trajectories
    #sys.exit(0)
    print "run MD"
    run_MD(syst,el,ao_list,e_list,sd_basis_list,params,label_list, Q_list, active_space)
    print "MD is done"

    t.stop();
    print "step 2 in main takes",t.show(),"sec"

    #return data, test_data
