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
from extract_qe import *
import include_mm

def construct_active_space(params):
##
# This function constructs the active space using:
# params["excitations"] - a list of user-provided excitations to include in the TD-SE basis
# params["nel"] - the number of electrons in the system - used to determine the index of HOMO
#

    active_space = []
    homo = params["nel"]/2 +  params["nel"] % 2  # the index of HOMO

    # Find the lowest orbital from which excitation occurs and find the highest orbital to where
    # electron is sent
    min_from = 0
    max_to = 0
    for ex in params["excitations"]:

        if min_from>ex.from_orbit[0]:
            min_from = ex.from_orbit[0]

        if max_to<ex.to_orbit[0]:
            max_to = ex.to_orbit[0]


    # This definition may be a bit excessive, but it is general and correct (at least for single
    # excitations)
    active_space = range(min_from + homo, max_to + homo + 1)

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
    nstates_init = len(params["excitations_init"])
    ninit = params["nconfig"]  
    SH_type = params["tsh_method"]
    # nspin = params["nspin"]  This parameter is used only in Libra-QE interface.
    uff = params["ff"]

    num_SH_traj = 1
    if SH_type >= 1: # calculate no SH probs.  
        num_SH_traj = params["num_SH_traj"]

    ntraj = nstates_init*ninit*num_SH_traj



    active_space = []
    mo_pool_alp, mo_pool_bet = None, None

    if params["interface"]=="QE":
        #ntraj = nstates*ninit*num_SH_traj
        pass
#        active_space = [5,6,7]  # For C2H4 
#    #********** active space is defined here *****************
    elif params["interface"]=="GAMESS":
        #ntraj = nstates_init*ninit*num_SH_traj
        for i in range(params["min_shift"],params["max_shift"]+1):
            active_space.append(i+params["HOMO"]+1) # Here MO order start from 1, not 0.
    #*********************************************************

    ################# Step 0: Use the initial file to create a working input file ###############

    if params["interface"]=="GAMESS": 
        os.system("cp %s %s" %(params["gms_inp0"], params["gms_inp"]))

    elif params["interface"]=="QE":
        for ex_st in xrange(nstates):
            os.system("cp x%i.scf.in x%i.scf_wrk.in" % (ex_st, ex_st))

    #### Step 1: Read initial input, run first calculation, and initialize the "global" variables ####
    #t = Timer()
    #t.start()
    # Initialize variables for a single trajectory first!
    label, Q, R, ao, tot_ene = None, [], None, [], None
    sd_basis = SDList()
    all_grads = []
    e = MATRIX(nstates,nstates)  # contains total energy of excited states!
    
    if params["interface"]=="GAMESS":

        params["gms_inp_templ"] = read_gms_inp_templ(params["gms_inp"])
        exe_gamess(params)
        label, Q, R, grads, E, c, ao, params["nel"] = gms_extract(params["gms_out"],params["excitations"],params["min_shift"],active_space,params["debug_gms_unpack"])
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

    elif params["interface"]=="QE":
        params["qe_inp_templ"] = []
        for ex_st in xrange(nstates): # modify here #

            params["qe_inp_templ"].append( read_qe_inp_templ("x%i.scf_wrk.in" % ex_st) )
            exe_espresso(ex_st)
            flag = 0

            # Basically, here we automatically determine the position of HOMO and will construct
            # the active space before actually using it
            if ex_st==0:
                tot_ene, params["norb"], params["nel"], params["nat"], params["alat"], icoord, iforce = qe_extract_info("x%i.scf.out" % ex_st, ex_st, flag)
                active_space = construct_active_space(params)

    ###########################################

            excitation = params["excitations"][ex_st]
            nspin = params["nspin"]
            nel = params["nel"]
            occ, occ_alp, occ_bet = excitation_to_qe_occ(params, excitation)
            status = -1

            while status != 0: #for i in xrange(5):
                write_qe_input_first("x%i.scf_wrk.in"%ex_st,occ,occ_alp,occ_bet,nspin)
                exe_espresso(ex_st)
                status = check_convergence("x%i.scf.out" % ex_st) # returns 0 if SCF converges, 1 if not converges
                if status == 0:
                    tot_ene, label, R, grads, mo_pool_alp, mo_pool_bet, params["norb"], params["nel"], params["nat"], params["alat"] = qe_extract("x%i.scf.out" % ex_st, active_space, ex_st, nspin, flag)

                else:
                if params["nspin"] == 2:
                    en_alp = qe_extract_eigenvalues("x%i.save/K00001/eigenval1.xml"%ex_st,nel)
                    en_bet = qe_extract_eigenvalues("x%i.save/K00001/eigenval2.xml"%ex_st,nel)
                    occ_alp = fermi_pop(en_alp,nel,params["nspin"],params["kT"])
                    occ_bet = fermi_pop(en_bet,nel,params["nspin"],params["kT"])
                if params["nspin"] == 1:
                    en_orb = qe_extract_eigenvalues("x%i.save/K00001/eigenval.xml"%ex_st,nel)
                    occ = fermi_pop(en_orb,nel,params["nspin"],params["kT"])


    ###########################################

            homo = params["nel"]/2 +  params["nel"] % 2
            #t.start()
            alp,bet = index_spin(params["excitations"][ex_st], active_space, homo)
            # use excitation object to create proper SD object for different excited state
            sd = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int(alp), Py2Cpp_int(bet) )

            #sd_basis2 = SDList()
            sd_basis.append(sd)
            #t.stop()
            #print "Time to cleade SD object",t.show(),"sec"
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


    ################## Step 2: Initialize molecular system and run MD part with TD-SE and SH####

    rnd = Random() # random number generator object
    syst = []
    syst_mm = []
    el = []

    Ttemp = 0.0
    if params["Nstart"] < 0 and params["MD_type"] == 1: # start NA-MD interacting with thermostat @ t=0
        Ttemp = params["Temperature"]

    # all excitations for each nuclear configuration
    for i in xrange(ninit):
        for i_ex in params["excitations_init"]:
            for itraj in xrange(num_SH_traj):
                df = 0 # debug flag
                # Here we use libra_py module!
                # Utilize the gradients on the ground (0) excited state
                x = init_system.init_system(label_list[i], R_list[i], grad_list[i][0], rnd, Ttemp, params["sigma_pos"], df, "elements.txt")

                # Add the connectivity - needed if we plan to use MM
                if params["is_MM"]:
                    if os.path.exists(params["ent_file"]):
                        LoadMolecule.Load_Molecule(params["U"], x, params["ent_file"], "pdb")
                    else:
                        print "is_MM is set to 1, which means you need to provide the \
                        a file containing the connectivity information for your MM system. \
                        As of now, such file is called = ", params["ent_file"], " but it cannot \
                        be found on your system. Pease check if it exists or set is_MM to 0. \
                        Exiting..."
                        sys.exit(0)

                syst.append(x)

                el.append(Electronic(nstates,i_ex))
    
    # set list of SH state trajectories
    print "run MD"
    run_MD(syst,el,ao_list,e_list,sd_basis_list,params,label_list, Q_list, active_space)
    print "MD is done"


