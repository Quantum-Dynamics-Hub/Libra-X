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

## \file md.py
# This module implements the functions which execute GAMESS, execute classical MD coupled to 
# Nose-Hoover thermostat, execute TD-SE and SH calculations, and set initial system.
#

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
from hamiltonian_vib import *
import print_results
#import print_results_qe


##############################################################

def init_files(params):
    ## 
    # This function initializes files.(make empty files)
    # \param[in] params Input data containing all manual settings and some extracted data.
    #                   Here, file prefix names will be used.  
    #
    # Used in md.py/run_MD 
    
    nconfig = params["nconfig"]
    nstates = len(params["excitations"])
    num_SH_traj = params["num_SH_traj"]

    # define prefixes
    traj_file_prefix = params["res"]+"md"
    ene_file_prefix = params["res"]+"ene"
    mu_file_prefix = params["res"]+"mu"
    se_pop_file_prefix = params["res"]+"se_pop"
    sh_pop_file_prefix = params["res"]+"sh_pop"
    se_pop_ex_file_prefix = params["res"]+"se_pop_ex"
    sh_pop_ex_file_prefix = params["res"]+"sh_pop_ex"

    for i in xrange(nconfig):
        for i_ex in xrange(nstates):
            index0 = "_"+str(i)+"_"+str(i_ex)

            se_pop_file = se_pop_file_prefix+index0+".txt"
            sh_pop_file = sh_pop_file_prefix+index0+".txt"
            fel = open(se_pop_file,"w"); fel.close();
            fel = open(sh_pop_file,"w"); fel.close();

            if params["print_aux_results"] == 1:
                for itraj in xrange(num_SH_traj):
                    index = index0+"_"+str(itraj)
                    ene_file = ene_file_prefix+index+".txt"
                    traj_file = traj_file_prefix+index+".xyz"
                    mu_file = mu_file_prefix+index+".txt"

                    fe = open(ene_file,"w"); fe.close();
                    ft = open(traj_file,"w"); ft.close();

                    if params["flag_ao"] == 1:
                        fm = open(mu_file,"w"); fm.close();

    for i_ex in xrange(nstates):

        se_pop_file = se_pop_ex_file_prefix+str(i_ex)+".txt"
        sh_pop_file = sh_pop_ex_file_prefix+str(i_ex)+".txt"
        fel = open(se_pop_file,"w"); fel.close();
        fel = open(sh_pop_file,"w"); fel.close();



def run_MD(syst,el,ao,E,sd_basis,params,label,Q, active_space):
    ##
    # This function handles a swarm of trajectories.
    # When NA-MD is utilized (by specifying the TSH method), we use the CPA with isotropic
    # velocity rescaling
    #
    # \param[in,out] syst list of  System objects that include atomic system information.
    # \param[in,out] el list of object containig electronic DOFs for the nuclear coordinate
    # given by syst. I have decided to go back a bit - one set of electronic DOF per set of
    # nuclear DOF. This is also needed when we do the velocity rescaling, even if we use the
    # ground state forces for propagation. This also brings a conceptual clarity
    # 
    # \param[in,out] ao   list pf Atomic orbital basis
    # \param[in,out] E    list of Molecular orbital energies
    # \param[in,out] sd_basis list of lists of MO-LCAO coefficients, such that 
    # sd_basis[i] is the list of CMATRIX objects representing SD for the "trajectory/initial condition/realization" i. Then sd_basis[i][j] corresponds to the determinant j of the initial condition i
    # \param[in,out] params list of input parameters from (gms/qe)_run.py , which will get some changes here.
    # \param[in] label    list of atomic labels e.g. H, He, Li, etc...
    # \param[in] Q        list of atomic charges
    # \param[in] active_space The list of indices (starting from 1) of the MOs to include in
    # calculations (and to read from the QE output files)


    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively and simulates excited electron dynamics with MF and SH way. 
    # It outputs MD trajectories, Energy evolutions, and SE and SH populations.
    #
    # Used in:  main.py/main

    rnd = Random()

    dt_nucl = params["dt_nucl"]
    el_mts = params["el_mts"] # multiple time stepping algorithm for electronic DOF propagation
    if el_mts < 1:
        print "Error in run_MD: el_mts must be positive integer"
        print "Value given = ", el_mts
        print "Exiting..."
        sys.exit(0)
    dt_elec = dt_nucl/float(el_mts)

    nconfig = params["nconfig"]
    flag_ao = params["flag_ao"]
    Nsnaps = params["Nsnaps"]
    Nsteps = params["Nsteps"]
    nstates = len(params["excitations"])
    print_coherences = params["print_coherences"]
    MD_type = params["MD_type"]
    SH_type = params["tsh_method"]

    # a flag for potential energy (Ehrenfest or SH)
    f_pot = 0 # Default: Ehrenfest
    if SH_type >= 1: # use SH potential
        f_pot = 1

    # TSH trajectories
    num_SH_traj = 1
    if SH_type >= 1: # use SH potential                                                                                                                   
        num_SH_traj = params["num_SH_traj"]

    #=============== Initialization =======================

    # Open and close energy and trajectory files - this will effectively
    # make them empty (to remove older info, in case we restart calculations)

    init_files(params)
    
    # prepare objects for MD
    ntraj = len(syst)
    nnucl = 3*syst[0].Number_of_atoms
    verbose = 0

    ham, ham_adi, d1ham_adi, ham_vib = init_ensembles.init_ext_hamiltonians(ntraj, nnucl, nstates, verbose)
    mol = init_ensembles.init_mols(syst, ntraj, nnucl, verbose)
    therm = init_ensembles.init_therms(ntraj, nnucl, params, verbose)


    # Initialize forces and Hamiltonians **********************************************
    #epot = data["tot_ene"]  # total energy from GAMESS which is the potential energy acting on nuclei
    #write_gms_inp(data, params, mol)
    #exe_gamess(params)
    #Grad, data, E_mol, D_mol, E_mol_red, D_mol_red = gamess_to_libra(params, ao, E, C, 0) # this will update AO and gradients
    #Hvib, D_SD = vibronic_hamiltonian(params,E_mol_red,D_mol_red,0) # create vibronic hamiltonian

    #sys.exit(0) # DEBUG!!!

    print "Starting propagation"

    #=============== Propagation =======================

    epot, ekin, etot, eext = 0.0, 0.0, 0.0, 0.0
    ens_sz = nconfig * nstates * num_SH_traj
    epot = [0.0]*ens_sz
    ekin = [0.0]*ens_sz
    etot = [0.0]*ens_sz
    eext = [0.0]*ens_sz
    mu = []
    for i in xrange(ens_sz):
        mu.append(MATRIX())


    for i in xrange(Nsnaps):   # number of printouts

#        tot_ene = []; mu = []; # initialize lists

        for j in xrange(Nsteps):   # number of integration steps per printout
            ij = i*Nsteps + j

            for iconf in xrange(nconfig):     # all initial nuclear configurations

                for i_ex in xrange(nstates):  # consider initial excitations to be on all the basis
                                              # states - this may be unnecessary for all cases, 
                                              # so we may want to make this part customizable

                    for itraj in xrange(num_SH_traj): # all stochastic SH realizations

                        cnt = iconf*nstates*num_SH_traj + i_ex*num_SH_traj + itraj

                        print "Initial geometry %i, initial excitation %i, tsh trajectory %i"%(iconf,i_ex,itraj)

                        # Electronic propagation: half-step
                        for k in xrange(el_mts):
                            el[cnt].propagate_electronic(0.5*dt_elec, ham[cnt])

                        # >>>>>>>>>>> Nuclear propagation starts <<<<<<<<<<<<
                        # Optional thermostat            
                        if MD_type == 1: # NVT-MD
                            for k in xrange(3*syst[cnt].Number_of_atoms):
                                mol[cnt].p[k] = mol[cnt].p[k] * therm[cnt].vel_scale(0.5*dt_nucl)

                        mol[cnt].propagate_p(0.5*dt_nucl)
                        mol[cnt].propagate_q(dt_nucl)


                        # ======= Compute forces and energies using external package ============
                        #tot_ene0 = 0.0
                        nac = CMATRIX()
                        all_grads = []
                        opt = 0 # default

                        if params["interface"]=="GAMESS":
                            opt = 0 # use 1-electron wavefunctions
                           
                            write_gms_inp(label[cnt], Q[cnt], params, mol[cnt])
                            exe_gamess(params)
                       
                            # update AO, MO, and gradients
                            E_SD, nac, sd_basis[cnt], all_grads, mu[cnt] = gamess_to_libra(params, ao[cnt], E[cnt], \
                                                                                                              sd_basis[cnt], active_space,str(ij) )
                            #tot_ene.append(tot_ene0); mu.append(mu0); # store total energy and dipole moment

                        elif params["interface"]=="QE":
                            opt = 1 # use true SD wavefunctions

                            # update MO and gradients
                            E_mol_red, nac, E[cnt], sd_basis[cnt], all_grads = qe_to_libra(params, E[cnt], sd_basis[cnt], label[cnt], mol[cnt], str(ij), active_space)
                            #tot_ene.append(E[cnt])

                        # ============== Common blocks ==================


                        # update forces
                        for k in xrange(syst[cnt].Number_of_atoms):
                            for st in xrange(nstates):
                                d1ham_adi[cnt][3*k+0].set(st,st,all_grads[st][k].x)
                                d1ham_adi[cnt][3*k+1].set(st,st,all_grads[st][k].y)
                                d1ham_adi[cnt][3*k+2].set(st,st,all_grads[st][k].z)

                        # Update the matrices that are bound to the Hamiltonian 
                        # Compose electronic and vibronic Hamiltonians
                        update_vibronic_hamiltonian(ham_adi[cnt], ham_vib[cnt], params, E_SD,nac, str(ij), opt)

           
                        # update potential energy
                        # according to new convention (yet to be implemented for GMS and need to
                        # check for QE - the Hamiltonians will contain the total energies of 
                        # excited stes, so no need for reference energy)
                        epot[cnt] = compute_forces(mol[cnt], el[cnt], ham[cnt], f_pot)  #  f_pot = 0 - Ehrenfest, 1 - TSH
                        ekin[cnt] = compute_kinetic_energy(mol[cnt])
                        etot[cnt] = epot[cnt] + ekin[cnt]
                        eext[cnt] = etot[cnt]
          
                        if MD_type == 1:
                            therm[cnt].propagate_nhc(dt_nucl, ekin[cnt], 0.0, 0.0)

                        mol[cnt].propagate_p(0.5*dt_nucl)

                        # optional thrmostat
                        if MD_type == 1: # NVT-MD
                            for k in xrange(3*syst[cnt].Number_of_atoms):
                                mol[cnt].p[k] = mol[cnt].p[k] * therm[cnt].vel_scale(0.5*dt_nucl)

                            eext[cnt] = eext[cnt] + therm[cnt].energy() 

                        # >>>>>>>>>>> Nuclear propagation ends <<<<<<<<<<<<

                        # Electronic propagation: half-step
                        for k in xrange(el_mts):
                            el[cnt].propagate_electronic(0.5*dt_elec, ham[cnt])

                        #******** end of itsh loop
                    #********* end of i_ex loop
                #********* end of iconf loop
            #***** End of TD-SE propagation for this step
                    
            ############ Add surface hopping ######################
            print "Before TSH"

            if SH_type>=1:
                tsh.surface_hopping_cpa2(mol, el, ham, rnd, params)

            ################### END of TSH ##########################
            print "Finished TSH"

        #************ end of j loop - all steps for this snap

        ################### Printing results ############################
        # print out SE and SH populations
        se_pop, sh_pop = print_results.pops_ave_TSH_traj(i,el,params)
        print_results.pops_ave_geometry(i,nstates,se_pop,sh_pop,params)

        # print auxiliary files: MD, Energy, and dipole moment trajectories
        if params["print_aux_results"]==1:
            print_results.print_ens_traj(i,mol,syst,mu,epot,ekin,etot,eext,params)

        print "       ********* %i snap ends ***********" % i
        print 



