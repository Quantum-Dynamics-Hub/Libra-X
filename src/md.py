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

from create_gamess_input import *
from gamess_to_libra import *
from vibronic_hamiltonian import *
from create_MD_objects import *
from surface_hopping import *
from print_results import *

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


##############################################################

def exe_gamess(params):
    ##
    # This is a function that call GAMESS execution on the compute node
    # \param[in] params Input data containing all manual settings and some extracted data.
    #
    # Used in main.py/main
    #         main.py/main/nve_MD

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

def run_MD(syst,el,ao,E,C,params,label,Q):
    ##
    # This function handles a SINGLE trajectory.
    # When NA-MD is utilized (by specifying the TSH method), we use the CPA with isotropic
    # velocity rescaling
    #
    # \param[in,out] syst a System object that includes atomic system information.
    # \param[in,out] el The the object containig electronic DOFs for the nuclear coordinate
    # given by syst. I have decided to go back a bit - one set of electronic DOF per set of
    # nuclear DOF. This is also needed when we do the velocity rescaling, even if we use the
    # ground state forces for propagation. This also brings a conceptual clarity
    # 
    # \param[in,out] ao   Atomic orbital basis
    # \param[in,out] E    Molecular orbital energies
    # \param[in,out] C    MO-LCAO coefficients
    # \param[in] label    atomic label e.g. H, He, Li, etc...
    # \param[in] Q        atomic charge

    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively and simulates excited electron dynamics with MF and SH way. 
    # It outputs MD trajectories, Energy evolutions, and SE and SH populations.
    #
    # Used in:  main.py/main
    
    dt_nucl = params["dt_nucl"]
    el_mts = params["el_mts"] # multiple time stepping algorithm for electronic DOF propagation
    if el_mts < 1:
        print "Error in run_MD: el_mts must be positive integer"
        print "Value given = ", el_mts
        print "Exiting..."
        sys.exit(0)
    dt_elec = dt_nucl/float(el_mts)

    nconfig = params["nconfig"]
    Nsnaps = params["Nsnaps"]
    Nsteps = params["Nsteps"]
    nstates = len(params["excitations"])
    print_coherences = params["print_coherences"]
    MD_type = params["MD_type"]
    SH_type = params["SH_type"]

    # a flag for potential energy (Ehrenfest or SH)
    f_pot = 0 # Default: Ehrenfest
    if SH_type >= 1: # use SH potential
        f_pot = 1

    #=============== Initialization =======================

    # Open and close energy and trajectory files - this will effectively
    # make them empty (to remove older info, in case we restart calculations)
    for i in xrange(nconfig):
        for i_ex in xrange(nstates):

            num_tmp = "_"+str(i)+"_"+str(i_ex)

            ene_file = params["ene_file_prefix"]+num_tmp+".txt"                                                                                    
            traj_file = params["traj_file_prefix"]+num_tmp+".xyz"                                                                                  
            mu_file = params["mu_file_prefix"]+num_tmp+".txt"                                                                                      
            se_pop_file = params["se_pop_file_prefix"]+num_tmp+".txt" 
            
            fe = open(ene_file,"w")
            fe.close()
            ft = open(traj_file,"w")
            ft.close()
            fm = open(mu_file,"w")
            fm.close()   
            fel = open(se_pop_file,"w")
            fel.close()

    # prepare objects for MD
    ham, ham_adi, d1ham_adi, ham_vib, mol, therm = md_objects(syst,nstates,params)

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
    SH_states = [] # stores SH states

    for i in xrange(Nsnaps):
        #syst.set_atomic_q(mol.q)
        #syst.print_xyz(params["traj_file"],i)       

        tot_ene = []; mu = []; # initialize lists

        for j in xrange(Nsteps):

            ij = i*Nsteps + j

            for iconf in xrange(nconfig):
                for i_ex in xrange(nstates):

                    cnt = iconf*nstates + i_ex

                    print "Initial geometry %i, initial excitation %i"%(iconf,i_ex)

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

                    # ======= Compute forces and energies using GAMESS ============
                    write_gms_inp(label[cnt], Q[cnt], params, mol[cnt])
                    exe_gamess(params)
                    
                    #sys.exit(0)
                    # update AO and gradients
                    tot_ene0, Grad, mu0, E_mol, D_mol, E_mol_red, D_mol_red = gamess_to_libra(params, ao[cnt], E[cnt], C[cnt], str(ij))
                    
                    tot_ene.append(tot_ene0); mu.append(mu0); # store total energy and dipole moment 
                    #========= Update the matrices that are bound to the Hamiltonian =========
                    #Compose electronic and vibronic Hamiltonians
                    update_vibronic_hamiltonian(ham_adi[cnt], ham_vib[cnt], params,E_mol_red,D_mol_red, str(ij))

                    print "Addresses of the ham matrices"
                    print "ham_adi = ", ham_adi[cnt]
                    print "ham_vib = ", ham_vib[cnt]
                    print "ham_adi "; ham_adi[cnt].show_matrix();
                    print "ham_vib "; ham_vib[cnt].show_matrix();
                    #sys.exit(0)

                    # update forces
                    for k in xrange(syst[cnt].Number_of_atoms):
                        for st in xrange(nstates):
                            d1ham_adi[cnt][3*k+0].set(st,st,Grad[k].x)
                            d1ham_adi[cnt][3*k+1].set(st,st,Grad[k].y)
                            d1ham_adi[cnt][3*k+2].set(st,st,Grad[k].z)
           
                    # update potential energy
                    epot = tot_ene0 + compute_forces(mol[cnt], el[cnt], ham[cnt], f_pot)  #  f_pot = 0 - Ehrenfest, 1 - TSH

                    print "epot= ", epot
                    #sys.exit(0)

                    ekin = compute_kinetic_energy(mol[cnt])
                    etot = epot + ekin
          
                    if MD_type == 1:
                        therm[cnt].propagate_nhc(dt_nucl, ekin, 0.0, 0.0)

                    mol[cnt].propagate_p(0.5*dt_nucl)

                    # optional thrmostat
                    if MD_type == 1: # NVT-MD
                        for k in xrange(3*syst[cnt].Number_of_atoms):
                            mol[cnt].p[k] = mol[cnt].p[k] * therm[cnt].vel_scale(0.5*dt_nucl)

                    # >>>>>>>>>>> Nuclear propagation ends <<<<<<<<<<<<


                    # Electronic propagation: half-step
                    for k in xrange(el_mts):
                        el[cnt].propagate_electronic(0.5*dt_elec, ham[cnt])

                    
                    ############ Add surface hopping ######################

                    print "Before TSH";# sys.exit(0)

                    if SH_type>=1:

                        mol[cnt], el[cnt], ham[cnt] = surface_hopping(mol[cnt],el[cnt],ham[cnt],params)

                    ################### END of TSH ##########################
         
                    print "Finished TSH";# sys.exit(0)

                    #********* end of i_ex loop

        #************ end of j loop

        # check d1ham_adi objects different forces
        #for l in xrange(len(ham)):
        #    print "%d th hamiltonian" %(l)
        #    for k in xrange(3*syst[l].Number_of_atoms):
        #        print "d1ham_adi[",k,"]= ", d1ham_adi[l][k].show_matrix()

        ################### Printing results ############################

        for iconf in xrange(nconfig):
            for i_ex in xrange(nstates):

                cnt = iconf*nstates + i_ex
                print_results(iconf,i_ex,i,mol[cnt],el[cnt],ham[cnt],syst[cnt],ao[cnt],therm[cnt],mu[cnt],tot_ene[cnt],f_pot,params)

        print "       ********* %i snap ends ***********" % i
        print 

    #print "SH_states = ", SH_states

    return SH_states


