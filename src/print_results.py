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
## \file print_results.py
# This module implements the functions which will print out the results of NA-MD calculation.

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def one_trajectory(i,iconf,i_ex,itsh,mol,el,ham,syst,ao,therm,mu,tot_ene,f_pot,params):
    # This function prints out the auxiliary results for only one trajectory.
    #
    # \param[in] i     number of snaps 
    # \param[in] iconf number of initial geometry
    # \param[in] i_ex  number of initial excitation
    # \param[in] itsh  number of TSH trajectory
    # \param[in] mol   Nuclear object                                                                                                                      
    # \param[in] el    Electronic object                                                                                                                    
    # \param[in] ham   hamiltonian object                                                                                                                   
    # \param[in] syst  System object                                                                                                                        
    # \param[in] ao    Atomic Orbital basis object                                                                                                          
    # \param[in] therm thermostat object                                                                                                                    
    # \param[in] mu    list of dipole moment                                                                                                                
    # \param[in] tot_ene total energy from GAMESS output                                                                                                    
    # \param[in] f_pot flag for potential energy : option 0 - Ehrenfest, 1 - TSH                                                                             
    # \param[in] params list of input parameters from run.py                                                                                                 
    #                                                                                                                                                        
    # Used in:  main.py/main/run_MD 

    kB = 3.166811429e-6 # Boltzmann constant in hartree unit                                                                                                 
    dt_nucl = params["dt_nucl"]
    MD_type = params["MD_type"]
    print_coherences = params["print_coherences"]
    nstates = len(params["excitations"])
    ij = i*params["Nsteps"]
    Ntsh = params["Ntsh"]

    # Re-compute energies, to print                                                                                                           
    epot = tot_ene + compute_potential_energy(mol, el, ham, f_pot)
    print "epot = ", epot
    
    ekin = compute_kinetic_energy(mol)
    etot = ekin + epot

    ebath = 0.0
    if MD_type == 1:
        ebath = therm.energy()
    eext = etot + ebath
    curr_T = 2.0*ekin/(3*syst.Number_of_atoms*kB)

    # set file name
    num_tmp = "_"+str(iconf)+"_"+str(i_ex)+"_"+str(itsh)
    ene_file = params["ene_file_prefix"]+num_tmp+".txt"
    traj_file = params["traj_file_prefix"]+num_tmp+".xyz"
    mu_file = params["mu_file_prefix"]+num_tmp+".txt"

    ##print 
    # Geometry
    syst.set_atomic_q(mol.q)
    syst.print_xyz(traj_file,i)

    # Energy
    fe = open(ene_file,"a")
    fe.write("t= %8.5f ekin= %8.5f  epot= %8.5f  etot= %8.5f  eext= %8.5f curr_T= %8.5f\n" % (ij*dt_nucl, ekin, epot, etot, eext,curr_T))
    fe.close()
        
    # Dipole moment components
    fm = open(mu_file,"a")
    line = "t= %8.5f " % (ij*dt_nucl)
    for k in xrange(len(ao)):
        line = line + " %8.5f %8.5f %8.5f " % (mu[0].get(k,k),mu[1].get(k,k),mu[2].get(k,k))
    line = line + "\n"
    fm.write(line)
    fm.close()


def auxiliary(i,mol,el,ham,syst,ao,therm,mu,tot_ene,f_pot,params):
    # This function prints out auxiliary results for a bunch of trajectories:
    # one trajctory contains the index of initial geometry, initial excitation state
    # (TSH trajectory, optionally)
    #
    # \param[in] i     number of snaps (ij=i*nsteps+j)
    # \param[in] mol   Nuclear objects
    # \param[in] el    Electronic objects
    # \param[in] ham   hamiltonian objects
    # \param[in] syst  System objects
    # \param[in] ao    Atomic Orbital basis objects  
    # \param[in] therm thermostat objects
    # \param[in] mu    list of dipole moments  
    # \param[in] tot_ene total energy from GAMESS outputs
    # \param[in] f_pot flag for potential energy : option 0 - Ehrenfest, 1 - TSH
    # \param[in] params list of input parameters from run.py
    #
    # Used in:  main.py/main/run_MD 

    nconfig = params["Ngeo"]
    nstates = len(params["excitations"])
    Ntsh = params["Ntsh"]

    for iconf in xrange(nconfig):
        for i_ex in xrange(nstates):
            for itsh in xrange(Ntsh):

                cnt = iconf*nstates*Ntsh + i_ex*Ntsh + itsh
                one_trajectory(i,iconf,i_ex,itsh,mol[cnt],el[cnt],ham[cnt],syst[cnt],ao[cnt],therm[cnt],mu[cnt],tot_ene[cnt],f_pot,params)


def pops_ave_TSH_traj(i,el,params):
    # This function prints out SE and SH populations averaged over TSH trajectories(only one trajectory without SH calculation).
    #                                                                                                                                                        
    # \param[in] i     number of snaps (ij=i*nsteps+j)                                                                                                       
    # \param[in] mol   Nuclear objects                                                                                                                       
    # \param[in] el    Electronic objects

    nconfig = params["Ngeo"]
    nstates = len(params["excitations"])
    Ntsh = params["Ntsh"]
    SH_type = params["tsh_method"]
    ij = i*params["Nsteps"]
    dt_nucl = params["dt_nucl"]
    print_coherences = params["print_coherences"]

    l_se_pop = []
    l_sh_pop = []

    for iconf in xrange(nconfig):
        for i_ex in xrange(nstates):

            #***********  SE population **************

            se_pop = []
            se_coh_re = []
            se_coh_im = []

            # average SE populations
            for st in xrange(nstates):
                p_tmp = 0.0
                for itsh in xrange(Ntsh):
                    cnt = iconf*nstates*Ntsh + i_ex*Ntsh + itsh
                    p_tmp += el[cnt].rho(st,st).real
                se_pop.append(p_tmp/float(Ntsh))
                l_se_pop.append(p_tmp/float(Ntsh)) # store SE pops

            # average coherences
            if print_coherences == 1:
                for st in xrange(nstates):
                    for st1 in xrange(st):
                        c1_tmp = 0.0; c2_tmp = 0.0
                        for itsh in xrange(Ntsh):
                            cnt = iconf*nstates*Ntsh + i_ex*Ntsh + itsh
                            c1_tmp += el[cnt].rho(st,st1).real
                            c2_tmp += el[cnt].rho(st,st1).imag                            
                        se_coh_re.append(c1_tmp/float(Ntsh))
                        se_coh_im.append(c2_tmp/float(Ntsh))

            # set file name
            num_tmp = "_"+str(iconf)+"_"+str(i_ex)
            se_pop_file = params["se_pop_file_prefix"]+num_tmp+".txt"

            # Populations
            fel = open(se_pop_file,"a")

            # Print time
            line_se = "t= %8.5f " % (ij*dt_nucl)

            # Print populations
            for st in xrange(nstates):
                line_se = line_se + " %8.5f " % se_pop[st]

            if print_coherences == 1:
                # Print coherences
                cnt = 0
                for st in xrange(nstates):
                    for st1 in xrange(st):
                        line_se = line_se + " %8.5f %8.5f " % (se_coh_re[cnt], se_coh_im[cnt])
                        cnt += 1

            line_se = line_se + "\n"

            fel.write(line_se)
            fel.close()

            #***********  SH population **************
            if SH_type>=1:

                # evaluate TSH probabilities
                tsh_probs = [0.0]*nstates
                for itsh in xrange(Ntsh): # count number of trajectories
                    cnt = iconf*nstates*Ntsh + i_ex*Ntsh + itsh
                    tsh_probs[el[cnt].istate] += 1.0
            
                for st in xrange(nstates):
                    tsh_probs[st] = tsh_probs[st]/float(Ntsh)
                    l_sh_pop.append(tsh_probs[st]) # store SH pops

                # set file name
                num_tmp = "_"+str(iconf)+"_"+str(i_ex)
                sh_pop_file = params["sh_pop_file_prefix"]+num_tmp+".txt"

                # Populations     
                fel1 = open(sh_pop_file,"a")

                # Print time      
                line_sh = "t= %8.5f " % (ij*dt_nucl)

                #Print populations
                for st in xrange(nstates):
                    line_sh = line_sh + " %8.5f " % (tsh_probs[st])

                line_sh = line_sh + "\n"

                fel1.write(line_sh)
                fel1.close()

    return l_se_pop, l_sh_pop


def pops_ave_geometry(i,se_pop,sh_pop,params):
    # This function prints out SE and SH populations averaged over initial geometry;
    # they have the index of initial excitation state.
    #                                                                                                                                                        
    # \param[in] i     number of snaps (ij=i*nsteps+j)                                                                                                       
    # \param[in] mol   Nuclear objects                                                                                                                       
    # \param[in] el    Electronic objects

    nconfig = params["Ngeo"]
    nstates = len(params["excitations"])
    SH_type = params["tsh_method"]
    ij = i*params["Nsteps"]
    dt_nucl = params["dt_nucl"]

    #print "length of se_pop is %d" %(len(se_pop))
    #print "length of sh_pop is %d" %(len(sh_pop))

    for i_ex in xrange(nstates):

        # set file name                                                                                                                               
        se_pop_file = params["se_pop_ex_file_prefix"]+str(i_ex)+".txt"
        sh_pop_file = params["sh_pop_ex_file_prefix"]+str(i_ex)+".txt"

        # Populations                                                                                                            
        fse = open(se_pop_file,"a")
        fsh = open(sh_pop_file,"a")

        # Print time                                                                                                                                  
        line_se = "t= %8.5f " % (ij*dt_nucl)
        line_sh = "t= %8.5f " % (ij*dt_nucl)

        for st in xrange(nstates):

            se_sum = 0.0; sh_sum = 0.0;
            for iconf in xrange(nconfig):
                cnt = iconf*nstates*nstates + i_ex*nstates + st

                se_sum += se_pop[cnt]
                sh_sum += sh_pop[cnt]

            #Print populations
            line_se = line_se + " %8.5f " % (se_sum)
            line_sh = line_sh + " %8.5f " % (sh_sum)

        line_se = line_se + "\n"
        line_sh = line_sh + "\n"

        fse.write(line_se); fse.close();
        fsh.write(line_sh); fsh.close();
