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

def one_trajectory(i,iconf,i_ex,itraj,mol,el,ham,syst,therm,tot_ene,f_pot,params):
    # This function prints out the auxiliary results for only one trajectory.
    #
    # \param[in] i          snap index 
    # \param[in] iconf      initial geometry index
    # \param[in] i_ex       initial excitation index
    # \param[in] itraj      TSH trajectory index
    # \param[in] mol[n]     contains the nuclear configuration of n-th trajectory
    # \param[in] el[n]      contains electronic DOFs of n-th trajectory
    # \param[in] ham[n]     hamiltonian in n-th trajectory
    # \param[in] syst[n]    System in n-th trajectory
    # \param[in] therm[n]   thermostat in n-th trajectory
    # \param[in] tot_ene[n] total energy from GAMESS output in n-th trajectory
    # \param[in] f_pot      flag for potential energy : option 0 - Ehrenfest, 1 - TSH 
    # \param[in] params     list of input parameters from run.py
    #                                                                                                                                                        
    # Used in:  print_results.py/auxiliary

    kB = 3.166811429e-6 # Boltzmann constant in hartree unit                                                                                                 
    dt_nucl = params["dt_nucl"]
    MD_type = params["MD_type"]
    print_coherences = params["print_coherences"]
    nstates = el.nstates
    ij = i*params["Nsteps"]

    # Re-compute energies, to print
    epot = tot_ene.get(i_ex,i_ex) + compute_potential_energy(mol, el, ham, f_pot)
    print "epot = ", epot
    
    ekin = compute_kinetic_energy(mol)
    etot = ekin + epot

    ebath = 0.0
    if MD_type == 1:
        ebath = therm.energy()
    eext = etot + ebath
    curr_T = 2.0*ekin/(3*syst.Number_of_atoms*kB)

    # set file name
    num_tmp = "_"+str(iconf)+"_"+str(i_ex)+"_"+str(itraj)
    ene_file = params["ene_file_prefix"]+num_tmp+".txt"
    traj_file = params["traj_file_prefix"]+num_tmp+".xyz"

    ##print 
    # Geometry
    syst.set_atomic_q(mol.q)
    syst.print_xyz(traj_file,i)

    # Energy
    fe = open(ene_file,"a")
    fe.write("t= %8.5f ekin= %8.5f  epot= %8.5f  etot= %8.5f  eext= %8.5f curr_T= %8.5f\n" % (ij*dt_nucl, ekin, epot, etot, eext,curr_T))
    fe.close()
        


def auxiliary(i,mol,el,ham,syst,therm,tot_ene,f_pot,params):
    # This function prints out auxiliary results for a bunch of trajectories:
    # one trajctory contains the index of initial geometry, initial excitation state
    # (TSH trajectory, optionally)
    #
    # \param[in] i       snap index
    # \param[in] mol     list of objects containing Nuclear DOF's whose size is nconfig*nstates*num_SH_traj. Lists below have the same size as this.
    # \param[in] el      list of objects containing Electronic DOF's 
    # \param[in] ham     list of objects containing Nuclear and Electornic DOF's (hamiltonian in this system)
    # \param[in] syst    list of objects containing atomic system information
    # \param[in] therm   list of objects containing thermostat variables
    # \param[in] tot_ene list of total energies from GAMESS outputs
    # \param[in] f_pot   flag for potential energy : option 0 - Ehrenfest, 1 - TSH
    # \param[in] params  list of input parameters from run.py
    #
    # Used in:  md.py/run_MD 

    nconfig = params["nconfig"]
    nstates = el[0].nstates
    num_SH_traj = params["num_SH_traj"]

    for iconf in xrange(nconfig):
        for i_ex in xrange(nstates):
            for itraj in xrange(num_SH_traj):

                cnt = iconf*nstates*num_SH_traj + i_ex*num_SH_traj + itraj
                one_trajectory(i,iconf,i_ex,itraj,mol[cnt],el[cnt],ham[cnt],syst[cnt],therm[cnt],tot_ene[cnt],f_pot,params)


