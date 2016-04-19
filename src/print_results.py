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
# This module implements the function which print out the results of NA-MD calculation.

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def print_results(iconf,i_ex,i,mol,el,ham,syst,ao,therm,mu,tot_ene,f_pot,params):
    ## find parameters below: 
    # \param[in] iconf number of initial geometry
    # \param[in] i_ex  number of initial excitation
    # \param[in] i    number of snaps (ij=i*nsteps+j)
    # \param[in] mol   Nuclear object
    # \param[in] el    Electronic object
    # \param[in] ham   hamiltonian object
    # \param[in] syst  System object
    # \param[in] ao    Atomic Orbital object  
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
    nstates = el.nstates
    ij = i*params["Nsteps"]

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
    num_tmp = "_"+str(iconf)+"_"+str(i_ex)
    ene_file = params["ene_file_prefix"]+num_tmp+".txt"
    traj_file = params["traj_file_prefix"]+num_tmp+".xyz"
    mu_file = params["mu_file_prefix"]+num_tmp+".txt"
    se_pop_file = params["se_pop_file_prefix"]+num_tmp+".txt"

    ## print 
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

    # Populations            
    fel = open(se_pop_file,"a")

    # Print time
    line_se = "t= %8.5f " % (ij*dt_nucl)

    # Print populations
    for st in xrange(nstates):
        line_se = line_se + " %8.5f " % el.rho(st,st).real

    if print_coherences == 1:
        # Print coherences
        for st in xrange(nstates):
            for st1 in xrange(st):
                line_se = line_se + " %8.5f %8.5f " % (el.rho(st,st1).real, el.rho(st,st1).imag)
             
    line_se = line_se + "\n"
                
    fel.write(line_se)
    fel.close()

