#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file md.py
# This module implements the functions which execute classical MD.
#
import os
import sys
import math


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *
from exe_espresso import*
from unpack_file import*
from create_qe_input import*
from export_wfc import *
from update_Hvib import *
##############################################################

def run_MD(label,syst,params,wfc):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in,out] syst System object that includes atomic system information.
    # \param[in,out] el   The list of the objects containig electronic DOFs for the nuclear coordinate
    #                     given by syst, but may correspond to differently-prepared coherent
    # wavefunctions (different superpositions or sampling over the wfc phase, initial excitations).
    # Under CPA, the propagation of several such variables corresponds to the same nuclear dynamics,
    # we really don't need to recompute electronic structure for each one, which can be used to 
    # accelerate the computations. Now, if you want to go beyond CPA - just use only one object in
    # the el list and run several copies of the run_MD function to average over initial conditions.
    # Also note that even under the CPA, we need to run this function several times - to sample
    # over initial nuclear distribution
    # \param[in,out] data Data extracted from QE output file, in the dictionary form.
    # \param[in,out] params Input data containing all manual settings and some extracted data.
    # \param[out] test_data  the output data for debugging, in the form of dictionary

    # This function executes classical MD in Libra and electronic structure calculation
    # in Quantum Espresso iteratively.
    #
    # Used in:  main.py/main
    # Open and close energy and trajectory files - this will effectively 
    # make them empty (to remove older info, in case we restart calculations)
    fe = open(params["ene_file"],"w")
    fe.close()
    ft = open(params["traj_file"],"w")
    ft.close()
    
    dt_nucl = params["dt_nucl"]
    Nsnaps = params["Nsnaps"]
    Nsteps = params["Nsteps"]
    n_el = params["nel"]
    n_mo = params["num_MO"]
    no_ex = len(params["excitations"])
    # Create a variable that will contain propagated nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    # Rydberg to Hartree conversion factor
    Ry_to_Ha = 0.5
    MD_type = params["MD_type"]
    kB = 3.1668e-6 # Boltzmann constant in a.u.

    # initialize Thermostat object
    if MD_type == 1: # NVT-MD
        print " Initialize thermostats......"
        therm = Thermostat({"nu_therm":params["nu_therm"], "NHC_size":params["NHC_size"], "Temperature":params["Temperature"], "thermostat_type":params["thermostat_type"]})
        therm.set_Nf_t(3*syst.Number_of_atoms)
        therm.set_Nf_r(0)
        therm.init_nhc()
    epot, ekin, etot, eext = 0.0, 0.0, 0.0, 0.0

    # Create an external Hamiltonian
    #ham = []
    ham_vib = CMATRIX(no_ex,no_ex)
    S_mat = CMATRIX(no_ex,no_ex)
    ham_adi = MATRIX(no_ex,no_ex)
    

    # Run actual calculations
    for ia in xrange(Nsnaps):

        for j in xrange(Nsteps):

            if MD_type == 1: # NVT-MD
                # velocity scaling
                for ka in xrange(3*syst.Number_of_atoms):
                    mol.p[ka] = mol.p[ka] * therm.vel_scale(0.5*dt_nucl)

            # >>>>>>>>>>> Nuclear propagation starts <<<<<<<<<<<<
            mol.propagate_p(0.5*dt_nucl)
            mol.propagate_q(dt_nucl) 

            params["epot1"] = 0.0
            # ======= Compute forces and energies using QUANTUM ESPRESSO ============
            # Running SCF calculation for different excited states, extracting their Energies and Forces
            for i in xrange(len(params["excitations"])):
                write_qe_input(params["qe_inp%i" %i],label,mol,params["excitations"][i],params)
                exe_espresso(params["qe_inp%i" % i], params["qe_out%i" % i] )
                wfc["coeff_%i"%i] = read_qe_wfc("x%i.export/wfc.1"%i, "Kpoint.1", n_el, n_mo)
                params["E%i" %i],label,R, params["Grad%i" %i] = unpack_file(params["qe_out%i" %i],params["qe_debug_print"],0)
                params["epot%i" %i] = Ry_to_Ha*params["E%i" %i]    # total energy from QE, the potential energy acting on nuclei
            epot = params["epot0"]
            epot_ex = params["epot1"]  #to print first excited state energy

            #Updating vibronic hamiltonian
            update_vibronic_hamiltonian(ham_adi, ham_vib, S_mat, wfc, params)

            # Old cofficient is now new coefficient matrix which is used to NACs calculation
            for i in xrange(len(params["excitations"])):
                wfc["coeff_old_%i"%i] = wfc["coeff_%i"%i]

            # Ry/au unit of Force in Quantum espresso
            # So, converting Rydberg to Hartree
            for k in xrange(syst.Number_of_atoms):
                mol.f[3*k]   = -1.0*params["Grad0"][k].x
                mol.f[3*k+1] = -1.0*params["Grad0"][k].y
                mol.f[3*k+2] = -1.0*params["Grad0"][k].z

            ekin = compute_kinetic_energy(mol)

            # Propagate Thermostat
            if MD_type == 1:
                therm.propagate_nhc(dt_nucl, ekin, 0.0, 0.0)

            mol.propagate_p(0.5*dt_nucl)
            # >>>>>>>>>>> Nuclear propagation ends <<<<<<<<<<<<
#            ekin = compute_kinetic_energy(mol)

            if MD_type == 1: # NVT-MD
                # velocity scaling
                for ka in xrange(3*syst.Number_of_atoms):
                    mol.p[ka] = mol.p[ka] * therm.vel_scale(0.5*dt_nucl)


            t = dt_nucl*(ia*Nsteps + j) # simulation time in a.u.
        ################### Printing results ############################
        # >>>>>>>>>>>>>> Compute energies <<<<<<<<<<<<<<<<<<<<<<<<<
        ekin = compute_kinetic_energy(mol)
        etot = ekin + epot

        ebath = 0.0
        if MD_type == 1:
            ebath = therm.energy()

        eext = etot + ebath
        curr_T = 2.0*ekin/(3*syst.Number_of_atoms*kB)

        ################### Printing results ############################

        # Energy
        fe = open(params["ene_file"],"a")
        fe.write("i= %3i ekin= %8.5f  epot= %8.5f  epot_ex= %8.5f etot= %8.5f  eext= %8.5f curr_T= %8.5f\n" % (ia, ekin, epot, epot_ex, etot, eext, curr_T)) 
        syst.set_atomic_q(mol.q)
        syst.print_xyz(params["traj_file"],ia)
        fe.close()
    # input test_data for debugging
    test_data = {}

    return test_data



