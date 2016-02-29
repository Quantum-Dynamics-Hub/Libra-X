#*********************************************************************************
#* Copyright (C) 2015 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file nve.py
# This module implements the functions which execute classical MD.
#

from create_gamess_input import *
from gamess_to_libra import *
from vibronic_hamiltonian import *

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_chemobjects_path"])
sys.path.insert(1,os.environ["libra_hamiltonian_path"])
sys.path.insert(1,os.environ["libra_dyn_path"])
#sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/hamiltonian_atomistic")

from libmmath import *
from libchemobjects import *
from libhamiltonian import *
from libdyn import *
#from libhamiltonian_atomistic import *
from LoadPT import * # Load_PT

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
    scr_dir = os.environ['SLURMTMPDIR']
    os.system("/usr/bin/time rungms.slurm %s 01 %s > %s" % (inp,nproc,out))

    # delete the files except input and output ones to do another GAMESS calculation.
    os.system("rm *.dat")              
    os.system("rm -r %s/*" %(scr_dir)) 


def run_MD(syst,el,ao,E,C,data,params):
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
    # \param[in,out] ao   Atomic orbital basis
    # \param[in,out] E    Molecular orbital energies
    # \param[in,out] C    MO-LCAO coefficients
    # \param[in,out] data Data extracted from GAMESS output file, in the dictionary form.
    # \param[in,out] params Input data containing all manual settings and some extracted data.

    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively.
    #
    # Used in:  main.py/main


    # Open and close energy and trajectory files - this will effectively 
    # make them empty (to remove older info, in case we restart calculations)
    fe = open(params["ene_file"],"w")
    fe.close()
    ft = open(params["traj_file"],"w")
    ft.close()
    
    dt_nucl = params["dt_nucl"]
    el_mts = params["el_mts"] # multiple time stepping algorithm for electronic DOF propagation
    if el_mts<1:
        print "Error in run_MD: el_mts must be positive integer"
        print "Value given = ", el_mts
        print "Exiting..."
        sys.exit(0)

    dt_elec = dt_nucl/float(el_mts)
    Nsnaps = params["Nsnaps"]
    Nsteps = params["Nsteps"]

    nstates = len(params["excitations"])

    print_coherences = params["print_coherences"]

    for k in xrange(nstates):
        tmp = params["se_pop_prefix"] + "se_pop_" + str(k)
        fel = open(tmp,"w")
        fel.close()

    # Create a variable that will contain propagated nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    # Debug printing
    if 0==1:
        for i in xrange(syst.Number_of_atoms):
            print "mass m=",mol.mass[3*i], mol.mass[3*i+1], mol.mass[3*i+2]
            print "coordinates q = ", mol.q[3*i], mol.q[3*i+1], mol.q[3*i+2]
            print "momenta p= ", mol.p[3*i], mol.p[3*i+1], mol.p[3*i+2]
            print "forces f= ",  mol.f[3*i], mol.f[3*i+1], mol.f[3*i+2]
            print "********************************************************"

    # Run actual calculations
    for i in xrange(Nsnaps):

        syst.set_atomic_q(mol.q)
        syst.print_xyz(params["traj_file"],i)

        for j in xrange(Nsteps):

            ij = i*Nsteps + j

            if ij > 0: # pass this function at t=0
                # Electronic propagation: half-step
                for k in xrange(el_mts):
                    for i_ex in range(0,nstates):  # loop over all initial excitations
                        propagate_electronic(0.5*dt_elec, el[i_ex], Hvib)

            # >>>>>>>>>>> Nuclear propagation starts <<<<<<<<<<<<

            mol.propagate_p(0.5*dt_nucl)
            mol.propagate_q(dt_nucl) 
          
            # ======= Compute forces and energies using GAMESS ============
            write_gms_inp(data, params, mol)
            exe_gamess(params)         

            Grad, data, E_mol, D = gamess_to_libra(params, ao, E, C, ij) # this will update AO and gradients
            Hvib = vibronic_hamiltonian(params,E_mol,D)

            epot = data["tot_ene"]         # total energy from GAMESS which is the potential energy acting on nuclei

            for k in xrange(syst.Number_of_atoms):
                mol.f[3*k]   = -Grad[k][0]
                mol.f[3*k+1] = -Grad[k][1]
                mol.f[3*k+2] = -Grad[k][2]

            mol.propagate_p(0.5*dt_nucl)

            # >>>>>>>>>>> Nuclear propagation ends <<<<<<<<<<<<

            # Electronic propagation: half-step
            for k in xrange(el_mts):
                for i_ex in range(0,nstates):  # loop over all initial excitations
                    propagate_electronic(0.5*dt_elec, el[i_ex], Hvib)


            ekin = compute_kinetic_energy(mol)
            t = dt_nucl*ij # simulation time in a.u.

        ################### Printing results ############################

        fe = open(params["ene_file"],"a")
        fe.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, ekin+epot)) 
        fe.close()
        
        for k in xrange(nstates):
            tmp = params["se_pop_prefix"] + "se_pop_" + str(k)
            fel = open(tmp,"a")

            # Print time
            line = "t= %8.5f " % t

            # Print populations
            for st in xrange(nstates):
                line = line + " %8.5f " % el[k].rho(st,st).real

            if print_coherences == 1:
                # Print coherences
                for st in xrange(nstates):
                    for st1 in xrange(st):
                        line = line + " %8.5f %8.5f " % (el[k].rho(st,st1).real, el[k].rho(st,st1).imag)
             
            line = line + "\n"

            fel.write(line)
            fel.close()

def init_system(data, g):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] data   The list of variables, containing atomic element names and coordinates
    # \param[in] g      The list of gradients on all atoms
    # This function returns System object which will be used in classical MD.
    #
    # Used in:  main.py/main

    # Create Universe and populate it
    U = Universe();   Load_PT(U, "elements.txt", 0)

    syst = System()

    sz = len(data["coor_atoms"])
    for i in xrange(sz):
        atom_dict = {} 
        atom_dict["Atom_element"] = data["l_atoms"][i]

        # warning: below we take coordinates in Angstroms, no need for conversion here - it will be
        # done inside
        atom_dict["Atom_cm_x"] = data["coor_atoms"][i][0]
        atom_dict["Atom_cm_y"] = data["coor_atoms"][i][1]
        atom_dict["Atom_cm_z"] = data["coor_atoms"][i][2]

        print "CREATE_ATOM ",atom_dict["Atom_element"]
        at = Atom(U, atom_dict)
        at.Atom_RB.rb_force = VECTOR(-g[i][0], -g[i][1], -g[i][2])

        syst.CREATE_ATOM(at)

    syst.show_atoms()
    print "Number of atoms in the system = ", syst.Number_of_atoms

    return syst



