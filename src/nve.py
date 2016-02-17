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
# Used in : main.py/


from create_gamess_input import *
from create_pdb_file import *
from exe_gamess import *
from gamess_to_libra import *

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
cwd= "/projects/academic/alexeyak/kosukesa/libracode-code/"
print "Current working directory", cwd
sys.path.insert(1,cwd+"_build/src/mmath")
sys.path.insert(1,cwd+"_build/src/chemobjects")
sys.path.insert(1,cwd+"_build/src/hamiltonian")
sys.path.insert(1,cwd+"_build/src/dyn")
#sys.path.insert(1,cwd+"/../../_build/src/hamiltonian/hamiltonian_atomistic")

print "\nTest 1: Importing the library and its content"
from libmmath import *
from libchemobjects import *
from libhamiltonian import *
from libdyn import *
#from cyghamiltonian_atomistic import *

from LoadPT import * # Load_PT
from LoadMolecule import * # Load_Molecule

##############################################################

def construct_systems(Grad):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] Grad : the forces extracted from gamess output file.
    # This function returns SYSTEM and NUCLEAR constructor 
    # which will be used in classical MD.
    #
    # Used in:  main.py/main/nve/nve_settings

    # Create Universe and populate it
    U = Universe()
    verbose = 0
    Load_PT(U, "elements.dat", verbose)


    verb = 0
    assign_rings = 0

    syst = System()
    Load_Molecule(U, syst, os.getcwd()+"/temp.pdb", "pdb_1")


    syst.show_atoms()
    #syst.show_fragments()
    #syst.show_molecules()
    print "Number of atoms in the system = ", syst.Number_of_atoms
    #atlst1 = range(1,syst.Number_of_atoms+1)

    #--------------------- Molecular dynamics ----------------------
    # Nuclear DOFs
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_mass(mol.mass)

    Natoms = len(Grad)

    for i in range(0,Natoms):
        mol.p[3*i], mol.p[3*i+1], mol.p[3*i+2] = 0.0, 0.0, 0.0
        mol.f[3*i], mol.f[3*i+1], mol.f[3*i+2] = -Grad[i][0], -Grad[i][1], -Grad[i][2]

    for i in range(0,Natoms):
        print "mass m=",mol.mass[3*i], mol.mass[3*i+1], mol.mass[3*i+2]
        print "coordinates q = ", mol.q[3*i], mol.q[3*i+1], mol.q[3*i+2]
        print "momenta p= ", mol.p[3*i], mol.p[3*i+1], mol.p[3*i+2]
        print "forces f= ",  mol.f[3*i], mol.f[3*i+1], mol.f[3*i+2]
        print "********************************************************"

    return syst, mol

def nve_MD(syst,mol,ao1,E1,C1,data,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] syst   : SYSTEM constructor including atomic system information.
    # \param[in] mol    : Nuclear constructor including atomic positions, momenta, forces etc...
    # \param[in] ao1    : atomic orbital basis at "t" old
    # \param[in] E1     : molecular energies at "t" old
    # \param[in] C1     : molecular coefficients at "t" old
    # \param[in] data   : data extracted from GAMESS output file, in the dictionary form.
    # \param[in] params : input data from a submit file, in the dictionary form.
    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively.
    #
    # Used in:  main.py/main/nve

    f = open("../out/_en_traj.txt","w")
    dt = params["dt_nuc"]
    job_tmp = params["GMS_JOB"] + "tmp"
    Nrounds = params["Nrounds"]
    Ncycles = params["Ncycles"]
    Natoms = syst.Number_of_atoms
    trajfile = "../out/_mol_traj.xyz"
    os.system("rm -r %s"%(trajfile))

    for i in range(0,Nrounds):

        syst.set_atomic_q(mol.q)
        syst.print_xyz(trajfile,i)

        for j in xrange(0,Ncycles):

            ij = i * Ncycles + j
            mol.propagate_p(0.5*dt)
            mol.propagate_q(dt)  #---> generate GAMESS input

            for k in range(0,Natoms):
                print k,"th atom"
                print "coordinates q = ", mol.q[3*k], mol.q[3*k+1], mol.q[3*k+2]
                print "momenta p= ", mol.p[3*k], mol.p[3*k+1], mol.p[3*k+2]
                print "forces f= ",  mol.f[3*k], mol.f[3*k+1], mol.f[3*k+2]
                print "********************************************************"

            # GAMESS iteration part
            write_GMS_INP(data["l_atoms"],data["l_charges"],params["l_gam_for"],mol,job_tmp ,params["GMS_DIR"])
            exe_gamess(params,job_tmp)                                              # call GAMESS with input file
            ao1, E1, C1, Grad, data = gamess_to_libra(params, ao1, E1, C1,job_tmp,ij)  # extract info from GAMESS
            epot = data["tot_ene"]                                                  # total energy from GAMESS

            for k in range(0,len(Grad)):
                mol.f[3*k], mol.f[3*k+1], mol.f[3*k+2] = -Grad[k][0], -Grad[k][1], -Grad[k][2]  # forces from GAMESS

            mol.propagate_p(0.5*dt)

            ekin = compute_kinetic_energy(mol)

        f.write("i= %3i ekin= %8.5f  epot= %8.5f  etot= %8.5f\n" % (i, ekin, epot, ekin+epot))
        
    f.close()

def nve_setting(l_atoms,coor_atoms,Grad):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_atoms     : a list of atoms.
    # \param[in] coor_atoms  : a list of atomic positions.
    # \param[in] Grad        : a list of foeces.
    # This function returns SYSTEM and NUCLEAR constructors which will be used in
    # classical MD.
    #
    # Used in:  main.py/main/nve

    create_pdb_file(l_atoms,coor_atoms)

    syst, mol = construct_systems(Grad)

    return syst, mol

def nve(data,params,ao1,E1,C1,Grad):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] data   : data extracted from GAMESS output file, in the dictionary form.
    # \param[in] params : input data from a submit file, in the dictionary form.
    # \param[in] ao1    : atomic orbital basis at "t" old
    # \param[in] E1     : molecular energies at "t" old
    # \param[in] C1     : molecular coefficients at "t" old
    # \param[in] Grad   : a list of foeces.
    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively ; it outputs the files for simulation of excited electron dynamics.
    #
    # Used in:  main.py/main/


    syst, mol = nve_setting(data["l_atoms"],data["coor_atoms"],Grad)

    nve_MD(syst,mol,ao1,E1,C1,data,params)
