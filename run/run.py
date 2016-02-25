import os
import sys
import math

# For Kosuke
libra_bin_path = "/projects/academic/alexeyak/kosukesa/libracode-code/_build/src"
libra_gamess_int_path = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/src"


# For Alexey
#libra_bin_path = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code/_build/src" # set the path name to the source files in libracode
#libra_gamess_int_path = "/user/alexeyak/Programming/libra-gamess_interface/src"

os.environ["src_path"] = libra_gamess_int_path
sys.path.insert(1,os.environ["src_path"]) # Path to the source code

from lazy import * # This module is originally PYXAID code from  https://sourceforge.net/p/pyxaid/code/ci/master/tree/

########## Setup all manual parameters here ####################

params = {}

params["gms_inp0"] = "H2O.inp"    # initial input file
params["gms_inp"] = "H2O_wrk.inp" # working input file 
params["gms_out"] = "H2O.out"  # output file
params["nproc"] = 1            # the number of processors
params["basis_option"] = 2 # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["dt_nucl"] = 20.0  # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
params["dt_ele"] = 1.0  # time step for electronic integration
params["Nsnaps"] = 1  # the number of MD rounds
params["Nsteps"] = 1  # the number of MD steps per snap
params["namdtime"] = 20  # Trajectory time, a.u. (<= dt_nucl * Nsnaps * Nsteps)

# For Kosuke
params["res"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/res/" # the directory where the energies and NACs files will be printed out

# For Alexey
#params["res"] = "/user/alexeyak/Programming/libra-gamess_interface/run/res/"

params["traj_file"] = params["res"]+"md.xyz"
params["ene_file"] = params["res"]+"ene.xyz"

# ***************************************************************
# Excited electronic states
# Define states:
# Example of indexing convention with Nmin = 5, HOMO = 5, Nmax = 8
# the orbitals indices are consistent with GAMESS indexing, which starts from 1
# [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] - all computed orbitals
# [ 1, 2, 3, 4, 5]                     - occupied orbitals
#                 [6, 7, 8, 9, 10, 11] - unoccupied orbitals
#              [5, 6, 7, 8]            - active space
 
 
# Excitations from HOMO to [LUMO,... Nmax]
Nmin = 1
HOMO = 4
LUMO = HOMO+1
Nmax = 6
 
# Initial condition excitation: I->J
I = 3
J = 5

# This is how I,J are mapped to index of the corresponding basis state
# see module lazy.py for more details. This is not the most general way
# of such mapping though.
ex_indx = 1 + (J-LUMO)*(HOMO+1 - Nmin) + (I-Nmin)
params["ex_indx"] = ex_indx
 
# Each entry of the list below is an initial condition. It is also a list
# but containing only 2 elements - first is the time step at which we start
# the dynamics, the second is the index of the excited state configuration
iconds = [0]

params["iconds"] = iconds
for ic in (0,len(iconds)):
    params["prop_ele_file"+str(ic)] = params["res"]+"prop_ele"+str(ic)+".xyz" 
 
# Set active space and the basis states
params["active_space"] = range(Nmin,Nmax+1)
 
# Generate basis states
GS = ground_state(Nmin,HOMO)  # ground state
SE = single_excitations(Nmin,Nmax,HOMO,1)  # single excitations
 
# Now combine the ground and singly excited states in one list
# In our convention, the GS configuration must be the first state in the
# list of the basis states.
params["states"] = []
params["states"].append(GS)
for se in SE:
    params["states"].append(se)

params["ex_num"] = len(params["states"])

print params

################################################################

from path_libra_lib import * # import module pathing the libra libralies 
path_libra_lib(libra_bin_path) # Path to the libra libraries

import main        # import main module of the libra-Gamess-interface code

main.main(params)  # run actual calculations

