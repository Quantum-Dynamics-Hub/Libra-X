# ******************************************************
# To confirm input parameters, see "run.manual".
# ******************************************************

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

########## Setup all manual parameters here ####################

params = {}

params["gms_inp0"] = "H2O.inp"    # initial input file
params["gms_inp"] = "H2O_wrk.inp" # working input file 
params["gms_out"] = "H2O.out"  # output file
params["nproc"] = 1            # the number of processors
params["basis_option"] = 2 # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["dt_nucl"] = 20.0  # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
params["el_mts"] = 1  # electronic time steps per one nuclear time step
params["Nsnaps"] = 1 # the number of MD rounds
params["Nsteps"] = 1  # the number of MD steps per snap

# For Kosuke
params["res"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/res/" # directory where the energies and trajectories files will be printed out
params["mo_ham"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/mo_ham/" # directory where MO basis vibronic hamiltonians will be printed out  
params["sd_ham"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/sd_ham/" # directory where SD basis vibronic hamiltonians will be printed out

# For Alexey
#params["res"] = "/user/alexeyak/Programming/libra-gamess_interface/run/res/"
#params["mo_ham"] = "/user/alexeyak/Programming/libra-gamess_interface/run/mo_ham/" # directory where MO basis vibronic hamiltonians will be printed out
#params["sd_ham"] = "/user/alexeyak/Programming/libra-gamess_interface/run/sd_ham/" # directory where SD basis vibronic hamiltonians will be printed out


params["traj_file"] = params["res"]+"md.xyz"
params["ene_file"] = params["res"]+"ene.dat"
params["mu_file"] = params["res"]+"mu.dat"
params["se_pop_prefix"] = "out/"  # where the results of the TD-SE calculations will be printed out 

params["print_coherences"] = 1 # a flag to compute and print electronic coherences (c^*_i * c_j) : option 0 -> no , 1 -> yes
params["debug_mu_output"] = 1 # print the debug info into standard output: transition dipole moment matrices
params["debug_densmat_output"] = 0 # print the debug info into standard output: density matrices, also including for the wavefunctions at different time steps
params["print_mo_ham"] = 1 # print full and reduced size MO basis vibronic Hamiltonian
params["debug_gms_unpack"] = 0 # print unpacked data from GAMESS

# ***************************************************************
# Excited electronic states

Nmin = 1   # lowest molecular orbital taken for creating excited states
HOMO = 4   # Highest Occupied Molecular Orbital : LUMO = HOMO + 1
Nmax = 6   # highest molecular orbital taken for creating excited states
spin = 0   # a flag to consider spin : option 0 -> no, 1 -> yes
flip = 0   # (if spin = 1,) a flag to consider spin-flip : option 0 -> no, 1 -> yes 

params["HOMO"] = HOMO

#print params

################################################################

from path_libra_lib import * # import path_libra_lib module 
path_libra_lib(libra_bin_path) # Path to the libra libraries

from create_states import *
# generate a list of "excitation" objects.
params["excitations"] = create_states(Nmin,HOMO,Nmax,spin,flip)

import main        # import main module of the libra-Gamess-interface code

data, test_data = main.main(params)  # run actual calculations
