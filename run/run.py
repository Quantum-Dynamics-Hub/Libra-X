# ******************************************************
# To confirm input parameters, see "input.manual".
# ******************************************************

import os
import sys
import math

user = 1 # 0 for Alexey, 1 for Kosuke, others should input the path they use

# input the paths of libra binary files and libra-gamess_interface source files. 

libra_bin_path = "" # set the path name to the source files in libracode
libra_gamess_int_path = "" # set the path name to the source files in libra-gamess_interface

if user==0:
    # For Alexey
    libra_bin_path = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code/_build/src" 
    libra_gamess_int_path = "/user/alexeyak/Programming/libra-gamess_interface/src" 
elif user==1:
    # For Kosuke
    libra_bin_path = "/projects/academic/alexeyak/kosukesa/libracode-code/_build/src"
    libra_gamess_int_path = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/src"

os.environ["src_path"] = libra_gamess_int_path
sys.path.insert(1,os.environ["src_path"]) # Path to the source code

########## Setup all manual parameters here ####################

params = {}

params["gms_inp0"] = "23waters.inp"    # initial input file of GAMESS
params["gms_inp"] = "23waters_wrk.inp" # working input file of GAMESS
params["gms_out"] = "23waters.out"     # output file of GAMESS
params["nproc"] = 2                    # the number of processors
params["basis_option"] = 2             # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["dt_nucl"] = 20.0               # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
params["el_mts"] = 1                   # electronic time steps per one nuclear time step
params["Nsnaps"] = 2                   # the number of MD rounds
params["Nsteps"] = 1                   # the number of MD steps per snap
params["nconfig"] = 2                  # the number of initial nuclei configurations

# Thermostat parameters for NVT MD (if MD_type=1)
params["MD_type"] = 1                       # option 1 -> NVT, otherwise -> NVE ; If this is 1, the parameters below should be selected.
params["nu_therm"] = 0.01                   # shows thermostat frequency
params["NHC_size"] = 3                      # the size of Nose-Hoover chains
params["Temperature"] = 300.0               # Target temperature in thermostat
params["thermostat_type"] = "Nose-Hoover"   # option : "Nose-Hoover" or "Nose-Poincare"

# Excited electronic states
# caution: start from 1, not 0.
Nmin = 92   # lowest molecular orbital taken for creating excited states
HOMO = 92   # Highest Occupied Molecular Orbital : LUMO = HOMO + 1
Nmax = 93   # highest molecular orbital taken for creating excited states
spin = 0    # a flag to consider spin : option 0 -> no, 1 -> yes
flip = 0    # (if spin = 1,) a flag to consider spin-flip : option 0 -> no, 1 -> yes

params["HOMO"] = HOMO

# Surface Hopping
params["SH_type"] = 1 # Surface Hopping type : option  1 -> FSSH, 2 -> GFSH , 3 -> MSSH
params["ntraj"] = 2   # number of excited states trajectories

# If you use boltzman factor, then params["use_boltz_factor"] = 1 and params["do_rescaling"] = 0
#            velocity rescaling, then params["use_boltz_factor"] = 0 and params["do_rescaling"] = 1 (This calculation will be expensive) 

params["use_boltz_factor"] = 0 # A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme: 0 -> no, 1-> yes
params["do_rescaling"] = 1     # The flag to turn on/off CPA: 0 - no velocity rescaling (CPA, no back-reaction)
params["do_reverse"] = 0       # The option that determines what to do if the hop was rejected because of the energy conservation(frustrated hop): 
                               # do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta (velocities) are inverted.

# select directories where the results will be printed out.
params["res"] = ""     # directory where the energies and trajectories files will be printed out
params["mo_ham"] = ""  # directory where MO basis vibronic hamiltonians will be printed out
params["sd_ham"] = ""  # directory where SD basis vibronic hamiltonians will be printed out

if user==0:
    # For Alexey
    params["res"] = "/user/alexeyak/Programming/libra-gamess_interface/run/res/"
    params["mo_ham"] = "/user/alexeyak/Programming/libra-gamess_interface/run/mo_ham/" 
    params["sd_ham"] = "/user/alexeyak/Programming/libra-gamess_interface/run/sd_ham/" 
elif user==1:
    # For Kosuke
    params["res"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/res/" 
    params["mo_ham"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/mo_ham/" 
    params["sd_ham"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/sd_ham/"

# output file
params["traj_file"] = params["res"]+"md" # containing MD trajectories
params["ene_file"] = params["res"]+"ene" # containing kinetic, potential, system, and thermostat-coupled system energies 
params["mu_file"] = params["res"]+"mu"   # containing dipole moment matrices
params["se_pop_prefix"] = "out/"         # where the results of the TD-SE calculations will be printed out 
params["sh_pop_prefix"] = "out/"         # where the results of the SH calculations will be printed out

# flags for debugging
params["print_coherences"] = 0              # compute and print electronic coherences (c^*_i * c_j) : option 0 -> no , 1 -> yes
params["print_sd_ham"] = 0                  # print SD basis vibronic Hamiltonian
params["print_mo_ham"] = 0                  # print full and reduced size MO basis vibronic Hamiltonian
params["print_SH_results_with_scaling"] = 1 # print MD, Energy, and dipole moment results of SH calculation with velocity rescaling  
params["debug_densmat_output"] = 0          # print the debug info into standard output: density matrices, also including for the wavefunctions at different time steps
params["debug_mu_output"] = 0               # print the debug info into standard output: transition dipole moment matrices
params["debug_gms_unpack"] = 0              # print the debug info into standard output: unpacked data from GAMESS
params["debug_ham_ex"] = 0                  # print the debug info into standard output: external hamiltonian matrices for SH calculation
params["debug_SH_cal"] = 0                  # print the debug info into standard output: hopping probabilities matrices and SH_states
params["check_hopping_probs"] = 1           # print the hopping probabilities if they are larger than 1.(To check whether dt_nucl is too large or not.)

# ***************************************************************

from path_libra_lib import * # import path_libra_lib module 
path_libra_lib(libra_bin_path) # Path to the libra libraries

from create_states import *
params["excitations"] = create_states(Nmin,HOMO,Nmax,spin,flip) # generate a list of "excitation" objects.

import main        # import main module of the libra-Gamess-interface code

data, test_data = main.main(params)  # run actual calculations
