# ******************************************************
# To confirm input parameters, see "input.manual".
# ******************************************************

import os
import sys
import math

user = 1 # 0 for Alexey, 1 for Kosuke, others should input the path they use
test = 0 # 0 for 1 water molecule; 1 for 23 water molecules

# input the paths of libra binary files and libra-gamess_interface source files. 

libra_bin_path = "" # set the path name to the source files in libracode
libra_gamess_int_path = "" # set the path name to the source files in libra-gamess_interface

if user==0:
    # For Alexey
    libra_bin_path = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code/_build/src" 
    libra_gamess_int_path = "/user/alexeyak/Programming/libra-gamess_interface/src" 
elif user==1:
    # For Kosuke
    libra_bin_path = "/home/e1667/install/libra-code/_build/src"
    libra_gamess_int_path = "/home/e1667/install/libra-gamess_interface/src"

os.environ["src_path"] = libra_gamess_int_path
sys.path.insert(1,os.environ["src_path"]) # Path to the source code

########## Setup all manual parameters here ####################

params = {}

## set variables in GAMESS
# Supposed we invoke "/usr/bin/time rungms gms_inp VERNO nproc > gms_out" in md.py/exe_gamess
# Variables such as SCR, USERSCR, GMSPATH, JOB, VERNO, and NCPUS in the rungms script will be defined later by the parameters below.
# (Supposed TARGET is already defined as sockets or mpi.)

params["GMSPATH"] = ""            # the directory including GAMESS binary files.
params["rungms"] = ""             # "rungms" name. On CCR @ UB, use "rungms.slurm".
params["gms_inp0"] = ""           # initial input file of GAMESS
params["gms_inp"] = ""            # working input file of GAMESS
params["gms_out"] = ""            # output file of GAMESS
params["nproc"] = 1               # the number of processors : default = 1
params["VERNO"] = ""              # Version No., e.g. 00, 01, etc....
params["scr_dir"] = ""            # scratch directory including GAMESS output files.

if user==0:
    # For Alexey (setting for CCR @ UB)
    params["GMSPATH"] = "" # GAMESS path is already taken.
    params["rungms"] = "rungms.slurm"
    params["VERNO"] = "01"
    params["scr_dir"] = os.environ['SLURMTMPDIR'] # slurm type

elif user==1:
    # For Kosuke
    params["GMSPATH"] = "/home/e1667/install/gamess"
    params["rungms"] =  params["GMSPATH"] + "/rungms"
    params["VERNO"] = "00"
    params["scr_dir"] = "/home/e1667/work_NAMD/gamess_scratch"

if test==0:
    params["gms_inp0"] = "H2O.inp"    # initial input file of GAMESS
    params["gms_inp"] = "H2O_wrk.inp" # working input file of GAMESS
    params["gms_out"] = "H2O.out"     # output file of GAMESS

elif test==1:
    params["gms_inp0"] = "23waters.inp"    # initial input file of GAMESS
    params["gms_inp"] = "23waters_wrk.inp" # working input file of GAMESS
    params["gms_out"] = "23waters.out"     # output file of GAMESS

params["basis_option"] = 2             # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["dt_nucl"] = 20.0               # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
params["el_mts"] = 1                   # electronic time steps per one nuclear time step
params["Nsnaps"] = 1                   # the number of MD rounds
params["Nsteps"] = 1                   # the number of MD steps per snap
params["nconfig"] = 1                  # the number of initial nuclear/velocity configurations

# Thermostat parameters for NVT MD (if MD_type=1)
params["MD_type"] = 0                       # option 1 -> NVT, otherwise -> NVE ; If this is 1, the parameters below should be selected.
params["nu_therm"] = 0.01                   # shows thermostat frequency
params["NHC_size"] = 3                      # the size of Nose-Hoover chains
params["Temperature"] = 300.0               # Target temperature in thermostat
params["sigma_pos"] = 0.01                  # Magnitude of random atomic displacements
params["thermostat_type"] = "Nose-Hoover"   # option : "Nose-Hoover" or "Nose-Poincare"

spin = 0    # a flag to consider spin : option 0 -> no, 1 -> yes
flip = 0    # (if spin = 1,) a flag to consider spin-flip : option 0 -> no, 1 -> yes

# Excited electronic states
# caution: start from 0
if test==0:
    params["HOMO"] = 3 # not 4 
elif test==1:
    params["HOMO"] = 91 # not 92

params["min_shift"] = -1  # HOMO-1, HOMO
params["max_shift"] = 1  # LUMO

# Surface Hopping
params["SH_type"] = 1 # Surface Hopping type : option  1 -> FSSH, 2 -> GFSH , 3 -> MSSH
params["ntraj"] = 1   # number of excited states trajectories

# If you use boltzman factor, then params["use_boltz_factor"] = 1 and params["do_rescaling"] = 0
#            velocity rescaling, then params["use_boltz_factor"] = 0 and params["do_rescaling"] = 1 (This calculation will be expensive) 

params["use_boltz_factor"] = 0 # A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme: 0 -> no, 1-> yes
params["do_rescaling"] = 1     # The flag to control velocity rescaling: 0 - no velocity rescaling, 1 - do rescaling
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
    cwd = os.getcwd()
    params["res"] =  cwd + "/res/" #; print "res is located on ",params["res"] ; 
    params["mo_ham"] =  cwd + "/mo_ham/" #; print "mo_ham is located on ",params["mo_ham"] ;
    params["sd_ham"] = cwd + "/sd_ham/" #; print "sd_ham is located on ",params["sd_ham"] ;

# output file
params["traj_file_prefix"] = params["res"]+"md" # containing MD trajectories
params["ene_file_prefix"] = params["res"]+"ene" # containing kinetic, potential, system, and thermostat-coupled system energies 
params["mu_file_prefix"] = params["res"]+"mu"   # containing dipole moment matrices
params["se_pop_file_prefix"] = "out/se_pop"         # where the results of the TD-SE calculations will be printed out
params["sh_pop_file_prefix"] = "out/sh_pop"         # where the results of the SH calculations will be printed out

# flags for debugging
params["print_coherences"] = 0              # compute and print electronic coherences (c^*_i * c_j) : option 0 -> no , 1 -> yes
params["print_sd_ham"] = 0                  # print SD basis vibronic Hamiltonian
params["print_mo_ham"] = 0                  # print full and reduced size MO basis vibronic Hamiltonian
params["print_SH_results_with_scaling"] = 1 # print MD, Energy, and dipole moment results of SH calculation with velocity rescaling  
params["debug_densmat_output"] = 0          # print the debug info into standard output: density matrices, also including for the wavefunctions at different time steps
params["debug_mu_output"] = 0               # print the debug info into standard output: transition dipole moment matrices
params["debug_gms_unpack"] = 1              # print the debug info into standard output: unpacked data from GAMESS
params["debug_ham_ex"] = 0                  # print the debug info into standard output: external hamiltonian matrices for SH calculation
params["debug_SH_cal"] = 0                  # print the debug info into standard output: hopping probabilities matrices and SH_states
params["check_hopping_probs"] = 1           # print the hopping probabilities if they are larger than 1.(To check whether dt_nucl is too large or not.)

# ***************************************************************

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from path_libra_lib import * # import path_libra_lib module 
path_libra_lib(libra_bin_path) # Path to the libra libraries

from create_states import *
#params["excitations"] = [ excitation(0,1,0,1), excitation(0,1,1,1), excitation(-1,1,1,1) ] 

HOMO = params["HOMO"]
Nmin = params["HOMO"] + params["min_shift"]
Nmax = params["HOMO"] + params["max_shift"]
params["excitations"] = create_states(Nmin,HOMO,Nmax,spin,flip) # generate a list of "excitation" objects.

import main        # import main module of the libra-Gamess-interface code

#data, test_data = main.main(params)  # run actual calculations
main.main(params)  # run actual calculations
