# ******************************************************
# To confirm input parameters, see "input.manual".
# ******************************************************

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

user = 1 # 0 for Alexey, 1 for Kosuke, and 2 for Ekadashi; others should input the path they use
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
    libra_gamess_int_path = "/home/e1667/dev/libra-gamess_interface/src"

elif user==2:
    # For Ekadashi
    libra_bin_path = "/projects/academic/alexeyak/ekadashi/libra-dev/libracode-code/_build/src"
    libra_gamess_int_path = "/projects/academic/alexeyak/ekadashi/devel/libra-gamess_interface/src"

os.environ["src_path"] = libra_gamess_int_path
sys.path.insert(1,os.environ["src_path"]) # Path to the source code

########## Setup all manual parameters here ####################

params = {}

# Of course, here we use GAMESS
params["interface"] = "GAMESS"

# GAMESS variables
# We invoke "/usr/bin/time rungms gms_inp VERNO nproc > gms_out" in x_to_libra_gms.py/exe_gamess
# Paths of SCR, USERSCR, GMSPATH in the rungms script will be defined by environmental variables later.
# (Supposed TARGET is already defined as sockets or mpi.)

params["GMSPATH"] = ""            # the directory including GAMESS binary files.
params["rungms"] = ""             # "rungms" name. On CCR @ UB, use "rungms.slurm".
params["gms_inp0"] = ""           # initial input file of GAMESS
params["gms_inp"] = ""            # working input file of GAMESS
params["gms_out"] = ""            # output file of GAMESS
params["nproc"] = 1               # the number of processors : default = 1
params["VERNO"] = ""              # Version No., e.g. 00, 01, etc....
params["scr_dir"] = ""            # scratch directory including GAMESS temporary files.This directory will be created and deleted every GAMESS calculation.
params["basis_option"] = 2        # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["ent_file"] = ""           # file including atomic coordinates and connectivity information for MM part 

if user==0 or user==2:
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
    params["scr_dir"] = "/home/e1667/work/scr"

if test==0:
    params["gms_inp0"] = "H2O.inp"    # initial input file of GAMESS
    params["gms_inp"] = "H2O_wrk.inp" # working input file of GAMESS
    params["gms_out"] = "H2O.out"     # output file of GAMESS
    params["ent_file"] = "H2O.ent"    # file including atomic coordinates and conncectivity information for MM part

elif test==1:
    params["gms_inp0"] = "23waters.inp"    # initial input file of GAMESS
    params["gms_inp"] = "23waters_wrk.inp" # working input file of GAMESS
    params["gms_out"] = "23waters.out"     # output file of GAMESS
    params["ent_file"] = "23waters.ent"    # file including atomic coordinates and conncectivity information for MM part

# MD variables

params["dt_nucl"] = 20.0                    # time step in a.u. for nuclear dynamics. 20 a.u. is close to 0.5 fsec.
params["Nsnaps"] = 5                        # the number of MD rounds
params["Nsteps"] = 1                        # the number of MD steps per snap
params["Ncool"]  = 3                        # the number of cooling rounds from t=0
params["Nstart"] = 6                        # the number of rounds for starting NA-MD
params["nconfig"] = 1                       # the number of initial nuclear/velocity geometry
params["flag_ao"] = 1                       # flag for atomic orbital basis : option 1 -> yes. otherwise -> no. Don't choose 1 when you use PM6: PM6 calculation doesn't output it at present.
params["MD_type"] = 1                       # option 1 -> NVT, otherwise -> NVE ; If this is 1, the parameters below should be selected.
params["nu_therm"] = 0.001                  # shows thermostat frequency
params["NHC_size"] = 5                      # the size of Nose-Hoover chains
params["Temperature"] = 300.0               # Target temperature in thermostat
params["thermostat_type"] = "Nose-Hoover"   # option : "Nose-Hoover" or "Nose-Poincare"
params["sigma_pos"] = 0.01                  # Magnitude of random atomic displacements 
params["f_vdw"] = 1                         # flag for including vdw(non-bonded) interaction : option 1 -> yes, otherwise -> no.

spin = 0    # a flag to consider spin : option 0 -> no, 1 -> yes
flip = 0    # (if spin = 1,) a flag to consider spin-flip : option 0 -> no, 1 -> yes

# Excited electron dynamics

if test==0:
    params["HOMO"] = 3 # In the case of H2O, which has 8 electrons, occupied orbitals are 0, 1, 2, and 3. 
elif test==1:
    params["HOMO"] = 91 # In the case of 23 waters, which has 184 electrons, occupied orbitals are 0, 1, ..., 90, and 91.

params["min_shift"] = -1               # e.g. -1 -> HOMO-1, HOMO
params["max_shift"] = 1                # e.g.  1 -> LUMO
params["el_mts"] = 1                   # electronic time steps per one nuclear time step
params["tsh_method"] = 1               # Surface Hopping type : option  1 -> FSSH, 2 -> GFSH , 3 -> MSSH
params["rep"] = 1                      # representation: 0 - diabatic, 1 - adiabatic
params["num_SH_traj"] = 1              # number of excited states trajectories per initial nuclei geometry and excited states
params["use_boltz_factor"] = 0         # A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme: 0 -> no, 1-> yes
params["do_rescaling"] = 0             # The flag to control velocity rescaling: 0 - no velocity rescaling, 1 - do rescaling
params["do_reverse"] = 0               # The option that determines what to do if the hop was rejected because of the energy conservation(frustrated hop): 
                                       # do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta (velocities) are inverted.
params["smat_inc"] = 0                 # 1 Including overlap matrix (S), 0 when overlap matrix (S) not included in el propagation

# select directories where the results will be printed out.
params["res"] = ""     # directory where the all results will be printed out
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

elif user==2:
    # For Ekadashi
    cwd = os.getcwd()
    params["res"] =  cwd + "/res/" #; print "res is located on ",params["res"] ; 
    params["mo_ham"] =  cwd + "/mo_ham/" #; print "mo_ham is located on ",params["mo_ham"] ;
    params["sd_ham"] = cwd + "/sd_ham/" #; print "sd_ham is located on ",params["sd_ham"] ;


# flags for debugging
params["print_aux_results"] = 1             # print auxiliary results ; a large amount of files(MD, Energy trajectories, etc..) will be printed out.
params["print_coherences"] = 1              # compute and print electronic coherences (c^*_i * c_j) : option 0 -> no , 1 -> yes
params["print_sd_ham"] = 1                  # print SD basis vibronic Hamiltonian
params["print_mo_ham"] = 1                  # print full and reduced size MO basis vibronic Hamiltonian
params["print_SH_results_with_scaling"] = 0 # print MD, Energy, and dipole moment results of SH calculation with velocity rescaling  
params["debug_densmat_output"] = 1          # print the debug info into standard output: density matrices, also including for the wavefunctions at different time steps
params["debug_mu_output"] = 0               # print the debug info into standard output: transition dipole moment matrices
params["debug_gms_unpack"] = 0              # print the debug info into standard output: unpacked data from GAMESS
#params["debug_ham_ex"] = 1                  # print the debug info into standard output: external hamiltonian matrices for SH calculation
params["print_tsh_probabilities"] = 0      # print the debug info into standard output: hopping probabilities matrices and SH_states
params["check_tsh_probabilities"] = 0      # print the hopping probabilities if they are larger than 1.(To check whether dt_nucl is too large or not.)

# ***************************************************************

from states import *

# create excitation list
params["excitations"] = [ excitation(0,1,0,1), excitation(0,1,1,1), excitation(-1,1,1,1) ] 
#params["excitations"] = [ excitation(0,1,0,1)]
params["excitations_init"] = [0]

# create thermostat
#params["therm"] = Thermostat({"thermostat_type":"Nose-Hoover","nu_therm":0.001,"Temperature":300.0,"NHC_size":5})

# create Universe
params["U"] = Universe(); LoadPT.Load_PT(params["U"], "elements.txt");

# Create force field                                                                                                                                 
params["uff"] = ForceField({"mb_functional":"LJ_Coulomb","R_vdw_on": 10.0,"R_vdw_off":15.0 })
LoadUFF.Load_UFF(params["uff"], "uff.d")

#HOMO = params["HOMO"]
#Nmin = params["HOMO"] + params["min_shift"]
#Nmax = params["HOMO"] + params["max_shift"]
#params["excitations"] = create_states(Nmin,HOMO,Nmax,spin,flip) # generate a list of "excitation" objects.

import main        # import main module of the libra-Gamess-interface code

#data, test_data = main.main(params)  # run actual calculations
main.main(params)  # run actual calculations
