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

user = 1 # 0 for Olga; others should input the path they use
test = 0 # 0 for 1 water molecule; 1 for 23 water molecules

# input the paths of libra binary files and libra-gamess_interface source files. 

libra_bin_path = "" # set the path name to the source files in libracode
libra_g09_int_path = "" # set the path name to the source files in libra-g09_interface

if user==0:
    # For Olga
    libra_bin_path = "/data/ob070/soft/libra-code/src"
    libra_g09_int_path = "/data/ob070/soft/Libra-X/src"

elif user==1:
    # For Ekadashi
    libra_bin_path = "/projects/academic/alexeyak/ekadashi/libracode-dev/libracode-code/_build/src"
    libra_x_path = "/gpfs/scratch/ekadashi/Libra-X/src"


os.environ["src_path"] = libra_x_path   # Path to the source code
sys.path.insert(1,os.environ["src_path"]) # Path to the source code

########## Setup all manual parameters here ####################

params = {}

# Of course, here we use G09
params["interface"] = "G09"

# G09 variables
# We invoke "run_g09a inp" in x_to_libra_g09.py/exe_g09

params["g09_inp0"] = ""           # initial input file of G09
params["g09_inp"] = ""            # working input file of G09
params["g09_out"] = ""            # output file of G09
params["nproc"] = 1               # the number of processors : default = 1
params["basis_option"] = 2        # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["ent_file"] = ""           # file including atomic coordinates and connectivity information for MM part 

params["mult"] = 1
params["charge"] = 0

if test==0:
    params["g09_inp0"] = "H2O_g09.inp"    # initial input file of G09
    params["g09_inp"] = "H2O_wrk_g09.inp" # working input file of G09
    params["g09_out"] = "H2O_g09.out"     # output file of G09
    params["ent_file"] = "H2O_g09.ent"    # file including atomic coordinates and conncectivity information for MM part

elif test==1:
    params["g09_inp0"] = "23waters.inp"    # initial input file of G09
    params["g09_inp"] = "23waters_wrk.inp" # working input file of G09
    params["g09_out"] = "23waters.out"     # output file of G09
    params["ent_file"] = "23waters.ent"    # file including atomic coordinates and conncectivity information for MM part

# MD variables

params["dt_nucl"] = 20.0                    # time step in a.u. for nuclear dynamics. 20 a.u. is close to 0.5 fsec.
params["Nsnaps"] = 5                        # the number of total MD snapshots
params["Nsteps"] = 1                        # the number of MD steps per 1 snapshot
params["Ncool"]  = 3                        # in the end of that many initial snapshots 
                                            # we will be cooling the system: resetting momenta to zero
                    # It is important to use a sufficiently large "Nsteps" variable, to make the
                    # annealing process more efficient. But, on the other had, if you are too far from
                    # equilibrium, make "Nsteps" smaller

params["Nstart"] = 6       # the printout cycle when we will initiate NA-MD and
                           # electronic dynamics with surface hoping
params["nconfig"] = 1                       # the number of initial nuclear/velocity geometry
params["flag_ao"] = 1                       # flag for atomic orbital basis : option 1 -> yes. otherwise -> no. Don't choose 1 when you use PM6: PM6 calculation doesn't output it at present.
params["MD_type"] = 0                       # option 1 -> NVT, otherwise -> NVE ; If this is 1, the parameters below should be selected.
params["sigma_pos"] = 0.01                  # Magnitude of random atomic displacements 
params["is_MM"] = 1                         # flag for including MM interaction : option 1 -> yes, otherwise -> no.
params["MM_fraction"] = 0.0              # For a QM/MM mixing: E_total = (1-f)*E(QM) + f*E(MM), same for forces!

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
    # For Olga
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
params["excitations"] = [ excitation(0,1,0,1), excitation(0,1,1,1) ] 
#params["excitations"] = [ excitation(0,1,0,1)]
params["excitations_init"] = [0]

# create thermostat
params["therm"] = Thermostat({"thermostat_type":"Nose-Hoover","nu_therm":0.001,"Temperature":300.0,"NHC_size":5})

params["debug_g09_unpack"] = 0              # print the debug info into standard output: unpacked data from GAMESS
params["non-orth"] = 0

# create Universe
params["U"] = Universe(); LoadPT.Load_PT(params["U"], "elements.txt");

# Create force field                                                                                                                                 
params["ff"] = ForceField({"mb_functional":"LJ_Coulomb","R_vdw_on": 10.0,"R_vdw_off":15.0 })
LoadUFF.Load_UFF(params["ff"], "uff.d")

#HOMO = params["HOMO"]
#Nmin = params["HOMO"] + params["min_shift"]
#Nmax = params["HOMO"] + params["max_shift"]
#params["excitations"] = create_states(Nmin,HOMO,Nmax,spin,flip) # generate a list of "excitation" objects.

import main        # import main module of the libra-Gamess-interface code
import defaults
defaults.set_defaults(params, "G09")
#data, test_data = main.main(params)  # run actual calculations
main.main(params)  # run actual calculations
