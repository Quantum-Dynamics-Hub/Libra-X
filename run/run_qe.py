import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *



user = 1  # 0 - Alexey, 1 - Ekadashi

################ System-specific settings ########################
if user==0:
    # For Alexey
    libra_bin_path = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code/_build/src" # set the path name to the source files in libracode
    libra_qe_int_path = "/user/alexeyak/Programming/libra-gamess_interface/src"
    #res_dir =  "/user/alexeyak/Programming/libra-qe_interface/run/res/"

elif user==1:
    # For Ekadashi
    libra_bin_path = "/projects/academic/alexeyak/ekadashi/libracode-dev/libracode-code/_build/src"
    libra_qe_int_path = "/projects/academic/alexeyak/ekadashi/devel/libra-gamess_interface/src"
    #res_dir = "/projects/academic/alexeyak/ekadashi/devel/libra-qe_interface/run/res/"


os.environ["src_path"] = libra_qe_int_path   # Path to the source code
sys.path.insert(1,os.environ["src_path"])    # Path to the source code

#from path_libra_lib import * # import path_libra_lib module
#path_libra_lib(libra_bin_path)               # Path to the libra libraries

import main # import main module of the libra-QE-interface code

########## Setup all manual parameters here ####################

params = {}
params["qe_debug_print"] = 0
params["nproc"] = 1              # the number of processors
params["dt_nucl"]=20.0  # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
params["Nsnaps"]=200      # the number of MD rounds
params["Nsteps"]=2      # the number of MD steps per snap
#params["res"]=res_dir   # the directory where the energies and NACs files will be printed out
#params["traj_file"] = params["res"]+"md.xyz"
#params["ene_file"] = params["res"]+"ene.dat"
#params["S_mat"] = params["res"]+"s_mat"
params["nspin"] = 1

params["nconfig"] = 1
params["el_mts"] = 1
params["tsh_method"] = 1     # Surface Hopping type : option  1 -> FSSH, 2 -> GFSH , 3 -> MSSH
params["num_SH_traj"] = 2
params["tsh_method"] = 1               # Surface Hopping type : option  1 -> FSSH, 2 -> GFSH , 3 -> MSSH
params["rep"] = 0                      # representation: 0 - diabatic, 1 - adiabatic
params["use_boltz_factor"] = 0         # A flag to select the Boltzmann scaling in lieu of hop rejection/velocity rescaling scheme: 0 -> no, 1-> yes
params["do_rescaling"] = 1             # The flag to control velocity rescaling: 0 - no velocity rescaling, 1 - do rescaling
params["do_reverse"] = 1               # The option that determines what to do if the hop was rejected because of the energy conservation(frustrated hop): 
                                       # do_reverse = 0 - nuclear momenta(velocities) stay unchanged; do_reverse = 1 - nuclear momenta (velocities) are inverted.


params["interface"] = "QE"   # "QE" for libra_qe_interface, "GAMESS" if libra_gamess_interface is used
params["MD_type"] = 0  # 1 NVT ensamble, 0 NVE ensamble
# Thermostat parameters
params["Temperature"] = 300.0
params["nu_therm"] = 0.01
params["NHC_size"] = 3
params["thermostat_type"] = "Nose-Hoover"
params["sigma_pos"] = 0.01  #Displace atomic position randomly

########### Now start actual calculations ###########################
#sys.path.insert(1,os.environ["libra_hamiltonian_path"] + "/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
#from libcontrol_parameters import *

#params["num_MO"] = 3  # number of MO basis used in constructing electronic wavefunction
params["excitations"] = [ excitation(0,1,0,1), excitation(0,1,1,1) ] 
params["HOMO"] = 0
params["min_shift"] = 0
params["max_shift"] = 1 
for i in range(0,len(params["excitations"])):
    params["qe_inp0%i" %i] = "x%i.scf.in" %i    # initial input file
    params["qe_inp%i" %i] = "x%i.scf_wrk.in" %i # working input file 
    params["qe_out%i" %i] = "x%i.scf.out" %i    # output file


# Flags and Debugging
params["flag_ao"] = 0   # flag for atomic orbital basis : option 1 -> yes. otherwise -> no.
params["print_coherences"] = 1
params["print_aux_results"] = 1
params["print_mo_ham"] = 0
params["print_sd_ham"] = 1
params["print_tsh_probabilities"] = 1
params["check_tsh_probabilities"] = 1 

# select directories where the results will be printed out.
params["res"] = ""     # directory where the all results will be printed out
params["mo_ham"] = ""  # directory where MO basis vibronic hamiltonians will be printed out
params["sd_ham"] = ""  # directory where SD basis vibronic hamiltonians will be printed out

if user==0:
    # For Alexey
    params["res"] = "/user/alexeyak/Programming/libra-gamess_interface/run/res/"
    params["mo_ham"] = "/user/alexeyak/Programming/libra-gamess_interface/run/mo_ham/" 
    params["sd_ham"] = "/user/alexeyak/Programming/libra-gamess_interface/run/sd_ham/" 

if user==1:
    # For Ekadashi
    cwd = os.getcwd()
    params["res"] = cwd+"/res/"
    params["mo_ham"] = cwd+"/mo_ham/"
    params["sd_ham"] = cwd+"/sd_ham/"


# output file
params["traj_file_prefix"] = params["res"]+"md"       # containing MD trajectories
params["ene_file_prefix"] = params["res"]+"ene"       # containing kinetic, potential, system, and thermostat-coupled system energies 
params["mu_file_prefix"] = params["res"]+"mu"         # containing dipole moment matrices
params["se_pop_file_prefix"] = params["res"]+"se_pop"           # containing the SE population (if velocity rescaling is applied, this is averaged over TSH trajectories). File name is defined as se_pop_"initial geometry"_"initial excitation"
params["sh_pop_file_prefix"] = params["res"]+"sh_pop"           # containing the SH population averaged over TSH trajectories. File name is defined in the SE way. 
params["se_pop_ex_file_prefix"] = params["res"]+"se_pop_ex"   # containing the SE population averaged over initial geometries. File name is se_pop_ex"initial excitation" 
params["sh_pop_ex_file_prefix"] = params["res"]+"sh_pop_ex"   # containing the SH population averaged over initial geometries. File name is defined in the SE way.

main.main(params)  # run actual calculations


