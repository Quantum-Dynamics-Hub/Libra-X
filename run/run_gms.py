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

user = 0 # 0 for Alexey, 1 for Kosuke, and 2 for Ekadashi; others should input the path they use
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


import defaults

########## Setup all manual parameters here ####################

params = {}

defaults.set_defaults(params, "GAMESS")


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
params["Nsnaps"] = 5                        # the number of total MD snapshots
params["Nsteps"] = 1                        # the number of MD steps per 1 snapshot
params["Ncool"]  = 3                        # in the end of that many initial snapshots 
                                            # we will be cooling the system: resetting momenta to zero
                                            # It is important to use a sufficiently large "Nsteps" variable, to make the
                                            # annealing process more efficient. But, on the other hand, if you are too far from
                                            # equilibrium, make "Nsteps" smaller
params["Nstart"] = 6                        # the printout cycle when we will initiate NA-MD and
                                            # electronic dynamics with surface hoping
params["nconfig"] = 1                       # the number of initial nuclear/velocity geometry
params["MD_type"] = 1                       # option 1 -> NVT, otherwise -> NVE ; If this is 1, the parameters below should be selected.
params["nu_therm"] = 0.001                  # shows thermostat frequency
params["NHC_size"] = 5                      # the size of Nose-Hoover chains
params["Temperature"] = 300.0               # Target temperature in thermostat
params["thermostat_type"] = "Nose-Hoover"   # option : "Nose-Hoover" or "Nose-Poincare"
params["sigma_pos"] = 0.01                  # Magnitude of random atomic displacements 
params["is_MM"] = 0                         # flag for including MM interaction : option 1 -> yes, otherwise -> no.
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
params["num_SH_traj"] = 1              # number of excited states trajectories per initial nuclei geometry and excited states
params["smat_inc"] = 0                 # 1 Including overlap matrix (S), 0 when overlap matrix (S) not included in el propagation



# ***************************************************************

from states import *

# create excitation list
params["excitations"] = [ excitation(0,1,0,1), excitation(0,1,1,1), excitation(-1,1,1,1) ] 
#params["excitations"] = [ excitation(0,1,0,1)]
params["excitations_init"] = [0]

# create thermostat
params["therm"] = Thermostat({"thermostat_type":"Nose-Hoover","nu_therm":0.001,"Temperature":300.0,"NHC_size":5})

# create Universe
params["U"] = Universe(); LoadPT.Load_PT(params["U"], "elements.txt");

# Create force field                                                                                                                                 
params["ff"] = ForceField({"mb_functional":"LJ_Coulomb","R_vdw_on": 10.0,"R_vdw_off":15.0 })
LoadUFF.Load_UFF(params["ff"], "uff.d")

params["ent_file"] = ""           # file including atomic coordinates and connectivity information for MM part

#HOMO = params["HOMO"]
#Nmin = params["HOMO"] + params["min_shift"]
#Nmax = params["HOMO"] + params["max_shift"]
#params["excitations"] = create_states(Nmin,HOMO,Nmax,spin,flip) # generate a list of "excitation" objects.

import main        # import main module of the libra-Gamess-interface code

#data, test_data = main.main(params)  # run actual calculations
main.main(params)  # run actual calculations
