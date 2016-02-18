import os
import sys
import math

# Path the the source code
sys.path.insert(1,"/user/alexeyak/Programming/libra-gamess_interface/src")

cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
print "Using the Libra installation at", cwd
sys.path.insert(1,cwd+"/_build/src/mmath")
sys.path.insert(1,cwd+"/_build/src/qchem")
sys.path.insert(1,cwd+"_build/src/dyn")
sys.path.insert(1,cwd+"_build/src/chemobjects")
sys.path.insert(1,cwd+"_build/src/hamiltonian")




########## Setup all manual parameters here ####################

params = {}

params["gms_inp0"] = "H2O.inp"    # initial input file
params["gms_inp"] = "H2O_wrk.inp" # working input file 
params["gms_out"] = "H2O.out"  # output file
params["nproc"] = 1             # the number of processors
params["basis_option"]=2 # ab initio or Semi-Empirical calculation?  Options: \"ab_initio\" = 1 , \"semi_empirical\" = 2
params["dt_nucl"]=20.0  # time step for nuclear dynamics  ex) 20 a.u. = 0.5 fsec
# the number of MD rounds
params["Nsnaps"]=5  # the number of MD rounds
params["Nsteps"]=1  # the number of MD steps per snap
params["res"]="/user/alexeyak/Programming/libra-gamess_interface/run/res/" # the directory where the energies and NACs files will be printed out
params["traj_file"] = params["res"]+"md.xyz"
params["ene_file"] = params["res"]+"ene.xyz"


################################################################

import main        # import main module of the libra-Gamess-interface code

main.main(params)  # run actual calculations

