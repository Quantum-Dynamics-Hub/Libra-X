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
params["Nsnaps"] = 1  # the number of MD rounds
params["Nsteps"] = 1  # the number of MD steps per snap

# For Kosuke
params["res"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/res/" # directory where the energies and trajectories files will be printed out
params["mo_ham"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/mo_ham/" # directory where MO basis vibronic hamiltonians will be printed out  
params["sd_ham"] = "/projects/academic/alexeyak/kosukesa/dev/libra-gamess_interface/run/sd_ham/" # directory where SD basis vibronic hamiltonians will be printed out


# For Alexey
#params["res"] = "/user/alexeyak/Programming/libra-gamess_interface/run/res/"

params["traj_file"] = params["res"]+"md.xyz"
params["ene_file"] = params["res"]+"ene.dat"
params["se_pop_prefix"] = "out/"  # where the results of the TD-SE calculations will be printed out 

params["print_coherences"] = 1 # a flag to compute and print electronic coherences (c^*_i * c_j) : option 0 -> no , 1 -> yes

# ***************************************************************
# Excited electronic states

Nmin = 2   # lowest molecular orbital taken for creating excited states
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

# generate a list of "excitation" objects manually.

excitations = []  #create_states(Nmin,HOMO,Nmax,spin,flip)
excitations.append(excitation(0,1,0,1)) # GS
excitations.append(excitation(0,1,1,1)) # SE0
excitations.append(excitation(0,1,2,1)) # SE1
excitations.append(excitation(-1,1,1,1)) # SE2
excitations.append(excitation(-1,1,2,1)) # SE3
excitations.append(excitation(-2,1,1,1)) # SE4
excitations.append(excitation(-2,1,2,1)) # SE5
params["excitations"] = excitations

#****************************************************************
#     function

sys.path.insert(1,os.environ["libra_hamiltonian_path"] + "/Hamiltonian_Atomistic/Hamiltonian_QM/Control_Parameters")
from libcontrol_parameters import *

def test1(test_data,Ngbf,Nmin,Nmax,_ex,HOMO):
    # This function checks if D_mol matrix is reduced to D_mol_red correctly
    # and if D_mol_red in MO basis corresponds to D_SD correctly.
    # It returns res=1 if the programs work correctly or res=0 otherwise.

    D_mol = test_data["D_mol"]
    D_mol_red = test_data["D_mol_red"]
    D_SD = test_data["D_SD"]

    res = 1
    # check if D_mol is reduced to  D_mol_red correctly.
    for i in range(Nmin,Nmax+1):
        for j in range(Nmin,Nmax+1):

            if not D_mol.get(i-1,j-1) == D_mol_red.get(i-Nmin,j-Nmin):
                print "(%i,%i) element in D_mol should be equal to (%i,%i) element in D_mol_red" %(i-1,j-1,i-Nmin,j-Nmin)
                print "(%i,%i) element in D_mol is %8.5f and (%i,%i) element in D_mol is %8.5f" %(i-1,j-1, D_mol.get(i,j),i-Nmin,i-Nmin,D_mol_red.get(i-Nmin,j-Nmin))
                res = 0; break;
        if res==0:
            break

    # check if D_mol_red in MO basis corresponds to D_SD in SD basiscorrectly.
    if res==1:
        i = 0
        for ex_i in _ex:
            ih = ex_i.from_orbit[0]
            ie = ex_i.to_orbit[0]
            j = 0
            for ex_j in _ex:
                jh = ex_j.from_orbit[0]
                je = ex_j.to_orbit[0]

                if ie == 0 and ih == 0 and not (je == 0 and jh == 0) : # GS -> EX
                    atmp = jh+HOMO-Nmin
                    btmp = je+HOMO-Nmin
                    if not -D_mol_red.get(atmp,btmp) == D_SD.get(i,j):
                        res = 0; break;

                elif not(ie == 0 and ih == 0) and je == 0 and jh == 0 : # EX -> GS
                    atmp = ie+HOMO-Nmin
                    btmp = ih+HOMO-Nmin                    
                    if not -D_mol_red.get(atmp,btmp) == D_SD.get(i,j):
                        res = 0; break;
                else: # EX -> EX
                    if ih == jh and not (ie == je) : # electron difference
                        atmp = ie+HOMO-Nmin
                        btmp = je+HOMO-Nmin
                        if not -D_mol_red.get(atmp,btmp) == D_SD.get(i,j):
                            res = 0; break;
                    if ie == je and not (ih == jh) : # hole difference
                        atmp = jh+HOMO-Nmin
                        btmp = ih+HOMO-Nmin
                        if not -D_mol_red.get(atmp,btmp) == D_SD.get(i,j):
                            res = 0; break;
                j += 1
            if res==0:
                print "(%i,%i) element in D_mol_red should be equal to (%i,%i) element in D_SD" %(atmp,btmp,i,j)
                print "(%i,%i) element in D_mol_red is %8.5f and (%i,%i) element in D_SD is %8.5f" %(atmp,btmp, D_mol_red.get(atmp,btmp),i,j,D_SD.get(i,j))
                break
            i += 1

    return res

#*****************************************************************

import main        # import main module of the libra-Gamess-interface code

data, test_data = main.main(params)  # run actual calculations

print "D_mol="
print test_data["D_mol"].show_matrix()
print "D_mol_red="
print test_data["D_mol_red"].show_matrix()
print "D_SD="
print test_data["D_SD"].show_matrix()

res = test1(test_data,data["Ngbf"],Nmin,Nmax,excitations,HOMO)

print "res = 1 : true , res = 0 : false"
print "res=",res

