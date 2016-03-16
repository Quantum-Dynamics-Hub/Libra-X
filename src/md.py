#*********************************************************************************
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file md.py
# This module implements the functions which execute GAMESS, execute classical MD coupled to 
# Nose-Hoover thermostat, execute TD-SE and SH calculations, and set initial system.
#

from create_gamess_input import *
from gamess_to_libra import *
from vibronic_hamiltonian import *

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
sys.path.insert(1,os.environ["libra_mmath_path"])
sys.path.insert(1,os.environ["libra_chemobjects_path"])
sys.path.insert(1,os.environ["libra_hamiltonian_path"])
sys.path.insert(1,os.environ["libra_dyn_path"])

from libmmath import *
from libchemobjects import *
from libhamiltonian import *
from libdyn import *
from LoadPT import * # Load_PT

##############################################################


def exe_gamess(params):
    ##
    # This is a function that call GAMESS execution on the compute node
    # \param[in] params Input data containing all manual settings and some extracted data.
    #
    # Used in main.py/main
    #         main.py/main/nve_MD

    inp = params["gms_inp"]
    out = params["gms_out"]
    nproc = params["nproc"]
    scr_dir = os.environ['SLURMTMPDIR']
    os.system("/usr/bin/time rungms.slurm %s 01 %s > %s" % (inp,nproc,out))

    # delete the files except input and output ones to do another GAMESS calculation.
    os.system("rm *.dat")              
    os.system("rm -r %s/*" %(scr_dir)) 

def run_MD(SYST,ao0,E0,C0,data0,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in,out] SYST a list of System objects that include atomic system information.
    # \param[in,out] el   The list of the objects containig electronic DOFs for the nuclear coordinate
    #                     given by syst, but may correspond to differently-prepared coherent
    # wavefunctions (different superpositions or sampling over the wfc phase, initial excitations).
    # Under CPA, the propagation of several such variables corresponds to the same nuclear dynamics,
    # we really don't need to recompute electronic structure for each one, which can be used to 
    # accelerate the computations. Now, if you want to go beyond CPA - just use only one object in
    # the el list and run several copies of the run_MD function to average over initial conditions.
    # Also note that even under the CPA, we need to run this function several times - to sample
    # over initial nuclear distribution
    # \param[in,out] ao0   Atomic orbital basis
    # \param[in,out] E0    Molecular orbital energies
    # \param[in,out] C0    MO-LCAO coefficients
    # \param[in,out] data0 Data extracted from GAMESS output file, in the dictionary form.
    # \param[in,out] params Input data containing all manual settings and some extracted data.
    # \param[out] test_data  the output data for debugging, in the form of dictionary

    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively and simulates excited electron dynamics with MF and SH way. 
    # It outputs MD trajectories, Energy evolutions, and SE and SH populations.
    #
    # Used in:  main.py/main

    dt_nucl = params["dt_nucl"]
    el_mts = params["el_mts"] # multiple time stepping algorithm for electronic DOF propagation
    dt_elec = dt_nucl/float(el_mts)
    Nsnaps = params["Nsnaps"]
    Nsteps = params["Nsteps"]
    nstates = len(params["excitations"])
    print_coherences = params["print_coherences"]
    MD_type = params["MD_type"]
    SH_type = params["SH_type"]
    do_reverse = params["do_reverse"]
    nconfig = params["nconfig"]
    ntraj = params["ntraj"]
    do_rescaling = params["do_rescaling"]
    use_boltz_factor = params["use_boltz_factor"]

    kB = 3.166811429e-6 # Boltzmann constant in a.u.

    # Open and close energy and trajectory files - this will effectively 
    # make them empty (to remove older info, in case we restart calculations)

    # MF results file (depending on only initial nuclei trajectories)
    for i in xrange(nconfig):        
        fe = open(params["ene_file"]+"%i_MF.dat"% i,"w")
        fe.close()
        ft = open(params["traj_file"]+"%i_MF.xyz"% i,"w")
        ft.close()
        fm = open(params["mu_file"]+"%i_MF.dat"% i,"w")
        fm.close()

    # SH results file (depending on only initial nuclei trajectories)
    if do_rescaling == 1 and use_boltz_factor == 0 and SH_type > 0 and SH_type < 4 :

        for i in xrange(nconfig):
            for j in xrange(nstates):
                for k in xrange(ntraj):

                    fe = open(params["ene_file"]+"%i_%i_%i_SH.dat"%(i,j,k),"w")
                    fe.close()
                    ft = open(params["traj_file"]+"%i_%i_%i_SH.xyz"%(i,j,k),"w")
                    ft.close()
                    fm = open(params["mu_file"]+"%i_%i_%i_SH.dat"%(i,j,k),"w")
                    fm.close()

    if el_mts < 1:
        print "Error in run_MD: el_mts must be positive integer"
        print "Value given = ", el_mts
        print "Exiting..."
        sys.exit(0)

    # initialize se_pop files
    for iconfig in xrange(nconfig):
        for k in xrange(nstates):
            tmp = params["se_pop_prefix"] + "se_pop_%i_%i"%(iconfig, k)
            fel = open(tmp,"w")
            fel.close()

    if SH_type > 0 and SH_type < 4: # FSSH=1 GSSH=2 MSSH=3

        # initialize sh_pop files
        for iconfig in xrange(nconfig):
            for k in xrange(nstates):
                tmp = params["sh_pop_prefix"] + "sh_pop_%i_%i"%(iconfig, k)
                fel = open(tmp,"w")
                fel.close()

        # external Hamiltonian to calculate SH populations
        syst = SYST[0] # extract syst
        ham_ex = Hamiltonian_Extern(nstates,3*syst.Number_of_atoms)  # (electronic DOF,nuclear DOF)
        ham_ex.set_rep(1)  # adiabatic
        ham_ex.set_adiabatic_opt(0)  # use the externally-computed adiabatic electronic Hamiltonian and derivatives
        ham_ex.set_vibronic_opt(0)  # use the externally-computed vibronic Hamiltonian and derivatives
        
        # bind actual matrices to external hamiltonian
        ham_adi = MATRIX(nstates,nstates);  ham_ex.bind_ham_adi(ham_adi); # bind adiabatic hamiltonian
        d1ham_adi = MATRIXList()
        for i in xrange(3*syst.Number_of_atoms):
            tmp = MATRIX(nstates,nstates)
            d1ham_adi.append(tmp)
        ham_ex.bind_d1ham_adi(d1ham_adi) # bind derivative of adiabatic hamiltonian
        ham_vib = CMATRIX(nstates,nstates);  ham_ex.bind_ham_vib(ham_vib); # bind vibronic hamiltonian

        # initial SH states
        sh_states0 = []
        for k in xrange(nstates):
            for itraj in xrange(ntraj):
                sh_states0.append(k)

        # "Random" object
        rnd = Random()

    print "set initial nuclear variables for each nuclei configuration"
    MOL0 = []
    for iconfig in xrange(nconfig): # select initial nuclei configuration
        syst = SYST[iconfig]
        mol = Nuclear(3*syst.Number_of_atoms)
        syst.extract_atomic_q(mol.q)
        syst.extract_atomic_p(mol.p)
        syst.extract_atomic_f(mol.f)
        syst.extract_atomic_mass(mol.mass)
        MOL0.append(mol0)

        if 0==1:
            print "mass m=",mol.mass[3*i], mol.mass[3*i+1], mol.mass[3*i+2]
            print "coordinates q = ", mol.q[3*i], mol.q[3*i+1], mol.q[3*i+2]
            print "momenta p= ", mol.p[3*i], mol.p[3*i+1], mol.p[3*i+2]
            print "forces f= ",  mol.f[3*i], mol.f[3*i+1], mol.f[3*i+2]
            print "********************************************************"

    print "set initial electronic variables for each excited state configuration"
    el0 = []
    for i_ex in xrange(nstates):  # loop over all initial excitations
        eltmp = Electronic(nstates,i_ex)
        el0.append(eltmp)

    #=============== Propagation =======================

    sh_nstates = 1
    sh_ntraj = 1
    if do_rescaling == 1 and use_boltz_factor == 0 : # each MD for each SH trajectory
        sh_nstates = nstates
        sh_ntraj = ntraj

    for iconfig in xrange(nconfig): # select initial nuclei configuration

        for for i_ex in xrange(sh_nstates):

            print "Initializing SH states"
            sh_states = sh_states0

            for itra in xrange(sh_ntraj): # select a trajectory of excited states propagation

                print "Initializing nuclear variables"
                mol = Nuclear(3*syst.Number_of_atoms)
                for at in xrange(3*syst.Number_of_atoms):
                    mol.q[at] = MOL0[iconfig].q[at]
                    mol.p[at] = MOL0[iconfig].p[at]
                    mol.f[at] = MOL0[iconfig].f[at]
                    mol.mass[at] = MOL0[iconfig].mass[at]

                print "Initializing electronic variables"
                el = []
                for i_ex in xrange(nstates):  # loop over all initial excitations
                    eltmp = Electronic(nstates,i_ex)
                    eltmp.q = el0[i_ex].q
                    eltmp.p = el0[i_ex].p
                    el.append(eltmp)

                # initialize Thermostat object
                if MD_type == 1: # NVT-MD
                    print " Initialize thermostats......"
                    therm = Thermostat({"nu_therm":params["nu_therm"], "NHC_size":params["NHC_size"], "Temperature":params["Temperature"], "thermostat_type":params["thermostat_type"]})
                    therm.set_Nf_t(3*SYST[0].Number_of_atoms)
                    therm.set_Nf_r(0)
                    therm.init_nhc()

                # set initial datas from GAMESS output 
                ao = []
                for i in range(0,len(ao0)):
                    ao.append(AO(ao0[i]))

                E = MATRIX(E0)
                C = MATRIX(C0)
                data = data0

                # set filename for each trajectories

                traj_file = params["traj_file"]+"%i_MF.xyz"% i
                ene_file = params["ene_file"]+"%i_MF.dat"% i
                mu_file = params["mu_file"]+"%i_MF.dat"% i

                for i in xrange(Nsnaps):

                    syst.set_atomic_q(mol.q)
                    syst.print_xyz(traj_file,i)

                    for j in xrange(Nsteps):

                        ij = i*Nsteps + j

                        if ij > 0: # pass this function at t=0
                            # Electronic propagation: half-step
                            for k in xrange(el_mts):
                                for i_ex in range(0,nstates):  # loop over all initial excitations
                                    propagate_electronic(0.5*dt_elec, el[i_ex], Hvib)

                # >>>>>>>>>>> Nuclear propagation starts <<<<<<<<<<<<

                if MD_type == 1: # NVT-MD
                    # velocity scaling
                    for k in xrange(3*syst.Number_of_atoms):
                        mol.p[k] = mol.p[k] * therm.vel_scale(0.5*dt_nucl)

                mol.propagate_p(0.5*dt_nucl)
                mol.propagate_q(dt_nucl)
          
                # ======= Compute forces and energies using GAMESS ============
                write_gms_inp(data, params, mol)
                exe_gamess(params)

                #************ should think of how to store the datas ************
                Grad, data, E_mol, D_mol, E_mol_red, D_mol_red = gamess_to_libra(params, ao, E, C, ij) # this will update AO and gradients
                Hvib, D_SD = vibronic_hamiltonian(params,E_mol_red,D_mol_red,ij) # create vibronic hamiltonian

                epot = data["tot_ene"]         # total energy from GAMESS which is the potential energy acting on nuclei

                for k in xrange(syst.Number_of_atoms):
                    mol.f[3*k]   = -Grad[k][0]
                    mol.f[3*k+1] = -Grad[k][1]
                    mol.f[3*k+2] = -Grad[k][2]

                #========== Propagate thermostat ==================
                ekin = compute_kinetic_energy(mol)

                if MD_type == 1:
                    therm.propagate_nhc(dt_nucl, ekin, 0.0, 0.0)

                mol.propagate_p(0.5*dt_nucl)

                if MD_type == 1: # NVT-MD
                    # velocity scaling
                    for k in xrange(3*syst.Number_of_atoms):
                        mol.p[k] = mol.p[k] * therm.vel_scale(0.5*dt_nucl)

                # >>>>>>>>>>> Nuclear propagation ends <<<<<<<<<<<<

                # Electronic propagation: half-step
                for k in xrange(el_mts):
                    for i_ex in range(0,nstates):  # loop over all initial excitations
                        propagate_electronic(0.5*dt_elec, el[i_ex], Hvib)

                t = dt_nucl*ij # simulation time in a.u.

                # >>>>>>>>>>>>>> SH propagation starts <<<<<<<<<<<<<<<<
                if SH_type > 0 and SH_type < 4: 

                    # update matrices
                    for l in xrange(nstates):
                        ham_adi.set(l,l, Hvib.get(l,l).real)
                        for m in xrange(nstates):
                            ham_vib.set(l,m, Hvib.get(l,m))

                    for l in xrange(3*syst.Number_of_atoms):
                        for m in xrange(nstates):
                            d1ham_adi[l].set(m,m,mol.f[l])

                    if params["debug_ham_ex"] == 1:
                        print "ham_adi =",ham_adi.show_matrix()
                        print "ham_vib =",ham_vib.show_matrix()
                        print "Hvib =",Hvib.show_matrix()
                        print "ham_ex.Hvib(0,0)= " , ham_ex.Hvib(0,0)
                        print "ham_ex.Hvib(0,1)= " , ham_ex.Hvib(0,1)
                        print "d1ham_adi[0] =",d1ham_adi[0].show_matrix()
                        print "mol.f[0] = %8.5f"% mol.f[0]

                    # select SH methods
                    for k in xrange(nstates):
                        g = MATRIX(nstates,nstates) # initialize a matrix of hopping probability
                        if SH_type == 1: # FSSH
                            compute_hopping_probabilities_fssh(mol, el[k], ham_ex, g, dt_nucl, use_boltz_factor, params["Temperature"])
                        if SH_type == 2: # GSSH
                            compute_hopping_probabilities_gssh(mol, el[k], ham_ex, g, dt_nucl, use_boltz_factor, params["Temperature"])
                        if SH_type == 3: # MSSH
                            compute_hopping_probabilities_mssh(mol, el[k], ham_ex, g, dt_nucl, use_boltz_factor, params["Temperature"])

                        # Hopping current states
                        for itraj in xrange(ntraj):
                            ksi = rnd.uniform(0.0,1.0) # generate random number for every trajectory
                            sh_states[k*ntraj+itraj]= hop(sh_states[k*ntraj+itraj], mol, ham_ex, ksi, g, do_rescaling, 1, do_reverse)

                        if params["debug_SH_calculations"] == 1:
                            print "%i state g matrix ="% k ;print g.show_matrix();
                            print "%i sh_states=" %k; print sh_states[k*ntraj:(k+1)*ntraj]

            # >>>>>>>>>>>>>> Compute energies <<<<<<<<<<<<<<<<<<<<<<<<<
            ekin = compute_kinetic_energy(mol)
            etot = ekin + epot

            ebath = 0.0
            if MD_type == 1:
                ebath = therm.energy()

            eext = etot + ebath
            curr_T = 2.0*ekin/(3*syst.Number_of_atoms*kB)

            # >>>>>>>>>>>>>> Compute SH populations <<<<<<<<<<<<<<<<<<<<

            sh_pops = []
            for k in xrange(nstates):
                sh_tmp = sh_states[k*ntraj:(k+1)*ntraj]
                for s in xrange(nstates):
                    sh_pops.append(float(sh_tmp.count(s))/float(ntraj))

            ################### Printing results ############################

            # Energy
            fe = open(ene_file,"a")
            fe.write("i= %5i ekin= %8.5f  epot= %8.5f  etot= %8.5f  eext= %8.5f curr_T= %8.5f\n" % (i, ekin, epot, etot, eext,curr_T))
            fe.close()
        
            # Dipole moment components
            fm = open(mu_file,"a")
            line = "t= %8.5f " % t
            for k in xrange(len(ao)):
                line = line + " %8.5f %8.5f %8.5f " % (data["mu_x"].get(k,k),data["mu_y"].get(k,k),data["mu_z"].get(k,k))
            line = line + "\n"
            fm.write(line)
            fm.close()

            # Populations
            for k in xrange(nstates):            
                tmp = params["se_pop_prefix"] + "se_pop_" + str(iconfig) +"_"+ str(k)
                fel_se = open(tmp,"a")

                tmp = params["sh_pop_prefix"] + "sh_pop_" + str(iconfig) +"_"+ str(k)
                fel_sh = open(tmp,"a")

                # Print time
                line_se = "t= %8.5f " % t
                line_sh = "t= %8.5f " % t

                # Print populations
                for st in xrange(nstates):
                    line_se = line_se + " %8.5f " % el[k].rho(st,st).real
                    line_sh = line_sh + " %8.5f " % sh_pops[k*nstates+st]

                if print_coherences == 1:
                    # Print coherences
                    for st in xrange(nstates):
                        for st1 in xrange(st):
                            line_se = line_se + " %8.5f %8.5f " % (el[k].rho(st,st1).real, el[k].rho(st,st1).imag)
             
                line_se = line_se + "\n"
                line_sh = line_sh + "\n"

                fel_se.write(line_se)
                fel_se.close()
                fel_sh.write(line_sh)
                fel_sh.close()

            print "********* No.%i snap ends" % i

        print "**************** Initial Nuclei Configuration No.%i finished**************" % iconfig

    # input test_data for debugging
    test_data = {}
    test_data["D_mol"] = D_mol
    test_data["D_mol_red"] = D_mol_red
    test_data["D_SD"] = D_SD

    return test_data

def init_system(data, g, T):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] data   The list of variables, containing atomic element names and coordinates
    # \param[in] g      The list of gradients on all atoms
    # \param[in] T      target temperature used to initialize momenta of atoms.
    # This function returns System object which will be used in classical MD.
    #
    # Used in:  main.py/main

    # Create Universe and populate it
    U = Universe();   Load_PT(U, "elements.txt", 0)

    syst = System()

    sz = len(data["coor_atoms"])
    for i in xrange(sz):
        atom_dict = {} 
        atom_dict["Atom_element"] = data["l_atoms"][i]

        # warning: below we take coordinates in Angstroms, no need for conversion here - it will be
        # done inside
        atom_dict["Atom_cm_x"] = data["coor_atoms"][i][0]
        atom_dict["Atom_cm_y"] = data["coor_atoms"][i][1]
        atom_dict["Atom_cm_z"] = data["coor_atoms"][i][2]

        print "CREATE_ATOM ",atom_dict["Atom_element"]
        at = Atom(U, atom_dict)
        at.Atom_RB.rb_force = VECTOR(-g[i][0], -g[i][1], -g[i][2])

        syst.CREATE_ATOM(at)

    syst.show_atoms()
    print "Number of atoms in the system = ", syst.Number_of_atoms

    # initialize momenta of the system where the temperature is T(K). 
    syst.init_atom_velocities(T)
    
    return syst
