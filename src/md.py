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
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


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

    scr_dir = params["scr_dir"]
    rungms = params["rungms"]
    VERNO = params["VERNO"]

    # set environmental variables for GAMESS execution
    os.environ["SCR"] = scr_dir
    os.environ["USERSCR"] = scr_dir
    os.environ["GMSPATH"] = params["GMSPATH"]

    #os.system("/usr/bin/time rungms.slurm %s 01 %s > %s" % (inp,nproc,out))
    print "/usr/bin/time %s %s %s %s > %s" % (rungms,inp,VERNO,nproc,out)
    os.system("/usr/bin/time %s %s %s %s > %s" % (rungms,inp,VERNO,nproc,out))
    sys.exit(0)

    # delete the files except input and output ones to do another GAMESS calculation.
    os.system("rm *.dat")              
    os.system("rm -r %s/*" %(scr_dir)) 

def run_MD(syst,el,ao0,E0,C0,params,label,Q):
    ##
    # This function handles a SINGLE trajectory.
    # When NA-MD is utilized (by specifying the TSH method), we use the CPA with isotropic
    # velocity rescaling
    #
    # \param[in,out] syst a System object that includes atomic system information.
    # \param[in,out] el The the object containig electronic DOFs for the nuclear coordinate
    # given by syst. I have decided to go back a bit - one set of electronic DOF per set of
    # nuclear DOF. This is also needed when we do the velocity rescaling, even if we use the
    # ground state forces for propagation. This also brings a conceptual clarity
    # 
    # \param[in,out] ao0   Atomic orbital basis
    # \param[in,out] E0    Molecular orbital energies
    # \param[in,out] C0    MO-LCAO coefficients
    # \param[in,out] data0 Data extracted from GAMESS output file, in the dictionary form.
    # \param[in,out] params Input data containing all manual settings and some extracted data.

    # This function executes classical MD in Libra and electronic structure calculation
    # in GAMESS iteratively and simulates excited electron dynamics with MF and SH way. 
    # It outputs MD trajectories, Energy evolutions, and SE and SH populations.
    #
    # Used in:  main.py/main
    
    dt_nucl = params["dt_nucl"]
    el_mts = params["el_mts"] # multiple time stepping algorithm for electronic DOF propagation
    if el_mts < 1:
        print "Error in run_MD: el_mts must be positive integer"
        print "Value given = ", el_mts
        print "Exiting..."
        sys.exit(0)
    dt_elec = dt_nucl/float(el_mts)

    Nsnaps = params["Nsnaps"]
    Nsteps = params["Nsteps"]
    nstates = len(params["excitations"])
    print_coherences = params["print_coherences"]
    MD_type = params["MD_type"]
    SH_type = params["SH_type"]
    kB = 3.166811429e-6 # Boltzmann constant in a.u.

    rnd = Random()

    # a flag for potential energy (Ehrenfest or SH)
    f_pot = 0 # Default: Ehrenfest
    if SH_type >= 1: # use SH potential
        f_pot = 1

    #=============== Initialization =======================

    # Open and close energy and trajectory files - this will effectively
    # make them empty (to remove older info, in case we restart calculations)
    fe = open(params["ene_file"],"w")
    fe.close()
    ft = open(params["traj_file"],"w")
    ft.close()
    fm = open(params["mu_file"],"w")
    fm.close()   
    fel = open(params["se_pop_file"],"w")
    fel.close()


    print "In run_MD. Creating Hamiltonians"

    # Create an External Hamiltonian object
    ham = Hamiltonian_Extern(nstates,3*syst.Number_of_atoms) # (electronic DOF,nuclear DOF)
    ham.set_rep(1)            # adiabatic representation
    ham.set_adiabatic_opt(0)  # use the externally-computed 
                                 # adiabatic electronic Hamiltonian and derivatives
    ham.set_vibronic_opt(0)   # use the externally-computed vibronic Hamiltonian and derivatives

        
    # bind actual matrices to external hamiltonian
    ham_adi = MATRIX(nstates,nstates);  ham.bind_ham_adi(ham_adi); # bind adiabatic hamiltonian
    d1ham_adi = MATRIXList()
    for k in xrange(3*syst.Number_of_atoms):
        tmp = MATRIX(nstates,nstates)
        d1ham_adi.append(tmp)
    ham.bind_d1ham_adi(d1ham_adi) # bind derivative of adiabatic hamiltonian
    ham_vib = CMATRIX(nstates,nstates);  ham.bind_ham_vib(ham_vib); # bind vibronic hamiltonian

    print "Addresses of the working matrices"
    print "ham_adi = ", ham_adi
    print "d1ham_adi = ", d1ham_adi
    for k in xrange(3*syst.Number_of_atoms):
        print "d1ham_adi[",k,"]= ", d1ham_adi[k]
    print "ham_vib = ", ham_vib

    print "Setting nuclear variables"

    # Initialize nuclear variables
    mol = Nuclear(3*syst.Number_of_atoms)
    syst.extract_atomic_q(mol.q)
    syst.extract_atomic_p(mol.p)
    syst.extract_atomic_f(mol.f)
    syst.extract_atomic_mass(mol.mass)

    # Initialize Thermostat object
    if MD_type == 1:
        therm = Thermostat({"nu_therm":params["nu_therm"], "NHC_size":params["NHC_size"], "Temperature":params["Temperature"], "thermostat_type":params["thermostat_type"]})
        therm.set_Nf_t(3*syst.Number_of_atoms)
        therm.set_Nf_r(0)
        therm.init_nhc()

    # set initial data from GAMESS output
    ao = []
    for i in range(0,len(ao0)):
        ao.append(AO(ao0[i]))

    E = MATRIX(E0)
    C = MATRIX(C0)
    #data = data0

    # Initialize forces and Hamiltonians
    #epot = data["tot_ene"]  # total energy from GAMESS which is the potential energy acting on nuclei
    #write_gms_inp(data, params, mol)
    #exe_gamess(params)
    #Grad, data, E_mol, D_mol, E_mol_red, D_mol_red = gamess_to_libra(params, ao, E, C, 0) # this will update AO and gradients
    #Hvib, D_SD = vibronic_hamiltonian(params,E_mol_red,D_mol_red,0) # create vibronic hamiltonian

    #sys.exit(0) # DEBUG!!!

    print "Starting propagation"

    #=============== Propagation =======================
    epot, ekin, etot, eext = 0.0, 0.0, 0.0, 0.0
    mu = []
    SH_states = [] # stores SH states

    for i in xrange(Nsnaps):
        syst.set_atomic_q(mol.q)
        syst.print_xyz(params["traj_file"],i)       

        for j in xrange(Nsteps):

            ij = i*Nsteps + j

            # Electronic propagation: half-step
            for k in xrange(el_mts):
                el.propagate_electronic(0.5*dt_elec, ham)


            # >>>>>>>>>>> Nuclear propagation starts <<<<<<<<<<<<
            # Optional thermostat            
            if MD_type == 1: # NVT-MD
                for k in xrange(3*syst.Number_of_atoms):
                    mol.p[k] = mol.p[k] * therm.vel_scale(0.5*dt_nucl)

            mol.propagate_p(0.5*dt_nucl)
            mol.propagate_q(dt_nucl)

            # ======= Compute forces and energies using GAMESS ============
            write_gms_inp(label, Q, params, mol)
            exe_gamess(params)

            # update AO and gradients
            tot_ene, Grad, mu, E_mol, D_mol, E_mol_red, D_mol_red = gamess_to_libra(params, ao, E, C, str(ij))

            #========= Update the matrices that are bound to the Hamiltonian =========
            # Compose electronic and vibronic Hamiltonians
            update_vibronic_hamiltonian(ham_adi, ham_vib, params,E_mol_red,D_mol_red, str(ij))

            print "Addresses of the ham matrices"
            print "ham_adi = ", ham_adi
            print "ham_vib = ", ham_vib
            print "ham_adi "; ham_adi.show_matrix();
            print "ham_vib "; ham_vib.show_matrix();

            for k in xrange(syst.Number_of_atoms):
                for st in xrange(nstates):
                    d1ham_adi[3*k+0].set(st,st,Grad[k].x)
                    d1ham_adi[3*k+1].set(st,st,Grad[k].y)
                    d1ham_adi[3*k+2].set(st,st,Grad[k].z)
           
            #epot = compute_forces(mol, el, ham, 1)  # 0 - Ehrenfest, 1 - TSH
            epot = tot_ene + compute_forces(mol, el, ham, f_pot)  #  f_pot = 0 - Ehrenfest, 1 - TSH

            print "epot= ", epot
            #sys.exit(0)

            ekin = compute_kinetic_energy(mol)
            etot = epot + ekin
          
            if MD_type == 1:
                therm.propagate_nhc(dt_nucl, ekin, 0.0, 0.0)

            mol.propagate_p(0.5*dt_nucl)

            # optional thrmostat
            if MD_type == 1: # NVT-MD
                for k in xrange(3*syst.Number_of_atoms):
                    mol.p[k] = mol.p[k] * therm.vel_scale(0.5*dt_nucl)

            # >>>>>>>>>>> Nuclear propagation ends <<<<<<<<<<<<


            # Electronic propagation: half-step
            for k in xrange(el_mts):
                el.propagate_electronic(0.5*dt_elec, ham)


            ############## Add surface hopping ######################

            print "Before TSH"

            if SH_type>=1:

                # Compute hopping probabilities
                g = MATRIX(nstates,nstates) # initialize a matrix of hopping probability
                use_boltz_factor = 0  # we don't need to use Boltzmann factor, since  
                                      # we are using velocity rescaling in the hopping procedure.
                                      # Although the rescaling doesn't account for the direction, but it
                                      # still accounts for energy partitioning between electronic and
                                      # nuclear DOFs

                if SH_type == 1: # FSSH
                    compute_hopping_probabilities_fssh(mol, el, ham, g, dt_nucl, use_boltz_factor, params["Temperature"])
                elif SH_type == 2: # GFSH
                    compute_hopping_probabilities_gfsh(mol, el, ham, g, dt_nucl, use_boltz_factor, params["Temperature"])
                elif SH_type == 3: # MSSH
                    compute_hopping_probabilities_mssh(mol, el, ham, g, dt_nucl, use_boltz_factor, params["Temperature"])

                #if params["debug_SH_cal"] = 1:
                    
                # Attempt to hop
                ksi = rnd.uniform(0.0,1.0) # generate random number for every trajectory  
                rep = 0 # velocity rescaling will be done based on the total energy conservation,
                        # no derivative couplings will be needed - we don't have them
                        # !!! This option makes do_rescaling and do_reverse not relevant - so
                        # we can set them to any value - they are not used
                do_rescaling = 1
                do_reverse = 1
                el.istate = hop(el.istate, mol, ham, ksi, g, do_rescaling, rep, do_reverse)

            ################### END of TSH ##########################
         
            print "Finished TSH"

            # Re-compute energies, to print
            epot = tot_ene + compute_potential_energy(mol, el, ham, f_pot) #  f_pot = 0 - Ehrenfest, 1 - TSH
            #epot = compute_potential_energy(mol, el, ham, 1)
            print "epot = ", epot

            ekin = compute_kinetic_energy(mol)
            etot = ekin + epot

            ebath = 0.0
            if MD_type == 1:
                ebath = therm.energy()
            eext = etot + ebath
            curr_T = 2.0*ekin/(3*syst.Number_of_atoms*kB)


        ################### Printing results ############################

        # Energy
        fe = open(params["ene_file"],"a")
        fe.write("t= %8.5f ekin= %8.5f  epot= %8.5f  etot= %8.5f  eext= %8.5f curr_T= %8.5f\n" % (ij*dt_nucl, ekin, epot, etot, eext,curr_T))
        fe.close()
        
        # Dipole moment components
        fm = open(params["mu_file"],"a")
        line = "t= %8.5f " % (ij*dt_nucl)
        for k in xrange(len(ao)):
            line = line + " %8.5f %8.5f %8.5f " % (mu[0].get(k,k),mu[1].get(k,k),mu[2].get(k,k))
        line = line + "\n"
        fm.write(line)
        fm.close()

        # Populations            
        fel = open(params["se_pop_file"],"a")

        # Print time
        line_se = "t= %8.5f " % (ij*dt_nucl)

        # Print populations
        for st in xrange(nstates):
            line_se = line_se + " %8.5f " % el.rho(st,st).real

        if print_coherences == 1:
            # Print coherences
            for st in xrange(nstates):
                for st1 in xrange(st):
                    line_se = line_se + " %8.5f %8.5f " % (el.rho(st,st1).real, el.rho(st,st1).imag)
             
        line_se = line_se + "\n"

        fel.write(line_se)
        fel.close()

        #**************** store states as a function of time **************
        if SH_type >= 1:
            SH_states.append(el.istate)


    print "       ********* %i snap ends ***********" % i
    print 

    test_data = {}

    print "SH_states = ", SH_states

    return SH_states


