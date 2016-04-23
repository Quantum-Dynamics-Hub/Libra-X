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
## \file create_MD_objects
# This module implements the function which prepares initial objects for NA-MD

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def md_objects(syst,nstates,params):
    ##
    # find the parameters below. 
    # \param[in] list of System objects
    # \param[in] nstates number of excitation states
    # \param[in] params input parameters from run.py
    #
    # Used in:  main.py/main/run_MD 

    print "In run_MD. Creating Hamiltonians"

    # Create an External Hamiltonian object
    ham = []
    for i in xrange(len(syst)):
        ham0 = Hamiltonian_Extern(nstates,3*syst[i].Number_of_atoms) # (electronic DOF,nuclear DOF)
        ham0.set_rep(1)            # adiabatic representation
        ham0.set_adiabatic_opt(0)  # use the externally-computed 
                                  # adiabatic electronic Hamiltonian and derivatives
        ham0.set_vibronic_opt(0)   # use the externally-computed vibronic Hamiltonian and derivatives
        ham.append(ham0)

    # bind actual matrices to external hamiltonian
    ham_adi = []; d1ham_adi = []; ham_vib = []
    for i in xrange(len(ham)):
        ham_adi0 = MATRIX(nstates,nstates)
        ham_adi.append(ham_adi0);  ham[i].bind_ham_adi(ham_adi[i]); # bind adiabatic hamiltonian
        d1ham_adi0= MATRIXList()
        for k in xrange(3*syst[i].Number_of_atoms):
            tmp = MATRIX(nstates,nstates)
            d1ham_adi0.append(tmp)
        d1ham_adi.append(d1ham_adi0); ham[i].bind_d1ham_adi(d1ham_adi[i]) # bind derivative of adiabatic hamiltonian
        ham_vib0 = CMATRIX(nstates,nstates)
        ham_vib.append(ham_vib0);  ham[i].bind_ham_vib(ham_vib[i]); # bind vibronic hamiltonian

    print "after bind hamiltonian";# sys.exit(0)

    print "Addresses of the working matrices"
    for i in xrange(len(ham)):
        print "%i th hamiltonian is" %(i)
        print "ham_adi = ", ham_adi[i]
        print "d1ham_adi = ", d1ham_adi[i]
        for k in xrange(3*syst[i].Number_of_atoms):
            print "d1ham_adi[",k,"]= ", d1ham_adi[i][k]
        print "ham_vib = ", ham_vib[i]

    #sys.exit(0)
    print "Setting nuclear variables"

    # Initialize nuclear variables
    mol = []
    for i in xrange(len(syst)):
        mol.append(Nuclear(3*syst[i].Number_of_atoms))
        syst[i].extract_atomic_q(mol[i].q)
        syst[i].extract_atomic_p(mol[i].p)
        syst[i].extract_atomic_f(mol[i].f)
        syst[i].extract_atomic_mass(mol[i].mass)

        print i,"mol=",mol[i]
        for k in xrange(syst[i].Number_of_atoms):
            print "mol[%d][%d]=%f,%f,%f"%(i,k,mol[i].q[3*k+0],mol[i].q[3*k+1],mol[i].q[3*k+2])  ;# sys.exit(0)

    # Initialize Thermostat object
    if params["MD_type"] == 1:
        therm = []
        for i in xrange(len(syst)):
            therm0 = Thermostat({"nu_therm":params["nu_therm"], "NHC_size":params["NHC_size"], "Temperature":params["Temperature"], "thermostat_type":params["thermostat_type"]})
            therm0.set_Nf_t(3*syst[i].Number_of_atoms)
            therm0.set_Nf_r(0)
            therm0.init_nhc()
            therm.append(therm0)

        print "therm =",therm #; sys.exit(0)

    return ham, ham_adi, d1ham_adi, ham_vib, mol, therm
