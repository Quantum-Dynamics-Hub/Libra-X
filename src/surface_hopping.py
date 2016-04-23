#*********************************************************************************                                                                          
#* Copyright (C) 2016 Kosuke Sato, Alexey V. Akimov                                                                                                         
#*                                                                                                                                                          
#* This file is distributed under the terms of the GNU General Public License                                                                               
#* as published by the Free Software Foundation, either version 2 of                                                                                        
#* the License, or (at your option) any later version.                                                                                                      
#* See the file LICENSE in the root directory of this distribution                                                                                          
#* or <http://www.gnu.org/licenses/>.                                                                                                                        
#***********************************************************************************
## \file surface_hopping.py 
# This module implements the functions which execute Surface Hopping algorithm and
# returns the mol, el, and ham objects

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def surface_hopping(mol,el,ham,params):
    ##  find parameters below
    # \param[in,out] mol Nuclear object
    # \param[in,out] el  Electronic object
    # \param[in,out] ham hamiltonian object
    # \param[in]     params input parameters from run.py
    #
    # Used in:  main.py/main/run_MD 

    dt_nucl = params["dt_nucl"]
    Temperature = params["Temperature"]
    kB = 3.166811429e-6 # Boltzmann constant in hartree unit    
    rnd = Random() # random number generator object
    nstates = el.nstates
    SH_type = params["SH_type"]

    #Compute hopping probabilities
    g = MATRIX(nstates,nstates) # initialize a matrix of hopping probability
    use_boltz_factor = 0        # we don't need to use Boltzmann factor, since  
                                # we are using velocity rescaling in the hopping procedure.
                                # Although the rescaling doesn't account for the direction, but it
                                # still accounts for energy partitioning between electronic and
                                # nuclear DOFs

    if SH_type == 1: # FSSH
        compute_hopping_probabilities_fssh(mol, el, ham, g, dt_nucl, use_boltz_factor, Temperature)
    elif SH_type == 2: # GFSH
        compute_hopping_probabilities_gfsh(mol, el, ham, g, dt_nucl, use_boltz_factor, Temperature)
    elif SH_type == 3: # MSSH
        compute_hopping_probabilities_mssh(mol, el, ham, g, dt_nucl, use_boltz_factor, Temperature)

    # output hopping probability
    if params["debug_SH_cal"] == 1:
        print "hopping_probability is"
        print g.show_matrix()

    # check elements of g matrix are less than 1 or not.
    if params["check_hopping_probs"] == 1:
        for st in xrange(nstates):
            for st1 in xrange(nstates):
                if g.get(st,st1) > 1:
                    print "g(%d,%d) is %f, larger than 1; better to decrease dt_nucl" %(st,st1,g.get(st,st1))

    #Attempt to hop
    ksi = rnd.uniform(0.0,1.0) # generate random number for every trajectory  
    rep = 0 # velocity rescaling will be done based on the total energy conservation,
            # no derivative couplings will be needed - we don't have them
            # !!! This option makes do_rescaling and do_reverse not relevant - so
            # we can set them to any value - they are not used
    do_rescaling = 1
    do_reverse = 1
    el.istate = hop(el.istate, mol, ham, ksi, g, do_rescaling, rep, do_reverse)

    return mol, el, ham
