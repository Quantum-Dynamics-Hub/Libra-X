#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Alexey V. Akimov
#* 
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

## \file spin_index.py 
# This module implements the functions that extract spin index of excitation object

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


def index_spin(params,active_space):
##
# This function creates alpha (up-spin) and beta (down spin) index list
# which are later used to create SD (Slater Determinant) data type
# \param[in] params Parameter dictionary, particularly, params["excitations"] and params["nel"]
#                   is used in this function
# \param[in] active_space List of active space orbital number
# \param[out] alp list of alpha spin index
# \param[out] bet list of beta spin index 
#
    print "in index_spin..."

    nstates = len(params["excitations"])
    # In case of non-spin-polarized calculation
    n_HOMO = params["nel"]/2 
    # Index of HOMO orbital in the active space, e.g., if active space is [4,5,6,7,8] and HOMO
    # orbital is 6, then the h_idx of HOMO in the active space is 2, new index in the active space is
    # [0,1,2,3,4]

    print "nel = ", params["nel"]
    print "n_HOMO = ", n_HOMO
    print "active_space = ", active_space

    h_idx = active_space.index(n_HOMO)
    for ex_st in xrange(nstates): # for each excited configuration
                                 
        idx = params["excitations"][ex_st]
        alp,bet = [],[]

        # use excitation object to create proper SD object for different excited state
        # Currently, excitation only from HOMO is included

        # Consider 2 electrons in the active space:
        # Ground state: 2 electrons are sitting on HOMO (HOMO-1 is not included into our active space)
        # Schematics: U - spin-up, D - spin-down, I - inactive electrons
        #
        #       SD0            SD1           SD2
        #
        #  3  ---------     ---------      ---------
        #  2  ---------     --   D --      -- U   --
        #  1  -- U D --     -- U   --      -- U   --
        #  0  -- I I --     -- I I --      -- I I --

        if idx.from_orbit[0] == idx.to_orbit[0]:
        # For SD0 or the ground state
            alp.append(idx.to_orbit[0] + h_idx)
            bet = alp
        elif idx.from_orbit[0] != idx.to_orbit[0]:
            if idx.from_spin[0] != idx.to_spin[0]:
            # For SD2 or excited state triplet
                bet.append(idx.from_orbit[0] + h_idx)
                bet.append(idx.to_orbit[0]+h_idx)
            else:
            # For SD1 or excited state singlet
                alp.append(idx.to_orbit[0] + h_idx)
                bet.append(idx.from_orbit[0] + h_idx)

    print "alp = ", alp
    print "bet = ", bet
    return alp, bet


