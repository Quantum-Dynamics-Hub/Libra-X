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
    nstates = len(params["excitations"])
    n_HOMO = params["nel"]/2
    h_idx = active_space.index(n_HOMO)
    for ex_st in xrange(nstates): # for each excited configuration
                                 
        idx = params["excitations"][ex_st]
        alp,bet = [],[]

        # use excitation object to create proper SD object for different excited state
        if idx.from_orbit[0] == idx.to_orbit[0]:
            alp.append(idx.to_orbit[0] + h_idx)
            bet = alp
            print "alp=",alp,"bet=",bet
        elif idx.from_orbit[0] != idx.to_orbit[0]:
            if idx.from_spin[0] != idx.to_spin[0]:
                bet.append(idx.from_orbit[0] + h_idx)
                bet.append(idx.to_orbit[0]+h_idx)
                print "alp=",alp,"bet=",bet

            else:
                alp.append(idx.to_orbit[0] + h_idx)
                bet.append(idx.from_orbit[0] + h_idx)
                print "alp=",alp,"bet=",bet

    return alp, bet


