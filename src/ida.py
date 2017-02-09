#*********************************************************************************
#* Copyright (C) 2017 Kosuke Sato, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/
## \file ida.py
# This module implements functions inducing ID-A.

import sys
import cmath
import math
import os

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *


def ida_py_modified(ele, old_st, new_st, E_old, E_new, T, ksi):
    ##
    # This function is a modified version of ida_py where coefficients are passed instead of "Electronic" objects.(see the next function) 
    # \param[in]     ele[in,out] An object containig electronic DOFs                                                          
    # \param[in]      old_st[in] an excitonic state before SH process
    # \param[in]      new_st[in] an excitonic state after SH process
    # \param[in]           T[in] a target Temperature 
    # \param[in]         ksi[in] a uniform random number in the range of (0.0, 1.0)
    #                                                                                                                                            
    # Used in md.py/run_MD

    kb = 3.166811429e-6  # Hartree/K
    res = old_st;  el = Electronic(ele)
    dE = (E_new - E_old);   boltz_f = 1.0

    if dE>0.0:
        argg = dE/(kb*T)
        if argg > 50.0:
            boltz_f = 0.0
        else:
            boltz_f = math.exp(-argg)

        if ksi<boltz_f:
            res = new_st  # accepted hop
            # Collapse the wavefunction to the new state                                                                                                     
            for st in xrange(el.nstates):
                el.q[st] = 0.0; el.p[st] = 0.0;
            el.q[new_st] = 1.0
        else:
            # Unsuccessful hop - collapse wfc back to the original state^M                                                                                   
            for st in xrange(el.nstates):
                el.q[st] = 0.0; el.p[st] = 0.0;
            el.q[old_st] = 1.0
    else:
        res = new_st

    el.istate = res # return res

    return el


def ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi, do_collapse):
    # original version from libra-code/src/developments/frangments-namd

    kb = 3.166811429e-6  # Hartree/K
    res = old_st;  C = CMATRIX(Coeff)
    dE = (E_new - E_old);   boltz_f = 1.0

    if dE>0.0:
        argg = dE/(kb*T)
        if argg > 50.0:
            boltz_f = 0.0
        else:
            boltz_f = math.exp(-argg)

        if ksi<boltz_f:
            res = new_st  # accepted hop
            
            # Collapse the wavefunction to the new state 
            if do_collapse:
                C *= 0.0; C.set(new_st, 1.0+0.0j)
        else:
            # Unsuccessful hop - collapse wfc back to the original state^M
            if do_collapse:
                C *= 0.0; C.set(old_st, 1.0+0.0j)
    else:
        res = new_st

    return res , C
