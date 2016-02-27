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
## \file reduce_matrix.py
# This module defines the functions that reduces matrix to be the same size as that of active space.
# and returns the reduced matrix.

import os
import sys
import math

sys.path.insert(1,os.environ["libra_mmath_path"])
from libmmath import *

def reduce_matrix(M_ori,active_space):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in,out]        M_ori   original matrix
    # \param[in]     active_space   the molecular orbitals taken in TD-SE calculation.
    # This function returns the reduced matrix "M_red".
    #
    # Used in: main.py/main/run_MD/gamess_to_libra

    # alined to python index
    Nmin = active_space[0] - 1  
    Nmax = active_space[-1] - 1  
    
    M_red = MATRIX(Nmax-Nmin+1,Nmax-Nmin+1)
    for i in range(Nmin,Nmax+1):
        for j in range(Nmin,Nmax+1):
            M_red.set(i-Nmin,j-Nmin,M_ori.get(i,j))


    return M_red
