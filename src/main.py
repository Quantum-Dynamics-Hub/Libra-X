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

# ******************************************************************************************
# This is a GAMESS/Libra "interface" program which communicates the GAMESS output data
# to Libra and vice versa.
# GAMESS is used for Semi-Empirical Calculation and Libra for Classical MD.
# ******************************************************************************************

from gamess_to_libra import *

import os
import sys
import math

# First, we add the location of the library to test to the PYTHON path
#cwd = "/projects/academic/alexeyak/alexeyak/libra-dev/libracode-code"
#print "Using the Libra installation at", cwd
#sys.path.insert(1,cwd+"/_build/src/mmath")
#sys.path.insert(1,cwd+"/_build/src/qchem")

#print "\nTest 1: Importing the library and its content"
#from libmmath import *
#from libqchem import *

#*********************************************************
#***************** input parameters **********************
#*********************************************************

params = {}

params["gamess_out1"] = "../gam_out/H2O_1.out"  # GAMESS output file for time "t"
#params["gamess_out2"] = "../gam_out/exam03_AM1_single.out"  # GAMESS output file for time "t+dt"    
params["gamess_out2"] = "../gam_out/H2O_2.out"              # GAMESS output file for time "t+dt"
params["basis_sets"] = 2      # ab initio or Semi-Empirical calculation? 
                              # Options: "ab_initio" = 1 , "semi_empirical" = 2
params["runtype"] = 1         # single point or optimization? 
                              # Options: single = 1 , "optimization" = 2
params["dt_nuc"] = 1.0        # time step for nuclear dynamics (in fsec)

print "params: ",params

# ************************************************************************* 
# extract parameters from gamess and communicate them to Libra.

gamess_to_libra(params)
