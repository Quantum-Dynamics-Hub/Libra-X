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


#from detect import *
#from detect1 import *
#from extract import *
#from ao_basis import *
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


#**********************************************************
#* This program extracts parameters from the GAMESS output
#* and communicate them to Libra modules.
#* Inside the modules, gradients acting on atoms are used 
#* for classical molecular dynamics
#* and eigenenergies and eigenfunctions are used 
#* for simulating excited electron dynamics. 
#**********************************************************

gamess_to_libra()
