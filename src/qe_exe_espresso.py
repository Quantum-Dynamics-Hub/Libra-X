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

## \file exe_espresso.py This file implements functions for executing calculations 
# with QE code
#

import os
import sys

def exe_espresso(inp, out):
##
# Function for executing calculations using Quantum Espresso
# once the calculations are finished, all the temporary data are
# deleted
# \param[in] inp The name of the input file
# \param[in] out The name of the output file
#
    a = inp.split('.')
    pref = a[0]
    inexp = pref+".exp.in"   #"x0.exp.in"
    outexp = pref+".exp.out" #"x0.exp.out"
    os.system("srun pw.x < %s > %s" % (inp,out))
    os.system("srun pw_export.x < %s > %s" % (inexp,outexp))

    # Delete scratch directory and unecessary files
    #os.system("rm *.dat *.wfc* *.igk* *.mix*")
    #os.system("rm -r *.save") # not sure if we  need to remove this directory
