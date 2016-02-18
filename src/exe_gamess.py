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

## \exe_gamess.py
# This module defines the function which executes GAMESS program in the "GAM_DIR" directory.
#
# Used in : main.py/initial_gamess_exe

import os
import sys
import math

def exe_gamess(params,job):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] params : input files from the submit file , in the directory form.
    # \param[in] job    : the name of GAMESS input file.
    # This function executes gamess program with the input file.
    #
    # Used in:  main.py/

    NPROCS = params["NPROCS"]
    GMS_DIR = params["GMS_DIR"]
    scratch = params["scratch"]
    SRC_DIR = params["SRC_DIR"]
    INP_DIR = params["INP_DIR"]

    os.chdir(GMS_DIR)

    os.system("/usr/bin/time rungms.slurm %s 01 %s > %s.out" % \
              (job,NPROCS,job))

    # copy FOCK MATRIX file to input directory
    #os.system("cp -r %s/%s.F18 %s/%s.fock"% (scratch,job,INP_DIR,job))
    # delete unnecessary files
    os.system("rm -r %s/%s.dat"% (INP_DIR,job))
    os.system("ls %s/"% (scratch))
    os.system("rm -r %s/%s.*"% (scratch,job)) 
    os.chdir(SRC_DIR)
    print "************GAMESS execution finished**************"
