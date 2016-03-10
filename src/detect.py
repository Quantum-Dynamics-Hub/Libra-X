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

## \file detect.py
# This module implements the functions that detect columns showing
# gradients , molecular energies, molecular orbitals, and atomic basis information
# written in gamess output file.
#
# Used in: main.py/main/nve/nve_MD/gamess_to_libra/unpack_file
#        : main.py/main/initial_gamess_exe/unpack_file

# **************************************************************
# This program detects the columns showing parameters.
# **************************************************************

import os
import sys
import math


def detect_columns(l_gam,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file. 
    # \param[in] params : The list which contains extracted data from l_gam file.
    # This function returns the data extracted from the file, in the form of dictionary
    #
    # Used in: main.py/main/nve/gamess_to_libra/unpack_file/detect
    #        : main.py/main/unpack_file/detect

    i = -1
    for line in l_gam:
        i += 1
        spline = line.split()

        # the number of electrons
        if len(spline) == 5 and spline[2] == "ELECTRONS":
            params["lele"] = i
            params["Nele"] = int(spline[4])

        # the number of occupied orbitals (alpha and beta)
        if len(spline) == 7 and spline[4] == "(ALPHA)":
            params["locc_alp"] = i
            params["Nocc_alp"] = int(spline[6])
        if len(spline) == 8 and spline[4] == "(BETA":
            params["locc_bet"] = i
            params["Nocc_bet"] = int(spline[7])

        # the number of cartesian gaussian basis functions
        if len(spline) == 8 and spline[5] == "FUNCTIONS":
            params["lgbf"] = i
            params["Ngbf"] = int(spline[7])
            #print params["Ngbf"]

        # the atomic basis sets
        if len(spline) == 3 and spline[1] == "BASIS" and spline[2] == "SET":
            params["ab_start"] = i + 7
        if len(spline) == 8 and spline[5] == "SHELLS":
            params["ab_end"] = i - 2

        #***********   single point calculation  ************

        # eigenvectors
        if len(spline) > 0 and spline[0] == "EIGENVECTORS":
            params["mo_start"] = i + 3
            
            params["mo_end"] = i + 1 + (params["Ngbf"] + 4) * int(math.ceil(params["Ngbf"]/5.0))

        # the coordinates of the atoms (in Bohr)
        if len(spline) == 4 and spline[2] == "COORDINATES" and spline[3] == "(BOHR)":
            params["coor_start"] = i + 2
        if len(spline) == 3 and spline[0] == "INTERNUCLEAR" and spline[1] == "DISTANCES":
            params["coor_end"] = i -2
            params["Natoms"] = params["coor_end"] - params["coor_start"] + 1

        # the gradients(in Hartree/Bohr)

        if len(spline) == 4 and spline[0] == "GRADIENT" and spline[3] == "ENERGY":
            params["grad_start"] = i + 4
            params["grad_end"] = params["grad_start"] + params["Natoms"] -1

        # total energy

        if len(spline) == 8 and spline[0] == "FINAL" and spline[2] == "ENERGY":
            params["ltot_ene"] = i
            params["tot_ene"] = float(spline[4])

        #***********   optimization   ***********************

            # molecular orbitals
            #if len(spline) == 2 and spline[0] == "MOLECULAR":
            #    params["mo_start"] = i + 3
            #if len(spline) > 0 and spline[0] == "PROPERTY":
            #    params["mo_end"] = i - 3

            # the coordinates of the atoms (in Angstrom)

            #if len(spline) == 6 and spline[0] == "COORDINATES" and spline[2] == "ALL":
            #    params["coor_start"] = i + 3
            #    params["coor_end"] = params["coor_start"] + params["Natoms"] -1

            # the gradients
            #if len(spline) == 2 and spline[0] == "GRADIENT" and spline[1] == "(HARTREE/BOHR)":
            #    params["grad_start"] = i + 4
            #    params["grad_end"] = params["grad_start"] + params["Natoms"] -1

            #else :
            #print "*********************************************************************"
            #print "********************* CAUTION ***************************************"
            #print "*********************************************************************"
            #print "**** run_type has an illegal value in construct ao_basis funtion ****"
            #print "***                                                               ***"
            #print "********************************************************************"
            #sys.exit()

def show_outputs(l_gam,params):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which includes extracted data from l_gam file.
    # This function shows the columns which includes the information.
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file/detect
    #        : main.py/main/unpack_file/detect

    print "******************************************"
    print "according to the",params["lele"]+1,"th column,"
    print l_gam[params["lele"]]
    print "Nele = ",params["Nele"]
    print "*******************************************"
    print "******************************************"
    print "according to the",params["locc_alp"]+1,"th column,"
    print l_gam[params["locc_alp"]]
    print "Nocc_alp = ",params["Nocc_alp"]
    print "according to the",params["locc_bet"]+1,"th column,"
    print l_gam[params["locc_bet"]]
    print "Nocc_bet = ",params["Nocc_bet"]
    print "*******************************************"
    print
    print "******************************************"
    print "according to the",params["lgbf"]+1,"th column,"
    print l_gam[params["lgbf"]]
    print "Ngbf = ",params["Ngbf"]
    print "*******************************************"
    print
    print "******************************************"
    print "ATOMIC BASIS SET is within",params["ab_start"]+1,"-",params["ab_end"]+1,"th lines."
    for l in range(params["ab_start"],params["ab_end"]+1):
        print l_gam[l]
    print "******************************************"
    print
    print "******************************************"
    print "MOLECULAR ORBITALS is within",params["mo_start"]+1,"-",params["mo_end"]+1,"th lines"
    for l in range(params["mo_start"],params["mo_end"]+1):
        print l_gam[l]
    print "******************************************"
    print
    print "******************************************"
    print "COORDINATES OF ATOMS (in Bohr) is within",params["coor_start"]+1,"-",params["coor_end"]+1,"th lines"
    for l in range(params["coor_start"],params["coor_end"]+1):
        print l_gam[l]
    print "And the number of atoms is ",params["Natoms"]
    print "******************************************"
    print
    print "******************************************"
    print "GRADIENT (in Hartree/Bohr) is within",params["grad_start"]+1,"-",params["grad_end"]+1,"th lines"
    for l in range(params["grad_start"],params["grad_end"]+1):
        print l_gam[l]
    print "******************************************"
    print
    print "******************************************"
    print "according to the",params["ltot_ene"]+1,"th column,"
    print l_gam[params["ltot_ene"]]
    print "total energy = ",params["tot_ene"]
    print "******************************************"
    print


def detect(l_gam,params,flag):
    ## This module detects the columns which includes the information
    #  of atomic basis sets, molecular energies , molecular orbitals,
    #  and gradients.
    # \param[in] l_gam : The list which contains the lines of the (GAMESS output) file.
    # \param[in] params : The list which includes extracted data from l_gam file.
    # \param[in] flag : a flag for debugging detect module
    # This function returns the data extracted from the file, in the form of dictionary
    #
    # Used in: main.py/main/nve_MD/gamess_to_libra/unpack_file
    #        : main.py/main/unpack_file

    detect_columns(l_gam,params)

    if flag == 1:
        show_outputs(l_gam,params)

    print "*********************************************"
    print "detect program ends"
    print "*********************************************"
    print 
