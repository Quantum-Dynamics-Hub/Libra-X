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


# **************************************************************
# This program detects the columns showing parameters.
# **************************************************************

import os
import sys
import math


def detect_columns(l_gam,params,runtype):
    #this program detects the columns showing the prameters like
    # coordinates of atoms, gradients acting on atoms,
    # atomic orbitals, eigenenergies and eigenvectors.
    i = -1
    for line in l_gam:
        i += 1
        spline = line.split()

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
        if runtype == 1: 
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

        #***********   optimization   ***********************
        elif runtype == 2: 
            # molecular orbitals
            if len(spline) == 2 and spline[0] == "MOLECULAR":
                params["mo_start"] = i + 3
            if len(spline) > 0 and spline[0] == "PROPERTY":
                params["mo_end"] = i - 3

            # the coordinates of the atoms (in Angstrom)

            if len(spline) == 6 and spline[0] == "COORDINATES" and spline[2] == "ALL":
                params["coor_start"] = i + 3
                params["coor_end"] = params["coor_start"] + params["Natoms"] -1

            # the gradients
            if len(spline) == 2 and spline[0] == "GRADIENT" and spline[1] == "(HARTREE/BOHR)":
                params["grad_start"] = i + 4
                params["grad_end"] = params["grad_start"] + params["Natoms"] -1

        else :
            print "*********************************************************************"
            print "********************* CAUTION ***************************************"
            print "*********************************************************************"
            print "**** run_type has an illegal value in construct ao_basis funtion ****"
            print "***                                                               ***"
            print "********************************************************************"
            sys.exit()

def show_outputs(l_gam,params):
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


def detect(l_gam,params,runtype):
    
    detect_columns(l_gam,params,runtype)

    show_outputs(l_gam,params)

    print "*********************************************"
    print "detect program ends"
    print "*********************************************"
    print 
