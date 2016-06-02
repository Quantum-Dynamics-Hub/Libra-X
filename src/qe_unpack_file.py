#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*^M
#*********************************************************************************/^M


## \file unpack_file.py This module implements functions for 
# extracting data from the QE output file 
#

import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


#sys.path.insert(1,os.environ["libra_mmath_path"])
#from libmmath import *


def extract_qe_coordinates(inp_str, alat, flag):
    ##
    # Extracts atomic labels, nuclear charges, and coordinates of all atoms
    # from the the list of input lines
    # each input line is assumed to have the format:
    #         1           C   tau(   1) = (  -5.0405786   2.3190577   0.0223876  )
    #         2           C   tau(   2) = (  -3.7211433   2.3070646  -0.0089072  )
    # ...

    # \param[in] inp_str  Strings containing the info for all atoms
    # \param[in] alat Lattice constant for conversion from Angstroms to Bohrs
    # \param[in] flag Controls printing options
    # label - returned list of atomic labels (strings)
    # R - returned list of nuclear coordinates (VECTOR objects)
    #

    label, R = [], []

    for a in inp_str: 
        spline = a.split() 

        # atom labels
        label.append(spline[1])

        # coordinates of atoms
        x = float(spline[6]) * alat
        y = float(spline[7]) * alat
        z = float(spline[8]) * alat
        r = VECTOR(x,y,z)        
        R.append(r)

    if flag == 1:
        print "label=", label
        print "coor_atoms="
        for r in R:
            print R.index(r), r, r.x, r.y, r.z

    # Returned R is in Bohr units
    return label, R


def extract_qe_gradients(inp_str,  flag):
    ##
    # Extracts atomic labels, nuclear charges, and coordinates of all atoms
    # from the the list of input lines
    # each input line is assumed to have the format:   
    # atom    1 type  1   force =     0.00679795    0.03542284    0.02666430
    # atom    2 type  1   force =    -0.00865430   -0.03138108   -0.02321979
    # ...

    # \param[in] inp_str  Strings containing the info for all atoms
    # \param[in] flag Controls printing options
    # grads - returned list of nuclear gradient (VECTOR objects)
    #

    grads = []
    Ry_to_Ha = 0.5

    for a in inp_str:
        spline = a.split()

        # forces acting on atoms   
        fx = Ry_to_Ha * float(spline[6]) 
        fy = Ry_to_Ha * float(spline[7])
        fz = Ry_to_Ha * float(spline[8])
        g = VECTOR(-fx,-fy,-fz)
        grads.append(g)

    if flag == 1:
        print "gradients="
        for g in grads:
            print grads.index(g), g, g.x, g.y, g.z

    # Gradients, in units Ha/Bohr
    return grads



#def unpack_file(filename,params, flag): 
def unpack_file(filename, flag,flag1): 
##
# Function for reading and extracting Quantum Espresso
# output. Extracted parameters are used in classical MD
# calculation using LIBRA in the next step.
# \param[in] filename The name of the QE output file which we unpack
# \param[out] params The dictionary containing control parameters
# \param[in] flag Controls the output: 0 - no additional printing, 1 - yes
#
    
    f_qe = open(filename, "r")
    A = f_qe.readlines()
    f_qe.close()
   

    alat = -1.0
    nel, norb, nat = -1, -1, -1
    icoord, iforce = -1, -1
    tot_ene = 0.0

    nlines = len(A)

    #for i in range(0,nlines):
    #    s = A[i].split()  #is this A or a????
    for a in A:
        s = a.split()  #is this A or a????
        # Lines where positions and forces start
        # example:
        #     site n.     atom                  positions (alat units)
        #         1           C   tau(   1) = (  -5.0405786   2.3190577   0.0223876  )
        #         2           C   tau(   2) = (  -3.7211433   2.3070646  -0.0089072  )
        # ...
        if len(s) > 0 and s[0] == "site"  and s[3] == "positions":
            icoord = A.index(a)           
        if len(s) > 0 and s[0] == "Forces" and s[1] == "acting":
            iforce = A.index(a)

        # Descriptors:
        # example:
        # number of electrons       =        12.00
        if len(s) > 0 and s[0] == "number" and s[2] == "electrons":
            nel = int(float(s[4]))

        # example:
        # number of Kohn-Sham states=           12  
        if len(s) > 0 and s[0] == "number" and s[2] == "Kohn-Sham":
            norb = int(float(s[4]))

        # example:
        # number of atoms/cell      =            6
        if len(s) > 0 and s[0] == "number" and s[2] == "atoms/cell":
            nat = int(float(s[4]))

        # example:
        # lattice parameter (alat)  =       1.8900  a.u.
        if len(s) > 0 and s[0] == "lattice" and s[1] == "parameter" and s[2] == "(alat)":
            alat = float(s[4])

        # example:
        # !    total energy              =     -27.62882078 Ry
        if len(s) > 0 and s[0] == "!" and s[1] == "total" and s[2] == "energy":
            tot_ene = float(s[4])


    if alat<0.0:
        print "Error in unpack_file\n"
        print "Lattice parameter is not found. Exiting...\n"
        sys.exit(0)
    if nel==-1:
        print "Error in unpack_file\n"
        print "The number of electronis is not found. Exiting...\n"
        sys.exit(0)
    if nat==-1:
        print "Error in unpack_file\n"
        print "The number of atoms is not found. Exiting...\n"
        sys.exit(0)
    if norb==-1:
        print "Error in unpack_file\n"
        print "The number of bands (orbitals) is not found. Exiting...\n"
        sys.exit(0)
    if icoord==-1:
        print "Error in unpack_file\n"
        print "Coordinates of atoms are not found. Exiting...\n"
        sys.exit(0)
    if iforce==-1:
        print "Error in unpack_file\n"
        print "Force of atoms are not found. Exiting...\n"
        sys.exit(0)
    


    # Reading atom names and xyz coordinates    
    label, R = extract_qe_coordinates(A[icoord+1:icoord+1+nat], alat, flag)

    # Get gradients
    grads = extract_qe_gradients(A[iforce+4:iforce+4+nat], flag)
    param = {}
    param["nel"] = nel
    param["norb"]= norb
    param["nat"] = nat
    param["alat"]= alat    
    #print params
    if flag1 == 1:
        return tot_ene, label, R, grads, norb, nel, nat,alat
    else:
        return tot_ene, label, R, grads


