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

#>>>>>>>>>>>>>>>> UNCOMMENT THE SECTION BELOW, if THERE IS A PROBLEM WITH PATH
#cwd = "/projects/academic/alexeyak/ekadashi/libracode-dev/libracode-code/_build"
#print "Current working directory", cwd
#sys.path.insert(1,cwd+"/src/mmath")
#sys.path.insert(1,cwd+"/src/context")
#print "\nTest 1: Importing the library and its content"
#from libmmath import *
#from libcontext import *
#<<<<<<<<<<<<<<<<<<<<<<<<<

def qe_extract_mo(filename, upper_tag, active_space):
##
# This function reads an ASCII/XML format file contaiing wavefunction
# and returns the coefficients of the plane wave that constitute the
# wavefunction
#
# \param[in] filename This is the name of the file we will be reading to construct MOs
# \param[in] upper_tag  Currently it is just "Kpoint.1"
# \param[in] active_space The list of indices for the MOs to consider in the calculations
#            In this case, the indices start with 1, not 0
#
# \param[in] n_mo Number of molecular orbital basis used to construct electronic wave function
# \param[in] n_el Number of electrons. e_el/2 is homo index here
#
#

#    homo = n_el/2

    ctx = Context(filename)
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()
    ctx.show_children(upper_tag)  # ("Kpoint.1") #

    ngw = int(float(ctx.get("Info/<xmlattr>/ngw","n")))
    nbnd = int(float(ctx.get("Info/<xmlattr>/nbnd","n")))
    print ngw, nbnd

    coeff = CMATRIX(ngw,n_mo)

    k = 0
#    for band in range(homo,homo+n_mo):   # [6,7,8]: #range(1,nbnd+1):
    for band in active_space:

        c = []
        print "band=",band
        all_coeff = ctx.get("Wfc."+str(band), "n").split(',')
        sz = len(all_coeff)

        for i in xrange(sz):
            a = all_coeff[i].split()
            for j in xrange(len(a)):
                c.append(a[j])

        sz = len(c)
        n = sz/2  # this should be equal to ngw - the number of plane waves 
                  # why /2 ? - because we first read all real and imaginary components
                  # for all planewaves - as two numbers
      
        # Now, we organize the previousely read real and imaginary parts into n complex numbers
        # these numbers are the coefficients with which all planewaves enter the expansion 
        # of the MO with index k (local index) - the one listed as "band" in the global indices
        # scheme. Note, "k" indexing now starts from 0, not 1
        for i in xrange(n):
            coeff.set(i, k, float(c[2*i]), float(c[2*i+1]))
        k = k+1

    nbnd = coeff.num_of_cols    # the number of MOs
    ngw = coeff.num_of_rows     # the number of the planewaves


    # The read MOs (KS orbitals) are not orthonormal, strictly-speaking, - becasue 
    # of the pseudopotentials. So we will normalize them, at least
    for i in xrange(n_mo):
        mo_i = coeff.col(n_mo)
        nrm = (mo_i.H() * mo_i).real
        nrm = (1.0/sqrt(nrm))

        for pw in xrange(ngw):
            coeff.set(pw,i,nrm*coeff.get(pw,i))

    # This returns normalized coefficients of MO
    return coeff

                                 

def qe_extract_coordinates(inp_str, alat, flag):
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


def qe_extract_gradients(inp_str,  flag):
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

def qe_extract(filename, flag, active_space, ex_st): 
##
# Function for reading and extracting Quantum Espresso
# output. Extracted parameters are used in classical MD
# calculation using LIBRA in the next step.
# \param[in] filename The name of the QE output file which we unpack
# \param[in] flag Controls the output: 0 - no additional printing, 1 - yes
# \param[in] active_space The list of indices (starting from 1) of the MOs to include in
# calculations (and to read from the QE output files)
# \param[in] ex_st The index of the currently computing electronic state. This index is
# also used in as a part of the corresponding input/output files
#
    Ry_to_Ha = 0.5
  
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
            tot_ene = Ry_to_Ha*float(s[4]) # so convert energy into atomic units


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
    label, R = qe_extract_coordinates(A[icoord+1:icoord+1+nat], alat, flag)

    # Get gradients
    grads = qe_extract_gradients(A[iforce+4:iforce+4+nat], flag)


    param = {}
    param["nel"] = nel
    param["norb"]= norb
    param["nat"] = nat
    param["alat"]= alat    
    #print params

    # Read the wavefunctions:
    MO = qe_extract_mo("x%i.export/wfc.1" % , "Kpoint.1", active_space)


    return tot_ene, label, R, grads, MO, norb, nel, nat, alat

