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

from libra_py import *

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
    #print "path=", ctx.get_path()
    #ctx.show_children(upper_tag)  # ("Kpoint.1") #

    ngw = int(float(ctx.get("Info/<xmlattr>/ngw","n")))
    nbnd = int(float(ctx.get("Info/<xmlattr>/nbnd","n")))
    #print ngw, nbnd

    n_mo = len(active_space)

    coeff = CMATRIX(ngw,n_mo)

    k = 0
#    for band in range(homo,homo+n_mo):   # [6,7,8]: #range(1,nbnd+1):
    for band in active_space:

        c = []
        #print "band=",band
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
        mo_i = coeff.col(i)
        nrm = (mo_i.H() * mo_i).get(0,0).real
        nrm = (1.0/math.sqrt(nrm))

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

#alp_d,bet_d = check_eig_deg()

def fermi_pop(e):

    #e = [-1.0, -0.5, -0.4] Energies of list of orbitals
    a = MATRIX(4,4)
    for i in range(0,4):
        for j in range(0,4):
            if i==j:
                a.set(i,j, e[i])
            else:
                a.set(i,j, 0.0)

    #a.show_matrix()
    Nel = 2.0 # Two electrons in 4 orbitals
    degen = 1.0 # One orbital can have 1 electrons, in this case of spin-polarization
    kT = 0.025
    etol = 0.0001
    #etol = 0.1
    Ef = fermi_energy(e, Nel, degen, kT, etol)

    bnds = order_bands(a)
    #print bnds
    #print "\n Test5: populate bands"
    Nocc = Nel
    pop_opt = 1
    occ_new = populate_bands(Nel, degen, kT, etol, pop_opt, bnds)
    #print occ
    #print occ[0][1]
    tot_elec = 0.0
    for i in xrange(4):
        occ_new[i][1] = "%4.2f"%occ_new[i][1]
        tot_elec = tot_elec + float(occ_new[i][1])
#        print occ[i][1], tot_elec
    #print occ[0][1]
#    print (2.0 - tot_elec)
#    print occ[3][1]
    occ_new[3][1] = "%4.2f"%(float(occ_new[3][1]) + (Nel - tot_elec))
#    print "%4.2f"%occ[3][1]
#######################################################
    if float(occ_new[3][1]) < 0.0 :
        occ_new[1][1] = float(occ_new[1][1]) + float(occ_new[3][1])
        occ_new[3][1] = 0.00
#######################################################

    return occ_new


#def get_active_mo_en(filename):
def qe_extract_eigenvalues(filename,nel):
#def extract_orb_energy(HOMO):
    f = open(filename,"r")
    #print "filename=",filename
    a = f.readlines()
    #f.close()
    na = len(a)
    HOMO = (nel/2 + nel%2) - 1
    en_alp = []
    for i in xrange(na):
        fa = a[i].split()
        if len(fa) > 0 and fa[0] =="<EIGENVALUES":
            eig_num = i + 1
            #print "Found Eigen values",eig_num

        #if len(fa) > 0 and fa[0] =="<OCCUPATIONS":
        #    occ_num = i + 1
        #    print "Found OCCUPATIONS", occ_num

    eig_val = []
    #occ_val = []
    for ia in range(eig_num,eig_num+HOMO+4):
        fa = a[ia].split()
        eig_val.append("%12.8f"%(float(fa[0])))
    en_alp.append(float(eig_val[HOMO-1]))
    en_alp.append(float(eig_val[HOMO]))
    en_alp.append(float(eig_val[HOMO+1]))
    en_alp.append(float(eig_val[HOMO+2]))
    #print "active space energies",en_alp

    return en_alp

def gen_new_occ(ex_st,nel):
#def extract_orb_energy(HOMO):

    cwd1 = os.getcwd()
    en_alp = get_active_mo_en(cwd1+"/x%i.save/K00001/eigenval1.xml"%ex_st,nel)
    en_bet = get_active_mo_en(cwd1+"/x%i.save/K00001/eigenval2.xml"%ex_st,nel)
    #print "cwd",cwd1

    # push HOMO-1, HOMO, LUMO and LUMO+1 energies to get fermi energies
    # and fermi population, so, another function needed.
    occ_alp_new = fermi_pop(en_alp)

    # Similarly, compute fermi population for Beta spin, as this is a spin-polarized calculation
    occ_bet_new = fermi_pop(en_bet)

    return occ_alp_new, occ_bet_new


def write_qe_input(filename, nel,norb, flag_a, flag_b):
    #norb = 16
    nHOMO = nel/2 + nel%2 # HOMO orbital index
    if norb%10 ==0:
        nl_spin_orb = norb/10      # Number of lines in spin orbital
    else:
        nl_spin_orb = norb/10 + 1  # For Ethylene, it is 2, 12/10 + 1

    #print "nHOMO=",nHOMO
    fa = filename.split('.')
    fb = fa[0]+".scf_wrk.in"
    f = open(fb,"r+")
    a = f.readlines()
    N = len(a)
    f.close()
    f = open(fb, "w")
    i_homo = -1 # initializing with random value
    i_homo_bet = -1 # initializing
    homo_idx = (nHOMO%10) -1
    for i in range(0,N):
        s = a[i].split()
        if len(s)>0 and s[0] =="OCCUPATIONS":
            i_alp = i+1
            #print "i_alp=",i_alp
            i_homo = i_alp + nHOMO/10
            i_bet = i + nl_spin_orb + 2
            i_homo_bet = i_bet + nHOMO/10
            #print "i_homo=",i_homo

        if i==i_homo:
            s[4],s[5],s[6],s[7] = flag_a[0][1],flag_a[1][1],flag_a[2][1],flag_a[3][1]
            a[i] = " ".join(str(x) for x in s)+'\n'
        if i==i_homo_bet:
            s[4],s[5],s[6],s[7] = flag_b[0][1],flag_b[1][1],flag_b[2][1],flag_b[3][1]
            a[i] = " ".join(str(x) for x in s)+'\n'

        f.write(a[i])
    f.close()

def extract_ene_force(filename):
    Ry_to_Ha = 0.5
    f_qe = open(filename, "r")
    A = f_qe.readlines()
    f_qe.close()
    iforce = -1
    tot_ene = 0.0
    nlines = len(A)
    for a in A:
        s = a.split() 
        # Lines forces start
        if len(s) > 0 and s[0] == "Forces" and s[1] == "acting":
            iforce = A.index(a)
        # !    total energy              =     -27.62882078 Ry
        if len(s) > 0 and s[0] == "!" and s[1] == "total" and s[2] == "energy":
            tot_ene = Ry_to_Ha*float(s[4]) # so convert energy into atomic units
    return tot_ene,iforce

def exe_espresso(i):
##
# Function for executing calculations using Quantum Espresso
# once the calculations are finished, all the temporary data are
# deleted
# \param[in] inp The name of the input file
# \param[in] out The name of the output file
#
    inp = "x%i.scf_wrk.in" % i # e.g. "x0.scf_wrk.in"
    out = "x%i.scf.out" % i    # e.g. "x0.scf.out"
    inexp = "x%i.exp.in" % i   # e.g. "x0.exp.in"
    outexp = "x%i.exp.out" % i # e.g "x0.exp.out"

    os.system("srun pw.x < %s > %s" % (inp,out))
    os.system("srun pw_export.x < %s > %s" % (inexp,outexp))

    # Delete scratch directory and unecessary files
    #os.system("rm *.dat *.wfc* *.igk* *.mix*")
    #os.system("rm -r *.save") # not sure if we  need to remove this directory



def robust_cal_extract(filename, ex_st, nel, flag_a, flag_b):
    write_qe_input(filename, nel,flag_a, flag_b)
    exe_espresso(ex_st)
    tot_ene,iforce = extract_ene_force(filename)
    return tot_ene, iforce
 
def check_convergence(filename):

    f_out = open(filename, "r")
    A = f_out.readlines()
    f_out.close()

    status = 1 # 0 is okay, non-zero is somethig else
    nlines = len(A)

    for a in A:
        s = a.split()
        # Line says "convergence has been achieved in -- iterations"
        if len(s) > 0 and s[0] == "convergence"  and s[3] == "achieved":
            status = 0
       # if len(s) > 0 and s[0] == "Forces" and s[1] == "acting":
       #     status = 0

    return status


def qe_extract_info(filename, ex_st, flag): 
##
# This function reads Quantum Espresso output and extracts 
# the descriptive data.
# \param[in] filename The name of the QE output file which we unpack
# \param[in] flag Controls the output: 0 - no additional printing, 1 - yes
#
    Ry_to_Ha = 0.5
  
    f_qe = open(filename, "r")
    A = f_qe.readlines()
    f_qe.close()
   

    alat = -1.0
    nel, norb, nat = -1, -1, -1
    icoord, iforce = -1, -1
    tot_ene = 0.0

    status = 0 # 0 is okay, non-zero is somethig else

    nlines = len(A)

    for a in A:
        s = a.split() 
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
        print "Error in unpack_file first time, Forces not found \n"
        #status = 1

        #for la in xrange(5):

####        a1,b1 = gen_new_occ(ex_st,nel)

        #    tot_ene, iforce = robust_cal_extract(filename, ex_st, nel, a1,b1)
        #    if iforce !=-1:
        #        break
        #alp_d,bet_d = check_eig_deg()
        #flag_a = 1
        #if iforce ==-1:
        #    print "Something new problem in SCF convergence"
        #    sys.exit(0)


    

    return tot_ene, norb, nel, nat, alat, icoord, iforce




def qe_extract(filename, active_space, ex_st, nspin, flag):
#def qe_extract(filename, active_space, ex_st, nspin, verbose_level):
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

    #f_qe = open(filename, "r")
    #A = f_qe.readlines()
    #f_qe.close()
    grads = []
    # Read the descriptive info
    tot_ene, norb, nel, nat, alat, icoord, iforce = qe_extract_info(filename, ex_st, flag)

    f_qe = open(filename, "r")
    A = f_qe.readlines()
    f_qe.close()

    # Read atom names and xyz coordinates    
    label, R = qe_extract_coordinates(A[icoord+1:icoord+1+nat], alat, flag)

    # Get gradients
    grads = qe_extract_gradients(A[iforce+4:iforce+4+nat], flag)

    MO_a, MO_b = None, None

    if nspin <= 1:
        # Read the wavefunctions:
        #MO_a = qe_extract_mo("x%i.export/wfc.1" % ex_st, "Kpoint.1", active_space)
        #MO_a = QE_methods.read_qe_wfc("x%i.export/wfc.1" % ex_st, "Kpoint.1", active_space)
        MO_a = QE_methods.read_qe_wfc("x%i.export/wfc.1" % ex_st, active_space,0)
        MO_b = CMATRIX(MO_a)

    if nspin == 2:
        # Read the wavefunctions:
        #MO_a = qe_extract_mo("x%i.export/wfc.1" % ex_st, "Kpoint.1", active_space)
        #MO_b = qe_extract_mo("x%i.export/wfc.2" % ex_st, "Kpoint.2", active_space)
        #MO_a = QE_methods.read_qe_wfc("x%i.export/wfc.1" % ex_st, "Kpoint.1", active_space)
        #MO_b = QE_methods.read_qe_wfc("x%i.export/wfc.2" % ex_st, "Kpoint.2", active_space)

        MO_a = QE_methods.read_qe_wfc("x%i.export/wfc.1" % ex_st, active_space,0)
        MO_b = QE_methods.read_qe_wfc("x%i.export/wfc.2" % ex_st, active_space,0)


    return tot_ene, label, R, grads, MO_a, MO_b, norb, nel, nat, alat


