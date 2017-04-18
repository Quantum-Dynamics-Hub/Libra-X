#*********************************************************************************  
#* Copyright (C) 2017 Kosuke Sato, Alexey V. Akimov 
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version. 
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>. 
#* 
#*********************************************************************************/
## \file eigenstates_order.py 
# This module implements functions that check if the order of eigenstates is changed 
# during nuclei propagation from "t" to "t+dt".
# As you know, orders of eigenstates are determined by energy order, not the MO's character.         
# Say, if states i and j are localized on B and A at "t+dt" while they were on A and B at "t", 
# we must commutate the indices of eigenstates. (phi_i(t+dt) <--> phi_j(t+dt)))

import os
import sys
import math
import copy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

def extract_indices(A):
    ## This function extracts the indices (i,j) where <phi_i(t)|phi_i(t+dt)> is not close to 1  
    #  ,but <phi_i(t)|phi_j(t+dt)> is close. 
    # param[in] A MATRIX object including density matrices <phi_i(t)|phi_j(t+dt)>.
    # returned values:
    # ind_old a list including indices where <phi_i(t)|phi_i(t+dt)> is not close to 1.
    # ind_new a list including indices where <phi_i(t)|phi_j(t+dt)> is close to 1.
    #
    # Used in x_to_libra_**.py

    # _thres = 0.9 # overlap threshold

    # extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    ind_old = []; ind_new = [];
    for i in xrange(A.num_of_rows):
        atmp = abs(A.get(i,i))
        for j in xrange(A.num_of_rows):
            if i != j and atmp < abs(A.get(i,j)):
                ind_old.append(i); ind_new.append(j);

    return ind_old, ind_new

def commutate_elements(ind_old,ind_new,E,sd):
    ## This function commutates elements of matrices following ind_old and ind_new indices.
    # param[in]     ind_old a list including indices where <phi_i(t)|phi_i(t+dt)> is not close to 1.
    # param[in]     ind_new a list including indices where <phi_i(t)|phi_j(t+dt)> is close to 1.
    # param[in,out]       E MATRIX object including eigenenergies in diagonal elements.
    # param[in,out]      sd SDlist including eigenvectors.
    #
    # Used in x_to_libra_**.py

    Etmp = []
    for i in ind_old:
        Etmp.append(E.get(i,i))

    # commutate eigenenergies
    for i in xrange(len(ind_old)):
        inew = ind_new[i]
        E.set(inew,inew,Etmp[i])

    # store sd objects
    sdtmp = []
    for i in ind_old:
        #sdtmp[i] = SD(sd[i])
        sdtmp.append(SD(sd[i]))

    # commutate indices of sd
    for i in xrange(len(ind_old)):
        inew = ind_new[i]
        sd[inew] = SD(sdtmp[i])


def test():
    ## This function tests the functions above.

    #************** "extract_indices" test starts *****************

    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = CMATRIX(4,4)
    a.set(0,0,1.0+0j);  a.set(1,1,1.0+0j);  a.set(2,2,1.0+0j);  a.set(3,3,1.0+0j); 
    print "'a' matrix is"
    a.show_matrix()

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = CMATRIX(4,4)
    b.set(0,0,1.0+0j);  b.set(1,2,1.0+0j);  b.set(2,3,1.0+0j);  b.set(3,1,1.0+0j);
    print "'b' matrix is"
    b.show_matrix()

    # non-identical (8x8) matrix
    # | 1 0 0 0 0 0 0 0 |
    # | 0 1 0 0 0 0 0 0 |
    # | 0 0 0 1 0 0 0 0 |
    # | 0 0 1 0 0 0 0 0 |
    # | 0 0 0 0 1 0 0 0 |
    # | 0 0 0 0 0 0 0 1 |
    # | 0 0 0 0 0 1 0 0 |
    # | 0 0 0 0 0 0 1 0 |

    c = CMATRIX(8,8)
    c.set(0,0,1.0+0j);  c.set(1,1,1.0+0j);  c.set(2,3,1.0+0j);  c.set(3,2,1.0+0j);
    c.set(4,4,1.0+0j);  c.set(5,7,1.0+0j);  c.set(6,5,1.0+0j);  c.set(7,6,1.0+0j);

    print "'c' matrix is"
    c.show_matrix()

    a_ind_old, a_ind_new = extract_indices(a)
    b_ind_old, b_ind_new = extract_indices(b)
    c_ind_old, c_ind_new = extract_indices(c)

    print "non-identical indices of 'a' are"
    print a_ind_old
    for i in xrange(len(a_ind_old)):
        print "(%i,%i) element has 1" % (a_ind_old[i],a_ind_new[i])

    print " non-identical indices of 'b' are"
    print b_ind_old
    for i in xrange(len(b_ind_old)):
        print "(%i,%i) element has 1" % (b_ind_old[i],b_ind_new[i])

    print " non-identical indices of 'c' are"
    print c_ind_old
    for i in xrange(len(c_ind_old)):
        print "(%i,%i) element has 1" % (c_ind_old[i],c_ind_new[i])

    #************* "extract_indices" test ends *****************

    #************* "commutate_elements" test starts ****************

    # | 0  0  0  0 |
    # | 0 11  0  0 |
    # | 0  0 22  0 |
    # | 0  0  0 33 |
    Eb = CMATRIX(4,4) # eigenenergy matrix (diagonal)
    Eb.set(0,0,0.0+0j); Eb.set(1,1,11.0+0j); Eb.set(2,2,22.0+0j); Eb.set(3,3,33.0+0j);
    print "'Eb' matrix is"
    Eb.show_matrix()

    mo_pool_alp = CMATRIX(4,4)
    mo_pool_bet = CMATRIX(4,4)

    for i in xrange(4):
        mo_pool_alp.set(i,i, 1.0, 0.0)
        mo_pool_bet.set(i,i, 1.0, 0.0)

    SD0 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1]), Py2Cpp_int([1]) ) # GS
    SD1 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1]), Py2Cpp_int([2]) ) # H->L, singlet
    SD2 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([1,2]), Py2Cpp_int([]) ) # H->L, triplet
    SD3 = SD(mo_pool_alp, mo_pool_bet, Py2Cpp_int([]), Py2Cpp_int([1,2]) ) # H->L, triplet

    sd1 = SDList(); sd1.append(SD0); sd1.append(SD1); sd1.append(SD2); sd1.append(SD3);
    sd2 = SDList(); sd2.append(SD0); sd2.append(SD2); sd2.append(SD3); sd2.append(SD1);

    ovlp11 = SD_overlap(sd1,sd1)
    ovlp12 = SD_overlap(sd1,sd2)

    print "overlap matrix of <sd1|sd1> is"
    SD_overlap(sd1,sd1).show_matrix()

    print "***** before commutation *****"
    print "'Eb' matrix is"
    Eb.show_matrix()
    print "overlap matrix of <sd1|sd2> is"
    SD_overlap(sd1,sd2).show_matrix()

    commutate_elements(b_ind_old,b_ind_new,Eb,sd2)

    print "***** after commutation *****"
    print "transferred 'Eb' matrix is"
    Eb.show_matrix()
    print "overlap matrix of <sd1|sd2> is"
    SD_overlap(sd1,sd2).show_matrix()

   #************* "commutate_elements" test ends ****************

#test()

    
