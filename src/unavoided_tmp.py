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
"""@package unavoided
 This module implements the functions to correct the state crossings (trivial, or 
 "unavoided") crossings.
 The order of the eigenstates is determined by the energies, not the orbital's character.         
 Say, if states i and j are localized on B and A, respectively, at "t+dt" while they were
 on A and B at "t", it is likely that this is just a manifestation of the trivial crossing, 
 not the actual transition. So one must correct the identity of the states.
"""

import os
import sys
import math
import copy
import unittest
import itertools

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *


    

def get_reordering(time_overlap):
    """ This function identifies which states have changed their identities via unavoided 
    (a.k.a. trivial) crossing when the system evolved in the interval from t to t+dt
    We start by looking at the time overlap matrix: <phi_i(t)|phi_i(t+dt)> 
    If no spurious state changes have happened during this time, the diagonal elements
    should be close to 1.0. If they are not - we locate to which state the transitions might 
    have happened.

    \param[in] time_overlap ( MATRIX ) the time overlap matrix, <phi_i(t)|phi_j(t+dt)>.
    Returns:
    perm - list of integers that describe the permutation. That is:
    perm[i] - is the index identifying the "older" state "i". Now, it may be labeled
    by some other index, j. 
    iperm[j] - if we were looking at a specific state labeled "j", iperm[j] this will be
    the index of that state in the present list of states (e.g. energies, etc.)

    """
    # define Random object
    rnd = Random()

    # extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    S = CMATRIX(time_overlap)   # just a temporary working object
    sz = time_overlap.num_of_rows
    perm = range(sz)  # original permutation

    # Before reordering matrices, we must find the indices pointing the same state.
    # If we get a list [0,1,1,3],
    # we find that new states 1 and 2 (starting from 0) will be assigned by old state 1.
    # In this case, they are also close to old state 2.
    # So, we determine randomly 2 transitions,
    # (1,2) -> (1,2) or (1,2) -> (2,1) with the same probabilities (both are 50%).
    # Once doing that, we no longer reorder the columns of the indices.

    pini = range(sz) # for storing initial indices 
    for col in xrange(sz):
        val = 0.0+0.0j
        # Find the max element in the given column "col"                                                                                                  
        [pini[col], val] = S.max_col_elt(col)

    print "pini is"; print pini
    # create a list "du" including duplicate numbers
    t = set()
    du = list(set([x for x in pini if x in t or t.add(x)])) 

    esc = [] # this list including indices for escaping from reordering.
    if len(du) > 0: # when duplication happens in pini
        print "duplicate numbers are"; print du
        for dnum in du:
            ksi = rnd.uniform(0.0, 1.0) # generate random number
            ind = [ i for i , x in enumerate(pini) if dnum == x] # extract indices related to duplication
            print "ksi is %f" % (ksi)
            print "indices concerning duplication are"; print ind
            esc.extend(ind) # extract elements of ind, not the list itself.
            
            ind_sz = len(ind)
            aperm = list(itertools.permutations(ind)) # generate permutation from the numbers of ind
            print aperm
            aperm_sz = len(aperm)
            itra = 0
            for i in xrange(aperm_sz):
                if i < aperm_sz*ksi and aperm_sz*ksi < i+1: # state assigning is done randomly
                    itra = i
            print "selected transition is"; print aperm[itra];

            # reindex the perm list
            for i in xrange(ind_sz): 
                ip = ind[i]
                perm[ip] = aperm[itra][i]
            #print aperm[itra][0], aperm[itra][1]

    print "esc list is"; print esc;

    '''     transition part ends     '''

    for col in xrange(sz):

        indx = -1
        val = 0.0+0.0j
        if all((col!=j for j in esc)): # do not reorder the indices of "esc" list
            while indx!=col and all((col!=j for j in esc)):

                # Find the max element in the given column "col"
                [indx, val] = S.max_col_elt(col)
            
                # Apply permutation (col, indx) to the present "perm" list
                tmp = perm[col]
                perm[col] = perm[indx]
                perm[indx] = tmp

                # Do the corresponding swap of the columns in the S matrix
                S.swap_cols(col,indx)

    return perm



def _test_setup():
    """Prepare the test matrices"""

    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = CMATRIX(4,4)
    a.set(0,0,1.0+0j);  a.set(1,1,1.0+0j); 
    a.set(2,2,1.0+0j);  a.set(3,3,1.0+0j); 
    #print "a = ";  a.show_matrix()

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = CMATRIX(4,4)
    b.set(0,0,1.0+0j);  b.set(1,2,1.0+0j); 
    b.set(2,3,1.0+0j);  b.set(3,1,1.0+0j);
    #print "b = ";  b.show_matrix()

    # non-identical (4x4) matrix
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 1 0 0 0 |
    # | 0 0 0 1 |
    c = CMATRIX(4,4)
    c.set(0,1,1.0+0j);  c.set(1,2,1.0+0j); 
    c.set(2,0,1.0+0j);  c.set(3,3,1.0+0j);
    #print "b = ";  b.show_matrix()

    # non-identical (8x8) matrix
    # | 1 0 0 0 0 0 0 0 |
    # | 0 1 0 0 0 0 0 0 |
    # | 0 0 0 1 0 0 0 0 |
    # | 0 0 1 0 0 0 0 0 |
    # | 0 0 0 0 1 0 0 0 |
    # | 0 0 0 0 0 0 0 1 |
    # | 0 0 0 0 0 1 0 0 |
    # | 0 0 0 0 0 0 1 0 |

    d = CMATRIX(8,8)
    d.set(0,0,1.0+0j);  d.set(1,1,1.0+0j);
    d.set(2,3,1.0+0j);  d.set(3,2,1.0+0j);
    d.set(4,4,1.0+0j);  d.set(5,7,1.0+0j);
    d.set(6,5,1.0+0j);  d.set(7,6,1.0+0j);

    # non-identical matrix (corresponding to doubly mixed states)
    # | 1    0    0 0 |
    # | 0 0.78 0.79 0 |
    # | 0 0.79 0.80 0 |
    # | 0    0    0 1 |

    e = CMATRIX(4,4)
    e.set(0,0,1.0+0j);
    e.set(1,1,0.78+0j); e.set(1,2,0.79+0j);
    e.set(2,1,0.79+0j); e.set(2,2,0.80+0j);
    e.set(3,3,1.0+0j);

    # non-identical matrix (corresponding to triply mixed states)
    # | 1    0    0    0 0 |
    # | 0 0.77 0.78 0.79 0 |
    # | 0 0.78 0.79 0.81 0 |
    # | 0 0.79 0.80 0.82 0 |
    # | 0    0    0    0 1 |

    f = CMATRIX(5,5)
    f.set(0,0,1.0+0j);
    f.set(1,1,0.77+0j); f.set(1,2,0.78+0j); f.set(1,3,0.79+0j);
    f.set(2,1,0.78+0j); f.set(2,2,0.79+0j); f.set(2,3,0.81+0j);
    f.set(3,1,0.79+0j); f.set(3,2,0.80+0j); f.set(3,3,0.82+0j);
    f.set(4,4,1.0+0j);

    # non-identical matrix (corresponding to very complicated condition)
    # | 1 0 0    0    0    0    0|
    # | 0 0 1    0    0    0    0|
    # | 0 1 0    0    0    0    0|
    # | 0 0 0 0.71 0.72    0    0|
    # | 0 0 0 0.74 0.78    0    0|
    # | 0 0 0    0    0 0.77 0.78|
    # | 0 0 0    0    0 0.79 0.81|
    
    g = CMATRIX(7,7)
    g.set(0,0,1.0+0j);
    g.set(1,2,1.0+0j);
    g.set(2,1,1.0+0j);
    g.set(3,3,0.71+0j); g.set(3,4,0.72+0j);
    g.set(4,3,0.74+0j); g.set(4,4,0.78+0j);
    g.set(5,5,0.77+0j); g.set(5,6,0.78+0j);
    g.set(6,5,0.79+0j); g.set(6,6,0.81+0j);

    return a, b, c, d, e, f, g


class TestUnavoided(unittest.TestCase):
    def test_reordering(self):
        """Tests the reordering algorithm"""
        a,b,c,d,e,f,g = _test_setup()

        perm_a = get_reordering(a)
        print "Input matrix "; a.show_matrix()
        print "Permutation = ", perm_a
        self.assertEqual(perm_a, [0,1,2,3])

        perm_b = get_reordering(b)
        print "Input matrix "; b.show_matrix()
        print "Permutation = ", perm_b
        self.assertEqual(perm_b, [0,2,3,1])

        perm_c = get_reordering(c)
        print "Input matrix "; c.show_matrix()
        print "Permutation = ", perm_c
        self.assertEqual(perm_c, [1,2,0,3])

        # Note: the 3 examples can't exactly output the permutation
        #       because some elements of it are determined randomly
        #       , that is, using random number.

        perm_e = get_reordering(e)
        print "Input matrix "; e.show_matrix()
        print "Permutation = ", perm_e

        perm_f = get_reordering(f)
        print "Input matrix "; f.show_matrix()
        print "Permutation = ", perm_f

        perm_g = get_reordering(g)
        print "Input matrix "; g.show_matrix()
        print "Permutation = ", perm_g

if __name__=='__main__':
    unittest.main()

