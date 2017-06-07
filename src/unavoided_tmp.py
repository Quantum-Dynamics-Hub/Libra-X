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
    perm - list of integers that describe all possible permutation. That is:
    perm[i] - contains all possible permutations
    (ex. [0,1,2,...,np-1,0,2,1,...,np-1,...],which can be divided into every "np".
    "np" is assumed to be more than 1 if we have mixed states during calculation.)
    In a permutation, say, [0,1,...,i-1,j,i+1,...], the index identifying the "older" state "i".
    Now, it may be labeled by some other index, j. 
    iperm[j] - if we were looking at a specific state labeled "j", iperm[j] this will be
    the index of that state in the present list of states (e.g. energies, etc.)

    """

    # extract the indices where <phi_i(t)|phi_i(t+dt)> is not close to 1. 
    S = CMATRIX(time_overlap)   # just a temporary working object
    sz = time_overlap.num_of_rows
    #perm = range(sz)  # original permutation
    perm_wrk = range(sz) # working permutation

    # Before getting indices for reordering, we must find ones pointing the same state.
    # If we get a list [0,1,1,3],
    # we find that new states 1 and 2 (starting from 0) will be assigned by old state 1.
    # In this case, they are also close to old state 2.
    # So, we have 2 permutations for transition: (1,2) or (2,1).

    perm_ini = range(sz) # for storing initial indices 
    for col in xrange(sz):
        val = 0.0+0.0j
        # Find the max element in the given column "col"                                                                   
        [perm_ini[col], val] = S.max_col_elt(col)

    #print "pini is"; print pini
    # create a list "du" including duplicate numbers
    t = set()
    du = list(set([x for x in perm_ini if x in t or t.add(x)])) 

    esc = [] # including indices of states for avoid reordering.
    perm_mix = [] # including groups of indices related to mixed states. 
    if len(du) > 0: # when duplication is found in pini
        #print "duplicate numbers are"; print du
        for dnum in du:
            ind = [ i for i , x in enumerate(perm_ini) if dnum == x] # extract indices related to duplication

            # extract the indices pointing the states related to duplication
            ltmp = []
            for i in xrange(len(perm_ini)):
                x = perm_ini[i]
                if x != dnum and x in ind: 
                    ltmp.append(i)
            ind.extend(ltmp)

            #print "mixed states are"; print ind
            esc.extend(ind) # extract elements of ind, not the list itself.
            
            perm_mix.append(ind) # extract groups of mixed states, say, [[0,1][2,3][4,5]]
            #print perm_mix

    #print "indices for mixed states are"; print esc;
    #print "groups of mixed states are"; print perm_mix

    #sys.exit(0)

    ''' reordering part starts '''
    for col in xrange(sz):

        indx = -1
        val = 0.0+0.0j
        cnt = 0
        if all((col!=x for x in esc)): # avoid indices of mixed states
            while indx!=col:

                # Find the max element in the given column "col"
                [indx, val] = S.max_col_elt(col)

                # Apply permutation (col, indx) to the present "perm" list
                tmp = perm_wrk[col]
                perm_wrk[col] = perm_wrk[indx]
                perm_wrk[indx] = tmp

                # Do the corresponding swap of the columns in the S matrix
                S.swap_cols(col,indx)

                # check if this loop ends or not.
                cnt+=1
                if cnt > sz:
                    print "***********************************************"
                    print "reordering counts reached the given threshold "
                    print "The original matrix is "; time_overlap.show_matrix();
                    print "row numbers indicating maximum value of each column are"; print perm_ini;
                    print "duplicate numbers are"; print du;
                    print "column numbers for escaping from reordering are"; print esc;
                    print "The reordered matrix is"; S.show_matrix();
                    print "The (col,indx) pair occuring errors is (%i,%i)" % (col,indx)
                    print "exitting..."; sys.exit(0)

    #print "After being reordered, the working permutation is"; print perm_wrk

    ''' Then, all posible permutations will be generated. '''

    Nperm = 1 # number of permutations
    for ip in perm_mix:
        Nperm *= math.factorial(len(ip)) # in the case [[1,2][3,4]], Nperm=2!x2!=4
    #print "number of permutations is %i" % Nperm

    perm_g = [] # contains groups of permutations, say, [[(1,2),(2,1)][(3,4),(4,3)]]
    for i in xrange(len(perm_mix)):
        perm_g.append(list(itertools.permutations(perm_mix[i])))
    #print "groups of permutations are "; print perm_g

    perm_com = [] # contains combined permutations of perm_g, say,
                  # [(1,2,3,4),(1,2,4,3),(2,1,3,4),(2,1,4,3)]
    if len(perm_g) == 1:
        for i in xrange(Nperm):
            perm_com.append(perm_g[0][i]) 
    if len(perm_g) > 1:
        for i in xrange(len(perm_mix)-1):
            ptmp = []
            for j in xrange(len(perm_g[i])):
                for k in xrange(len(perm_g[i+1])):
                    ptmp.append(perm_g[i][j]+perm_g[i+1][k])
            perm_g[i+1] = ptmp
        perm_com = perm_g[len(perm_mix)-1]
    #print "combined permutations are"; print perm_com

    perm = perm_wrk*Nperm # generate perm_wrk "Nperm" times,
    # say, [0,1,2,3,4,0,1,2,4,3,0,2,1,3,4,0,2,1,4,3]
    for i in xrange(len(perm_com)):
        for j in xrange(len(perm_com[0])):
                perm[i*sz+perm_com[0][j]] = perm_com[i][j]

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

    # non-identical matrix (corresponding to 3 groups of doubly mixed states)
    # | 1    0    0    0    0    0    0|
    # | 0 0.73 0.75    0    0    0    0|
    # | 0 0.74 0.76    0    0    0    0|
    # | 0    0    0 0.71 0.72    0    0|
    # | 0    0    0 0.74 0.78    0    0|
    # | 0    0    0    0    0 0.77 0.78|
    # | 0    0    0    0    0 0.79 0.81|
    
    g = CMATRIX(7,7)
    g.set(0,0,1.0+0j);
    g.set(1,1,0.73+0j); g.set(1,2,0.75+0j);
    g.set(2,1,0.74+0j); g.set(2,2,0.76+0j);
    g.set(3,3,0.71+0j); g.set(3,4,0.72+0j);
    g.set(4,3,0.74+0j); g.set(4,4,0.78+0j);
    g.set(5,5,0.77+0j); g.set(5,6,0.78+0j);
    g.set(6,5,0.79+0j); g.set(6,6,0.81+0j);

    # non-identical matrix (corresponding to triply mixed states)
    # | 1    0    0    0 0 |
    # | 0 0.77 0.78 0.65 0 |
    # | 0 0.76 0.75 0.89 0 |
    # | 0 0.65 0.54 0.77 0 |
    # | 0    0    0    0 1 | 

    h = CMATRIX(5,5)
    h.set(0,0,1.0+0j);
    h.set(1,1,0.77+0j); h.set(1,2,0.78+0j); h.set(1,3,0.65+0j);
    h.set(2,1,0.76+0j); h.set(2,2,0.75+0j); h.set(2,3,0.89+0j);
    h.set(3,1,0.65+0j); h.set(3,2,0.54+0j); h.set(3,3,0.77+0j);
    h.set(4,4,1.00+0j);

    # non-identical matrix (corresponding to triply mixed states)
    # | 1    0    0    0 0 | 
    # | 0 0.65 0.78 0.77 0 |
    # | 0 0.89 0.75 0.76 0 |
    # | 0 0.77 0.54 0.65 0 |
    # | 0    0    0    0 1 |

    h1 = CMATRIX(5,5)
    h1.set(0,0,1.0+0j);
    h1.set(1,3,0.77+0j); h1.set(1,2,0.78+0j); h1.set(1,1,0.65+0j);
    h1.set(2,3,0.76+0j); h1.set(2,2,0.75+0j); h1.set(2,1,0.89+0j);
    h1.set(3,3,0.65+0j); h1.set(3,2,0.54+0j); h1.set(3,1,0.77+0j);
    h1.set(4,4,1.00+0j);

    return a, b, c, d, e, f, g, h, h1


class TestUnavoided(unittest.TestCase):
    def test_reordering(self):
        """Tests the reordering algorithm"""
        a,b,c,d,e,f,g,h, h1 = _test_setup()

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

        perm_e = get_reordering(e)
        print "Input matrix "; e.show_matrix()
        print "Permutation = ", perm_e
        self.assertEqual(perm_e, [0,1,2,3,0,2,1,3])

        perm_f = get_reordering(f)
        print "Input matrix "; f.show_matrix()
        print "Permutation = ", perm_f
        perm_f_eq = [0, 1, 2, 3, 4, 0, 1, 3, 2, 4, 0, 2, 1, 3, 4, 0, 2, 3, 1, 4, 0, 3, 1, 2, 4, 0, 3, 2, 1, 4]
        self.assertEqual(perm_f, perm_f_eq)

        perm_g = get_reordering(g)
        print "Input matrix "; g.show_matrix()
        print "Permutation = ", perm_g
        perm_g_eq = [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 6, 5, 0, 1, 2, 4, 3, 5, 6, 0, 1, 2, 4, 3, 6, 5,\
                         0, 2, 1, 3, 4, 5, 6, 0, 2, 1, 3, 4, 6, 5, 0, 2, 1, 4, 3, 5, 6, 0, 2, 1, 4, 3, 6, 5]
        self.assertEqual(perm_g, perm_g_eq)

        perm_h = get_reordering(h)
        print "Input matrix "; h.show_matrix()
        print "Permutation = ", perm_h
        perm_h_eq = [0, 1, 2, 3, 4, 0, 1, 3, 2, 4, 0, 2, 1, 3, 4, 0, 2, 3, 1, 4, 0, 3, 1, 2, 4, 0, 3, 2, 1, 4]
        self.assertEqual(perm_h, perm_h_eq)

        perm_h1 = get_reordering(h1)
        print "Input matrix "; h1.show_matrix()
        print "Permutation = ", perm_h1
        perm_h1_eq = [0, 1, 2, 3, 4, 0, 3, 2, 1, 4, 0, 1, 3, 2, 4, 0, 2, 3, 1, 4, 0, 3, 1, 2, 4, 0, 2, 1, 3, 4]
        self.assertEqual(perm_h1, perm_h1_eq)

if __name__=='__main__':
    unittest.main()

