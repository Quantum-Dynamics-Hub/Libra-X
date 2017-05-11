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
## \file reorder_matrices.py 
# This module implements a function that reorders elements of density matrices 
# or energy ones according to a permutation "perm". 
# Also it includes unittest on the function with some examples.

import os
import sys
import math
import copy
import unittest

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *

from libra_py import *

def reorder(p,A,E):
    # param[in]     p a list of permutation indices
    # param[in,out] A MATRIX object including density matrix <phi_i(t)|phi_j(t+dt)>
    # param[in,out] E MATRIX object including energies in diagonal elements

    perm = p[:] # create a working list

    for col in xrange(len(perm)):
        indx=perm[col]
        if indx!=col:
            A.swap_cols(indx,col)
            E.swap_cols(indx,col); E.swap_rows(indx,col)

            indx2 = perm[indx]
            if indx2!=col:
                perm[col] = col
                perm[indx2] = indx
            else:
                perm[col] = col
                perm[indx] = indx
            #print "col is %i and indx is %i" % (col,indx)
            #print "temporal perm_c is"; print perm_c
            #print "temporal matrix is"; c.show_matrix()

def _prepare_density_matrices():
    
    # identical
    # | 1 0 0 0 |
    # | 0 1 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    a = CMATRIX(4,4)
    a.set(0,0,1.0+0j);  a.set(1,1,1.0+0j);  a.set(2,2,1.0+0j);  a.set(3,3,1.0+0j); 

    # non-identical (4x4) matrix
    # | 1 0 0 0 |
    # | 0 0 1 0 |
    # | 0 0 0 1 |
    # | 0 1 0 0 |
    b = CMATRIX(4,4)
    b.set(0,0,1.0+0j);  b.set(1,2,1.0+0j);  b.set(2,3,1.0+0j);  b.set(3,1,1.0+0j);

    # non-identical (6x6) matrix
    # | 0 0 1 0 0 0|
    # | 0 0 0 1 0 0|
    # | 1 0 0 0 0 0|
    # | 0 0 0 0 0 1|
    # | 0 0 0 0 1 0|
    # | 0 1 0 0 0 0|
    c = CMATRIX(6,6)
    c.set(0,2,1.0+0j);  c.set(1,3,1.0+0j);  c.set(2,0,1.0+0j);  c.set(3,5,1.0+0j);
    c.set(4,4,1.0+0j);  c.set(5,1,1.0+0j);

    return a,b,c

def _prepare_energy_matrices():

    # | 0  0  0  0 |
    # | 0 11  0  0 |
    # | 0  0 22  0 |
    # | 0  0  0 33 |
    Ea = CMATRIX(4,4)
    Ea.set(0,0,0.0+0j); Ea.set(1,1,11.0+0j); Ea.set(2,2,22.0+0j); Ea.set(3,3,33.0+0j);
    
    Eb = CMATRIX(Ea)

    # | 0  0  0  0  0  0 |
    # | 0 11  0  0  0  0 |
    # | 0  0 22  0  0  0 |
    # | 0  0  0 33  0  0 |
    # | 0  0  0  0 44  0 |
    # | 0  0  0  0  0 55 |

    Ec = CMATRIX(6,6) # eigenenergy matrix (diagonal)
    Ec.set(0,0,0.0+0j); Ec.set(1,1,11.0+0j); Ec.set(2,2,22.0+0j); Ec.set(3,3,33.0+0j);
    Ec.set(4,4,44.0+0j); Ec.set(5,5,55.0+0j);

    return Ea,Eb,Ec

class Test_unavoided(unittest.TestCase):
    def test_reorder(self):

        a,b,c = _prepare_density_matrices()
        Ea,Eb,Ec = _prepare_energy_matrices()

        p4 = range(a.num_of_rows) 
        p6 = range(c.num_of_rows) 

        print "p4 is ",p4
        print "p6 is ",p6

        '''extract indices for reordering '''
        perm_a = unavoided.get_reordering(a)
        print "Input density matrix a"; a.show_matrix()
        print "Permutation a = ", perm_a
        self.assertEqual(perm_a, [0,1,2,3])

        perm_b = unavoided.get_reordering(b)
        print "Input density matrix b"; b.show_matrix()
        print "Input energy matrix Eb" ; Eb.show_matrix()
        print "Permutation b = ", perm_b
        self.assertEqual(perm_b, [0,2,3,1])

        perm_c = unavoided.get_reordering(c)
        print "Input density matrix c"; c.show_matrix()
        print "Input energy matrix Ec"; Ec.show_matrix()
        print "Permutation c = ", perm_c
        self.assertEqual(perm_c, [2,3,0,5,4,1])

        ''' permutation according to "perm" list '''

        if p4 != perm_a:
            reorder(perm_a,a,Ea); print "Matrix a is reordered"
        print "Output density matrix a"; a.show_matrix()
        print "Output energy matrix Ea"; Ea.show_matrix()
        print "Permutation a = ", perm_a

        if p4 != perm_b:
            reorder(perm_b,b,Eb); print "Matrix b is reordered"
        print "Output density matrix b"; b.show_matrix()
        print "Output energy matrix Eb"; Eb.show_matrix()
        print "Permutation b = ", perm_b

        if p6 != perm_c:
            reorder(perm_c,c,Ec); print "Matrix c is reordered"
        print "Output density matrix c"; c.show_matrix()
        print "Output energy matrix Ec"; Ec.show_matrix()
        print "Permutation c = ", perm_c


if __name__=='__main__':
    unittest.main()
    
