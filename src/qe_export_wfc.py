#*********************************************************************************
#* Copyright (C) 2016 Ekadashi Pradhan, Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

##
# \file text_context.py This file demonstrates to read QE wavfunctions using
# context class of the libra package. We also show how to use the created object 
#
import os
import sys
import math
import numpy

# Fisrt, we add the location of the library to test to the PYTHON path
#print os.getcwd()
cwd = "/projects/academic/alexeyak/ekadashi/libracode-dev/libracode-code/_build"
print "Current working directory", cwd
sys.path.insert(1,cwd+"/src/mmath")
sys.path.insert(1,cwd+"/src/context")


print "\nTest 1: Importing the library and its content"
from libmmath import *
from libcontext import *

def read_qe_wfc(filename, upper_tag, n_el, n_mo):
##
# This function reads an ASCII/XML format file contaiing wavefunction
# and returns the coefficients of the plane wave that constitute the
# wavefunction
#
# \param[in] filename This is the name of the file we will be reading to construct a wavefunction
# \param[in] upper_tag Kpoint.1
# \param[in] n_mo Number of molecular orbital basis used to construct electronic wave function
# \param[in] n_el Number of electrons. e_el/2 is homo index here
# 
#
    homo = n_el/2
    ctx = Context(filename)
    ctx.set_path_separator("/")
    print "path=", ctx.get_path()
    ctx.show_children(upper_tag)  # ("Kpoint.1") #

    ngw = int(float(ctx.get("Info/<xmlattr>/ngw","n")))
    nbnd = int(float(ctx.get("Info/<xmlattr>/nbnd","n")))
    print ngw, nbnd
    coeff = CMATRIX(ngw,n_mo)
    k = 0
    for band in range(homo,homo+n_mo):   # [6,7,8]: #range(1,nbnd+1):

        c = []
        print "band=",band
        all_coeff = ctx.get("Wfc."+str(band), "n").split(',')
        sz = len(all_coeff)

        for i in xrange(sz):
            a = all_coeff[i].split()
            for j in xrange(len(a)):
                c.append(a[j])
        sz = len(c)
        n = sz/2  # this should be equal to ngw
        #n = sz

        for i in xrange(n):
            coeff.set(i, k, float(c[2*i]), float(c[2*i+1]))
        k = k+1
    nbnd = coeff.num_of_cols
    ngw = coeff.num_of_rows
    ovlp  = CMATRIX(nbnd, nbnd)
    ovlp  = coeff.H() * coeff

    for i in xrange(ngw):
        for j in xrange(n_mo):
            coeff.set(i,j,(1/math.sqrt((ovlp.get(j,j)).real))*coeff.get(i,j))

    # This returns normalized coefficients of MO
    return coeff
