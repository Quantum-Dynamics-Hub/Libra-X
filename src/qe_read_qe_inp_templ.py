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

## \file read_qe_inp_templ.py This file implements functions for reading the template
# input file for QE 

#
def read_qe_inp_templ(inp_filename):
##
# Reading and storing template input file for QE calculations. The input file is essentially a
# normal input file, but we store only the constant section (control option), not the 
# coordinates. The latter will be updated at each iteration using the propagated objects
#
# \param[in] inp_filename The name of the initial input file, which will serve as a template
#

    f = open(inp_filename,"r")
    templ = f.readlines()
    f.close()

    for a in templ:
        s = a.split()
        if len(s) > 0 and s[0] == "celldm(1)" and s[1] == "=":
            sa = s[2].split(',')
            cell_dm = float(sa[0])
            break

    # Find the line preceeding the actual atomic coordinates
    for a in templ:
        s = a.split()
        if len(s) > 0 and s[0] == "ATOMIC_POSITIONS":
            ikeep = templ.index(a)
            break

    N = len(templ)
    # Blank space for the atomic positions
    templ[ikeep+1:N] = []
    for i in xrange(ikeep+1):
        print templ[i]
    

    return  templ

