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

## \file create_pdb_file.py
# This module defines the function which creates a .pdb file to be read in "Load_Molecule" module.
#
# Used in: main.py/nve/nve_setting

import os
import sys
import math

#def read_connect_txt():

def create_pdb_file(l_atoms,coor_atoms):
    ##
    # Finds the keywords and their patterns and extracts the parameters
    # \param[in] l_atoms : a list of atoms.
    # \param[in] coor_atoms : a list of atomic positions.
    # This function outputs a .pdb file to be read in "Load_Molecule" module.
    # 
    # Used in:  main.py/nve/nve_setting

    B_to_A = 0.529177208 # Bohr to Angstrom
    f = open("temp.pdb","w")
    Natoms = len(l_atoms)
    for i in range(0,Natoms):

        f.write("HETATM      %d     %s    %d    %8.5f  %8.5f  %8.5f   %s  \n" % \
                (i+1, l_atoms[i], 0, B_to_A*coor_atoms[i][0], B_to_A*coor_atoms[i][1], B_to_A*coor_atoms[i][2], l_atoms[i]))

#f.write("%8.5f   %8.5f %8.5f %8.5f   %8.5f %8.5f %8.5f \n" % (x, ss, pxpx, dz2dz2, g_ss, g_pxpx, g_dz2dz2) )
