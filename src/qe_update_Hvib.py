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

import os
import sys
import math
import numpy

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *



#def find_det(A0,A1,n_mo):
def overlap_sd(MO1, MO2):
##
# \brief Compute the overlap of two determinants SD1 and SD2 
# \param[in] MO1 A first set of orbitals - the CMATRIX object of dimensions: Npw x Norb
# \param[in] MO2 A first set of orbitals - the CMATRIX object of dimensions: Npw x Norb
# Here, Norb - the number of MOs in the set, should be the same for each object
# Npw - the number of basis functions (plane waves) used to represent the MO
#
# Each determinant is given in terms of the set of MOs: MO1, MO2, such that:
# SD1 = det(MO1[0], MO1[1], ... , MO1[Norb-1])
# SD2 = det(MO2[0], MO2[1], ... , MO2[Norb-1])
# The sets of MOs are obtained from different calculations (e.g. different times or
# different occupation schemes)
#
# First the overlap matrix is constructed by A0.H()*A1
# Then the determinant of the matrix is calculated using linalg.det module of numpy library
# return determinant, detaa

    if MO1.num_of_rows != MO2.num_of_rows:
        print "The vertical dimensions of the two matrices do not match"
        sys.exit(0)
    
    if MO1.num_of_cols != MO2.num_of_cols:
        print "The horizontal dimensions of the two matrices do not match"
        sys.exit(0)

    Norb = MO1.num_of_cols  # the number of 1-el orbitals in each determinant

    ovlp = CMATRIX(Norb, Norb)
    ovlp = MO1.H()*MO2     # the overlap of the 1-el MOs

    det = []
    for i in xrange(Norb):    
        row_i = []
        for j in xrange(Norb):   
            row_i.append(ovlp.get(i,j))
        det.append(row_i)

    OVLP = numpy.linalg.det(det)/FACTORIAL(Norb)

    return OVLP



#def find_nac(A0,A1,Ao0,Ao1,dt,n_mo):
def compute_nac_sd(MO_old, MO_cur, dt):
##
# \param[in] MO_old A list of MO sets, each defining SD for a given electronic state: at time t-dt
# \param[in] MO_cur A list of MO sets, each defining SD for a given electronic state: at time t
# \param[in] dt nuclear time step
# Returned result: NAC non-adiabatic coupling matrix

# Although the sets MO1 and MO2 are not mutually-orthogonal, so there would be a dS/dt contribution,
# we compute only the Hermitian part, since the non-Hermitian will cancel out in the solving TD-SE

    nstates = len(MO_cur)  # the number of electronic states
    NAC = CMATRIX(nstates,nstates)

    for i in xrange(nstates):
        for j in xrange(nstates):
            s_01 = overlap_sd(MO_old[i],MO_cur[j])   # <SD_i(t-dt)|SD_j(t)>
            s_10 = overlap_sd(MO_cur[i],MO_old[j])   # <SD_i(t)|SD_j(t-dt)>
            NAC.set(i,j,(s_01 - s_10)/(2.0*dt))

    return NAC


#def update_H_vib(D_mat,E_mat,no_ex):
def compute_Hvib(H_el,NAC):
##
# Compute the vibronic Hamiltonian
# \param[in] H_el Electronic Hamiltonian matrix (diagonal) - of type MATRIX: nstates x nstates
# \param[in] NAC Nonadiabatic couplings matrix - of type CMATRIX: nstates x nstates
# Here, nstates - is the number of excited states included into consideration
# Assume atomic units: hbar = 1

    nstates = H_el.num_of_cols

    H_vib = CMATRIX(nstates,nstates)

    for i in xrange(nstates):
            H_vib.set(i,i,H_el.get(i,i))

    for i in xrange(nstates):
        for j in xrange(nstates):
            if j != i:
                H_vib.set(i,j,(-1.0j+0.0)*NAC.get(i,j))

    return H_vib



##======================== All the functions below are not needed ============


def update_E_mat(params):
##
# This is later termed as ham_adi
    no_ex = len(params["excitations"])
    E_mat = MATRIX(no_ex,no_ex)

    for i in xrange(no_ex):
        E_mat.set(i,i,params["epot%i"%i])  # <------- populate H_el matrix directly, not via params[]

    return E_mat


def update_S_matrix(wfc,no_ex,n_mo):
##
# This updates the overlap matrix
    S_mat = CMATRIX(no_ex,no_ex)
    for i in xrange(no_ex):
        for j in xrange(no_ex):
            S_mat.set(i,j,1.0)
            if j != i :
                detaa = find_det(wfc["coeff_%i"%i],wfc["coeff_%i"%j],n_mo)
                S_mat.set(i,j,detaa)

    return S_mat


def update_D_matrix(wfc,dt,no_ex,n_mo):
##
# Updates NACs
    D_mat = CMATRIX(no_ex,no_ex)  #no_ex being number of excited states
    for i in xrange(no_ex):
        for j in xrange(no_ex):
            if j != i :
                nac = find_nac(wfc["coeff_%i"%i],wfc["coeff_%i"%j],wfc["coeff_old_%i"%i],wfc["coeff_old_%i"%j],dt,n_mo)
                if i < j:
                    D_mat.set(i,j,-1.0*nac)
                else:
                    D_mat.set(i,j,nac)

    return D_mat
                                                 

def update_vibronic_hamiltonian(ham_adi, ham_vib, S_mat, wfc, params):
    ##
    # \param[out] ham_adi Electronic (adiabatic) Hamiltonian (MATRIX), diagonal energies
    # \param[out] ham_vib Vibronic Hamiltonian (CMATRIX)
    # \param[out] S_mat Overlap matrix (CMATRIX)
    # \param[in] params  contains the dictionary of the input parameters
    # \param[in] wfc dictionary containing wave function coefficients of MO basis
    #
    # Used in: md.py/run_MD
    # Update vibronic hamiltonian, S_matrix, ham_adi

    dt = params["dt_nucl"]
    no_ex = len(params["excitations"])
    n_mo = params["num_MO"]

    # Update S_mat the overlap matrix (CMATRIX)
    S_mat = update_S_matrix(wfc, no_ex, n_mo)

    # Update NAC matrix  (CMATRIX)
    D_mat = update_D_matrix(wfc,dt,no_ex,n_mo)

    # Update adiabatic hamiltonian (ham_adi)
    ham_adi = update_E_mat(params)

    #Update Vibronic hamiltonian  (CMATRIX)
    ham_vib = update_H_vib(D_mat,ham_adi,no_ex)

