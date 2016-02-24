#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

def ground_state(Nmin,Nmax):
    st = []
    for i in range(Nmin,Nmax+1):
        st.append(i)
        st.append(-i)
    return ["GS",st]

def single_excitations(Nmin,Nmax,HOMO,nspin):
    # form all single excitations from [Nmin,HOMO] to [LUMO,Nmax]
    # where LUMO = HOMO + 1
    # nspin - defines if we want all spin orientations(2) or only one(1)

    st = []  # states
    LUMO = HOMO + 1

    cnt = 0
    for j in range(LUMO,Nmax+1):
        for i in range(Nmin,HOMO+1):
            ex1  = [] # i->j
            ex2 = []
            for v in range(Nmin,HOMO+1):
                if v==i: 
                    ex1.append(j)
                    ex1.append(-i)
                    ex2.append(-j)
                    ex2.append(i)
                else:
                    ex1.append(v)
                    ex1.append(-v)
                    ex2.append(v)
                    ex2.append(-v)

            if nspin>=1:
                st.append(["SE"+str(cnt),ex1])
                cnt = cnt + 1
            if nspin>=2:
                st.append(["SE"+str(cnt),ex2])
                cnt = cnt + 1

    return st
