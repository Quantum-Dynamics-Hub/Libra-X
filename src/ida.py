

def ida_py(Coeff, old_st, new_st, E_old, E_new, T, ksi, do_collapse):

    kb = 3.166811429e-6  # Hartree/K
    res = old_st;  C = CMATRIX(Coeff)
    dE = (E_new - E_old);   boltz_f = 1.0

    if dE>0.0:
        argg = dE/(kb*T)
        if argg > 50.0:
            boltz_f = 0.0
        else:
            boltz_f = math.exp(-argg)

        if ksi<boltz_f:
            res = new_st  # accepted hop
            
            # Collapse the wavefunction to the new state 
            if do_collapse:
                C *= 0.0; C.set(new_st, 1.0+0.0j)
        else:
            # Unsuccessful hop - collapse wfc back to the original state^M
            if do_collapse:
                C *= 0.0; C.set(old_st, 1.0+0.0j)
    else:
        res = new_st

    return res , C
