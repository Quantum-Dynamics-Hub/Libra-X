def file_names(params):
    l_atoms = params["l_atoms"]
    atom_spec = params["atom_spec"]
    basis_type = params["basis_type"]
    basis_expo = params["basis_expo"]
    basis_coef = params["basis_coef"]
    nGTO = params["nGTO"]
    Ngbf = params["Ngbf"]

    expo_1 = []
    expo_2 = []
    expo_3 = []

    coef_1s = []
    coef_2s = []
    coef_2p = []
    coef_3s = []
    coef_3p = []
    coef_3d = []

    for la in l_atoms: # all atoms
        for j in range(0,len(atom_spec)): # specify the kind of the atom
            if la == atom_spec[j]:
                i = j
        expo_1tmp = []
        expo_2tmp = []
        #expo_3tmp = []
        coef_1stmp = []
        coef_2stmp = []
        coef_2ptmp = []
        #coef_3stmp = []
        #coef_3ptmp = []
        #coef_3dtmp = []
        for j in range(0,len(basis_type[i])): # basis number of atoms
            b_tmp = basis_type[i][j]
            if b_tmp == "S":
                expo_1tmp.append(basis_expo[i][j])
                coef_1stmp.append(basis_coef[i][j][0])
            elif b_tmp == "L":
                expo_2tmp.append(basis_expo[i][j])
                coef_2stmp.append(basis_coef[i][j][0])
                coef_2ptmp.append(basis_coef[i][j][1])
            #elif b_tmp == "D":
            #    expo_3tmp.append(basis_expo[i][j])
            #    coef_3dtmp.append(basis_coef[i][j][0])
            #    print "you inputed illegal character (or D), so exit"
            #    sys.exit
            # f orbitals are not taken into account, so should add them.
        expo_1.append(expo_1tmp)
        expo_2.append(expo_2tmp)
        coef_1s.append(coef_1stmp)
        coef_2s.append(coef_2stmp)
        coef_2p.append(coef_2ptmp)

    print "expo_1=",expo_1
    print "expo_2=",expo_2
    print "coef_1s=",coef_1s
    print "coef_2s=",coef_2s
    print "coef_2p=",coef_2p

    params["expo_1"] = expo_1
    params["expo_2"] = expo_2
    #params["expo_3"] = expo_3
    params["coef_1s"] = coef_1s
    params["coef_2s"] = coef_2s
    params["coef_2p"] = coef_2p

    return

def add_PrimitiveG(params):
    
    l_atoms = params["l_atoms"]
    coor_atoms = params["coor_atoms"]
    expo_1 = params["expo_1"]
    expo_2 = params["expo_2"]
    coef_1s =  params["coef_1s"]
    coef_2s =  params["coef_2s"]
    coef_2p =  params["coef_2p"]
    aoa = params["aoa"]
    priG =  params["priG"]
    vec =  params["vec"]
    orb_name = params["orb_name"]

    ksi = 1.3

    # define atomic orbitals using gaussian basis
    for i in range(0,len(l_atoms)): # all atoms
        vec.x, vec.y, vec.z = coor_atoms[i][0], coor_atoms[i][1], coor_atoms[i][2]
        #vec.x, vec.y, vec.z = 0 ,0, 0
        priG.set_position(vec)
        for j in range(0,len(aoa[i])):
            c_tmp = aoa[i][j]
            o_tmp = orb_name[i][j]
            if o_tmp == "1s":
                expo_tmp = expo_1[i]
                coef_tmp = coef_1s[i]
                nx = 0; ny = 0 ; nz = 0
            elif o_tmp == "2s":
                expo_tmp = expo_2[i]
                coef_tmp = coef_2s[i]
                nx = 0; ny = 0 ; nz = 0
            elif o_tmp[0:2] == "2p":
                expo_tmp = expo_2[i]
                coef_tmp = coef_2p[i]
                if o_tmp[2] == "x":
                    nx = 1 ; ny = 0 ; nz = 0
                elif o_tmp[2] == "y":
                    nx = 0 ; ny = 1 ; nz = 0
                elif o_tmp[2] == "z":
                    nx = 0 ; ny = 0 ; nz = 1
            priG.set_x_exp(nx)
            priG.set_y_exp(ny)
            priG.set_z_exp(nz)
            for k in range(0,len(expo_tmp)):
                priG.set_alpha(ksi*ksi*expo_tmp[k])
                c_tmp.add_primitive(coef_tmp[k],priG)
            aoa[i][j] = c_tmp

    params["aoa"] = aoa

    return

