def atomic_basis_set(l_gam,params):

    # detect atoms
    ab_start = params["ab_start"]
    ab_end = params["ab_end"]
    l_atom_spec = []
    atom_spec = []
    for i in range(ab_start,ab_end+1):
        spline = l_gam[i].split()
        if len(spline) == 1:
            l_atom_spec.append(i)
            atom_spec.append(spline[0])
    params["atom_spec"] = atom_spec

    # produce the lists for exponent and coefficient
    #o = ["s","p","d"]
    #for s in atom_spec:
    #    for i in range(1,4):
    #        expo_ao = "expo_" + s + str(i)
    #        expo_ao = []
    #        for l in o:                
    #            coef_ao = "coef_" + s + str(i) + l
                #print coef_ao
    #            coef_ao = []

    # extract atomic basis sets
    basis_type = []
    basis_expo = []
    basis_coef = []
    for i in range(0,len(atom_spec)):
        stmp = atom_spec[i]
        type_tmp = []
        expo_tmp = []
        coef_tmp = []
        if i < len(atom_spec) - 1:
            tmp_start = l_atom_spec[i] + 2
            tmp_end = l_atom_spec[i+1] - 2
        else:
            tmp_start = l_atom_spec[i] + 2
            tmp_end = ab_end
        for j in range(tmp_start,tmp_end+1):
            spline = l_gam[j].split()
            coef_tmp1 = []
            if len(spline) > 1:
                type_tmp.append(spline[1])
                expo_tmp.append(float(spline[3]))
                if spline[1] == "L":
                    coef_tmp1.append(float(spline[4]))
                    coef_tmp1.append(float(spline[5]))
                elif spline[1] == "S":
                    coef_tmp1.append(float(spline[4]))
                coef_tmp.append(coef_tmp1)
        basis_type.append(type_tmp)
        basis_expo.append(expo_tmp)
        basis_coef.append(coef_tmp)

    params["basis_type"] = basis_type
    params["basis_expo"] = basis_expo
    params["basis_coef"] = basis_coef

    print "params=",params
    return

def molecular_orbitals(l_gam,params):

    import math  # import should always be in the top of the module

    mo_start = params["mo_start"]
    mo_end = params["mo_end"]
    Ngbf = params["Ngbf"]
    stat_span = 4 + Ngbf
    #num_orb_colm = math.ceil(stat_s)

    mol_ene = []
    l_tmp = []
    mol_coef = []

    for i in range(mo_start,mo_end+1):
        spline = l_gam[i].split()
        di = i - mo_start
        # molecular energy
        if di % stat_span == 1 :
            l_tmp.append(i+2)
            for j in range(0,len(spline)):
                mol_ene.append(spline[j])

    # molecular coefficients
    for i in l_tmp:
        for j in range(i,i+Ngbf):
            coef_tmp = []
            spline = l_gam[j].split()
            if len(mol_coef) < Ngbf:
                for k in range(4,len(spline)):
                    coef_tmp.append(float(spline[k]))
                mol_coef.append(coef_tmp)
            else:
                for k in range(4,len(spline)):
                    mol_coef[j-i].append(float(spline[k]))


    params["mol_ene"] = mol_ene
    params["mol_coef"] = mol_coef
    print "mol_ene=",params["mol_ene"]
    print "mol_coef=",params["mol_coef"]

def coordinates_of_atoms(l_gam,params):

    coor_start = params["coor_start"]
    coor_end = params["coor_end"]

    l_atoms = []
    coor_atoms = []
    for i in range(coor_start,coor_end+1):
        spline = l_gam[i].split()
        # explicit atoms
        l_atoms.append(spline[0])
        # coordinates of atoms
        coor_tmp = []
        for j in range(2,5):
            coor_tmp.append(float(spline[j]))
        coor_atoms.append(coor_tmp)

    params["l_atoms"] = l_atoms
    params["coor_atoms"] = coor_atoms
    print "l_atoms=",params["l_atoms"]
    print "coor_atoms=",params["coor_atoms"]

def gradient(l_gam,params):

    grad_start = params["grad_start"]
    grad_end = params["grad_end"]

    gradient = []
    for i in range(grad_start,grad_end+1):
        spline = l_gam[i].split()
        grad_tmp = []
        for j in range(3,6):
            grad_tmp.append(float(spline[j]))
        gradient.append(grad_tmp)

    params["gradient"] = gradient
    print "gradient=",params["gradient"]
