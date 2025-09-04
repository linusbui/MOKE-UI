'''
Handles all basic for computing the MO Tensors K, G and H 
for different crystals and orientations
'''


from tensor_utilities import *


# ==================================================================================================================
# Basic Tensor creation
# ==================================================================================================================

# creates MO Tensor K in voigt Matrix notation, all onsager relations, all sympy symbols from the tensor and the name
def general_onsager_K():
    i_ = [1, 2, 3, 2, 3, 1, 3, 1, 2]
    j_ = [1, 2, 3, 3, 1, 2, 2, 3, 1]
    
    symbols = []
    eqs = []
    resmat = []

    for e in range(9):
        row = []

        i = i_[e]
        j = j_[e]
        for k in range(3):
            symb = Symbol('K' + str(i) + str(j) + str(k+1))
            symbols.append(symb)
            row.append(symb)
            if i == j:
                eqs.append(symb)
            if i > j:
                eqs.append(symb + Symbol('K' + str(j) + str(i) + str(k+1)))
            eqs.append(symb + Symbol('K' + str(j) + str(i) + str(k+1)))
        resmat.append(row)

    return Matrix(resmat), eqs, symbols, 'K'


# creates MO Tensor G in voigt Matrix notation, all onsager relations, all sympy symbols from the tensor and the name
def general_onsager_G():
    i_ = [1, 2, 3, 2, 3, 1, 3, 1, 2]
    j_ = [1, 2, 3, 3, 1, 2, 2, 3, 1]

    symbols = []
    eqs = []
    resmat = []

    for e in range(9):
        row = []

        i = i_[e]
        j = j_[e]
        for e_ in range(9):
            k = i_[e_]
            l = j_[e_]
            symb = Symbol('G' + str(i) + str(j) + str(k) + str(l))
            symbols.append(symb)
            if i == j:
                eqs.append(symb - Symbol('G' + str(i) + str(j) + str(l) + str(k)))
            if i > j:
                eqs.append(symb - Symbol('G' + str(j) + str(i) + str(k) + str(l)))
                eqs.append(Symbol('G' + str(j) + str(i) + str(k) + str(l)) - Symbol('G' + str(i) + str(j) + str(l) + str(k)))
                eqs.append(Symbol('G' + str(i) + str(j) + str(l) + str(k)) - Symbol('G' + str(j) + str(i) + str(k) + str(l)))
            row.append(symb)

        resmat.append(row)

    return Matrix(resmat), eqs, symbols, 'G'


# creates MO Tensor H in voigt Matrix notation, all onsager relations, all sympy symbols from the tensor and the name
def general_onsager_H():
    i_ = [1, 2, 3, 2, 3, 1, 3, 1, 2]
    j_ = [1, 2, 3, 3, 1, 2, 2, 3, 1]

    symbols = []
    eqs = []
    resmat = []

    for e in range(9):
        row = []

        i = i_[e]
        j = j_[e]
        for k in range(3):
            for l in range(3):
                for m in range(3):
                    symb = Symbol('H' + str(i) + str(j) + str(k+1) + str(l+1) + str(m+1))
                    symbols.append(symb)
                    row.append(symb)
                    if i == j:
                        eqs.append(symb)
                    if i > j:
                        # positives
                        eqs.append(symb - Symbol('H' + str(i) + str(j) + str(k+1) + str(m+1) + str(l+1)))

                        eqs.append(Symbol('H' + str(i) + str(j) + str(k+1) + str(m+1) + str(l+1)) - 
                                   Symbol('H' + str(i) + str(j) + str(l+1) + str(k+1) + str(m+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(l+1) + str(k+1) + str(m+1)) -
                                   Symbol('H' + str(i) + str(j) + str(l+1) + str(m+1) + str(k+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(l+1) + str(m+1) + str(k+1)) -
                                   Symbol('H' + str(i) + str(j) + str(m+1) + str(k+1) + str(l+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(m+1) + str(k+1) + str(l+1)) - 
                                   Symbol('H' + str(i) + str(j) + str(m+1) + str(l+1) + str(k+1)))
                        # negatives
                        eqs.append(symb + Symbol('H' + str(j) + str(i) + str(k+1) + str(l+1) + str(m+1)))

                        eqs.append(Symbol('H' + str(i) + str(j) + str(k+1) + str(m+1) + str(l+1)) + 
                                   Symbol('H' + str(j) + str(i) + str(k+1) + str(m+1) + str(l+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(l+1) + str(k+1) + str(m+1)) + 
                                   Symbol('H' + str(j) + str(i) + str(l+1) + str(k+1) + str(m+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(l+1) + str(m+1) + str(k+1)) + 
                                   Symbol('H' + str(j) + str(i) + str(l+1) + str(m+1) + str(k+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(m+1) + str(k+1) + str(l+1)) + 
                                   Symbol('H' + str(j) + str(i) + str(m+1) + str(k+1) + str(l+1)))
                        
                        eqs.append(Symbol('H' + str(i) + str(j) + str(m+1) + str(l+1) + str(k+1)) + 
                                   Symbol('H' + str(j) + str(i) + str(m+1) + str(l+1) + str(k+1)))
        resmat.append(row)
    
    return Matrix(resmat), eqs, symbols, 'H'


# ==================================================================================================================
# Tensor transformations for presentation
# Used to make nice printable matrix representations
# ==================================================================================================================

# Takes unreduced MO Tensor K and removes trivial elements (due to onsager relation)
# fixes sign conventions aswell
# returns reduced Tensor K
def fix_K(tens: Matrix):
    res = []
    for i in range(6):
        row = []
        for j in range(3):
            row.append(0)
        res.append(row)
    
    for i in range(6):
        for j in range(3):
            symb = tens[i + 3, j]
            symb_str = str(symb)
            if symb_str[0] == '-':
                res[i][j] = Symbol(symb_str[1:])
                res[i+3][j] = -Symbol(symb_str[1:])
    return Matrix(res)


# fixes sign convention for H tensor -> just flips all signs according to onsager relation
# could use same code for K tensor
# onsager relation in voigt notation: H[i, j] = -H[i+3, j]
def fix_H(tens: Matrix):
    height, width = tens.shape
    res_tens = tens.copy()
    
    for i in range(height - 3):
        for j in range(width):
                upper_symb_str = str(tens[i, j])
                if upper_symb_str[0] == '-':
                    res_tens[i, j] = -res_tens[i, j]
                    res_tens[i+3, j] = -res_tens[i, j]

    return res_tens


# Takes unreduced MO Tensor G and simplifies trivial elements (due to onsager relation)
# returns reduced Tensor G
def prettify_G(tens: Matrix):
    res = []
    for i in range(6):
        row = []
        for j in range(6):
            row.append(0)
        res.append(row)

    for i in range(6):
        for j in range(6):
            symb = tens[i, j]
            symb_ = tens[i, j + 3]
            if j > 2:
                if not symb == symb_:
                    print('WIERD')
                    continue
                res[i][j] = 2*symb
            else:
                res[i][j] = symb
    
    return Matrix(res)


# takes full H tensor and reduces it to entries with unique magnetization combination
def prettify_H(tens: Matrix):
    res = []
    for i in range(6):
        row = []
        for j in range(10):
            row.append(0)
        res.append(row)
    
    for i in range(3):
        for j in range(3):
            # entries where i == j are trivial: exclude them
            # cuts off first three rows of tensor
            if i == j:
                continue
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        i_, j_ = indices_to_matrix_H([i, j, k, l, m])
                        k_ = magnet_comb_to_index_H(k, l, m)
                        entr = tens[i_, j_]
                        res[i_ - 3][k_] += entr
    return Matrix(res)


# takes tensor in voigt notation and name of the tensor (K, G or H)
# assumes that only simple tensor entries are present, no additions
# returns Martix with ennumerated entries, and set of those 
def substitute(tensor: Matrix, ident: str):
    height, width = tensor.shape

    # collect distinct entries
    sep_strs = []
    for i in range(height):
        for j in range(width):
            sep = sep_symbols(tensor[i, j], ident)[0]
            entr_str = str(sep[1])
            if not entr_str in sep_strs and not entr_str == '0':
                sep_strs.append(entr_str)

    #print(sep_strs)
    symbs = set()
    res = []
    for i in range(height):
        row = []
        for j in range(width):
            sep = sep_symbols(tensor[i, j], ident)[0]
            if sep[1] == 0:
                row.append(0)
            else:
                ind = sep_strs.index(str(sep[1]))
                symb = Symbol(ident + str(ind + 1))
                symbs.add(symb)
                row.append(sep[0] * symb)
        res.append(row)
    return Matrix(res), symbs


# ==================================================================================================================
# Tensor transformations for calculations
# Used to fix sign conventions or introduce new parameters for simplifications
# ==================================================================================================================


# rewrites G using delta G and delta gamma
# usage of delta gamma is optional and disabled by default, can be enabled with use_dgamma = True
# FOR CUBIC / TETRAGONAL G12 = G21, MAYBE NOT FOR SOME OTHERS
def rewrite_G(G: Matrix, use_dgamma: bool = False):
    G11 = G[0, 0]
    G12 = G[0, 1]
    G21 = G[0, 2]
    G44 = G[8, 8]

    d_gamma = 0
    if not G12 == G21 and use_dgamma:
        d_gamma = Symbol('Delta_Gamma')

    G12_new = G11 - 2*G44 - Symbol('Delta_G') + d_gamma
    G21_new = G11 - 2*G44 - Symbol('Delta_G') - d_gamma

    G_new = G.copy()
    G_new = G_new.subs(G12, G12_new)
    G_new = G_new.subs(G21, G21_new)
    return G_new


# introduce delta gamma?
# ONLY WORKS FOR CUBIC FOR NOW, PROBABLY
def rewrite_H(H: Matrix):
    # convention for cubic crystals
    #pprint(H, num_columns = 300)
    H123 = H[3, 0]
    H125 = H[3, 4]

    d_H = Symbol('Delta_H')
    H_temp = Symbol('temp')

    H123_new = d_H + 3*H_temp
    #pprint(H123_new)
    H125_new = (H123 - d_H)/3

    H_new = H.copy()
    H_new = H_new.subs(H123, H123_new)
    H_new = H_new.subs(H125, H125_new)
    H_new = H_new.subs(H_temp, H125)
    #pprint(H_new, num_columns = 300)
    return H_new


# ==================================================================================================================
# Tensor simplification
# Used to calculate reduced Tensors based on crystal structure symmetries
# ==================================================================================================================


# perm indices from 0-2
def rotate_entry(perm: list, factors: list, expr: Symbol, ident: str):
    decomp = sep_symbols(expr, ident)

    res = 0
    for term in decomp:
        term_fac = term[0]
        symb = term[1]

        symstr = str(symb)
        symnums = symstr[1:]

        ind = []
        for c in symnums:
            ind.append(int(c)-1)

        # not nice
        if ident == 'K':
            i_ = perm[ind[0]]
            j_ = perm[ind[1]]
            k_ = perm[ind[2]]
            f_i = factors[ind[0]]
            f_j = factors[ind[1]]
            f_k = factors[ind[2]]

            new_term = 0
            for i in range(len(i_)):
                for j in range(len(j_)):
                    for k in range(len(k_)):
                        fac = f_i[i]*f_j[j]*f_k[k]
                        new_num_str = str(i_[i]+1) + str(j_[j]+1) + str(k_[k]+1)
                        symb = Symbol(ident + new_num_str)
                        #print(fac * symb, fac, symb)
                        new_term += fac * Symbol(ident + new_num_str)
            res += term_fac*new_term

        if ident == 'G':
            i_ = perm[ind[0]]
            j_ = perm[ind[1]]
            k_ = perm[ind[2]]
            l_ = perm[ind[3]]
            f_i = factors[ind[0]]
            f_j = factors[ind[1]]
            f_k = factors[ind[2]]
            f_l = factors[ind[3]]

            new_term = 0
            for i in range(len(i_)):
                for j in range(len(j_)):
                    for k in range(len(k_)):
                        for l in range(len(l_)):
                            fac = f_i[i]*f_j[j]*f_k[k]*f_l[l]
                            new_num_str = str(i_[i]+1) + str(j_[j]+1) + str(k_[k]+1) + str(l_[l]+1)
                            symb = Symbol(ident + new_num_str)
                            #print(fac * symb, fac, symb)
                            new_term += fac * Symbol(ident + new_num_str)
            res += term_fac*new_term
        
        if ident == 'H':
            i_ = perm[ind[0]]
            j_ = perm[ind[1]]
            k_ = perm[ind[2]]
            l_ = perm[ind[3]]
            m_ = perm[ind[4]]
            f_i = factors[ind[0]]
            f_j = factors[ind[1]]
            f_k = factors[ind[2]]
            f_l = factors[ind[3]]
            f_m = factors[ind[4]]

            new_term = 0
            for i in range(len(i_)):
                for j in range(len(j_)):
                    for k in range(len(k_)):
                        for l in range(len(l_)):
                            for m in range(len(m_)):
                                fac = f_i[i]*f_j[j]*f_k[k]*f_l[l]*f_m[m]
                                new_num_str = str(i_[i]+1) + str(j_[j]+1) + str(k_[k]+1) + str(l_[l]+1) + (str(m_[m]+1))
                                symb = Symbol(ident + new_num_str)
                                #print(fac * symb, fac, symb)
                                new_term += fac * Symbol(ident + new_num_str)
            res += term_fac*new_term
    return res


# tensor is in matrix notation
# returns all rotated tensors
def rotate_tensor(rot: Rotmat, tensor: Matrix, ident: str):
    res = [tensor]

    for _ in range(rot.n - 1):
    # initialize array for rotated tensor
    # maybe np.zeroes?
        last_rotated = res[len(res) - 1]
        tens_entr = []
        height, width = last_rotated.shape
        for i in range(height):
            row = []
            for j in range(width):
                row.append(0)
            tens_entr.append(row)

        # create index transform. and factors from rotmat
        rotmat = rot.mat
        new_indices = []
        factors = []
        r_height, r_width = rotmat.shape
        for i in range(r_height):
            ind = []
            fac = []
            for j in range(r_width):
                entry = rotmat[i, j]
                if entry == 0:
                    continue
                ind.append(j)
                fac.append(entry)
            new_indices.append(ind)
            factors.append(fac)
        
        #print(new_indices, factors)

        for i in range(height):
            for j in range(width):
                symb = last_rotated[i, j]
                new_symb = rotate_entry(new_indices, factors, symb, ident)
                tens_entr[i][j] = expand(new_symb)

        res.append(Matrix(tens_entr))

    return res[1:]


# Assume given tensor is in 'base form' and entries only have one symbol
def simplify_tensor(rotmats: list, tensor: Matrix, symbols: list, onsager: list, ident: str):
    rotations = [tensor]
    for rotmat in rotmats:
        rot = rotate_tensor(rotmat, tensor, ident)
        for r in rot:
            rotations.append(r)
    
    #print(rotations)
    # construct equations from rotated tensor entries
    eqs = []
    for i in range(len(rotations)-1):
        t = rotations[i]
        t_ = rotations[i+1]
        diff = t - t_
        height, width = diff.shape
        for i in range(height):
            for j in range(width):
                eqs.append(diff[i, j])
    solv = linsolve(eqs + onsager, symbols).args[0]
    #print(solv)

    res = []
    ind = 0
    height, width = tensor.shape
    for i in range(height):
        row = []
        for j in range(width):
            row.append(nsimplify(solv[ind]))
            ind += 1
        res.append(row)

    return Matrix(res)


# ==================================================================================================================
# Tensor rotation
# Used to calculate reduced Tensors based on crystal structure symmetries
# ==================================================================================================================


# takes full tensor G and crystal orientation (e.g '001')
# returns rotated tensor
def rotate_G(G: Matrix, dir: str):
    a = Symbol('alpha')
    R = Matrix([[0, 0, 0], 
                [0, 0, 0], 
                [0, 0, 0]])
    Rz = Matrix([[cos(a), sin(a), 0], 
                 [-sin(a), cos(a), 0], 
                 [0, 0, 1]])

    if not dir == '001' and not dir == '011' and not dir == '111':
        return
    if dir == '001':
        R = Rz
    elif dir == '011':
        R = Rz * Matrix([[1, 0, 0], 
                         [0, sqrt(2)/2, -sqrt(2)/2], 
                         [0, sqrt(2)/2, sqrt(2)/2]])
    elif dir == '111':
        R = Rz * Matrix([[1/2 + 1/(2*sqrt(3)), -1/2 + 1/(2*sqrt(3)), -1/sqrt(3)], 
                         [-1/2 + 1/(2*sqrt(3)), 1/2 + 1/(2*sqrt(3)), -1/(sqrt(3))], 
                         [1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]])
    
    G_res = G.copy()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    sum = 0
                    for n in range(3):
                        for o in range(3):
                            for p in range(3):
                                for q in range(3):
                                    n_, o_ = indices_to_matrix_G([n, o, p, q])
                                    sum += nsimplify(R[i, n]*R[j, o]*R[k, p]*R[l, q]*G[n_, o_])
                    i_, j_ = indices_to_matrix_G([i, j, k, l])
                    G_res[i_, j_] = custom_simp_1(trigsimp(sum))

    return G_res


# takes full H tensor and crystal orientation (e.g '001')
# returns rotated tensor
def rotate_H(H: Matrix, dir: str):
    a = Symbol('alpha')
    R = Matrix([[0, 0, 0], 
                [0, 0, 0], 
                [0, 0, 0]])
    Rz = Matrix([[cos(a), sin(a), 0], 
                 [-sin(a), cos(a), 0], 
                 [0, 0, 1]])

    if not dir == '001' and not dir == '011' and not dir == '111':
        return
    if dir == '001':
        R = Rz
    elif dir == '011':
        R = Rz * Matrix([[1, 0, 0], 
                         [0, sqrt(2)/2, -sqrt(2)/2], 
                         [0, sqrt(2)/2, sqrt(2)/2]])
    elif dir == '111':
        R = Rz * Matrix([[1/2 + 1/(2*sqrt(3)), -1/2 + 1/(2*sqrt(3)), -1/sqrt(3)], 
                         [-1/2 + 1/(2*sqrt(3)), 1/2 + 1/(2*sqrt(3)), -1/(sqrt(3))], 
                         [1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]])
    
    H_res = H.copy()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        sum = 0
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
                                    for q in range(3):
                                        for r in range(3):
                                            n_, o_ = indices_to_matrix_H([n, o, p, q, r])
                                            sum += nsimplify(R[i, n]*R[j, o]*R[k, p]*R[l, q]*R[m, r]*H[n_, o_])
                        i_, j_ = indices_to_matrix_H([i, j, k, l, m])
                        H_res[i_, j_] = custom_simp_1(trigsimp(sum))
    
    return H_res