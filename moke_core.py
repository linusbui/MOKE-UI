'''
Brings all modules together and creates MO / eps tensor from given crystal structures and orientations
'''


from tensor_utilities import *
from eps_tensor import *
from mo_tensor import *


# ==================================================================================================================
# Full Tensor creation routines
# ==================================================================================================================
def core_routine(crystal: str, orient: str, magn_rep: str, polar: str):
    rotmats = string_to_rotmat(crystal)
    M = magnetization_rep(magn_rep)

    K_basic, K_rot = mo_tensor_K(rotmats)
    G_basic, G_rot = mo_tensor_G(rotmats, orient, False, M)
    H_basic, H_rot = mo_tensor_H(rotmats, orient, M)

    lat_k = latex(K_basic) + latex_magn_k(M)
    lat_g = latex(G_basic) + latex_magn_g(M)
    lat_h = latex(H_basic) + latex_magn_h(M)

    pprint(G_rot.shape)

    full_eps = eps_from_K(K_rot, M) + eps_from_G(G_rot, M) + eps_from_H(H_rot, M)

    kerr = kerr_angles(full_eps, M, polar)

    eight = eightdir(magnetization_polar(kerr, M))


    latex_eps = latex_eps_comps(full_eps)
    lat_kerr = latex_kerr_angles(kerr)
    lat_eight = latex_eight_dir(eight)

    #pprint(kerr, num_columns = 300)
    #print('====================')
    #pprint(intro_complex(kerr), num_columns = 300)

    return [lat_k, lat_g, lat_h], latex_eps, lat_kerr, lat_eight


# returns latex mo tensors for given settings
def core_mo(crystal: str, orient: str, magn_rep: str):
    rotmats = string_to_rotmat(crystal)
    M = magnetization_rep(magn_rep)

    K_show, K_rot = mo_tensor_K(rotmats)
    G_show, G_rot = mo_tensor_G(rotmats, orient, False, M)
    H_show, H_rot = mo_tensor_H(rotmats, orient, M)

    return [str(K_show), str(G_show), str(H_show)], [str(K_rot), str(G_rot), str(H_rot)]


# returns latex eps tensor
def core_eps(mo: list, magn_rep: str):
    M = magnetization_rep(magn_rep)
    full_eps = eps_from_K(mo[0], M) + eps_from_G(mo[1], M) + eps_from_H(mo[2], M)
    return full_eps


def core_kerr(full_eps: Matrix, magn_rep: str, polar: str):
    M = magnetization_rep(magn_rep)
    kerr = kerr_angles(full_eps, M, polar)
    return kerr


def core_eight(kerr: Symbol, magn_rep: str):
    M = magnetization_rep(magn_rep)
    eight = eightdir(magnetization_polar(kerr, M))
    return eight


def mo_tensor_K(crystal: list):
    tens, symm_eqs, symb, ident = general_onsager_K()
    simpl = simplify_tensor(crystal, tens, symb, symm_eqs, ident)
    fix = fix_K(simpl)
    subst, symbols = substitute(fix, ident)

    return fix, subst


def mo_tensor_G(crystal: list, orient: str, use_dgamma: bool, M: list):
    tens, symm_eqs, symb, ident = general_onsager_G()
    simpl = simplify_tensor(crystal, tens, symb, symm_eqs, ident)
    subst, symbols = substitute(simpl, ident)
    subst_rew = rewrite_G(subst, use_dgamma)
    subst_rot = rotate_G(subst_rew, orient)

    return prettify_G(subst), subst_rot


def mo_tensor_H(crystal: list, orient: str, M: list):
    tens, symm_eqs, symb, ident = general_onsager_H()
    simpl = simplify_tensor(crystal, tens, symb, symm_eqs, ident)
    fix = fix_H(simpl)
    subst, symbols = substitute(fix, ident)
    subst_rew = rewrite_H(subst)
    #pprint(prettify_H(subst_rew), num_columns = 300)
    subst_rot = rotate_H(subst_rew, orient)

    return prettify_H(subst), subst_rot

# ==================================================================================================================
# Rotation Matrices for different crystal structures
# ==================================================================================================================

c2z_m = Matrix([[-1, 0, 0], 
              [0, -1, 0],
              [0, 0, 1]])
C2z = Rotmat(c2z_m, 2)

c2x_m = Matrix([[1, 0, 0],
              [0, -1, 0],
              [0, 0, -1]])
C2x = Rotmat(c2x_m, 2)

c2y_m = Matrix([[-1, 0, 0],
                [0, 1, 0],
                [0, 0, -1]])
C2y = Rotmat(c2y_m, 2)

c3z_m = Matrix([[-1/2, -sqrt(3)/2, 0],
              [sqrt(3)/2, -1/2, 0],
              [0, 0, 1]])
C3z = Rotmat(c3z_m, 3)

c2a_m = Matrix([[0, 1, 0],
                [1, 0, 0],
                [0, 0, -1]])
C2a = Rotmat(c2a_m, 2)

c3_m = Matrix([[0, 0, 1],
               [1, 0, 0],
               [0, 1, 0]])
C3 = Rotmat(c3_m, 3)

c4z_m = Matrix([[0, 1, 0],
              [-1, 0, 0],
              [0, 0, 1]])
C4z = Rotmat(c4z_m, 4)

c6z_m = Matrix([[1/2, sqrt(3)/2, 0],
              [-sqrt(3)/2, 1/2, 0],
              [0, 0, 1]])
C6z = Rotmat(c6z_m, 6)


monoclinic = [C2z]
orthorhombic = [C2z, C2x, C2y]
tetragonal_1 = [C4z]
tetragonal_2 = [C4z, C2x]
trigonal_1 = [C3z]
trigonal_2 = [C3z, C2x]
hexagonal_1 = [C6z]
hexagonal_2 = [C6z, C2x]
cubic = [C2x, C2y, C2z, C3]

# does that even make a difference?
# C2x?
tet_g_1 = [C2z, C4z]
tet_g_2 = [C2z, C4z, C2a]

hex_g_1 = [C3z, C2z]
hex_g_2 = [C3z, C2z, C2x]

cubic_g_2 = [C2x, C2y, C2z, C3, C2a]


#M = [Symbol('M_T'), Symbol('M_L'), Symbol('M_P')]


# ==================================================================================================================
# UI Compatibility
# ==================================================================================================================


# gets string representations of crystal structure from UI and returns corresponding rotmats
def string_to_rotmat(crystal: str):
    rotmats = []
    match crystal:
        case 'monoclinic':
            rotmats = [C2z]
        case 'orthorhombic':
            rotmats = [C2z, C2x, C2y]
        case 'tetragonal':
            rotmats = [C2z, C4z]
        case 'tetragonal2':
            rotmats = [C2z, C4z, C2a]
        case 'trigonal':
            rotmats = [C3z]
        case 'trigonal2':
            rotmats = [C3z, C2x]
        case 'hexagonal':
            rotmats = [C3z, C2z]
        case 'hexagonal2':
            rotmats = [C3z, C2z, C2x]
        case 'cubic':
            rotmats = [C2x, C2y, C2z, C3]
        case 'cubic2':
            rotmats = [C2x, C2y, C2z, C3, C2a]
    
    return rotmats


def magnetization_rep(magn_rep: str):
    M = []
    match magn_rep:
        case 'letters':
            M = [Symbol('M_T'), Symbol('M_L'), Symbol('M_P')]
        case 'numbers':
            M = [Symbol('M1'), Symbol('M2'), Symbol('M3')]
    
    return Matrix(M)


def magnetization_type(magn_type: str):
    M = []
    match magn_type:
        case 'spherical':
            t = Symbol('theta')
            p = Symbol('phi')
            M = [sin(t)*cos(p), sin(p)*sin(p), cos(p)]
        case 'polar':
            p = Symbol('phi')
            M = [cos(p), sin(p), 0]
    
    return Matrix(M)


# returns latex representation of magnetization vector for printing K
def latex_magn_k(M):
    m = Matrix(M)
    return latex(m)


def latex_magn_g(M):
    i_ = [0, 1, 2, 1, 2, 0]
    j_ = [0, 1, 2, 2, 0, 1]

    res = []
    for e in range(6):
        res.append(M[i_[e]]*M[j_[e]])
    m = Matrix(res)
    return latex(m)


def latex_magn_h(M):
    i_ = [0, 1, 2, 1, 2, 0, 2, 0, 1]
    j_ = [0, 1, 2, 2, 0, 1, 1, 2, 0]

    res = []
    for e in range(9):
        res.append(M[i_[e]]*M[j_[e]]**2)
    res.append(M[0]*M[1]*M[2])
    m = Matrix(res)
    return latex(m)


def latex_eps_comps(eps: Matrix):
    res = []
    for i in range(3):
        res.append(latex(Symbol('epsilon' + str(i+1) + str(i+1))) + ' = ' + latex(eps[i, i]))
    res.append(latex(Symbol('epsilon' + str(1) + str(2))) + ' = ' + latex(eps[0, 1]))
    res.append(latex(Symbol('epsilon' + str(2) + str(1))) + ' = ' + latex(eps[1, 0]))

    res.append(latex(Symbol('epsilon' + str(1) + str(3))) + ' = ' + latex(eps[0, 2]))
    res.append(latex(Symbol('epsilon' + str(3) + str(1))) + ' = ' + latex(eps[2, 0]))

    res.append(latex(Symbol('epsilon' + str(2) + str(3))) + ' = ' + latex(eps[1, 2]))
    res.append(latex(Symbol('epsilon' + str(3) + str(2))) + ' = ' + latex(eps[1, 2]))
    res.append(latex(eps[2, 1]))

    return res


def latex_kerr_angles(kerr):
    return latex(Symbol('Phi_s')) + '=' + latex(kerr)

def latex_eight_dir(eightdir: list):
    return [r'\Phi_{M_L, M_L^3}=' + latex(eightdir[0]),  r'\Phi_{M_T^3}=' + latex(eightdir[1]),
            r'\Phi_{M_LM_T}=' + latex(eightdir[2]),
            r'\Phi_{M_T^2-M_L^2}=' + latex(eightdir[3])]