'''
Handles all tasks for computing the epsilon tensor from the MO Tensors K, G H
for different crystal structures and orientations
'''


from tensor_utilities import *


# ==================================================================================================================
# Epsilon Tensor creation
# all work for MO Tensors of arbitrary rotation
# ==================================================================================================================


def eps_from_G(G: Matrix, M: list):
    res = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(0)
        res.append(row)
    
    for i in range(3):
        for j in range(3):
            sum = 0
            for k in range(3):
                for l in range(3):
                    i_, j_ = indices_to_matrix_G([i, j, k, l])
                    sum += G[i_, j_]*M[k]*M[l]
            res[i][j] = custom_simp_1(trigsimp(sum))

    return Matrix(res)


def eps_from_H(H: Matrix, M: list):
    res = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(0)
        res.append(row)

    for i in range(3):
        for j in range(3):
            sum = 0
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        i_, j_ = indices_to_matrix_H([i, j, k, l, m])
                        sum += H[i_, j_]*M[k]*M[l]*M[m]
            res[i][j] = custom_simp_1(trigsimp(sum))
    
    return Matrix(res)


# takes prettified tensor in matrix notation
def eps_from_K(K: Matrix, M: list):
    i_ = [2, 3, 1, 3, 1, 2]
    j_ = [3, 1, 2, 2, 3, 1]

    res = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(0)
        res.append(row)

    for e in range(6):
        sum = 0
        i = i_[e] - 1
        j = j_[e] - 1
        for k in range(3):
            sum += K[e, k] * M[k]
        res[i][j] = nsimplify(sum)
    
    return Matrix(res)


# ==================================================================================================================
# Introduction of magnetization
# ==================================================================================================================


# replaces arbitrary magnetization components M_T, M_L, M_P in given expression with actual values
# Assumes magnetization components have shape M_T, M_L, M_P or M1, M2, M3
def magnetization_dir(expr, M_old: list, M_new: list):
    for k in range(3):
        expr = expr.subs(M_old[k], M_new[k])
    return expr


# introduces magnetization in form of polar coordinates
def magnetization_polar(expr: Matrix, M: list):
    phi = Symbol('phi')
    return magnetization_dir(expr, M, [cos(phi), sin(phi), 0])


def magnetization_spherical(expr: Matrix, M: list):
    phi = Symbol('phi')
    thet = Symbol('theta')
    return magnetization_dir(expr, M, [sin(thet) * cos(phi), sin(thet) * sin(phi), cos(thet)])


# ==================================================================================================================
# Kerr angles
# ==================================================================================================================
def kerr_angles(eps: Matrix, M: list, polar: str):
    eps_d = Symbol('epsilon_d')
    A = [Symbol('A_s'), Symbol('A_p')]
    B = [Symbol('B_s'), Symbol('B_p')]
    if polar == 's':
        return custom_simp_1(A[0]*(eps[1, 0] - custom_simp_2(eps[1, 2]*eps[2, 0], M)/eps_d) + B[0]*eps[2, 0])
    if polar == 'p':
        return custom_simp_1(-A[1]*(eps[0, 1] - custom_simp_2(eps[2, 1]*eps[0, 2], M)/eps_d) + B[1]*eps[0, 2])


# eight directional method for kerr angles
# either for s or p polarized, depending on whats given
def eightdir(kerr: Symbol):
    ml = nsimplify(kerr.subs('phi', 90*pi/180) - kerr.subs('phi', 270*pi/180))/2
    mt = nsimplify(kerr.subs('phi', 0) - kerr.subs('phi', pi))/2
    mlmt = nsimplify(kerr.subs('phi', 45*pi/180) + kerr.subs('phi', 225*pi/180) - kerr.subs('phi', 135*pi/180) - kerr.subs('phi', 315*pi/180))/2
    mtmlsq = nsimplify(kerr.subs('phi', 0) + kerr.subs('phi', pi) - kerr.subs('phi', 90*pi/180) - kerr.subs('phi', 270*pi/180))/2

    return [custom_simp_1(simplify(ml)), custom_simp_1(simplify(mt)), custom_simp_1(simplify(mlmt)), custom_simp_1(simplify(mtmlsq))]