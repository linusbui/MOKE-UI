from sympy import *
from sympy.simplify.fu import TR5, TR7, TR8
from dataclasses import dataclass
from re import *


# ==================================================================================================================
# Tensor Utilities
# ==================================================================================================================


@dataclass
class Rotmat:
    mat: Matrix
    n: int


# for full tensor G
# takes i, j, k, l in {0..2} and returns i, j
def indices_to_matrix_G(ind: list):
    i_ = [0, 1, 2, 1, 2, 0, 2, 0, 1]
    j_ = [0, 1, 2, 2, 0, 1, 1, 2, 0]

    i = 0
    j = 0
    for e in range(9):
        if ind[0] == i_[e] and ind[1] == j_[e]: i = e
    for e in range(9):
        if ind[2] == i_[e] and ind[3] == j_[e]: j = e
    return i, j


# takes list of tensor indices [i, j, k, l, m] in {0..2}
# returns i and j of full tensor H in voigt notation
# uses the same way that H is generated
def indices_to_matrix_H(ind: list):
    i_ = [0, 1, 2, 1, 2, 0, 2, 0, 1]
    j_ = [0, 1, 2, 2, 0, 1, 1, 2, 0]

    i_ind = 0
    j_ind = 0

    for e in range(9):
        if ind[0] == i_[e] and ind[1] == j_[e]: i_ind = e

    found = False
    for i in range(3):
        if found: break
        for j in range(3):
            if found: break
            for k in range(3):
                if ind[2] == i and ind[3] == j and ind[4] == k:
                    found = True
                    break
                j_ind += 1

    return i_ind, j_ind


# maybe rewrite to use [i, j, k] like indice_to_matrix
def magnet_comb_to_index_H(i, j, k):
    if i == j and j == k:
        return i 
    if (i == 1 and j == 2 and k == 2) or (i == 2 and j == 1 and k == 2) or (i == 2 and j == 2 and k == 1):
        return 3
    if (i == 2 and j == 0 and k == 0) or (i == 0 and j == 2 and k == 0) or (i == 0 and j == 0 and k == 2):
        return 4
    if (i == 0 and j == 1 and k == 1) or (i == 1 and j == 0 and k == 1) or (i == 1 and j == 1 and k == 0):
        return 5
    if (i == 2 and j == 1 and k == 1) or (i == 1 and j == 2 and k == 1) or (i == 1 and j == 1 and k == 2):
        return 6
    if (i == 0 and j == 2 and k == 2) or (i == 2 and j == 0 and k == 2) or (i == 2 and j == 2 and k == 0):
        return 7
    if (i == 1 and j == 0 and k == 0) or (i == 0 and j == 1 and k == 0) or (i == 0 and j == 0 and k == 1):
        return 8
    if (i == 0 and j == 1 and k == 2) or (i == 1 and j == 0 and k == 2) or (i == 0 and j == 2 and k == 1) or (i == 2 and j == 1 and k == 0) or (i == 1 and j == 2 and k == 0) or (i == 2 and j == 0 and k == 1):
        return 9
    print(i, j, k)


# seperates given sympy expression into all of its terms
def sep_symbols(expr, symb):
    terms = expr.args
    # WORKAROUND - Assumes no constants present
    if expr.func == Mul:
        return [terms]
    if len(terms) == 0:
        return [(1, expr)]
    res = []
    for i in range(len(terms)):
        exp = terms[i]
        args = exp.args
        if len(args) == 0:
            res.append((1, exp))
        else:
            fac = 1
            #entry = Symbol('Ä')             #needed?
            for t in args:
                #print(t)
                if str(t)[0] == symb:
                    entry = t
                    #print(entry)
                    continue
                else:
                    fac *= t
            #print(fac, entry)
            res.append((fac, entry))
    return res


# ==================================================================================================================
# Custom sympy simplifications
# ==================================================================================================================


# first reduces sin**3, sin**4, sin**5, sin**6 
# then simplifies powers of sin/cos with tr8 and simplifies
def custom_simp_1(expr: Symbol):
    a = Symbol('alpha')
    with evaluate(False):
        expr = expr.subs(sin(a)**2, 1/2*(1 - cos(2*a)))
        expr = expr.subs(cos(a)**2, 1/2*(1 + cos(2*a)))
        expr = expr.subs(sin(a)**3, 1/4*(3*sin(a) - sin(3*a)))
        expr = expr.subs(sin(a)**4, 1/4*(1 - 2*cos(2*a) + 1/2*(1 + cos(4*a))))
        expr = expr.subs(sin(a)**5, 1/16*(10*sin(a) - 5*sin(3*a) + sin(5*a)))
        expr = expr.subs(sin(a)**6, 1/8*(1 - 3*cos(2*a) + 3/2*(1 + cos(4*a)) - 1/4*(3*cos(2*a) + cos(6*a))))
        expr = expr.subs(cos(a)**4, 1/4*(1 + 2*cos(2*a) + 1/2*(1 + cos(4*a))))
        expr = expr.subs(cos(a)**5, 1/16*(10*cos(a) + 5*cos(3*a) + cos(5*a)))
    return (TR8(nsimplify(expr)))


# used to neglect terms of order 4 or higher in kerr angle cross-term
# M is list op Symbols
# check for addition?
def custom_simp_2(expr: Symbol, M: list):
    # string representation of M
    Mstr = [str(m) for m in M]

    resexpr = 0

    # rewrite into polynomial in M
    expr = expr.expand()
    parts = expr.args
    for p in parts:
        pstr = str(p)
        mcount = 0
        for m in Mstr:
            inds = [match.start() for match in finditer(m, pstr)]             # indices of all occurences of m in pstr
            #print(m, inds)
            offset = len(m) - 1                                               # offset chosen such that ind + offset is last char of M occurence
            # check for exponent
            for ind in inds:
                ind_off = ind + offset
                mcount += 1
                if ind_off + 3 < len(pstr) and pstr[ind_off + 1] == '*' and pstr[ind_off + 2] == '*':
                    mcount += int(pstr[ind_off + 3]) - 1
        if mcount > 3:
            continue
        resexpr += p

    # nsimplify?
    return nsimplify(resexpr)


# ==================================================================================================================
# Visualization
# ==================================================================================================================


# prints given epsilon tensor
def print_epsilon(eps: Matrix):
    for i in range(3):
        for j in range(3):
            symb = Symbol('epsilon' + str(i+1) + str(j+1))
            pprint(symb)
            pprint(eps[i, j], num_columns = 1000)
            #print(str(eps[i, j]))
            print('------------------------------')


# replace all free variables with complex numbers
def intro_complex(expr: Symbol):
    symbols = expr.free_symbols
    for symb in symbols:
        symbstr = str(symbols)
        expr = expr.subs(symb, Symbol(symbstr + '_R') + I*Symbol(symbstr + '_I'))
    return expr