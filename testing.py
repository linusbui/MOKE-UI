from tensor_utilities import *

A = Symbol('A')
B = Symbol('B')
a = Symbol('alpha')

expr = 3*A*(-2*cos(2*a) + cos(4*a)/2 + 3/2)/8
#pprint((simplify(expr)))
pprint((trigsimp(expr)))