from scipy.linalg import block_diag 
from numpy import array, arange, zeros, ones, identity, vstack, hstack, append, square
from definitions import OptimError, __EPS


## objective is radius
def __objective(A,b):
	return lambda(s): square(A * s - b).min()
	
## objective is radius
def __constraint(A,b):
	return {'type': "ineq", 'fun' : lambda(s): ((A * s - b) < 0).all()}

## Compute the chebyshev ball
# of a given polytope A x <= b
# \param A inequality matrix
# \param b inequality vector
# \return center and radius of the ball
def chebyball(A, b, init_guess = array([0,0,0]), verbose = False):
	objective = __objective(A,b)
	constraints = _constraint(A,b)
	res = minimize(objective, init_guess, constraints=cons, options={'disp': verbose, 'ftol': 1e-03, 'maxiter' : 200})
	print res
	
	

A_t = identity(2)
b_t = ones(2)
