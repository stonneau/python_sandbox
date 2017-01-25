from scipy.linalg import block_diag, norm
from numpy import array, arange, zeros, ones, identity, vstack, hstack, append, sqrt, square
from definitions import __EPS
import math

## objective is  similar to constraint.
# costs functions are chosen from a objective factory, with
# user-defined weights.
# the factory automatically builds the cost function
# and returns it

def __normsquared(v):
	#~ return norm(v)
	return square(norm(v))

def __sum_over(var_name, weight):
	return lambda variables: weight * reduce(lambda a, b :a + b, [ __normsquared(s) for s in variables[var_name] ])

## ("dL")
#  Minimize angular momentum variation
#  \return  cone matrices for each phase
def min_dL(param, weight):
	return __sum_over('dL', weight)
	
## ("ddL")
#  Minimize angular momentum variation
#  \return  cone matrices for each phase
def min_ddL(param, weight):
	return __sum_over('ddL', weight)
	
## ("ddc")
#  Minimize angular momentum variation
#  \return  cone matrices for each phase
def min_ddc(param, weight):
	return __sum_over('ddc', weight)
	
## ("dddc")
#  Minimize jerk
#  \return  cone matrices for each phase
def min_dddc(param, weight):
	return __sum_over('dddc', weight)
	
def __line(line_points, weight):
	def fun(variables):		
		c = variables["c"]
		return weight * reduce(lambda a, b :a + b, [ __normsquared(line_points[i]-c[i]) for i in range(len(c)) ])
	return fun
	
## ("c")
#  Make the COM as much as possible follow a straight line
#  \return  cone matrices for each phase
def line(params, weight):
	x_init = array(params["x_init"][0:3]);
	x_end  = array(params["x_end"][0:3]);
	v_dir  = x_end - x_init; 
	#v_dir = v_dir / norm(v_dir)
	#director vector is scaled such that last dt equals to xend
	t_max =  params["t_init_phases"][-1]; dt = params["dt"]
	v_dir = v_dir /	t_max
	line_points = [x_init + v_dir * t_i for t_i in arange(dt, t_max + dt-__EPS, dt)]
	#line between start and end point
	return __line(line_points, weight)
	
## ("c", "dc")
#  constrains end com position and velocities to be equal to x_end
#  \param param requires "c_end", "t_init_phases"
#  \return 
def end_reached(param, weight):
	x_end  = array(param["x_end"][0:3])
	v_end  = array(param["x_end"][3:6])
	return lambda variables : weight * (__normsquared(variables["c"][-1] - x_end) + 0.5 * __normsquared(variables["dc"][-1] - v_end))

# constraint alpha value to be 1
def alpha_one(param,weight):
	print "objective alpha = ", param['alpha']
	return lambda variables : weight * (math.fabs(1.-variables["alpha"]))
	
__objective_factory = { 
	'min_ddL' : min_ddL,
	'min_dL' : min_dL,
	'min_ddc': min_ddc,
	'min_dddc': min_dddc,
	'line': line,
	'end_reached' : end_reached,
	'alpha_one' : alpha_one}

## From a user selected set of cost functions, generate objective functions
#  \param objective list of couple (string, value) describing the 
#  selected cost function and the weight given to it
#  \param params require at least 'simulate'
#  \return the cost function
def init_objective(objective, params):
	#retrieve all constraints
	simulate = params['simulate']
	cost_functions = [__objective_factory[name](params, weight) for name, weight in objective]
	if(not cost_functions):
		cost_functions = [lambda var: 0]
	def objective_fun (u):
		variables = simulate(u)
		return reduce(lambda a, b :a + b, [cost(variables) for cost in cost_functions])
	return objective_fun
	
