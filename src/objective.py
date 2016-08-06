from scipy.linalg import block_diag, norm
from numpy import array, arange, zeros, ones, identity, vstack, hstack, append, sqrt

## objective is  similar to constraint.
# costs functions are chosen from a objective factory, with
# user-defined weights.
# the factory automatically builds the cost function
# and returns it


def __sum_over(var_name, weight):
	return lambda variables: weight * norm(reduce(lambda a, b :a + b, variables[var_name]))

## ("dL")
#  Minimize angular momentum variation
#  \return  cone matrices for each phase
def min_dL(param, weight):
	return __sum_over('dL', weight)
	
## ("ddc")
#  Minimize angular momentum variation
#  \return  cone matrices for each phase
def min_ddc(param, weight):
	return __sum_over('ddc', weight)
	
## ("x")
#  constrains end com position and velocities to be equal to x_end
#  \param param requires "c_end", "t_init_phases"
#  \return 
def end_reached(param, weight):
	x_end  = array(param["x_end"])
	return lambda variables : weight * norm(variables["x"][-1] - x_end)
	
__objective_factory = { 
	'min_dL' : min_dL,
	'min_ddc': min_ddc,
	'end_reached' : end_reached}

## From a user selected set of cost functions, generate objective functions
#  \param objective list of couple (string, value) describing the 
#  selected cost function and the weight given to it
#  \param params require at least 'simulate'
#  \return the cost function
def init_objective(objective, params):
	#retrieve all constraints
	simulate = params['simulate']
	cost_functions = [__objective_factory[name](params, weight) for name, weight in objective]
	def objective_fun (u):
		variables = simulate(u)
		return reduce(lambda a, b :a + b, [cost(variables) for cost in cost_functions])
	return objective_fun
	
