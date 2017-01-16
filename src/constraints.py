from scipy.linalg import block_diag 
from numpy import array, arange, zeros, ones, identity, vstack, hstack, append
from definitions import OptimError, __EPS

## constraint can be of kind equality ("eq"),
# inequality ("ineq"). 
# They can apply to either the control ("u")
# or the state ("x")
# Each constraint returns  2 arguments: a premultiplication
# matrix A, and a vector b, such that the tested constraint
# is given by A * var -b, where var is the selected variable


## ("ineq","w")
#  Defines the cone constraints for the COM acceleration
#  and angular momentum. The active cone constraint for a 
#  phase is given by the switching time phase. Thus the 
#  resulting matrix is simply a diagonal stacking of the cone
#  active at a given frame
#  \param param requires "cones", "phases", "t_init_phases" and "dt"
#  \return cone matrices for each phase
def cones_constraint(param):
	cones  = param["cones"]()	
	phases = param["t_init_phases"]
	dt = param["dt"]
	indexes = [index for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1]-__EPS,dt))]
	for i in range(0,len(indexes)-1):
		if not (indexes[i] == indexes[i+1]):
			indexes[i] = indexes[i+1] #state before flight constrained to be in flight cone
			break
	A = block_diag(*[cones[index] for _, index in enumerate(indexes)])
	#~ A = block_diag(*[cones[index] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1]-__EPS,dt))])
	b = zeros(A.shape[0])
	b = -1 * ones(A.shape[0])
	return A, b 
	
## ("ineq","x")
#  Defines the kinematic constraints for the COM position
#  depending on the contact positions. The active com constraint for a 
#  phase is given by the switching time phase. Thus the 
#  resulting matrix is simply a diagonal stacking of the constraint
#  active at a given frame
#  \param param requires "COMCons", "t_init_phases" and "dt"
#  \return cone matrices for each phase
def com_kinematic_constraint(param):
	comCons  = param["COMCons"]
	if(comCons == None):
		raise OptimError("com_kinematic_constraint activated but no constraints given")
	consA = []; consb = [];
	for _, c in enumerate(comCons):
		A = c[0]; b = c[1]
		a = zeros([A.shape[0],A.shape[1]+3])
		a [:,:-3] = A
		consA.append(a)
		consb.append(b)
		
	phases = param["t_init_phases"]
	dt = param["dt"]
	A = block_diag(*[consA[index] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1]-__EPS,dt))])
	b = hstack([consb[index] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1]-__EPS,dt))])
	return A, b 
	
#~ ## ("eq","x_end")
#~ #  constrains end com position and velocities to be equal to x_end
#~ #  \param param requires "c_end", "t_init_phases"
#~ #  \return 
def end_reached_constraint(param):
	x_end  = param["x_end"]
	#~ x_size = len(x_end)
	#~ phases = param["t_init_phases"]
	#~ dt = param["dt"]
	#~ zero_matrix = zeros((x_size,x_size))
	#~ A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [identity(x_size)])
	A = 10000 * identity(3)	
	b = 10000 * array(x_end[0:3])
	return A, b
	
#~ ## ("eq","c")
#~ #  constrains end com position and velocities to be equal to x_end
#~ #  \param param requires "c_end", "t_init_phases"
#~ #  \return 
def waypoint_reached_constraint(param, step, value):
	return lambda(c): c[step] - value
	
#~ ## ("eq","dc_end")
#~ #  constrains end com position and velocities to be equal to x_end
#~ #  \param param requires "c_end", "t_init_phases"
#~ #  \return 
def end_speed_constraint(param):
	x_end  = param["x_end"]
	#~ x_size = len(x_end)
	#~ phases = param["t_init_phases"]
	#~ dt = param["dt"]
	#~ zero_matrix = zeros((x_size,x_size))
	#~ A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [identity(x_size)])
	A = identity(3)	
	b = x_end[3:6]
	return A, b
	
#~ ## ("eq","w")
#~ #  constrains end com position and velocities to be equal to x_end
#~ #  \param param requires "c_end", "t_init_phases"
#~ #  \return 
def end_null_acceleration_constraint(param):
	x_size = 6
	phases = param["t_init_phases"]
	dt = param["dt"]
	x_end = zeros((x_size))
	zero_matrix = zeros((x_size,x_size))
	A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt-__EPS,dt))] + [identity(x_size)])	
	b = append([zeros(A.shape[0]-x_size)], [x_end])
	return A, b
	
def __make_id_half(x_size):
	res = identity(x_size)
	for i in range(x_size / 2, x_size):
		res[i,i] = 0
	return res


__constraint_factory = { 
	'end_reached_constraint'			: {'type': 'eq'  , 'var' : 'c_end', 'fun': end_reached_constraint},
	'end_speed_constraint'				: {'type': 'eq'  , 'var' : 'dc_end', 'fun': end_speed_constraint},
	'cones_constraint' 		        	: {'type': 'ineq', 'var' : 'w', 'fun': cones_constraint},
	'com_kinematic_constraint' 		   	: {'type': 'ineq', 'var' : 'x', 'fun': com_kinematic_constraint},
	'end_null_acceleration_constraint'  : {'type': 'eq'  , 'var' : 'w', 'fun': end_null_acceleration_constraint}}

__parametric_constraint_factory = {
	'waypoint_reached_constraint'		: {'type': 'eq'  , 'var' : 'c', 'fun': waypoint_reached_constraint},
}

def __filter_cons(factory, cond_name, cond_value):
	return [c for c in factory if c[cond_name] == cond_value]

from numpy import where

#initialize list of constraints, stack them in a single matrix
#and create the constraint function
def __stack_filter_cons(factory, cons_type, var_type, params):
	mat_and_vector= [c['fun'](params) for c in factory if c['var'] == var_type]
	if mat_and_vector:
		# stack matrices
		A = vstack([c[0] for c in mat_and_vector])
		b = array([el for els in mat_and_vector for el in els[1] ])
		def fun(var):
			var_of_interest = params['simulate'](var)[var_type]
			#~ if(var_type == "w"):
				#~ print "debug"
				#~ print "raa", len(where(-(A.dot(var_of_interest) - b)<=0))
				#~ print "end debug"
			return -(A.dot(var_of_interest) - b)
		return {'type': cons_type, 'fun' : fun}
		
#initialize list of constraints, stack them in a single matrix
#and create the constraint function
def __gen_parametric_cons(param_constraint, params, args):
	fun = param_constraint['fun'](params, *args)
	var_type = param_constraint['var']
	cons_type = param_constraint['type']
	def fun2(var):
		var_of_interest = params['simulate'](var)[var_type]
		return fun(var_of_interest)
	return {'type': cons_type, 'fun' : fun2}


## From a user selected set of constraints, initialize the constraints array to pass to a minimization
#  problem.
#  \param constraints list of strings describing the selected constraints
#  \param params parameter dictionnary obtained from a call to init_problem
#  \return a dictionnary of constraints to pass to the minimization problem
def init_constraints(constraints, params):
	#retrieve all constraints type
	c_types = ['eq', 'ineq'];
	v_types = set([]); [v_types.add(__constraint_factory[name]['var' ]) for name in __constraint_factory]	
	#retrieve all constraints
	cons = [__constraint_factory[name] for name in set(constraints)]
	# combine all possible constraint type and variable type => [['eq', 'w'], ['eq', 'x'], ..., ['ineq', 'w'] ...]
	# and create constraint method for each of those
	constraint_list = [__stack_filter_cons([c for c in cons if c['type'] == c_type], c_type, v_type, params) 
								for c_type in  c_types for v_type in v_types]
	return tuple([c for c in constraint_list if c != None])
	
## From a user selected set of constraints, initialize the constraints array to pass to a minimization
#  problem. These constraints can be parametrized by extra arguments associated with 
# their definition
#  \param constraints list of tuple ('constraint_name', tuple(arg0, arg1...) ) describing the selected constraints and parameters
#  \param params parameter dictionnary obtained from a call to init_problem
#  \return a dictionnary of constraints to pass to the minimization problem
def init_parametric_constraints(constraints, params):
	constraint_list = [__gen_parametric_cons(__parametric_constraint_factory[cons_name], params, args) for cons_name, args in constraints]	
	return tuple([c for c in constraint_list if c != None])
	
		
