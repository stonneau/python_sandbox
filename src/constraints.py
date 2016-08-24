from scipy.linalg import block_diag 
from numpy import array, arange, zeros, ones, identity, vstack, hstack, append

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
	A = block_diag(*[cones[index] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1],dt))])
	b = zeros(A.shape[0])
	return A, b 
	
## ("ineq","x")
#  Defines the cone constraints for the COM acceleration
#  and angular momentum. The active cone constraint for a 
#  phase is given by the switching time phase. Thus the 
#  resulting matrix is simply a diagonal stacking of the cone
#  active at a given frame
#  \param param requires "cones", "phases", "t_init_phases" and "dt"
#  \return cone matrices for each phase
def com_kinematic_constraint(param):
	comCons  = param["COMCons"]
	if(comCons == None):
		print "Error: com_kinematic_constraint activated but no constraints given "
	phases = param["t_init_phases"]
	dt = param["dt"]
	A = block_diag(*[comCons[index][0] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1],dt))])
	b = hstack([comCons[index][1] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1],dt))])
	return A, b 
	
#~ ## ("eq","x")
#~ #  constrains end com position and velocities to be equal to x_end
#~ #  \param param requires "c_end", "t_init_phases"
#~ #  \return 
def end_reached_constraint(param):
	x_end  = param["x_end"]
	x_size = len(x_end)
	phases = param["t_init_phases"]
	dt = param["dt"]
	zero_matrix = zeros((x_size,x_size))
	A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [identity(x_size)])	
	b = append([zeros(A.shape[0]-x_size)], [x_end])
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
	A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [identity(x_size)])	
	b = append([zeros(A.shape[0]-x_size)], [x_end])
	return A, b
	
## ("ineq","x")
#  constrains end com position and velocities to be almost equal to x_end
#  \param param requires "c_end", "t_init_phases"
#  \return 
def end_reached_constraint_plus(param):
	epsilon = 0.00000001
	x_end  = (array(param["x_end"]) + epsilon).tolist()
	x_size = len(x_end)
	phases = param["t_init_phases"]
	dt = param["dt"]
	zero_matrix = zeros((x_size,x_size))
	A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [identity(x_size)])	
	b = append([zeros(A.shape[0]-x_size)], [x_end])
	return A, b
	
def end_reached_constraint_minus(param):
	epsilon = 0.00000001
	x_end  = (epsilon -array(param["x_end"])).tolist()
	x_size = len(x_end)
	phases = param["t_init_phases"]
	dt = param["dt"]
	zero_matrix = zeros((x_size,x_size))
	A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [-identity(x_size)])	
	b = append([zeros(A.shape[0]-x_size)], [x_end])
	return A, b	

__constraint_factory = { 
	'end_reached_constraint'			: {'type': 'eq'  , 'var' : 'x', 'fun': end_reached_constraint},
	'end_reached_constraint_plus'		: {'type': 'ineq', 'var' : 'x', 'fun': end_reached_constraint_plus},
	'end_reached_constraint_minus'		: {'type': 'ineq', 'var' : 'x', 'fun': end_reached_constraint_minus},
	'cones_constraint' 		        	: {'type': 'ineq', 'var' : 'w', 'fun': cones_constraint},
	'com_kinematic_constraint' 		   	: {'type': 'ineq', 'var' : 'c', 'fun': com_kinematic_constraint},
	'end_null_acceleration_constraint'  : {'type': 'eq'  , 'var' : 'w', 'fun': end_null_acceleration_constraint}}

def __filter_cons(factory, cond_name, cond_value):
	return [c for c in factory if c[cond_name] == cond_value]

#initialize list of constraints, stack them in a single matrix
#and create the constraint function
def __stack_filter_cons(factory, cons_type, var_type, params):
	mat_and_vector= [c['fun'](params) for c in factory if c['var'] == var_type]
	if mat_and_vector:
		# stack matrices
		A = vstack([c[0] for c in mat_and_vector])
		b = array([c[1] for c in mat_and_vector]).flatten()
		def fun(var):
			var_of_interest = params['simulate'](var)[var_type]
			return A.dot(var_of_interest) - b
		return {'type': cons_type, 'fun' : fun}
		

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
	cons = [__constraint_factory[name] for name in constraints]
	# combine all possible constraint type and variable type => [['eq', 'w'], ['eq', 'x'], ..., ['ineq', 'w'] ...]
	# and create constraint method for each of those
	constraint_list = [__stack_filter_cons([c for c in cons if c['type'] == c_type], c_type, v_type, params) 
								for c_type in  c_types for v_type in v_types]
	return tuple([c for c in constraint_list if c != None])
	
		
