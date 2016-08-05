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
	
## ("eq","x")
#  constrains end com position and velocities to be equal to x_end
#  \param param requires "c_end", "t_init_phases"
#  \return 
def end_reached_constraint(param):
	x_end  = param["x_end"]
	x_size = len(x_end)
	phases = param["t_init_phases"]
	dt = param["dt"]
	zero_matrix = zeros((x_size,x_size))
	A = block_diag(*[zero_matrix for _ in (arange(phases[0],phases[-1]-dt,dt))] + [identity(x_size)])	
	b = append([zeros(A.shape[0]-x_size)], [x_end])
	return A, b
	
__constraint_factory = { 
	'end_reached_constraint': {'type': 'eq'  , 'var' : 'x',  'fun': end_reached_constraint},
	'cones_constraint' 		: {'type': 'ineq', 'var' : 'w',  'fun': cones_constraint		}}

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
		
__cons_type = ['eq','ineq']
__var_type = ['w','x']

## From a user selected set of constraints, initialize the constraints array to pass to a minimization
#  problem.
#  \param constraints list of strings describing the selected constraints
#  \param params parameter dictionnary obtained from a call to init_problem
#  \return a dictionnary of constraints to pass to the minimization problem
def init_constraints(constraints, params):
	#retrieve all constraints
	cons = [__constraint_factory[name] for name in constraints]
	# combine all possible constraint type and variable type => [['eq', 'w'], ['eq', 'x'], ..., ['ineq', 'w'] ...]
	# and create constraint method for each of those
	constraint_list = [__stack_filter_cons([c for c in cons if c['type'] == c_type], c_type, v_type, params) 
								for c_type in  __cons_type for v_type in __var_type]
	return tuple([c for c in constraint_list if c != None])
	
		
def test_cones_constraint():
	param = {'cones' : lambda : [[[1,2],[3,4]],[[5,6],[7,8]]],
	 'dt' : 0.5,
	 't_init_phases' : [0,1,4]}
	A, b = cones_constraint(param)
	print A
	assert(b.shape==(16,))
	assert((b==zeros(16)).all())
	assert(A.shape == (16,16))
	for i in range(0,4,2):
		assert A[i][i] == 1
		assert A[i][i+1] == 2
		assert A[i+1][i] == 3
		assert A[i+1][i+1] == 4
	for i in range(4,16,2):
		assert A[i][i] == 5
		assert A[i][i+1] == 6
		assert A[i+1][i] == 7
		assert A[i+1][i+1] == 8
	print "test exited normally" 
		
def test_end_reached_constraint():
	param = {'dt' : 0.5, "x_end" : [1,2,3,4,5,6], 't_init_phases' : [0,1.5]}
	A, b = end_reached_constraint(param)
	test_vec = append(zeros(12),param["x_end"])
	assert(((A.dot(test_vec) - b) == zeros(18)).all())
	print "test exited normally" 
	

def test_init_constraints():
	params = {'cones' : lambda : [[[1,1,1,1,1,1],[1,1,1,1,1,1]],[[1,1,1,1,1,1],[1,1,1,1,1,1]]],
	 'dt' : 0.5,
	 't_init_phases' : [0,1.5],
	 "x_end" : [1,2,3,4,5,6],
	 "simulate" : lambda x: { 'x' : array([i+1 for _ in arange(0,1.5,0.5) for i in range(1,7)]),
							  'w' : array([1 for _ in arange(0,1.5,0.5) for i in range(6)])}}
	 
	cons = init_constraints(["end_reached_constraint","end_reached_constraint",'cones_constraint','cones_constraint'], params)
	assert(len(cons)==2)
	espected_eq_x = array(zeros(12).tolist() + ones(6).tolist() + zeros(12).tolist() + ones(6).tolist());
	expected_ineq_w = 6*(ones(12));
	assert((cons[0]['fun']([1])==espected_eq_x).all()) 
	assert((cons[1]['fun']([1])==expected_ineq_w).all())
	print "test exited normally" 
