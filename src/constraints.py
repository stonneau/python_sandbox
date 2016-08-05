from scipy.linalg import block_diag 
from numpy import array, arange, zeros, ones, identity, vstack, append

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
	#~ return block_diag(*[cones[cone_id] for cone_id, phase in enumerate(phases) for _ in (range(phase[0],phase[1]))])
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
		
def test_end_reached_constraint():
	param = {'dt' : 0.5, "x_end" : [1,2,3,4,5,6], 't_init_phases' : [0,1.5]}
	A, b = end_reached_constraint(param)
	test_vec = append(zeros(12),param["x_end"])
	assert(((A.dot(test_vec) - b) == zeros(18)).all())
