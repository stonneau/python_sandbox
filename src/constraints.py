from scipy.linalg import block_diag 
from numpy import array
from numpy import arange

## constraint can be of kind equality ("eq"),
# inequality ("ineq"). 
# They can apply to either control ("u")
# or the state ("x")
# TODO should I put state and control in the same matrix ?


## ("ineq","u")
#  Defines the cone constraints for the COM acceleration
#  and angular momentum. The active cone constraint for a 
#  phase is given by the switching time phase. Thus the 
#  resulting matrix is simply a diagonal stacking of the cone
#  active at a given frame
#  \param param requires "cones", "phases", and "dt"
#  \return size of configuration
def cones_constraint(param):
	cones  = param["cones"]()	
	phases = param["t_init_phases"]
	dt = param["dt"]
	#~ return block_diag(*[cones[cone_id] for cone_id, phase in enumerate(phases) for _ in (range(phase[0],phase[1]))])
	return block_diag(*[cones[index] for index, _ in enumerate(phases[:-1]) for _ in (arange(phases[index],phases[index+1],dt))])
	
def test_cones_constraint():
	param = {'cones' : lambda : [[['A','B'],['C','D']],[['E','F'],['G','H']]],
	 'dt' : 0.5,
	 't_init_phases' : [0,1,4]}
	cones = cones_constraint(param)
	print cones
	assert(cones.shape == (16,16))
	for i in range(0,4,2):
		assert cones[i][i] == "A"
		assert cones[i][i+1] == "B"
		assert cones[i+1][i] == "C"
		assert cones[i+1][i+1] == "D"
	for i in range(4,16,2):
		assert cones[i][i] == "E"
		assert cones[i][i+1] == "F"
		assert cones[i+1][i] == "G"
		assert cones[i+1][i+1] == "H"
