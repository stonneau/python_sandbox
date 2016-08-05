from scipy.linalg import block_diag 
from numpy import array, hstack
	
__lastcomputed = None

## Given an initial 6d X (c, c_p), and a command vector u
# computes the state c, c_p at each step of the simulation
#  \param param requires "x_init", and "dt"
#  \return dictionnary of variables.
def create_integrate(param):
	init_c   = array(param["x_init"][0:3]) #init position
	init_c_p = array(param["x_init"][3:6]) #init velocity
	init_l_p = array([0,0,0])  #init angular momentum
	dt = param["dt"] #time step
	def res_fun(u):
		# return computed variables if already done at this step
		global __lastcomputed
		if __lastcomputed != None and (__lastcomputed['u'] == u):
			return __lastcomputed
		u_n = [array(u[i:i+3]) for i in range(0, len(u), 3)]
		# else init variables and integrate forward in time
		res_c   = [init_c]; 
		res_c_p = [init_c_p];
		res_l_p = [init_l_p];
		[res_c_p.append(res_c_p[-1] + dt *u_i                            ) for u_i in u_n                  ];
		[  res_c.append(res_c  [-1] + dt *0.5*(res_c_p[i] + res_c_p[i-1])) for i   in range(1,len(res_c_p))];
		__lastcomputed = { 'c' : res_c[1:], 'c_d': res_c_p[1:], 'x' : hstack([res_c[1:],res_c_p[1:]]) , 'w':  u_n, 'u':  u}
		return __lastcomputed
	return res_fun


def test_integrate():
	integrate = create_integrate({ "x_init" : [1,2,3,0,0,0], "dt" : 0.5})
	u = [1,0,0,
		 1,0,0,
		 0,0,1]
	integration = integrate(u)
	assert((integration['c'][-1] == array([ 2., 2., 3.125])).all())
	assert(len(integration['x']) == len(integration['c']) == len(integration['c_d']) == len(integration['w']))
	
