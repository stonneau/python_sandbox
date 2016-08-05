from scipy.linalg import block_diag 
from numpy import array, hstack, cross
	
__lastcomputed = None

## Given an initial 6d X (c, dc), and a command vector u
# computes the state c, dc at each step of the simulation
#  \param param requires "x_init", and "dt"
#  \return dictionnary of variables.
def create_simulation(param):
	init_c   = array(param["x_init"][0:3]) #init position
	init_dc = array(param["x_init"][3:6]) #init velocity
	init_dL = array([0,0,0])  #init angular momentum
	g_vec = array([0,0,-param["g"]])
	m = param["mass"]
	dt = param["dt"] #time step
	# u is the control [ddc[0],...,ddc[n],L_p[0],L_p[n]]
	def res_fun(u):
		# return computed variables if already done at this step
		global __lastcomputed
		if __lastcomputed != None and (__lastcomputed['u'] == u):
			return __lastcomputed
		# else init variables and integrate forward in time
		ddc = [array(u[i:i+3]) for i in range(0, len(u)/2, 3)]
		dL 	= [array(u[i:i+3]) for i in range(len(u)/2, len(u), 3)]
		c  	= [init_c]; 
		dc 	= [init_dc];		
		# dc
		# ambiguity comes from the fact that control vector u has one less entry than the state vector x
		# euler integration for velocity gives:
		# dc[i+1] = ddc[i+1] *dt + dc[i].
		# but is effectively implemented as dc[i+1] = ddc[i] *dt + dc[i] since ddc is offset
		[dc.append(dc[-1] + dt *ddc_i ) for ddc_i in ddc ];		
		#c
		[ c.append(c[-1] + dt *0.5*(dc[i] + dc[i-1])) for i in range(1,len(dc))];		
		#remove init values from c and dc
		c  =  c[1:]
		dc = dc[1:]
		#w
		y = [m * (ddc_i - g_vec) for ddc_i in ddc]
		w = [y[i].tolist() + (cross(c[i], y[i]) + dL[i]).tolist() for i,_ in enumerate(c)]
		x = [c[i].tolist() + dc[i].tolist() for i,_ in enumerate(c)]
		__lastcomputed = { 'c' : c, 'dc': c, 'ddc' : ddc, 'x' : array(x) , 'w':  array(w), 'u':  u, 'dL' : dL}
		return __lastcomputed
	return res_fun


def test_create_simulation():
	simulation = create_simulation({ "x_init" : [1,2,3,0,0,0], "dt" : 0.5, "g" : 10, "mass" : 10})
	u = [1,0,0, #ddc
		 1,0,0,
		 0,0,1,
		 1,0,0, #dL
		 1,0,0,
		 0,0,1]
	simulation = simulation(u)
	assert((simulation['c'][-1] == array([ 2., 2., 3.125])).all())
	assert(len(simulation['x']) == len(simulation['c']) == len(simulation['dc']) == len(simulation['w']))
	print "test exited normally" 
		
