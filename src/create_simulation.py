from scipy.linalg import block_diag 
from numpy import array, hstack, cross, arange, copy
from numpy.linalg import norm
from definitions import __EPS
	

#~ __lastcomputed = []
#~ __last_u = []
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
		#~ global __lastcomputed
		#~ global __last_u
		# return computed variables if already done at this step
		#~ for i, u_i in enumerate(__last_u):
			#~ if (u_i == u).all():
				#~ return __lastcomputed[i]
		# else init variables and integrate forward in time
		ddc = [array(u[i:i+3]) for i in range(0, len(u)/2, 3)]
		dL 	= [array(u[i:i+3]) for i in range(len(u)/2, len(u), 3)]
		c  	= [init_c]; 
		dc 	= [init_dc];		
		dddc= [(ddc[i+1] - ddc[i]) / dt for i in range(len(ddc)-1)];	
		ddL= [(dL[i+1] - dL[i]) / dt for i in range(len(dL)-1)];		
		# dc
		# ambiguity comes from the fact that control vector u has one less entry than the state vector x
		# euler integration for velocity gives:
		# dc[i+1] = ddc[i+1] *dt + dc[i].
		# but is effectively implemented as dc[i+1] = ddc[i] *dt + dc[i] since ddc is offset
		[dc.append(dc[-1] + dt *ddc_i ) for ddc_i in ddc ];		
		#c
		#~ [ c.append(c[-1] + dt *0.5*(dc[i] + dc[i-1])) for i in range(1,len(dc))];		
		[ c.append(c[-1] + dt *(dc[i]) + dt *dt * 0.5 * ddc[i]) for i in range(0,len(ddc))];		
		#remove init values from c and dc
		c  =  c[1:]
		dc = dc[1:]
		#w
		y = [m * (ddc_i - g_vec) for ddc_i in ddc]
		#~ w = [y[i].tolist() + (cross(c[i], y[i]) + dL[i]).tolist() for i,_ in enumerate(c)]
		w = [y[i].tolist() + (cross(c[i], y[i])).tolist() for i,_ in enumerate(c)] #dL is 0
		# WARNING ! in robustequilibrium lib, the GIWC computation is reversed. Therefore it is required
		# to test equilibrium with -w rather than w ...  changing this while robust-equilibrium lib is not updated
		x = [c[i].tolist() + dc[i].tolist() for i,_ in enumerate(c)]
		result = { "x_init" : param["x_init"], 'c' : c, 'c_end' : c[-1], 'dc_end' : dc[-1], 'dc': dc, 'ddc' : ddc, 'dddc' : dddc, 
		'x' : array(x).flatten() , 'w':  -array(w).flatten(), 'u':  u, 'dL' : dL, 'ddL' : ddL}		
		#~ __last_u.append(u)
		#~ __lastcomputed.append(result)
		#~ if len(__last_u) > 2* len(u):
			#~ __last_u = __last_u[len(u):]
			#~ __lastcomputed = __lastcomputed[len(u):]
		return result
	return res_fun


		
