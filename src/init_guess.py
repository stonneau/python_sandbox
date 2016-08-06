# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 5

@author: stonneau
"""



from numpy import array, arange
from numpy.linalg import norm
from math import floor

def normalize(v):
    n = norm(v)
    if n==0: 
       return v
    return v/n

## really naive initial guess
# straight line from start to goal
# accelerate constantly half the way, then decelerate constantly 
# by the same amount. For this solve for ddc such that
# at half time, the robot is mid-way.
# angular momentum variation set to 0 all along
# assume velocities at start and goal are zeros
# assume velocities at start and goal are zeros
# \param requires x_init, x_end, dt and t_init_phases
# \return velocities at start and goal are zeros
def initial_guess_naive(params):
	phases = params["t_init_phases"]; x_init = params["x_init"];
	x_end  = params["x_end"]; dt = params["dt"];
	c_0 = array(x_init[0:3]); c_end = array(x_end[0:3])
	c_mid_min_c0 = (c_end - c_0) / 2.; # mid-point
	t_mid = phases[-1]/2.; # mid-time
 	ddc = 2*(c_mid_min_c0 / t_mid**2); 
 	#now find closest x*dt < t_mid, with x integer
 	t_mid = floor(t_mid /dt)*dt
 	#now just apply it for half the time, then negatively the rest
 	acc = [el for _ in arange(0,t_mid,dt) for el in ddc] + [-el for _ in arange(t_mid,phases[-1],dt) for el in ddc ]
 	# now add 0 dL
 	return acc + [0 for _ in range(0,3) for _ in arange(0,phases[-1],dt)]
 	


def test_initial_guess_naive():
	from create_simulation import create_simulation
	params = { 'dt'          : 0.5,
			  't_init_phases': [0,1,2],
			  'x_init'       : [0,0,0,0,0,0],
			  'x_end'        : [1,0,1,0,0,0],
			  'g' 			 : 10,
			  'mass' 		 : 10 }
	  
	res = initial_guess_naive(params)
	assert(len(res)==24)

	sim =create_simulation(params)
	res = sim(res)['c']
	assert( (res[3][0:3]==array([1,0,1])).all() )
	print "test_initial_guess_naive exited normally" 
		  
