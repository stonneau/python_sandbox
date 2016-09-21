# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 5

@author: stonneau
"""

from definitions import __EPS

from numpy import array, arange, zeros, identity
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
 	acc = [el for _ in arange(0,t_mid-__EPS,dt) for el in ddc] + [-el for _ in arange(t_mid,phases[-1]-__EPS,dt) for el in ddc ]
 	# now add 0 dL
 	res= acc + [0 for _ in range(0,len(acc))]
	assert len(res)%2 == 0, "in initial_guess_naive, u is not par"		
 	return res
 	
def __barycenter(points):
	x = zeros(3)
	for point in points:
		x = x + point
	return x / len(points)
 	
## initial guess with 2 lines, the middle point lying at the centroid
# of the support polygon
# straight line from start to goal
# accelerate constantly half the way, then decelerate constantly 
# by the same amount. For this solve for ddc such that
# at half time, the robot is mid-way.
# angular momentum variation set to 0 all along
# assume velocities at start and goal are zeros
# assume velocities at start and goal are zeros
# \param requires x_init, x_end, dt and t_init_phases
# \return velocities at start and goal are zeros
def initial_guess_support(params):
	phases = params["t_init_phases"]; x_init = params["x_init"];
	x_end  = params["x_end"]; dt = params["dt"];
	barycenter = __barycenter(params["p"][1])
	barycenter[2] = (x_init[2] + x_end[2]) /2.
	c_0 = array(x_init[0:3]); c_end = array(x_end[0:3])
	c_mid_1 = (barycenter - c_0) / 2.; # mid-point
	c_mid_2 = (c_end - barycenter) / 2.; # mid-point
	t_mid = phases[-1]/2.; # mid-time
 	ddc1 = 2*(c_mid_1 / t_mid**2); 
 	ddc2 = 2*(c_mid_2 / t_mid**2); 
 	#now find closest x*dt < t_mid, with x integer
 	t_mid = floor(t_mid /dt)*dt
 	#now just apply it for half the time, then negatively the rest
 	acc = [el for _ in arange(0,t_mid-__EPS,dt) for el in ddc1] + [-el for _ in arange(t_mid,phases[-1]-__EPS,dt) for el in ddc2 ]
 	# now add 0 dL
 	#~ res = acc + [0 for _ in range(0,3) for _ in arange(0,phases[-1]-__EPS,dt)] 	
 	res = acc + [0 for _ in range(0,len(acc))]
 	assert len(res)%2 == 0, "in initial_guess_support, u is not par"
 	return res
 	
## initial guess
# which tries to compute two straight lines, with one middle
# point in the middle cone.
# accelerate constantly half the way, then decelerate constantly 
# by the same amount. For this solve for ddc such that
# at half time, the robot is mid-way.
# angular momentum variation set to 0 all along
# assume velocities at start and goal are zeros
# assume velocities at start and goal are zeros
# \param requires x_init, x_end, dt and t_init_phases
# \return velocities at start and goal are zeros
def initial_guess_naive_noise(params):
	phases = params["t_init_phases"]; x_init = params["x_init"];
	x_end  = params["x_end"]; dt = params["dt"];
	c_0 = array(x_init[0:3]); c_end = array(x_end[0:3])
	c_mid_min_c0 = (c_end - c_0) / 2.; # mid-point
	t_mid = phases[-1]/2. + 4 * dt; # mid-time
 	ddc = 2*(c_mid_min_c0 / t_mid**2); 
 	#now find closest x*dt < t_mid, with x integer
 	t_mid = floor(t_mid /dt)*dt
 	#now just apply it for half the time, then negatively the rest
 	acc = [el for _ in arange(0,t_mid-__EPS,dt) for el in ddc] + [-el for _ in arange(t_mid,phases[-1]-__EPS,dt) for el in ddc ]
 	# now add 0 dL
 	#~ res = acc + [0 for _ in range(0,3) for _ in arange(0,phases[-1]-__EPS,dt)]
 	res = acc + [0 for _ in range(0,len(acc))]
 	assert len(res)%2 == 0, "in initial_guess_naive_noise, u is not par"
 	return res
 	
## initial guess
# that simply compensates for the gravity
# at half time, the robot is mid-way.
# angular momentum variation set to 0 all along
# assume velocities at start and goal are zeros
# assume velocities at start and goal are zeros
# \param requires x_init, x_end, dt and t_init_phases
# \return velocities at start and goal are zeros
def initial_guess_naive_gravity_compensation(params):
	#compute naive initial guess
	#~ naive_guess = initial_guess_naive(params)
	gdt = params["dt"] * params["g"];
	phases = params["t_init_phases"];
	ddc = [0,0,gdt]
	res = [el for _ in arange(0,phases[-1]-__EPS,dt) for el in ddc]
	res += [0 for _ in range(0,len(acc))]
 	assert len(res)%2 == 0, "in initial_guess_naive_gravity_compensation, u is not par"
 	return res
