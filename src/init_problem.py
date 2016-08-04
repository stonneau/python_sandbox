# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 4

@author: stonneau
"""

from compute_CWC import compute_CWC
from numpy import array

## optional parameters are only computed if necessary when requested
def __def_access_param_method(param, keyword, fun):
	private_var = "__"+keyword
	def retfun():
		if not param.has_key(private_var):
			param[private_var] = [f() for f in fun] 
		return param[private_var]
	return retfun

## 
#  The optimization problem consists in finding a trajectory
#  between two or three contact phases (two if no contact are broken, then recreated).
#  To initialize it, for each phase the inputs expected
#  are the list of contact points and normals, the COM positions at start and
#  end phases, the duration of each phase, and finally discretization step.
#  which formulation is adopted: either COM acceleration is the variable, or forces are the variables
#  \param p list of 3d contact positions
#  \param N list of 3d contact normals
#  \param c_input: list of exactly two 3D vectors, indicating start and end positions of the com 
#  \param t_phases list of ending times for each contact phases. For instance if first phase last 1 second, and 
#  second lasts 0.5 second, t_phases=[1,1.5] 
#  \param dt: step size
#  \return a dictionnary param, used throughout the optimization problem to set
#  up constraints and cost functions
def init_problem(p, N, c_input, t_end_phases, dt, mu =0.5, mass = 75, g = 9.81):
	param = {"p"      		: [array(phase) for phase in p], 
	         "N" 	  		: [array(phase) for phase in N],
	         "c_init" 		: c_input[0],
	         "c_end" 		: c_input[1],
	         "dt" 			: dt,
	         "mass" 		: mass,
	         "mu" 			: mu,
	         "g" 			: g,
	         #t_phases in updated to each starting phase time and one final phase
	         "t_init_phases": [0] +[t_end_phases[i] + dt for i in range(len(t_end_phases)-1)] + [t_end_phases[-1]] }
	         
	#defining cone method compute all cones on first call, otherwise return hidden variable __cones	
	def f(phase ): return lambda: compute_CWC(param["p"][phase], param["N"][phase], param)
	param['cones'] = __def_access_param_method(param, "cones", [f(phase) for phase in range(len(p))])
														  
	return param
	
def test():
	from test_vars import p, N, t_end_phases
	param =  init_problem(p, N, ["com_s", "com_g"], t_end_phases,0.1)
#check phases
	assert(param["t_init_phases"]) == [0, 1.1, 1.5]
#check cones
	cones = param["cones"]()
	assert(len(cones)) == 2
	assert(cones[0].shape[0]) == 16
	assert(cones[0].shape[1]) == 6
	#~ return param
