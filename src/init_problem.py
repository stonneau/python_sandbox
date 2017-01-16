# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 4

@author: stonneau
"""

from compute_CWC import compute_CWC
from numpy import array
from create_simulation  import create_simulation

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
#  which formulation is adopted: either COM acceleration is the variable, or forces are the variables
#  \param p list of 3d contact positions
#  \param N list of 3d contact normals
#  \param x_input: list of exactly two 6D vectors, indicating start and end positions and velocities of the com 
#  \param t_phases list of ending times for each contact phases. For instance if first phase last 1 second, and 
#  second lasts 0.5 second, t_phases=[1,1.5] 
#  \param cones: None by default. If cones is a list of cones, they will be used
#   as CWC, instead of computing them
#  \param dt: step size
#  \return a dictionnary param, used throughout the optimization problem to set
#  up constraints and cost functions
def init_problem(p, N, x_input, t_end_phases, dt, cones = None, COMConstraints = None, mu =0.5, mass = 75, g = 9.81, simplify_cones = True):
	param = {"p"      		: [array(phase) for phase in p], 
	         "N" 	  		: [array(phase) for phase in N],
	         "x_init" 		: x_input[0],
	         "x_end" 		: x_input[1],
	         "dt" 			: dt,
	         "mass" 		: mass,
	         "mu" 			: mu,
	         "g" 			: g,
	         "COMCons" 		: COMConstraints,
	         #t_phases in updated to each starting phase time and one final phase
	         #~ "t_init_phases": [dt] +[t_end_phases[i] + dt for i in range(len(t_end_phases)-1)] + [t_end_phases[-1]] }
	         "t_init_phases": [0] + [t_end_phases[i] for i in range(len(t_end_phases)-1)] + [t_end_phases[-1]] }
	         
	#defining cone method compute all cones on first call, otherwise return hidden variable __cones	
	if cones == None:
		def f(phase ): return lambda: compute_CWC(param["p"][phase], param["N"][phase], simplify_cones, param)
		param['cones'] = __def_access_param_method(param, "cones", [f(phase) for phase in range(len(p))])
	else:
		assert(len(cones)==len(p))
		def f(phase ): return lambda: cones[phase]
		param['cones'] = __def_access_param_method(param, "cones", [f(phase) for phase in range(len(p))])
	param['simulate'] = create_simulation(param)
														  
	return param
