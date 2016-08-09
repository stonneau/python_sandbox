# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 4

@author: stonneau
"""

from init_problem import init_problem
from constraints  import init_constraints
from objective    import init_objective
from init_guess   import initial_guess_naive

from scipy.optimize import minimize

#formulation of the optimization problem
def cone_optimization(p, N, x_input, t_end_phases, dt, cones =None, mu =0.5, mass = 75, g = 9.81, simplify_cones = True):
	params = init_problem(p, N, x_input, t_end_phases, dt, cones, mu, mass, g, simplify_cones)
	cons = init_constraints(['cones_constraint'], params)
	objective = init_objective([["end_reached", 10]],params)
	init_guess = initial_guess_naive(params)
	res = minimize(objective, init_guess, constraints=cons, method='SLSQP', options={'disp': True})
	var_final = params['simulate'](res['x'])
	return var_final, params

