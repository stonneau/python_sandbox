# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 4

@author: stonneau
"""

from init_problem import init_problem
from constraints  import init_constraints
from objective    import init_objective
from init_guess   import initial_guess_naive, initial_guess_naive_noise, initial_guess_naive_gravity_compensation

from scipy.optimize import minimize

#formulation of the optimization problem
def cone_optimization(p, N, x_input, t_end_phases, dt, cones =None, COMConstraints = None, mu =0.5, mass = 75, g = 9.81, simplify_cones = True, verbose=False):
	params = init_problem(p, N, x_input, t_end_phases, dt, cones, COMConstraints, mu, mass, g, simplify_cones)
	cons = init_constraints(['cones_constraint', 'end_reached_constraint_plus', 'end_reached_constraint_minus'], params)
	objective = init_objective([["line", 10]],params)
	init_guess = initial_guess_naive(params)
	res = minimize(objective, init_guess, constraints=cons, method='SLSQP', options={'disp': verbose, 'ftol': 1e-03})
	var_final = params['simulate'](res['x'])
	return var_final, params

