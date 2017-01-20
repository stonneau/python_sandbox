# -*- coding: utf-8 -*-
"""
Created on Thurs Aug 4

@author: stonneau
"""

from init_problem import init_problem
from constraints  import init_constraints, init_parametric_constraints
from objective    import init_objective
from init_guess   import initial_guess_naive, initial_guess_naive_noise, initial_guess_naive_gravity_compensation, initial_guess_support, initial_guess_support_velocity
from definitions  import OptimError


from scipy.optimize import minimize

#formulation of the optimization problem
def cone_optimization(p, N, x_input, t_end_phases, dt, cones =None, COMConstraints = None, mu =0.5, mass = 75, g = 9.81, simplify_cones = True, verbose=False,
constraint_set=['cones_constraint', 'end_reached_constraint','end_speed_constraint'], parametric_constraint_set = [], score_treshold = 50000, initial_guess = []):
	params = init_problem(p, N, x_input, t_end_phases, dt, cones, COMConstraints, mu, mass, g, simplify_cones)
	#print "params :"
	#print params
	#print "p ="
	#print params["p"]
	#print "p[1] = "
	#print params["p"][1]
	cons = init_constraints(constraint_set, params) + init_parametric_constraints(parametric_constraint_set, params)
	objective = init_objective([["min_dddc", 0.001], ["min_ddc", 1],["alpha_one",1]],params)
	#~ objective = init_objective([["min_dddc", 0.001], ["min_ddc", 0.001], ["min_ddL", 40]],params)
	res= {'success' : False}
	if(len(initial_guess)>0):
			initial_guess.append(1) # alpha variable
			print "try with velocity_initial_guess : ",initial_guess
			res = minimize(objective, initial_guess, constraints=cons, method='SLSQP', options={'disp': verbose, 'ftol': 1e-06, 'maxiter' : 500})
	if (res ['success'] == False or res ['fun'] > score_treshold):
		if(verbose):
			print "error in minimization, trying with naive heuristic"
		init_guess = initial_guess_naive(params)
		res = minimize(objective, init_guess, constraints=cons, method='SLSQP', options={'disp': verbose, 'ftol': 1e-06, 'maxiter' : 500})
	if (res ['success'] == False or res ['fun'] > score_treshold):
		if(verbose):
			print "error in minimization, trying with support heuristic"
		objective = init_objective([["min_dddc", 0.001], ["min_ddc", 1],["alpha_one",1]],params)
		init_guess = initial_guess_support(params)
		res = minimize(objective, init_guess, constraints=cons, method='SLSQP', options={'disp': verbose, 'ftol': 1e-06, 'maxiter' : 500})
	if res ['success'] == False:
		raise OptimError("OPTIMIZATION FAILED EVERY TIME")
	var_final = params['simulate'](res['x'])
	print "final c_end", var_final["c_end"]
	print "final dc_end", var_final["dc_end"]
	params['alpha'] = res['x'][-1]
	print "final alpha = ", params['alpha']
	#print "final X = ", res['x'][:-1]
	return var_final, params

