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
def cone_optimization(p, N, x_input, t_end_phases, dt, mu =0.5, mass = 75, g = 9.81, simplify_cones = True):
	params = init_problem(p, N, x_input, t_end_phases, dt, mu, mass, g, simplify_cones)
	cons = init_constraints(['cones_constraint'], params)
	objective = init_objective([["end_reached", 1]],params)
	init_guess = initial_guess_naive(params)
	return minimize(objective, init_guess, constraints=cons, method='SLSQP', options={'disp': True}), params


#~ def test():
from test_vars import *
res, params = cone_optimization(p_3, N_3, [x_init_3, x_end_3], t_end_phases_3, dt, mu, mass, g)

var_final = params['simulate'](res['x'])

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 100
points = var_final['x']
xs = [points[i] for i in range(0,len(points),6)]
ys = [points[i] for i in range(1,len(points),6)]
zs = [points[i] for i in range(2,len(points),6)]
ax.scatter(xs, ys, zs, c='b')

colors = ["r", "b", "g"]
#print contact points of first phase
for id_c, points in enumerate(p_3):
	xs = [point[0] for point in points]
	ys = [point[1] for point in points]
	zs = [point[2] for point in points]
	ax.scatter(xs, ys, zs, c=colors[id_c])
	
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
