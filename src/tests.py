"""
Created on Thurs Aug 6

@author: stonneau
"""

from compute_CWC       import *
from constraints  	   import *
from objective         import *
from init_guess        import *
from init_problem      import *
from main              import *
from test_vars         import *
from create_simulation import *

from scipy.linalg import block_diag 
from numpy import array, arange, zeros, ones, identity, vstack, hstack, append
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#~ *****************************
#~ *  CWC
#~ *****************************
def test_compute_CWC():
	from test_vars import p, N, mass, g, mu
	return compute_CWC(array(p[1]),array(N[1]), True, {"mass" : mass, "g": g, "mu":mu })


#~ *****************************
#~ *  simulation
#~ *****************************
def test_create_simulation():
	simulation = create_simulation({ "x_init" : [1,2,3,0,0,0], "dt" : 0.5, "g" : 10, "mass" : 10})
	u = array([1,0,0, #ddc
		 1,0,0,
		 0,0,1,
		 1,0,0, #dL
		 1,0,0,
		 0,0,1])
	simulation = simulation(u)
	assert((simulation['c'][-1] == array([ 2., 2., 3.125])).all())
	assert(len(simulation['x']) == len(simulation['w']))
	assert(len(simulation['c']) == len(simulation['dc']) == len(simulation['ddc'])  == len(simulation['dL']))
	print "test exited normally" 

#~ *****************************
#~ *  Problem
#~ *****************************
def test_init_problem():
	from test_vars import p, N, t_end_phases
	param =  init_problem(p, N, ["com_s", "com_g"], t_end_phases,0.1)
#check phases
	assert(param["t_init_phases"]) == [0, 1, 1.5]
#check cones
	cones = param["cones"]()
	assert(len(cones)) == 2
	assert(cones[0].shape[0]) == 16
	assert(cones[0].shape[1]) == 6
	print "test exited normally" 

#~ *****************************
#~ *  constraints
#~ *****************************
def test_cones_constraint():
	param = {'cones' : lambda : [array([[1,2],[3,4]]),array([[5,6],[7,8]])],
	 'dt' : 0.5,
	 't_init_phases' : [0,1,4]}
	A, b = cones_constraint(param)
	assert(b.shape==(16,))
	assert((b==zeros(16)).all())
	assert(A.shape == (16,16))
	for i in range(0,4,2):
		assert A[i][i] == 1
		assert A[i][i+1] == 2
		assert A[i+1][i] == 3
		assert A[i+1][i+1] == 4
	for i in range(4,16,2):
		assert A[i][i] == 5
		assert A[i][i+1] == 6
		assert A[i+1][i] == 7
		assert A[i+1][i+1] == 8
	print "test exited normally" 
	
#~ *****************************
#~ *  constraints
#~ *****************************
def test_com_kinematic_constraint():
	param = {'COMCons' : [[identity(3),array([0,0,0])],[identity(3), array([0,0,1])]],
	 'dt' : 0.5,
	 't_init_phases' : [0,1,4]}
	A, b = com_kinematic_constraint(param)
	assert(b.shape==(24,))
	assert(A.shape == (24,24))
	test = zeros(24); test[8]=1;
	assert((A.dot(test) - b) <= 0).all()
	test[9] = 1 
	assert(not((A.dot(test) - b) <= 0).all())
	print "test exited normally" 
		
def test_end_reached_constraint():
	param = {'dt' : 0.5, "x_end" : [1,2,3,4,5,6], 't_init_phases' : [0,1.5]}
	A, b = end_reached_constraint(param)
	test_vec = append(zeros(12),param["x_end"])
	assert(((A.dot(test_vec) - b) == zeros(18)).all())
	print "test exited normally" 

def test_init_constraints():
	params = {'cones' : lambda : [array([[1,1,1,1,1,1],[1,1,1,1,1,1]]) for _ in range(2)],
	 'dt' : 0.5,
	 't_init_phases' : [0,1.5],
	 "x_end" : [1,2,3,4,5,6],
	 "simulate" : lambda x: { 'x' : array([i+1 for _ in arange(0,1.5,0.5) for i in range(1,7)]),
							  'w' : array([1 for _ in arange(0,1.5,0.5) for i in range(6)])}}
	 
	cons = init_constraints(["end_reached_constraint","end_reached_constraint",'cones_constraint','cones_constraint'], params)
	assert(len(cons)==2)
	espected_eq_x = array(zeros(12).tolist() + ones(6).tolist() + zeros(12).tolist() + ones(6).tolist());
	expected_ineq_w = 6*(ones(12));
	assert((cons[0]['fun']([1])==espected_eq_x).all()) 
	assert((cons[1]['fun']([1])==expected_ineq_w).all())
	print "test exited normally" 

#~ *****************************
#~ *  initial guess
#~ *****************************
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
	

#~ *****************************
#~ *  objective
#~ *****************************
def test_init_objective():
	variables = {'dL' : [array([1,0,0]) for _ in range(4)], 
				 'ddc': [array([0.5,0,0]) for i in range (4)], 
				 'x'  :  array([ [0 for _ in range (6)] for _ in range(4)]),
				 'c'  :  array([ [0 for _ in range (3)] for _ in range(4)]),
				 'dc' :  array([ [0 for _ in range (3)] for _ in range(4)])}

	params = {'x_end' : [i for i in range(6)], 'simulate' : lambda(_): variables}
	# test each cost individually
	objective = init_objective([["min_dL", 1]], params)
	assert(objective(variables)==4.)

	objective = init_objective([["min_dL", 0.5]], params)
	assert(objective(variables)==2.)

	objective = init_objective([["min_ddc", 1]], params)
	assert(objective(variables)==2.)

	objective = init_objective([["min_ddc", 2]], params)
	assert(objective(variables)==4.)

	objective = init_objective([["end_reached", 1]], params)
	assert(objective(variables)==norm(array(params['x_end'][0:3])) + 0.5 * norm(array(params['x_end'][3:6])))

	objective = init_objective([["end_reached", 2]], params)
	assert(objective(variables)==2* (norm(array(params['x_end'][0:3])) + 0.5 * norm(array(params['x_end'][3:6]))))

	objective = init_objective([["min_dL", 0.5],["min_ddc", 2],["end_reached", 1]], params)
	assert(objective(variables)==2 + 4 + norm(array(params['x_end'][0:3])) + 0.5 * norm(array(params['x_end'][3:6])))

#~ *****************************
#~ *  integration test
#~ *****************************
def test_optimize():

	#~ print len(p_hyq)
	#~ print len(N_hyq)

	var_final, params = cone_optimization(p_3, N_3, [x_init_3, x_end_3], t_end_phases_3, dt, None, None, mu, mass, g)
	#~ var_final, params = cone_optimization(p, N, [x_init, x_end], t_end_phases, dt, mu, mass, g)
	#~ p = [[p_hyq_j[i] for i in range(0,len(p_hyq_j),4)] for p_hyq_j in p_hyq]
	#~ N = [[N_hyq_j[i] for i in range(0,len(N_hyq_j),4)] for N_hyq_j in N_hyq]
	#~ print N
	#~ print len(N[0])
	#~ print len(N_hyq[0])
	#~ res, params = cone_optimization(p, N, [x_init_hyq, x_end_hyq], t_end_phases_3, dt, mu_hyq, mass, g)


	#~ var_final = params['simulate'](res['x'])

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


#comment tests you do not want to run
tests_run = {
	'test_compute_CWC' : test_compute_CWC,
	'test_create_simulation' : test_create_simulation,
	'test_init_problem' : test_init_problem,
	'test_cones_constraint' : test_cones_constraint,
	'test_end_reached_constraint' : test_end_reached_constraint,
	'test_init_constraints' : test_init_constraints,
	'test_initial_guess_naive' : test_initial_guess_naive,
	'test_init_objective' : test_init_objective,
	'test_optimize' : test_optimize,
	'test_end_reached_minus_plus_constraint' : test_end_reached_minus_plus_constraint,
	'test_com_kinematic_constraint' : test_com_kinematic_constraint
}

def run_tests():
	for name, fun in tests_run.iteritems():
		print "running tests: " + name
		fun()		
		
run_tests()
