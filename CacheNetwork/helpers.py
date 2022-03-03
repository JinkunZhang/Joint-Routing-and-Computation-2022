from cvxopt import solvers
from cvxopt import matrix
import numpy as np
from scipy.stats import rv_discrete

def projectToSimplex(d,cap):
    keys, vals = zip( *[ (key,d[key]) for key in d] )

    #print "Projecting:",d

    n = len(vals)
    q = -matrix(vals)
    P = matrix(np.eye(n))

    G = matrix(np.concatenate( (np.eye(n), -np.eye(n), np.ones((1,n)) ) ))
    
    h = matrix( n*[1.0] + n*[0.0] +[cap]  )   
     
    solvers.options['show_progress'] = False
    res = solvers.qp(P,q,G,h)


   
    sol = res['x']
    return dict(zip(keys,sol)), res

def constructDistribution(d,cap):
    epsilon = 1.e-3

    #Remove very small values, rescale the rest
    dd = dict( (key,d[key]) for key in d if d[key]>epsilon  )
    keys, vals = zip( *[ (key,d[key]) for key in dd] )
    ss = sum(vals)
    vals = [ val/ss*cap for val in vals]
    dd = dict(zip(keys,vals))

    intvals = [int(np.round(x/epsilon)) for x in vals]
    intdist = int(1/epsilon)
    intdd = dict(zip( keys,intvals  ))

    s = {}
    t = {}
    taus = [] 
    sumsofar = 0
    for item in keys:
	s[item] = sumsofar
	t[item] = sumsofar + intdd[item]
	taus.append(t[item] % intdist )
	sumsofar = t[item]

    #print s,t,taus
    taus = sorted(set(taus))
    #print taus

    if intdist not in taus:
	taus.append(intdist)

    placements = {}
    prob = {}

    for i in range(len(taus)-1):
	x = []
	t_low = taus[i]
	t_up = taus[i+1]

	diff = t_up -t_low

        for ell in range(int(cap)):
    	   lower = ell*intdist + t_low
           upper = ell*intdist + t_up
	   for item in keys:
		#print lower,upper,' inside ', s[item],t[item], '?',
		if lower>=s[item] and upper<=t[item]:
		    x.append(item)
		#    print ' yes'
		#else: print ' no'
	prob[i] = 1.*diff/intdist 
        placements[i] = x

    totsum = np.sum(prob.values())
    if not np.allclose(totsum,1):
	for i in prob:
	    prob[i] = 1.*prob[i]/totsum
	# round to 1
	

    #    print "Placements ",placements,"with prob",prob,"summing to",np.sum(prob.values())
 
    
    return placements,prob, rv_discrete(values=(prob.keys(),prob.values()))	 

