#! /usr/bin/env python
'''
    Cache Utilities

'''

from abc import ABCMeta,abstractmethod
from Priority import PriorityStruct
from EWMA import EWMA
import random
import numpy as np
from helpers import projectToSimplex,constructDistribution

def addDictionaries(d1,d2):
	"""
	Add the elements of two dictionaries. Missing entries treated as zero.
	"""
	lsum =  [( key, d1[key]+d2[key]) for key in d1 if key in d2]
	lsum +=  [( key, d1[key]) for key in d1 if key not in d2]
	lsum +=  [( key, d2[key]) for key in d2 if key not in d1]
	return dict(lsum)

def scaleDictionary(d,s):
	return dict( [ (key,s*d[key]) for key in d  ]  )

def dictL1(d):
	return sum( (np.abs(x) for x in d.values()) )


def dictL2(d):
	return np.sqrt(sum( (x*x for x in d.values()) ))


class Cache:
    """ A generic cache object

    An abstract cache object is determined by a capacity (integer) and an id, used to identify the cache. The cache contents are represented as a set

    Attributes:
    	_id : the unique cache id
        _capacity: the cache size
        _contents: the set containing the cache contents
    """

    __metaclass__ = ABCMeta

    def __init__(self,capacity, _id):
    	""" Initialize an empty cache with given capacity and id.
	"""
        self._id = _id
	self._capacity = capacity
	self._contents = set([])

    def __str__(self):
	"""String representation.
	"""
        return `self._id`+":"+str(self._contents)     

    # Storage Operations 
    def __len__(self):
	return len(self._contents)

    def is_empty(self):
	return len(self) == 0

    def __contains__(self,item):
	return item in self._contents
    
    def __iter__(self):
	return self._contents.__iter__()

    @abstractmethod
    def add(self,item):
	pass


    @abstractmethod
    def remove(self,item): 
	pass

    @abstractmethod
    def clear(self):
	pass


class PriorityCache(Cache): 
    ''' A priority cache. It can be used to implement LRU (priority=time accessed) and LFU (priority = number of accesses).  '''
    def  __init__(self,capacity, _id ):
	Cache.__init__(self,capacity,_id)
	self._contents = PriorityStruct()

    def __contains__(self,item):
	return item in self._contents

    def add(self,item,priority=float("-inf")):
	    if len(self._contents) < self._capacity or item in self._contents :
	       self._contents.insert(item,priority)
	       return None
	    else:
	       return self._contents.insert_pop(item,priority) #warning: item inserted will be popped if it is the lowest priority element

    def add_pop_if_needed(self,item,priority=float("-inf")):
	    if len(self) < self._capacity or item in self :
	       self._contents.insert(item,priority)
	       return None
	    else:
	       return self._contents.pop_insert(item,priority) #warning: item inserted will be added even if it is the lowest priority element


    def remove(self,item):
	return self._contents.remove(item)
    
    def clear(self):
	self._contents= PriorityStruct()
  
    def priority(self,item):
	return self._contents.priority(item) 



class EWMACache(Cache):
    """ An EWMA cache object """

    

    def __init__(self,capacity, _id,beta):
    	""" Initialize an empty cache with given capacity and id.
	"""
        self._id = _id
	self._capacity = capacity
        self._beta = beta
        self._ewma = EWMA(self._beta)

    def __str__(self):
	"""String representation.
	"""
        return `self._id`+'('+`self._capacity` +'):'+ str(self._ewma.topK(self._capacity)) +'EWMA: '+ str(self._ewma)    

    def __len__(self):
	return min(len(self.ewma),self._capacity)


    def __contains__(self,item):
	return item in self._ewma.topK(self._capacity)

    def __iter__(self):
	return self._ewma.topK(self._capacity).__iter__()

    def add(self,item,value,time): #THIS DOES NOT REALLY ADD AN ITEM, IT UPDATES ITS EWMA ESTIMATE
        self._ewma.add(item,value,time)
       

    def remove(self,item): 
	return self._ewma.remove(item)

    def clear(self):
	del self._ewma
        self.ewma = EWMA()

class LMinimalCache(Cache):
    """ A cache that minimizes the relaxation L. """

    def __init__(self,capacity,_id,gamma = 0.1,expon = 0.5,interpolate = False):
    	""" Initialize an empty cache with given capacity and id.
	"""
        self._id = _id
	self._capacity = capacity
	self._gamma = gamma
	self._expon = expon
        self._state = {}
	self._past_states = {}
	self._grad = {}
	self._interpolate = interpolate
	self._last_shuffle_time = 0.0
	self._cache = set()
	self._k = 0

    def __str__(self):
	"""String representation.
	"""
        return `self._id`+'('+`self._capacity` +'):'+ str(self._cache) 

    def __len__(self):
	return len(self._cache)


    def __contains__(self,item):
	return item in self._cache

    def __iter__(self):
	return self._cache.__iter__()

    
    def add(self,item): 
        self._cache.add(item)

    def remove(self,item): 
	self._cache.remove(item)
	pass

    def clear(self):
	self.cache.clear()

    def interpolate(self,states,k):
	 
        if k<=3:
	     return states[k]
	sum_so_far = 0.0
	dict_so_far = {}

	for  ell in range(int(np.ceil(k*0.5)  ),k):
	     gain = self._gamma*(ell** (-self._expon))
	     dict_so_far = addDictionaries(dict_so_far,scaleDictionary(self._past_states[ell],gain))
	     sum_so_far += gain

	return scaleDictionary(dict_so_far,1.0/sum_so_far)     


    def updateGradient(self,item,value):
	if item in self._grad:
	    self._grad[item] += value
	else:
	    self._grad[item] = value
	#print self._grad[item]
	
    def shuffle(self,time):
	T = time - self._last_shuffle_time
	#print "Time passed", T
	k = self._k+1

	#print "About to update state:", self._state,
	#print "with grad",self._grad
	if T>0.0:
	    self._grad = scaleDictionary(self._grad,1.0/T)
	
	if dictL2(self._grad) > 100*self._capacity:
		rescale = scaleDictionary(self._grad,1.0/dictL2(self._grad))
	else:
		rescale = self._grad

	for item in rescale:
	    if item in self._state:
		self._state[item] += self._gamma * (k** (-self._expon)) *rescale[item]
 	    else: 
		self._state[item] = self._gamma * (k** (-self._expon)) *rescale[item]
	
	if self._state == {} or self._state == None:
	    return
	#print "After grad", self._state, 'now will project'

	# Project. Solution maybe approximate, may have negative values, etc.	
	self._state, res = projectToSimplex(self._state, self._capacity)
	#print "After projection ",self._state,"now will clean up"

	#cleanup, make kosher
	#print "Cleaning up",
	self._state = dict( [(x,min(self._state[x],1.0)) for x in self._state if self._state[x]>1.e-18])
        
        #print "After cleanup:",self._state
	
        
	self._past_states[k] = dict(self._state)

	if self._interpolate:
	    sliding_average = self.interpolate(self._past_states,k)
	else:
	    sliding_average = dict(self._state)

	#print k, dictL1(self._state),dictL1(sliding_average),dictL2(self._grad)

	#if dictL1(sliding_average)>self._capacity:
	#     print self._past_states
        #print "state:",self._state
        #print "sliding:",sliding_average
        #print "past:",self._past_states 

        placements, probs, distr = constructDistribution(sliding_average,self._capacity)
		
	u = random.random()
	sumsofar = 0.0

	for key in probs:
	    sumsofar += probs[key]
	    if sumsofar >= u:
		break
	#key = distr.rvs(size=1)[0]
	#print key

	self._cache = set(placements[key])
        #print self._cache, placements, probs
	
	self._k = k
	self._grad = {}
	self._last_shuffle_time = time

    def state(self,item):
	
	if self._interpolate:
	    if self._k >= 3:
		sliding_average = self.interpolate(self._past_states,self._k)
	    else:
		sliding_average = self._state
	else:
	    sliding_average = self._state

	if item in sliding_average:
	    return sliding_average[item]
	else:
	    return 0.0


if __name__=="__main__":
    cache = LMinimalCache(4,0)
    cache.updateGradient(2,2.0)
    cache.updateGradient(20,2.0)
    cache.updateGradient(21,2.0)
    cache.updateGradient(21,3.0)
    cache.updateGradient(1,3.0)
    cache.updateGradient(1,30.0)
    cache.updateGradient(12,30.0)
    cache.updateGradient(112,30.0)
    cache.updateGradient(11,3.0)
    cache.shuffle(1.0)

    for i in cache:
	print i

