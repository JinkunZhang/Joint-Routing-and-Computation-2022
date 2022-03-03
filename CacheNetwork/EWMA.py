#! /usr/bin/env python
'''
    EWMA Counter

'''

import numpy as np
from operator import itemgetter

class TimeException(Exception):
    def __init__(self, previous_time, current_time):
        self.pt = previous_time
	self.ct = current_time
	 
    def __str__(self):
	return 'Current time '+str(self.ct)+' occurs before previous time '+str(self.pt)


class EWMA:
    """ An EWMA data structure.  """

    def __init__(self,beta):
	self.beta = beta
	self.time = 0
	self.v = {}	

    def update(self,time):
	elapsed_time = time-self.time
	if elapsed_time<0:
	   raise TimeException(self.time,time)
	self.v = dict([ (key, value*np.exp(-self.beta*elapsed_time)) for key,value in self.v.iteritems() ])
	self.time = time

    def __len__(self):
	return len(self.v)

    def __iter__(self):
	return self.v.__iter__()

    def __getitem__(self,key):
	return self.v[key]

    def __setitem__(self,key,value):
	self.v[key] = value

    def __delitem__(self,key):
	del self.v[key]

    def time(self):
	return self.time

    def add(self,item,value,time): 
	self.update(time)
	if item in self.v:
	    self.v[item] += self.beta*value
	else:
	    self.v[item] = self.beta*value

    def sorted_itemvalues(self):
	return sorted( self.v.items(), key = itemgetter(1,0), reverse=True)

    def topK(self,k):
	topkeys= [key for (key,value) in self.sorted_itemvalues()]
	if k<len(topkeys):
	    return topkeys[:k]
	else:
	    return topkeys

    def remove(self,item):
	if item in self.v:
	   val = self.v[item]
	   del self.v[item]
           return (item,val)
	else:
	   return None

    def __str__(self):
	return str(self.time)+': '+','.join ([ str(key)+':'+str(value) for key,value in  self.sorted_itemvalues()])

if __name__ == "__main__":
    counter = EWMA(2.0)
    counter.add(1,5.0,1)
    print counter
    counter.update(1.5)
    print counter
    counter.add(2,4,2)
    print counter
    counter.add(1,3,3)
    print counter
    counter.add(4,10,3.5)
    print counter
    print map(str,counter.topK(2))
    counter.update(4)
    print counter
    print map(str,counter.topK(2))
    counter.add('b',13,5)
    print counter
    print map(str,counter.topK(2))
    counter.update(2)
