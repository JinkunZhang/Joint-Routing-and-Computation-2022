#! /usr/bin/env python
from bintrees import AVLTree

class PriorityStruct(object):
   """ A priority data structure built over AVL trees.  

   This structure allows for fast insertion, deletion, membership, priority, and extrema access. Implementation uses an AVL tree for fast insertion, deletion, and extrema access, and a dictionary for fast membership testing and priority computation.  

   """
 
   def __init__(self):
       	""" Initialize an empty priority data structure."""
	self._counter = 0L
       	self._T = AVLTree()
       	self._data = {}

   def __len__(self):
	""" Returns the total number of items stored in the structure."""
	return len(self._data)

   def is_empty(self):
	"""True if empty."""
       	return len(self._data)==0
 
   def clear(self):
	"""Delete all items in data structure."""
       	self._data = {}
       	self._T = AVLTree()

   def __str__(self):
	"""Maps to string, printing tuples of (item,priority)."""
	return  "PriorityStruct("+str( [ (item, priority[0]) for item, priority in self._data.iteritems() ] ) + ")"


   def __iter__(self):
	return self._data.__iter__()
	
   def __contains__(self,item):
       """Membership testing. """
       return item in self._data


   def insert(self, item, priority = float("-inf")):
	""" Add an item to the data structure.  

	If the item is not present, it is added with the specified priority.  An internal integer field is added to the priority to break ties between items with the same priority, and allow unique identification of a node in the AVL tree (used, e.g., for node removal).  If the item is already present, insert simply updates its priority.
	"""
       	if item not in self._data: 	
          self._T[(priority,self._counter)] = item
	  self._data[item] = (priority,self._counter)
	  self._counter += 1L
   	else:  
	  (old_priority,old_id) = self._data[item]
	  del self._T[(old_priority,old_id)]
	  self._T[(priority,old_id)] = item
	  self._data[item] = (priority,old_id)	  


   def remove(self,item):
        """Remove an item from the data structure.  
	
	The return value is a tuple (priority,item) containing the removed item's priority followed by the item. If the item is not present, None is returned.
	"""
	if item not in self._data:
	   return None
	else:
	   key = self._data[item]
	   del self._T[key]
	   del self._data[item] 
	   return (key[0],item)
	
   def min_priority(self):
	"""Return item of lowest priority """
	( (priority,_id), item ) = self._T.min_item()
	return (priority,item)

   def max_priority(self):
	"""Return item of highest priority. """
	( (priority,_id), item ) = self._T.max_item()
	return (priority,item)


   def pop_insert(self,item,priority = float("-inf")):
	"""Remove the item with lowest priority, and then insert new item. """
        (low_priority,low_item) = self.min_priority()
	self.remove(low_item)
	self.insert(item,priority)
	return (low_priority,low_item)

   def insert_pop(self,item,priority= float("-inf")):
	"""Insert item, and then remove the item with lowest priority. """
       	self.insert(item,priority)
	(low_priority,low_item) = self.min_priority()
	self.remove(low_item)
	return (low_priority,low_item)

   def priority(self, item):
	"""Return an item's priority """
	return self._data[item][0]
	

