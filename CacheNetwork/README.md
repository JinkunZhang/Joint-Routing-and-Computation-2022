A Caching Network Simulator
===========================

Simulator of experiments presented in ["Adaptive Caching Networks with Optimality Guarantees"](http://ieeexplore.ieee.org/document/8279650/), S. Ioannidis and E. Yeh, IEEE/ACM Transactions of Networking, 2018. Please cite this paper if you intend to use this code for your research.

Usage
-----
An example execution is as follows:

	python CacheNetwork.py outputfile --graph_type balanced_tree --graph_degree 2 --cache_type LRU --graph_size 60 --min_capacity 3 --max_capacity 3 --query_nodes 20 --catalog_size 300 --demand_size 1000 --max_weight 100 --time 5000


This simulates an experiment over a balanced_tree topology, where each node has 2 children, using path replication with LRU caches. The graph size contains 60 nodes of capacity 3 (excluding capacity for designated source content), storing items from a catalog of 30 items. Only 20 of these nodes generate queries. There are 1000 types of requests (items and paths followed). The maximum weight of an edge is 1000, and the simulation runs for 5000 time units.

More execution examples can be found in file `execution_examples`. 

Dependencies
------------
Please use `networkx` version 1.11. The code is not compatible with later versions.


Simulator Overview
------------------


A cache network comprises a weighted graph and a list of demands. Each node in the graph is associated with a cache of finite capacity.
NetworkCaches must support a message receive operation, that determines how they handle incoming messages.

The cache network handles messaging using `simpy` stores and processes. In partiqular, each cache, edge and demand is associated with a 
Store object, that receives, stores, and processes messages from simpy processes.

In more detail:
* Each demand is associated with two processes, one that generates new queries, and one that monitors and logs completed queries (existing only for logging purposes)
* Each cache/node is associated with a process that receives messages, and processes them, and produces new messages to be routed, e.g., towards neigboring edges
* Each edge is associated with a process that receives messages to be routed over the edge, and delivers them to the appropriate target node.
During "delivery", messages are (a) delayed, according to configuration parameters, and (b) statistics about them are logged (e.g., number of hops, etc.)
     
Finally, a global monitoring process computes the social welfare at poisson time intervals.



Command-line arguments
----------------------
Several program parameters can be controlled from the command line.


	usage: CacheNetwork.py [-h] [--max_capacity MAX_CAPACITY]
				       [--min_capacity MIN_CAPACITY] [--max_weight MAX_WEIGHT]
				       [--min_weight MIN_WEIGHT] [--max_rate MAX_RATE]
				       [--min_rate MIN_RATE] [--time TIME] [--warmup WARMUP]
				       [--catalog_size CATALOG_SIZE]
				       [--demand_size DEMAND_SIZE]
				       [--demand_change_rate DEMAND_CHANGE_RATE]
				       [--demand_distribution {powerlaw,uniform}]
				       [--powerlaw_exp POWERLAW_EXP]
				       [--query_nodes QUERY_NODES]
				       [--graph_type {erdos_renyi,balanced_tree,hypercube,cicular_ladder,cycle,grid_2d,lollipop,expander,hypercube,star,barabasi_albert,watts_strogatz,regular,powerlaw_tree,small_world,geant,abilene,dtelekom,servicenetwork}]
				       [--graph_size GRAPH_SIZE] [--graph_degree GRAPH_DEGREE]
				       [--graph_p GRAPH_P] [--random_seed RANDOM_SEED]
				       [--debug_level {INFO,DEBUG,WARNING,ERROR}]
				       [--cache_type {LRU,FIFO,LFU,RR,EWMAGRAD,LMIN}]
				       [--query_message_length QUERY_MESSAGE_LENGTH]
				       [--response_message_length RESPONSE_MESSAGE_LENGTH]
				       [--monitoring_rate MONITORING_RATE]
				       [--interpolate INTERPOLATE] [--beta BETA]
				       [--gamma GAMMA] [--expon EXPON] [--T T]
				       outputfile

		Simulate a Network of Caches

		positional arguments:
		  outputfile            Output file

		optional arguments:
		  -h, --help            show this help message and exit
		  --max_capacity MAX_CAPACITY
					Maximum capacity per cache (default: 2)
		  --min_capacity MIN_CAPACITY
					Minimum capacity per cache (default: 2)
		  --max_weight MAX_WEIGHT
					Maximum edge weight (default: 100.0)
		  --min_weight MIN_WEIGHT
					Minimum edge weight (default: 1.0)
		  --max_rate MAX_RATE   Maximum demand rate (default: 1.0)
		  --min_rate MIN_RATE   Minimum demand rate (default: 1.0)
		  --time TIME           Total simulation duration (default: 1000.0)
		  --warmup WARMUP       Warmup time until measurements start (default: 0.0)
		  --catalog_size CATALOG_SIZE
					Catalog size (default: 100)
		  --demand_size DEMAND_SIZE
					Demand size (default: 1000)

		  --demand_change_rate DEMAND_CHANGE_RATE
		                        Demand change rate (default: 0.0)
		  --demand_distribution {powerlaw,uniform}
					Demand distribution (default: powerlaw)
		  --powerlaw_exp POWERLAW_EXP
					Power law exponent, used in demand distribution
					(default: 1.2)
		  --query_nodes QUERY_NODES
					Number of nodes generating queries (default: 100)
		  --graph_type {erdos_renyi,balanced_tree,hypercube,cicular_ladder,cycle,grid_2d,lollipop,expander,hypercube,star,barabasi_albert,watts_strogatz,regular,powerlaw_tree,small_world,geant,abilene,dtelekom,servicenetwork}
					Graph type (default: erdos_renyi)
		  --graph_size GRAPH_SIZE
					Network size (default: 100)
		  --graph_degree GRAPH_DEGREE
					Degree. Used by balanced_tree, regular,
					barabasi_albert, watts_strogatz (default: 4)
		  --graph_p GRAPH_P     Probability, used in erdos_renyi, watts_strogatz
					(default: 0.1)
		  --random_seed RANDOM_SEED
					Random seed (default: 123456789)
		  --debug_level {INFO,DEBUG,WARNING,ERROR}
					Debug Level (default: INFO)
		  --cache_type {LRU,FIFO,LFU,RR,EWMAGRAD,LMIN}
					Networked Cache type (default: LRU)
		  --query_message_length QUERY_MESSAGE_LENGTH
					Query message length (default: 0.0)
		  --response_message_length RESPONSE_MESSAGE_LENGTH
					Response message length (default: 0.0)
		  --monitoring_rate MONITORING_RATE
					Monitoring rate (default: 1.0)
		  --interpolate INTERPOLATE
					Interpolate past states, used by LMIN (default: False)
		  --beta BETA           beta used in EWMA (default: 1.0)
		  --gamma GAMMA         gamma used in LMIN (default: 0.1)
		  --expon EXPON         exponent used in LMIN (default: 0.5)
		  --T T                 Suffling period used in LMIN (default: 5.0)


Using `pydoc`
------------

Type `pydoc CacheNetwork` to see documentation for the `CacheNetwork` class; same applies to other `.py` files.

