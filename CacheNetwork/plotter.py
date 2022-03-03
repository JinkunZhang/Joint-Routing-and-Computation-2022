import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
import argparse
import colormaps as cmaps

methods = ['LRU','LFU','FIFO','RR','EWMAGRAD','LMIN']
pretty_methods = ['LRU','LFU','FIFO','RR','GRD','PGA']
cmap = dict(zip(methods,pretty_methods))
all_graphs = ['cycle','lollipop','grid_2d','balanced_tree','hypercube','expander','erdos_renyi','regular','watts_strogatz','small_world','barabasi_albert','cicular_ladder','star', 'powerlaw_tree','geant','abilene','dtelekom' ] 
threshold = 1000.0

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

colors =['b', 'g', 'r', 'c' ,'m' ,'y' ,'k' ,'w']


ecg_avg = dict(  (  (g,{}) for g in all_graphs   )  ) 


if __name__=='__main__':
    caches = set()
    graphs = set()

    parser = argparse.ArgumentParser(description='Process output files.')
    parser.add_argument('filenames', metavar='N', type=str, nargs='+',
                   help='pickled file to be processed')
    parser.add_argument('--legend', dest='legend', action='store_true')
    parser.add_argument('--no-legend', dest='legend', action='store_false')
    parser.set_defaults(legend=True)
    parser.add_argument('--hatch', dest='legend', action='store_true')
    parser.add_argument('--no-hatch', dest='legend', action='store_false')
    parser.set_defaults(hatch=True)



    myargs = parser.parse_args()
    seen_so_far = set()

    

    for filename in myargs.filenames:
 
	print 'Processing...',
	try:
	    with open(filename,'rb') as f:
	    	args,construct_stats,optimal_stats,demand_stats,node_stats,network_stats = pickle.load(f)
	except (IOError, OSError ) as e:
	    print e
	    continue

        
    	cache = cmap[args.cache_type]
    	graph = args.graph_type	
	
	if cache == 'PGA':
	    cache += str(int(args.T))
            if args.interpolate:
		print 'Interpolated, skipping'
		continue
		cache +='I'
	
	_id = graph+'_'+cache

	while _id in seen_so_far:
	    _id +='_'
 	seen_so_far.add(_id)

	print graph,cache,
 
    	F = optimal_stats['F']
    	L = optimal_stats['L']

    	times1 = sorted(network_stats['fun'].keys())[1:]
    	times2 = sorted(network_stats['demand'].keys())[1:]
    	ecgs = [ network_stats['fun'][t][0] for t in times1 ]
        tot = network_stats['fun'][times1[0]][3]  
   	tags = [ tot-network_stats['demand'][x]['weight'] for x in times2 ]
    
   
    	fig =plt.figure(figsize=(4*1.61803398875,4))
    	ax = fig.add_subplot(1,1,1)
    	names= ['F(Y*)','ECG','TACG'] 
    	forms = ['k:','rs-','b^-']
    	optline, =ax.plot(times1,F*np.ones(len(times1)),forms[0])
    	ecline, =ax.plot(times1,ecgs,forms[1])
    	taline, =ax.plot(times2,tags,forms[2])
    	ax.set_ylabel('Gain')
    	ax.set_xlabel('Time')
	ax.set_ylim(0,100)
    	plt.legend([optline,ecline,taline],names)
    	ax.set_title(graph.replace('_','-')+' '+cache )
    	fig.savefig(_id+'.pdf',bbox_inches='tight')
    	plt.close(fig)

	kosher = [ ecgs[i]  for i in range(len(ecgs))  if times1[i]>threshold ]
	
	ecg_avg[graph][cache] = sum( kosher )/ len(kosher) / F
 
	print 	ecg_avg[graph][cache]
        caches.add(cache)
	graphs.add(graph)    

    graphs = [  x for x in all_graphs if x in graphs  ]
    caches = ['LRU', 'LFU','FIFO', 'RR', 'GRD', 'PGA1','PGA10','PGA20'] 
    hatches = ['////', '/', '\\', '\\\\', '-', '--', '+', '']
    hatchmap = dict(zip(caches,hatches))
    ecg_list = {}
    for cache in caches:
	ecg_list[cache]  =[ 0.0 if cache not in ecg_avg[graph] else ecg_avg[graph][cache]    for graph in graphs]

    fig =plt.figure(figsize = (3*(len(graphs)+1), 3*1.0))
    ax = fig.add_subplot(1,1,1)
    ind = 1.0*np.arange(len(graphs))
    width = 0.9/len(caches)
    print width
    rects = []
    newind = ind
    i = 0
    for cache in caches:
	print newind
        if myargs.hatch:
		rects += ax.bar(newind,ecg_list[cache],width,color = cmaps.plasma(0.5+ 0.5* i / len(caches) ),hatch=hatchmap[cache],linewidth=2,label = cache)
        else:
		rects += ax.bar(newind,ecg_list[cache],width,color = cmaps.plasma(0.5+ 0.5* i / len(caches) ),linewidth=1,label = cache)
        newind += width 
	i += 1
    ax.set_xticks(ind-width*len(caches)/2 )
    ax.set_xticklabels([g.replace("_","-") for g in  graphs])
    ax.grid(False)
    if myargs.hatch:
	hatchstring = "_withhatch"
    else:
	hatchstring = ""
    if myargs.legend:
	#lgd=ax.legend( [rects[i*len(graphs)]  for i in range(len(caches)) ], caches,loc = 'center left',bbox_to_anchor =(1.0,0.5))
    	lgd=ax.legend( [rects[i*len(graphs)]  for i in range(len(caches)) ], caches,loc=3, bbox_to_anchor=(0., 1.02, 1., .102),mode='expand',ncol=len(caches),borderaxespad=0.)
    	#plt.tight_layout()
    	fig.savefig('barplot%s%s.pdf' % ('_'.join(graphs),hatchstring),bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
    	fig.savefig('barplot%s%s.pdf' % ('_'.join(graphs),hatchstring), bbox_inches='tight')
    plt.close(fig)


