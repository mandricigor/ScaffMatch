import logging
import networkx as nx
import pickle




class Matcher(object):
    
    _settings = None
    
    def __init__(self, graph):
        self._IGORgraph = graph
    
    def set_settings(self, settings):
        self._settings = settings
        
    def add_one_two(self, matching_graph):
        nodes = set([node[:-2] for node in matching_graph.nodes()])
        for node in nodes:
            matching_graph.add_edge(node + "_1", node + "_2", weight=0)
        return matching_graph
    
    
    def multikeysort(self, items, columns):
        from operator import itemgetter
        comparers = [ ((itemgetter(col[1:].strip()), -1) if col.startswith('-') else (itemgetter(col.strip()), 1)) for col in columns]
        def comparer(left, right):
            for fn, mult in comparers:
                result = cmp(fn(left), fn(right))
                if result:
                    return mult * result
            else:
                return 0
        return sorted(items, cmp=comparer)
    
    
        
    def greedy_matching(self, graph, nr_edges=0):
        """This can be used for extra-large genomes"""
        from Queue import PriorityQueue
        pq = PriorityQueue()
        for x, y in graph.edges():
            pq.put((-graph.edge[x][y]["weight"], (x, y)))
        matchings = {}
        while not pq.empty() or len(matchings) < nr_edges:
            el = pq.get()
            w, el = el
            if not matchings.get(el[0]) and not matchings.get(el[1]):
                matchings[el[0]] = el[1]
                matchings[el[1]] = el[0]
        return matchings
        




    def match(self):
        logger = self._settings.get("logger")
        fragsize = self._settings.get("ins_size")
        bundle_threshold = int(self._settings.get("bundle_threshold"))
        matching_type = self._settings.get("matching")
        if matching_type == "max_weight" or matching_type == "backbone":
            matching_function = nx.max_weight_matching
        elif matching_type == "greedy":
            matching_function = self.greedy_matching
        else:
            raise Exception("Unrecognized matching heuristic type: %s" % matching_type)
        matching_graph = self._IGORgraph # original graph

        ourgraph = matching_graph.copy() # 1st copy of the graph
        ourgraph2 = ourgraph.copy()

        nodes = set() # node set
        repeat_nodes = set() # the set of nodes that are deemed to have repeats
        
        # repeats
        for node in ourgraph.nodes():
            if ourgraph.node[node]["cov"] > self._settings.get("mean_cov") + 2.5 * self._settings.get("disp_cov"):
                repeat_nodes.add(node[:-2]) # we are storing here the nodes that are deemed to be repeats

        bad_edges = []
        for x, y in ourgraph.edges():
            if x in repeat_nodes or y in repeat_nodes:
                bad_edges.append((x, y))

        for x, y in bad_edges:
            if ourgraph.has_edge(x, y):
                ourgraph.remove_edge(x, y)

        # remove repeat nodes temporarily
        for node in repeat_nodes:
            if ourgraph.has_node(node + "_1"):
                ourgraph.remove_node(node + "_1")
            if ourgraph.has_node(node + "_2"):
                ourgraph.remove_node(node + "_2")

        # EDGES TO REMOVE DUE TO NOT ENOUGH WEIGHT
        edges_to_remove = []
        for x, y in ourgraph.edges():
            dist = ourgraph.edge[x][y]["dist"]
            weight = ourgraph.edge[x][y]["weight"]
            if weight <= bundle_threshold: # this is the place where we skip edges that are suspicious due to low weight
                edges_to_remove.append((x, y))
        for x, y in edges_to_remove: # and now boom! remove them at all (for now at all, who knows...)
            ourgraph.remove_edge(x, y)



        goodourgraph = ourgraph.copy() # this snapshot will be used later, this is the so-called CLEAN GRAPH, without repeats

        for i in range(1):
            matchings = matching_function(ourgraph) # MAX WEIGHT MATCHING IS BETTER, GUYS!
            price = 0
            for x, y in matchings.iteritems():
                if x < y:
                    price += ourgraph.edge[x][y]["weight"]

            newgraph = nx.Graph() # we need this graph for building the chains!
            for node in ourgraph.nodes():
                newgraph.add_node(node)

            newgraph = self.add_one_two(newgraph) # this graph now contains "double" edges between the strands of the same node
            for x, y in matchings.iteritems():
                count = ourgraph.edge[x][y]['weight']
                if not newgraph.has_edge(x, y):
                    newgraph.add_edge(x, y, weight=count) # populate newgraph with edges from original graph

            """At this point we have a graph which contains matched edges, but DOES NOT contain double edges"""

            mines = []

            chains = [] # now do chains    
            """Now, we split the newgraph into connected components. Each connected component represents
            either a chain, or a cycle."""
            for comp in nx.connected_components(newgraph):
                edges = nx.subgraph(newgraph, comp).edges() # connected component's edges
                gr = nx.Graph()
                for x, y in edges:
                    gr.add_edge(x, y)
                ends = [x for x in gr.nodes() if gr.degree(x) == 1]
                if len(ends) == 2: # this means that we have a chain, that's good, guys!
                    start, end = ends
                    chain = nx.shortest_path(gr, source=start, target=end) # shortest path is the chain we are searching for
                    chains.append(chain)
                else: # we have a cycle
                    edges2 = [(x, y) for (x, y) in edges if x[:-2] != y[:-2]]
                    minimal_edge = edges2[0]
                    for edge in edges2[1:]:
                        if ourgraph.edge[edge[0]][edge[1]]['weight'] < ourgraph.edge[minimal_edge[0]][minimal_edge[1]]['weight']:
                            minimal_edge = edge
                    mines.append(minimal_edge)
                    gr.remove_edge(*minimal_edge)
                    ends = [x for x in gr.nodes() if gr.degree(x) == 1]
                    start, end = ends
                    chain = nx.shortest_path(gr, source=start, target=end)
                    chains.append(chain)
            for mine1, mine2 in mines:
                if ourgraph.has_edge(mine1, mine2):
		    ourgraph.remove_edge(mine1, mine2)


        """And now, we have chains matched by matching"""

        final_graph = nx.DiGraph()

        # it is possible to introduce a check here
        # check how the chains were constructed
        # for example, if a node is "jumpable", 
        # check whether there exists a link
        # between its neighbor nodes


        # start filling out the graph!!!!!!!!!!
 
        NEXT = {}
        PREV = {}
        for chain in chains:
            if len(chain) < 4:
                continue
            chain = list(reversed(chain))
            if len(chain) == 2:
                final_graph.add_node(chain[0][:-2], orien=False)
            elif len(chain) >= 2:
                nodes = [chain[i][:-2] for i in range(len(chain)) if i % 2 == 0]
                chain_no_ones = []
                for i in range(len(chain) - 1):
                    if chain[i][:-2] == chain[i + 1][:-2]:
                        continue
                    chain_no_ones.append((chain[i], chain[i + 1]))
                orients = []
                for i in range(len(chain_no_ones)):
                    chain = chain_no_ones[i]
                    x, y = chain
                    ori = (x[-1], y[-1])
                    if i == 0:
                        if ori == ("1", "2"):
                            final_graph.add_node(x[:-2], orien=False)
                            final_graph.add_node(y[:-2], orien=False)
                            orients.extend([False, False])
                        elif ori == ("2", "1"):
                            final_graph.add_node(x[:-2], orien=True)
                            final_graph.add_node(y[:-2], orien=True)
                            orients.extend([True, True])
                        elif ori == ("1", "1"):
                            final_graph.add_node(x[:-2], orien=False)
                            final_graph.add_node(y[:-2], orien=True)
                            orients.extend([False, True])
                        elif ori == ("2", "2"):
                            final_graph.add_node(x[:-2], orien=True)
                            final_graph.add_node(y[:-2], orien=False)
                            orients.extend([True, False])
                        pair = tuple(sorted([x, y]))
                        distance = matching_graph.edge[pair[0]][pair[1]]["dist"]
                        final_graph.add_edge(x[:-2], y[:-2], dist=max(0, distance))
                        NEXT[x[:-2]] = y[:-2]
                        PREV[y[:-2]] = x[:-2]
                    else:
                        if ori == ("1", "2"):
                            if orients[-1] == False:
                                final_graph.add_node(y[:-2], orien=False) # done
                                orients.append(False)
                            else:
                                final_graph.add_node(y[:-2], orien=False) # this is very strange!!!!!!!!!1
                                orients.append(False)
                        elif ori == ("2", "1"):
                            if orients[-1] == True:
                                final_graph.add_node(y[:-2], orien=True) # done
                                orients.append(True)
                            else:
                                final_graph.add_node(y[:-2], orien=True) # hz ego znaet
                                orients.append(True)
                        elif ori == ("1", "1"):
                            if orients[-1] == False:
                                final_graph.add_node(y[:-2], orien=True) # done
                                orients.append(True)
                            else:
                                final_graph.add_node(y[:-2], orien=True) # Igor, check
                                orients.append(True)
                        elif ori == ("2", "2"):
                            if orients[-1] == True:
                                final_graph.add_node(y[:-2], orien=False) # done
                                orients.append(False)
                            else:
                                final_graph.add_node(y[:-2], orien=False) # Igor, check
                                orients.append(False)
                        pair = tuple(sorted([x, y]))
                        distance = matching_graph.edge[pair[0]][pair[1]]["dist"]
                        final_graph.add_edge(x[:-2], y[:-2], dist=max(0, distance))
                        NEXT[x[:-2]] = y[:-2]
                        PREV[y[:-2]] = x[:-2]



        deci_removed_contigs = set([x[:-2] for x in goodourgraph.nodes() if x[:-2] not in final_graph.nodes()])
        deci_removed_contigs = list(deci_removed_contigs)
        removed_contigs = [x for x in deci_removed_contigs] + list(repeat_nodes)
        #print len(removed_contigs), "THIS MANY NOT MATCHED"   
	#print len(final_graph.nodes()), "THIS MANY IN THE FINAL GRAPH"
        #print len(final_graph.edges()), "THIS MANY EDGES IN THE FINAL GRAPH"

        # first round of construction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if matching_type != "backbone":
	    for i in range(1):
		SLOTS = {}
		for x, y in final_graph.edges():
		    SLOTS[(x, y)] = []
		for contig in removed_contigs:
		    if contig  in repeat_nodes:
			ggoodourgraph = ourgraph2
		    else:
			ggoodourgraph = goodourgraph
		    ones = {}
		    for x, y in ggoodourgraph.edge[contig + "_1"].iteritems():
		    #for x, y in ourgraph2.edge[contig + "_1"].iteritems():
			if x[:-2] in final_graph.nodes():
			    ones[x] = y
		    twos = {}
		    for x, y in ggoodourgraph.edge[contig + "_2"].iteritems():
		    #for x, y in ourgraph2.edge[contig + "_2"].iteritems():
			if x[:-2] in final_graph.nodes():
			    twos[x] = y
		    reject = False
		    slot_array = []
		    for her in ones:
			slot_dict = {}
			strand = int(her[-1])
			support = ourgraph2.edge[contig + "_1"][her]["weight"]
			if her[:-2] not in final_graph.nodes():
			    continue
			else:
			    ori = final_graph.node[her[:-2]]["orien"]
			if ori == True and strand == 2:
			    newori = True
			    direction = "->"
			elif ori == True and strand == 1:
			    newori = False
			    direction = "<-"
			elif ori == False and strand == 2:
			    newori = False
			    direction = "<-"
			elif ori == False and strand == 1:
			    newori = True
			    direction = "->"
	
	
			if direction == "->":
			    try:
				slot = (her[:-2], NEXT[her[:-2]])
				distance = final_graph.edge[her[:-2]][NEXT[her[:-2]]]["dist"]
				dist = ones[her]["dist"]
				aa = matching_graph.node[her[:-2] + "_1"]["width"] > fragsize and matching_graph.node[NEXT[her[:-2]] + "_1"]["width"] > fragsize
				if aa:
				    pass
				else:
				    if dist > distance:
					slot = (NEXT[her[:-2]], NEXT[NEXT[her[:-2]]])
			    except Exception:
				reject = True
				break
			elif direction == "<-":
			    try:
				slot = (PREV[her[:-2]], her[:-2])
				distance = final_graph.edge[PREV[her[:-2]]][her[:-2]]["dist"]
				dist = ones[her]["dist"]
				aa= matching_graph.node[her[:-2] + "_1"]["width"] > fragsize and matching_graph.node[PREV[her[:-2]] + "_1"]["width"] > fragsize
				if aa:
				    pass
				else:
				    if dist > distance:#+ matching_graph.node[contig + "_1"]["width"] > distance + 900:
					slot = (PREV[PREV[her[:-2]]], PREV[her[:-2]])
				dist = max(0, distance - dist)
			    except Exception:
				reject = True
				break
	
			slot_dict["slot"] = slot
			slot_dict["orien"] = newori
			slot_dict["support"] = support
			slot_dict["dist"] = dist
	
			este = False
			for i in range(len(slot_array)):
			    x = slot_array[i]
			    if x["slot"] == slot and x["orien"] == newori:
				x["support"] += support
				este = True
				slot_array[i] = x
	
			if not este:
			    slot_array.append(slot_dict)
	
		    if reject == True:
			continue
	
		    for her in twos:
			slot_dict = {}
			strand = int(her[-1])
			#support = goodourgraph.edge[contig + "_2"][her]["weight"]
			support = ourgraph2.edge[contig + "_2"][her]["weight"]
	
			if her[:-2] not in final_graph.nodes():
			    continue
			else:
			    ori = final_graph.node[her[:-2]]["orien"]
	
			if ori == True and strand == 2: 
			    newori = False
			    direction = "->"
			elif ori == True and strand == 1:
			    newori = True
			    direction = "<-"
			elif ori == False and strand == 2:
			    newori = True
			    direction = "<-"
			elif ori == False and strand == 1:
			    newori = False
			    direction = "->"
	
	
			if direction == "->":
			    try:
				slot = (her[:-2], NEXT[her[:-2]])
				distance = final_graph.edge[her[:-2]][NEXT[her[:-2]]]["dist"]
				dist = twos[her]["dist"]
				aa = matching_graph.node[her[:-2] + "_1"]["width"] > fragsize and matching_graph.node[NEXT[her[:-2]] + "_1"]["width"] > fragsize
				if aa:
				    pass
				else:
				    if dist > distance:#+ matching_graph.node[contig + "_1"]["width"] > distance + 900:
					slot = (NEXT[her[:-2]], NEXT[NEXT[her[:-2]]])
			    except Exception:
				break
			elif direction == "<-":
			    try:
				slot = (PREV[her[:-2]], her[:-2])
				distance = final_graph.edge[PREV[her[:-2]]][her[:-2]]["dist"]
				dist = twos[her]["dist"]
				aa = matching_graph.node[her[:-2] + "_1"]["width"] > fragsize and matching_graph.node[PREV[her[:-2]] + "_1"]["width"] > fragsize
				if aa:
				    pass
				else:
				    if dist > distance:#+ matching_graph.node[contig + "_1"]["width"] > distance + 900:
					slot = (PREV[PREV[her[:-2]]], PREV[her[:-2]])
				dist = max(0, distance - dist)
			    except Exception:
				break
	
			slot_dict["slot"] = slot
			slot_dict["orien"] = newori
			slot_dict["support"] = support
			slot_dict["dist"] = dist
	
	
			este = False
			for i in range(len(slot_array)):
			    x = slot_array[i]
			    if x["slot"] == slot and x["orien"] == newori:
				x["support"] += support
				este = True
				slot_array[i] = x
	
			if not este:
			    slot_array.append(slot_dict)
	
	
	
		    if reject == True:
			continue
	
		    sorted_slot = self.multikeysort(slot_array, ["-support"])
		    if not sorted_slot:
			continue
		    best_slot = sorted_slot[0]
	
		    if len(sorted_slot) > 1:
			second_best_slot = sorted_slot[1]
			if best_slot["support"] == second_best_slot["support"]:
			    continue
	
		    if best_slot["support"] < 0:
			continue
	
	
		    slot = best_slot["slot"]
		    del best_slot["slot"]
		    best_slot["contig"] = contig
		    try:
			SLOTS[slot].append(best_slot)
		    except Exception:
			#print slot, contig
                        pass
	
	
		to_be_inserted = []
	
		for x, y in SLOTS:
		    if len(SLOTS[(x, y)]) == 0:
			continue
		    elif len(SLOTS[(x, y)]) == 1:
			slots = self.multikeysort(SLOTS[(x, y)], ["-support"])
			contig = slots[0]
			node = contig["contig"]
			orien = contig["orien"]
			dist = contig["dist"]
			distance = final_graph.edge[x][y]["dist"]
			final_graph.add_node(node, orien=orien)
			final_graph.add_edge(x, node, dist=dist)
			final_graph.add_edge(node, y, dist=distance)
			final_graph.remove_edge(x, y)
			to_be_inserted.append(node)
		    else:
			slots = self.multikeysort(SLOTS[(x, y)], ["dist"])
			slots = self.multikeysort(SLOTS[(x, y)], ["dist"])
			for sl in slots:
			    final_graph.add_node(sl["contig"], orien=sl["orien"])
			final_graph.add_edge(x, slots[0]["contig"], dist=100)
			final_graph.add_edge(slots[-1]["contig"], y, dist=100)
			for i in range(len(slots) - 1):
			    final_graph.add_edge(slots[i]["contig"], slots[i + 1]["contig"], dist=100)
 
    	#print len(final_graph.nodes()), "THIS MANY IN THE FINAL GRAPH"
        #print len(final_graph.edges()), "THIS MANY EDGES IN THE FINAL GRAPH"

        for node in final_graph.nodes():
            if "orien" in final_graph.node[node].keys():
                pass
            else:
                final_graph.node[node]["orien"] = False
            final_graph.node[node]['width'] = matching_graph.node[node + "_1"]["width"]

        for x in matching_graph.nodes():
            if not final_graph.has_node(x[:-2]):
                final_graph.add_node(x[:-2], orien=False, width=matching_graph.node[x[:-2] + "_1"]["width"])
           
        wdir = self._settings.get("scaff_dir")
 
        #print len(final_graph.edges()), "EDGES"
        #print len(final_graph.nodes()), 'NODES'
        a = set()
        for x, y in final_graph.edges():
            if x < y:
                a.add((x, y))
            else:
                a.add((x, y))
        #print len(a)

        nx.write_gpickle(final_graph, wdir + "/final_graph.cpickle")
    
