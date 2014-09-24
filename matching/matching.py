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
            matching_graph.add_edge(node + "_1", node + "_2", weight=-1000)
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
        from Queue import PriorityQueue
        pq = PriorityQueue()
        for x, y in graph.edges():
            pq.put((-graph.edge[x][y]["weight"], (x, y)))
        matchings = {}
        while not pq.empty() or len(matchings) < nr_edges:
            el = pq.get()
            w, el = el
            if matchings.get(el[0], "vafli") == "vafli" and matchings.get(el[1], "vafli") == "vafli":
                matchings[el[0]] = el[1]
                matchings[el[1]] = el[0]
        return matchings
        




    def match(self):
        fragsize = self._settings.get("ins_size")
        matching_graph = self._IGORgraph
        ourgraph = matching_graph.copy()   
        removed_nodes = {}
        nodes = set()
        repeat_nodes = set()
        for node in ourgraph.nodes():
            if ourgraph.node[node]["cov"] > self._settings.get("mean_cov") + 2 * self._settings.get("disp_cov"):
                print "REPEAT:", node
                repeat_nodes.add(node[:-2])
       
        edges_to_remove = []
        for x, y in ourgraph.edges():
            dist = ourgraph.edge[x][y]["dist"]
            weight = ourgraph.edge[x][y]["weight"]
            if weight < 1:
                edges_to_remove.append((x, y))


        for x, y in edges_to_remove:
            ourgraph.remove_edge(x, y)

        #EXPERIMENT
        edges_to_delete = []
        for x, y in ourgraph.edges():
            if x.split(":")[0] != y.split(":")[0]:
                edges_to_delete.append((x, y))
        for x, y in edges_to_delete:
            if ourgraph.has_edge(x, y):
                ourgraph.remove_edge(x, y)




        for node in repeat_nodes:
            if ourgraph.has_node(node + "_1"):
                ourgraph.remove_node(node + "_1")
            if ourgraph.has_node(node + "_2"):
                ourgraph.remove_node(node + "_2")

        bad_edges = []
        for x, y in ourgraph.edges():
            if x in repeat_nodes or y in repeat_nodes:
                bad_edges.append((x, y))

        for x, y in bad_edges:
            if ourgraph.has_edge(x, y):
                ourgraph.remove_edge(x, y)

        goodourgraph = ourgraph.copy() # this snapshot will be used later
        for node in removed_nodes: # remove these nodes temporarily from the graph
            ourgraph.remove_node(node)
	matchings = self.greedy_matching(ourgraph)
        #matchings = nx.max_weight_matching(ourgraph, maxcardinality=True)
        newgraph = nx.Graph() # we need this graph for building the chains!
        for node in ourgraph.nodes():
            newgraph.add_node(node)
        newgraph = self.add_one_two(newgraph)
        for x, y in matchings.iteritems():
            count = ourgraph.edge[x][y]['weight']
            if not newgraph.has_edge(x, y):
                newgraph.add_edge(x, y, weight=count) # populate newgraph with edges from original graph
        chains = [] # now do chains    
        for comp in nx.connected_components(newgraph):
            edges = nx.subgraph(newgraph, comp).edges()
            gr = nx.Graph()
            for x, y in edges:
                gr.add_edge(x, y)
            ends = [x for x in gr.nodes() if gr.degree(x) == 1]
            if len(ends) == 2:
                start, end = ends
                chain = nx.shortest_path(gr, source=start, target=end)
                chains.append(chain)
            else: # we have a cycle
                edges2 = [(x, y) for (x, y) in edges if x[:-2] != y[:-2]]
                minimal_edge = edges2[0]
                for edge in edges2[1:]:
                    if ourgraph.edge[edge[0]][edge[1]]['weight'] < ourgraph.edge[minimal_edge[0]][minimal_edge[1]]['weight']:
                        minimal_edge = edge
                gr.remove_edge(*minimal_edge)
                ends = [x for x in gr.nodes() if gr.degree(x) == 1]
                print ends
                start, end = ends
                chain = nx.shortest_path(gr, source=start, target=end)
                chains.append(chain)
        ######################################################### join some of the chains now ##############
        mbchains = []
        remove_from_chains = []
        for i in range(len(chains)):
            ch = chains[i]
            if len(ch) <= 2:
                mbchains.extend(ch)
                remove_from_chains.append(ch)
            else:
                mbchains.extend([ch[0], ch[-1]])
    
        chains = [x for x in chains if x not in remove_from_chains]
        mbchains = mbchains + removed_nodes.keys()
    
    
        ourgraph2 = goodourgraph.copy()
        nodes_to_remove = []
        for node in goodourgraph.nodes():
            if node not in mbchains:
                nodes_to_remove.append(node)
        for node in nodes_to_remove:
            ourgraph2.remove_node(node)
    
        ourgraph2 = self.add_one_two(ourgraph2)
    
    
        #matchings = nx.max_weight_matching(ourgraph2, maxcardinality=True)
        matchings = self.greedy_matching(ourgraph2)
    
        newgraph = nx.Graph()
    
    
        for node in ourgraph2.nodes():
            newgraph.add_node(node)
    
    
        newgraph = self.add_one_two(newgraph)
    
        for x, y in matchings.iteritems():
            if x[:-2] == y[:-2]:
                continue
            count = goodourgraph.edge[x][y]['weight']
            if not newgraph.has_edge(x, y):
                newgraph.add_edge(x, y, count=count)
        chains2 = []
        for comp in nx.connected_components(newgraph):
            edges = nx.subgraph(newgraph, comp).edges()
            gr = nx.Graph()
            for x, y in edges:
                gr.add_edge(x, y)
            ends = [x for x in gr.nodes() if gr.degree(x) == 1]
            if len(ends) == 2:
                start, end = ends
                chain = nx.shortest_path(gr, source=start, target=end)
                chains2.append(chain)
            else: # we have a cycle
                edges2 = [(x, y) for (x, y) in edges if x[:-2] != y[:-2]]
                minimal_edge = edges2[0]
                for edge in edges2[1:]:
                    if ourgraph2.edge[edge[0]][edge[1]]['weight'] < ourgraph2.edge[minimal_edge[0]][minimal_edge[1]]['weight']:
                        minimal_edge = edge
                gr.remove_edge(*minimal_edge)
                ends = [x for x in gr.nodes() if gr.degree(x) == 1]
                start, end = ends
                chain = nx.shortest_path(gr, source=start, target=end)
                chains2.append(chain)
        real_chains = []
        for chain in chains2:
            chain_ends = []
            for ch in chains:
                chain_ends.extend([ch[0][:-2], ch[-1][:-2]])
            if chain[0][:-2] not in chain_ends and chain[-1][:-2] not in chain_ends:
                real_chains.append(chain) # this is a separate chain
            else:
                chain1, chain2 = None, None
                for ch in chains:
                    if chain[0] in ch:
                        chain1 = ch
                    if chain[-1] in ch:
                        chain2 = ch
                if chain1 != chain2 and chain1 != None and chain2 != None:
                    chains = [x for x in chains if x not in (chain1, chain2)]
                    if chain[0][:-2] == chain1[-1][:-2] and chain[-1][:-2] == chain2[0][:-2]:
                        print "deci bine"
                        newchain = chain1[:-2] + chain + chain2[2:]
                    elif chain[0][:-2] == chain1[0][:-2] and chain[-1][:-2] == chain2[-1][:-2]:
                        print "deci iaca bine"
                        newchain = list(reversed(chain1))[:-2] + chain + list(reversed(chain2))[2:]
                    elif chain[0][:-2] == chain1[0][:-2] and chain[-1][:-2] == chain2[0][:-2]:
                        print "deci super"
                        newchain = list(reversed(chain1))[:-2] + chain + chain2[2:]
                    elif chain[0][:-2] == chain1[-1][:-2] and chain[-1][:-2] == chain2[-1][:-2]:
                        print "deci super puper"
                        newchain = chain1[:-2] + chain + list(reversed(chain2))[2:]
                    chains.append(newchain)
    
        for ch in chains:
            real_chains.append(ch)
        # now filter just those chains that are of length > 3
        chains = []
        for ch in real_chains:
            if len(ch) > 4:
                chains.append(ch)
        final_graph = nx.DiGraph()
    
        NEXT = {}
        PREV = {}
        for chain in chains:
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
                        print "ADDING EDGE:", x[:-2], y[:-2]
                        NEXT[x[:-2]] = y[:-2]
                        PREV[y[:-2]] = x[:-2]
                    else:
                        if ori == ("1", "2"):
                            if orients[-1] == False:
                                final_graph.add_node(y[:-2], orien=False) # done
                                orients.append(False)
                            else:
                                final_graph.add_node(y[:-2], orien=False)
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
                        print "ADDING EDGE:", x[:-2], y[:-2]
                        NEXT[x[:-2]] = y[:-2]
                        PREV[y[:-2]] = x[:-2]
    
    
        deci_removed_contigs = set([x[:-2] for x in goodourgraph.nodes() if x[:-2] not in final_graph.nodes()])
        deci_removed_contigs = list(deci_removed_contigs)
        inserted_contigs = []
        removed_contigs = [x for x in deci_removed_contigs if x not in inserted_contigs and x not in repeat_nodes]
        print len(removed_contigs), "LEN OF REMOVED CONTIGS"   
 
        # first round of construction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for i in range(1):
            SLOTS = {}
            for x, y in final_graph.edges():
                SLOTS[(x, y)] = []
            for contig in removed_contigs:
                ones = goodourgraph.edge[contig + "_1"]
                twos = goodourgraph.edge[contig + "_2"]
    
                reject = False
    
                slot_array = []
    
    
                for her in ones:
    
                    slot_dict = {}
                    strand = int(her[-1])
                    support = goodourgraph.edge[contig + "_1"][her]["weight"]
    
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
                                #reject = True
                                #break
                                ##if dist > distance + matching_graph.node[NEXT[her[:-2]] + "_1"]["width"]:
                                    ##slot = (NEXT[her[:-2]], NEXT[NEXT[her[:-2]]])
                                if dist > distance:#matching_graph.node[contig + "_1"]["width"] > distance + 900:
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
                                #reject = True
                                #break
                                ##if dist > distance + matching_graph.node[PREV[her[:-2]] + "_1"]["width"]:
    
                                    ##slot = (PREV[PREV[her[:-2]]], PREV[her[:-2]])
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
                    support = goodourgraph.edge[contig + "_2"][her]["weight"]
    
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
                                #reject = True
                                #break
                                if dist > distance:#+ matching_graph.node[contig + "_1"]["width"] > distance + 900:
                                    slot = (NEXT[her[:-2]], NEXT[NEXT[her[:-2]]])
                                #if dist > distance + matching_graph.node[NEXT[her[:-2]] + "_1"]["width"]:
                                    #slot = (NEXT[her[:-2]], NEXT[NEXT[her[:-2]]]
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
                                #reject = True
                                #break
                                #if dist > distance + matching_graph.node[PREV[her[:-2]] + "_1"]["width"]:
                                    #slot = (PREV[PREV[her[:-2]]], PREV[her[:-2]])
    
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
                    print slot, contig, final_graph.edge["NC_010079:120:2112133:2114005"]
    
                inserted_contigs.append(contig)
    
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
                    if "NC_010063" in node:
                        print "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS"
                        continue
                    final_graph.add_node(node, orien=orien)
                    final_graph.add_edge(x, node, dist=dist)
                    final_graph.add_edge(node, y, dist=distance)
                    final_graph.remove_edge(x, y)
                    to_be_inserted.append(node)
                else:
                    slots = self.multikeysort(SLOTS[(x, y)], ["dist"])
                    hui = False
                    for sl in slots:
                        if "NC_010063" in sl:
                            print "SSSSSSSSSSSSJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
                            hui = True
                            break
                    if hui == True:
                        continue
                    slots = self.multikeysort(SLOTS[(x, y)], ["dist"])
                    for sl in slots:
                        final_graph.add_node(sl["contig"], orien=sl["orien"])
                    final_graph.add_edge(x, slots[0]["contig"], dist=100)
                    final_graph.add_edge(slots[-1]["contig"], y, dist=100)
                    for i in range(len(slots) - 1):
                        final_graph.add_edge(slots[i]["contig"], slots[i + 1]["contig"], dist=100)
    
    
            removed_contigs = [x for x in removed_contigs if x not in inserted_contigs]
            inserted_contigs = []
    
    
        for node in final_graph.nodes():
            if "orien" in final_graph.node[node].keys():
                pass
            else:
                final_graph.node[node]["orien"] = False
            final_graph.node[node]['width'] = matching_graph.node[node + "_1"]["width"]

        for x in repeat_nodes:
            print "ADDING:", x
            final_graph.add_node(x, orien=False, width=matching_graph.node[x + "_1"]["width"])
           
        print len(final_graph.nodes())
        wdir = self._settings.get("scaff_dir")
        nx.write_gpickle(final_graph, wdir + "/final_graph.cpickle")
    
    

