'''
code for dealing with data I/O
'''
import networkx as nx
import numpy as np
import subprocess
import logging
import sys
import string

agp_dt = np.dtype([\
    ('scaf_name', 'S255'),\
    ('scaf_start', np.long),\
    ('scaf_stop', np.long),\
    ('scaf_idx', np.long),\
    ('comp_type', 'S50'),\
    ('comp_name', 'S255'),\
    ('comp_start', np.long),\
    ('comp_stop', np.long),\
    ('comp_orien', np.long),\
    ('comp_linkage', np.long),\
])


def to_agp(G, agp_file):
    ''' writes graph to an AGP file '''

    # grab some info.
    nsz = G.number_of_nodes()
    esz = G.number_of_edges()
    tsz = nsz + esz

    # create AGP array.
    agp_arr = np.zeros(tsz, dtype=agp_dt)

    # locate connected components.
    components = nx.weakly_connected_components(G)

    # loop over components.
    idx = 0
    agp_idx = 0
    for component in components:

        # create subgraph.
        subg = nx.DiGraph()
        subg.add_nodes_from(component)

        # add edges.
        for edge in G.edges(component):

            # get distance.
            #dist = G[edge[0]][edge[1]]['means'][G[edge[0]][edge[1]]['state']]
            dist = G[edge[0]][edge[1]]['dist']

            # add the edge.
            subg.add_edge(edge[0],edge[1],dist=dist)

        # compute topological order.
        path = nx.topological_sort(subg)

        # validate path.
        pn = path[0]
        for i in range(1,len(path)):
            if subg.has_edge(pn, path[i]) == False:
                logging.error("Not a linear graph.")
                sys.exit()
            pn = path[i]


        # setup scaffold vars.
        scaf_idx = 1
        scaf_start = 1
        scaf_stop = G.node[path[0]]['width']

        # add first node.
        agp_arr[agp_idx]['scaf_name'] = "scaf_%i" % idx
        agp_arr[agp_idx]['scaf_start'] = scaf_start
        agp_arr[agp_idx]['scaf_stop'] = scaf_stop
        agp_arr[agp_idx]['scaf_idx'] = scaf_idx
        agp_arr[agp_idx]['comp_type'] = "W"
        agp_arr[agp_idx]['comp_name'] = path[0]
        agp_arr[agp_idx]['comp_start'] = 1
        agp_arr[agp_idx]['comp_stop'] = G.node[path[0]]['width']
        agp_arr[agp_idx]['comp_orien'] = G.node[path[0]]['orien']
        agp_arr[agp_idx]['comp_linkage'] = 0

        # increment pointers.
        agp_idx += 1
        pn = path[0]
        scaf_idx += 1
        scaf_start = scaf_stop + 1

        # fill AGP.
        for i in range(1,len(path)):

            # defualt gapsize to 7.
            gapsz = subg[pn][path[i]]['dist']

            #print gapsz

            # get stop location.
            scaf_stop = scaf_start + gapsz

            # add gap info.
            agp_arr[agp_idx]['scaf_name'] = "scaf_%i" % idx
            agp_arr[agp_idx]['scaf_start'] = scaf_start
            agp_arr[agp_idx]['scaf_stop'] = scaf_stop
            agp_arr[agp_idx]['scaf_idx'] = scaf_idx
            agp_arr[agp_idx]['comp_type'] = "N"
            agp_arr[agp_idx]['comp_name'] = "fragment"
            agp_arr[agp_idx]['comp_start'] = 1
            agp_arr[agp_idx]['comp_stop'] = gapsz
            agp_arr[agp_idx]['comp_linkage'] = 0

            # increment pointers.
            agp_idx += 1
            scaf_idx += 1
            scaf_start = scaf_stop + 1

            # set stop.
            scaf_stop = scaf_start + G.node[path[i]]['width']

            # add nodes.
            agp_arr[agp_idx]['scaf_name'] = "scaf_%i" % idx
            agp_arr[agp_idx]['scaf_start'] = scaf_start
            agp_arr[agp_idx]['scaf_stop'] = scaf_stop
            agp_arr[agp_idx]['scaf_idx'] = scaf_idx
            agp_arr[agp_idx]['comp_type'] = "W"
            agp_arr[agp_idx]['comp_name'] = path[i]
            agp_arr[agp_idx]['comp_start'] = 1
            agp_arr[agp_idx]['comp_stop'] = G.node[path[i]]['width']
            agp_arr[agp_idx]['comp_orien'] = G.node[path[i]]['orien']

            # increment pointers.
            pn = path[i]
            agp_idx += 1
            scaf_idx += 1
            scaf_start = scaf_stop + 1


        # increment scaffold idx
        idx += 1

    # save to disk
    save_agps(agp_file, agp_arr)

def draw_graph(graph, name):
    ''' draw the graph using neato'''

    # write dot file.
    dot = "./tmp.dot"
    nx.write_dot(graph, dot)

    # convert to picture using neato.
    pic = "/home/jlindsay/transfer/%s.png" % name
    subprocess.call(["neato", "-Tpng", "-o", pic, dot])

    # remove dot.
    subprocess.call(["rm","-f",dot])

def write_agp(graph_file, agp_file, runtime=None):
    """ translates oriented graph to agp file
    Parameters
    ----------
    paths.gap_file       : file
    args.agp_file           : file
    args.run_time           : float
    """

    # load the oriented graph.
    SG = nx.read_gpickle(graph_file)

    # ensure node degree is low.
    deg_list = [len(SG.neighbors(x)) for x in SG.nodes()]
    if max(deg_list) > 2:
        logging.error("is not a path")
        sys.exit(1)

    # ensure its a DAG.
    if nx.is_directed_acyclic_graph(SG) == False:
        logging.error("not a DAG?")
        sys.exit(1)

    # ensure positive gap size.
    for p,q in SG.edges():
        if SG[p][q]['dist'] < 10:
            SG[p][q]['dist'] = 10

    # save it to disk.
    to_agp(SG, agp_file)

    # append the timing.
    fin = open(agp_file, "rb")
    lines = fin.readlines()
    fin.close()
    if runtime != None:
        txt = "RUNTIME:\t%f\n" % (runtime)
        lines.append(txt)
    
    fout = open(agp_file, "wb")
    fout.write(''.join(lines))
    fout.close()

def save_agps(agp_file, agp):
    ''' saves agp to disk.'''

    # write to file.
    fout = open(agp_file, "w")

    # write each entry.
    z = len(agp_dt.names)
    for i in range(agp.size):

        # sanity skip.
        if agp[i]['scaf_name'] == "":
            continue

        # format result.
        tmp = agp[i]
        if tmp['comp_type'] == "W":
            # get orientation.
            if tmp["comp_orien"] == 0:
                o = "+"
            else:
                o = "-"

            # write contig.
            txt = str(tmp['scaf_name']) + "\t"
            txt += str(tmp['scaf_start']) + "\t"
            txt += str(tmp['scaf_stop']) + "\t"
            txt += str(tmp['scaf_idx']) + "\t"
            txt += str(tmp['comp_type']) + "\t"
            txt += str(tmp['comp_name']) + "\t"
            txt += str(tmp['comp_start']) + "\t"
            txt += str(tmp['comp_stop']) + "\t"
            txt += o + "\n"

        else:
            # get linkage.
            if tmp['comp_linkage'] == 0:
                o = "no"
            else:
                o = "yes"

            # write gap.
            txt = str(tmp['scaf_name']) + "\t"
            txt += str(tmp['scaf_start']) + "\t"
            txt += str(tmp['scaf_stop']) + "\t"
            txt += str(tmp['scaf_idx']) + "\t"
            txt += str(tmp['comp_type']) + "\t"
            txt += str(tmp['comp_stop'] - tmp['comp_start']) + "\t"
            txt += str(tmp['comp_name']) + "\t"
            txt += o + "\n"

        fout.write(txt)

    # close file.
    fout.close()

def load_agp(fpath):
    ''' read agp file into array.'''

    # read in agp.
    fin = open(fpath, "rb")
    lines = fin.readlines()
    fin.close()

    # count number of lines minus comments.
    cnt = 0
    for line in lines:
        if line[0] != "#" and len(line) != 0 and line.strip().split()[0] != "RUNTIME:":
            cnt += 1

    # instantiate array.
    agp_edges = np.zeros(cnt, dtype=agp_dt)

    # parse agp.
    idx = 0
    for line in lines:
        # tokenize.
        if line[0] == "#": continue
        tmp = line.strip().split()
        if len(tmp) == 0: continue
        if tmp[0] == "RUNTIME:": continue

        # get general tokenize.
        agp_edges[idx]['scaf_name'] = tmp[0]
        agp_edges[idx]['scaf_start'] = int(float(tmp[1]))
        agp_edges[idx]['scaf_stop'] = int(float(tmp[2]))
        agp_edges[idx]['scaf_idx'] = int(float(tmp[3]))
        agp_edges[idx]['comp_type'] = tmp[4]

        # contig.
        if tmp[4] == "W":
            # get parts.
            agp_edges[idx]['comp_name'] = tmp[5]
            agp_edges[idx]['comp_start'] = int(tmp[6])
            agp_edges[idx]['comp_stop'] = int(tmp[7])
            if tmp[8] == "+":
                agp_edges[idx]['comp_orien'] = 0
            else:
                agp_edges[idx]['comp_orien'] = 1

        else:

            # save entry.
            agp_edges[idx]['comp_name'] = tmp[6]
            agp_edges[idx]['comp_start'] = 1
            agp_edges[idx]['comp_stop'] = int(tmp[5])
            if tmp[7] != "yes":
                agp_edges[idx]['comp_linkage'] = 0
            else:
                agp_edges[idx]['comp_linkage'] = 1


        # update index.
        idx += 1

    # shirnk array.
    agp_edges.resize(idx)

    return agp_edges


def gen_sam_pair(fpath1, fpath2, flip_1, flip_2):
    ''' pulls info from paired sam file '''
    fin1 = open(fpath1)
    fin2 = open(fpath2)
    for line1 in fin1:
        line2 = fin2.readline()

        # sanity.
        if line1[0] == "@":
            continue

        # tokenize.
        tok1 = line1.strip().split()
        tok2 = line2.strip().split()

        # simplify.
        if tok1[1] == '0':
            orien1 = 0
        else:
            orien1 = 1

        if tok2[1] == '0':
            orien2 = 0
        else:
            orien2 = 1

        if flip_1 == True:
            orien1 = 1 - orien1
        if flip_2 == True:
            orien2 = 1 - orien2

        ctg1 = tok1[2]
        ctg2 = tok2[2]

        left1 = int(tok1[3])
        left2 = int(tok2[3])

        right1 = left1 + len(tok1[9])
        right2 = left2 + len(tok2[9])

        state = misc.determine_state(ctg1, ctg2, orien1, orien2)

        # yield the info.
        yield ctg1, left1, right1, ctg2, left2, right2, state, tok1, tok2

    fin1.close()
    fin2.close()

def write_gap_info(paths, args):
        
    # load the graphs.
    EG = nx.read_gpickle(paths.edge_file)
    DG = nx.read_gpickle(paths.order_file)
    
    # openoutput.
    fout = open("%s/gap_info.txt" % paths.work_dir, "w")
    
    # loop over each edge.
    for p,q in DG.edges():
        for i in EG[p][q]:
            e = EG[p][q][i]
            
            # write it.
            tmp = [
                p, q, DG[p][q]['state'], e['left1'], e['right1'], e['left2'], e['right2'],\
                e['state'], e['dist'], e['std_dev']
            ]
            fout.write(' '.join([str(x) for x in tmp]) + '\n')
        
    # done.
    fout.close()

def write_text(graph_file, node_file, bundle_file):
    ''' writes bundle graph to txt '''

    # load the bundle graph.
    BG = nx.read_gpickle(graph_file)

    # write out nodes.
    with open(node_file, "wb") as fout:
        for n in BG.nodes():
            fout.write("%s\t%i\t%s\n" % (n, BG.node[n]['width'], BG.node[n]['seq']))

    with open(bundle_file, "wb") as fout:
        for p, q in BG.edges():
            fout.write("%s\t%s\t%i\t%i\t%i\t%i\n" % (p, q, BG[p][q]['bcnts'][0], BG[p][q]['bcnts'][1], BG[p][q]['bcnts'][2], BG[p][q]['bcnts'][3]))


def load_fasta(file_path):
    ''' loads fasta file into dictionary'''

    # read file into memory.
    fin = open(file_path)
    lines = fin.readlines()
    fin.close()

    # build dictionary.
    data = dict()
    seq = ""
    for line in lines:

        # Skip blanks.
        if len(line) < 2: continue
        if line[0] == "#": continue

        # remove blanks.
        line = line.strip()

        # Check for ids.
        if line.count(">") > 0:

            # Check if ending seq.
            if len(seq) > 0:

                # save.
                data[head] = seq

            # reset head.
            head = line.replace(">","")
            seq = ""

            # skip to next line.
            continue

        # Filter chars.
        seq += line

    # save the last one.
    data[head] = seq

    # return dictionary.
    return data

def write_fasta(ctg_file, agp_file, scf_file):
    """ translates AGP and contigs to fasta file"""
    
    # load the contigs.
    contigs = load_fasta(ctg_file)
    
    for z in contigs.keys():
        a = z.split()[0]
        contigs[a] = contigs[z]
    
    # load the agp.
    agp = load_agp(agp_file)
        
    # convert to scaf fasta.
    scf_fasta = dict()
    for i in range(agp.shape[0]):
        
        # boot strap.
        if agp[i]['scaf_name'] not in scf_fasta:
            scf_fasta[agp[i]['scaf_name']] = ""
            
        # find the length.
        length = agp[i]['comp_stop'] - agp[i]['comp_start']
        
        # choose how to handle.
        if agp[i]['comp_type'] == 'N':
            # gap
            if length < 0:
                # minimum gap of 10
                #scf_fasta[agp[i]['scaf_name']] = scf_fasta[agp[i]['scaf_name']][0:length]
                scf_fasta[agp[i]['scaf_name']] = 'N' * 10
            else:
                # add gap to end.
                scf_fasta[agp[i]['scaf_name']] += 'N' * length
        else:
            # contig.
            comp_name = agp[i]['comp_name']
            if comp_name not in contigs:
                comp_name = comp_name.split()[0]
                if comp_name not in contigs:
                    print "whut name is not here"
                    print comp_name
                    sys.exit()
            seq = contigs[comp_name].upper()
            if agp[i]['comp_orien'] != 0:
                seq = seq.translate(string.maketrans("ATCGN", "TAGCN"))[::-1]
            scf_fasta[agp[i]['scaf_name']] += seq    
                
    # add any missing contigs back in.
    hits = set(list(agp[:]['comp_name']))
    if 'fragment' in hits:
        hits.remove('fragment')
    
    # comment it out
    #for name in contigs:
    #    if name not in hits:
    #        scf_fasta[name] = contigs[name]
                                
    # write out scaffolds to file.
    with open(scf_file, 'wb') as fout:
        for name in scf_fasta:
            if len(scf_fasta[name]) < 1:
                print "WHUT EMPTY SEQ"
                sys.exit()
            fout.write('>%s\n%s\n' % (name, scf_fasta[name].strip("N")))
