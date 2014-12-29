
from random import random
from collections import deque, Counter
import networkx as nx
import fasta.io as io
from math import sqrt, log
import numpy
import operator
import os.path

class Link(object):
    def __init__(self, **args):
        for x, y in args.iteritems():
            self.__dict__[x] = y

    def __str__(self):
        return "distance: %s" % self.dist
    

class GraphConstructor(object):
    
    _settings = None
    
    def __init__(self):
        self._covs = {}
        self._hits = {}
        self._sizes = {}
        self._dist = {}
        self._IGORgraph = nx.Graph() # the main scaffolding graph
        
    def set_settings(self, settings):
        self._settings = settings
        
        
    def _read_in_chunks(self, file_object, chunk_size=10000000000):
        """Lazy function (generator) to read a file piece by piece.
        Default chunk size: 1k."""
        while True:
            data = file_object.read(chunk_size)
            if not data:
                break
            yield data
        

    
    def scaffolding_graph(self):
        """Read mapping files chunk by chunk and feed the lines
        corresponding to valid read pairs to graph construction"""
        self._contigs_coverage()
        libraries = self._settings.get("libraries")
        mapped_ = self._settings.get("mapped_file")
        unmapped_ = self._settings.get("unmapped_file")
        wdir = self._settings.get("scaff_dir") 

   
        curlen1 = curlen2 = 0
  
        for lib_id in libraries.keys(): # stub for multiple libraries
            lib = libraries[lib_id]
            sam1, sam2 = lib["sam1"], lib["sam2"]
            with open(os.path.join(wdir, mapped_), 'w') as mapped, open(os.path.join(wdir, unmapped_), 'w') as unmapped, \
                    open(sam1) as xx, open(sam2) as yy, open(os.path.join(wdir, "multiple.txt"), "w") as mtp, open(os.path.join(wdir, "same.txt"), "w") as same:
                fileiter1 = self._read_in_chunks(xx)
                fileiter2 = self._read_in_chunks(yy)
    
                leftover1 = leftover2 = ""
    
                while True:
                    try: # grab a bit from input
                        output1 = fileiter1.next()
                    except StopIteration:
                        output1 = ""
            
                    try: # grab a bit from input
                        output2 = fileiter2.next()
                    except StopIteration:
                        output2 = ""
                    
                    output1 = leftover1 + output1
                    output2 = leftover2 + output2
                    
                    if curlen1 == len(output1) and curlen2 == len(output2):
                        break
                    else:
                        curlen1 = len(output1) 
                        curlen2 = len(output2)
                    
                    if output1 == output2 == "":
                        break
                    
                    leftoverpos1 = output1.rfind("\n")
                    base1 = output1[:leftoverpos1 + 1]
                    leftover1 = output1[leftoverpos1 + 1:]
            
                    leftoverpos2 = output2.rfind("\n")
                    base2 = output2[:leftoverpos2 + 1]
                    leftover2 = output2[leftoverpos2 + 1:]
            
                    base1, base2 = deque(base1.split("\n")[:-1]), deque(base2.split("\n")[:-1])
                    
                    reads1 = map(lambda x: x.split()[0], base1)
                    reads2 = map(lambda x: x.split()[0], base2)
                    
                    reads11 = []
                    curread, curcounter = "", 0
                    for read in reads1:
                        if read == curread:
                            curcounter += 1
                        else:
                            reads11.append((curread, curcounter))
                            curread, curcounter = read, 1
                    else:
                        reads11.append((curread, curcounter))
                    reads11 = deque(reads11[1:-1])
                    
                    reads22 = []
                    curread, curcounter = "", 0
                    for read in reads2:
                        if read == curread:
                            curcounter += 1
                        else:
                            reads22.append((curread, curcounter))
                            curread, curcounter = read, 1
                    else:
                        reads22.append((curread, curcounter))
                    reads22 = deque(reads22[1:-1])
        
                    while reads11 and reads22:
                        r1, c1 = reads11.popleft()
                        r2, c2 = reads22.popleft()
    
                        if c1 > 1 or c2 > 1:
                            """skip these reads:
                            they are multiple"""
                            for i in range(c1):
                                mtp.write("%s\n" % base1.popleft())
                            for i in range(c2):
                                mtp.write("%s\n" % base2.popleft())
                        else:
                            line1, line2 = base1.popleft().split(), base2.popleft().split()

			    if line1[2] != "*":
				first = self._covs.get(line1[2], {"len": 0, "count": 0})
				first["len"] += len(line1[9])
                                first["count"] += 1
				self._covs[line1[2]] = first
			    if line2[2] != "*":
				second = self._covs.get(line2[2], {"len": 0, "count": 0})
				second["len"] += len(line2[9])
                                second["count"] += 1
                                self._covs[line2[2]] = second

                            if line1[2] == "*" or line2[2] == "*": # we use this with bowtie2
                                unmapped.write(self._format_sam_line(line1, line2))
                            else:
                                if line1[2] != line2[2]: # discard those who has mapped to the same contig
				    mismatches1, mismatches2 = line1[13].split(":")[-1], line2[13].split(":")[-1]
				    if int(mismatches1) > 1000 and int(mismatches2) > 1000:
				        continue
				    else:
                                        mapped.write(self._format_sam_line(line1, line2))
                                        self._paired_read_to_graph(line1, line2, lib_id)
                                else:
                                    same.write(self._format_sam_line(line1, line2))
                    leftover1 = "\n".join(base1) + "\n" + leftover1
                    leftover2 = "\n".join(base2) + "\n" + leftover2 
        self._build_graph()       
        for x in self._covs.keys():
            width = self._IGORgraph.node[x + "_1"]["width"]
	    self._IGORgraph.node[x + "_1"]["cov"] = self._covs[x]["len"] * 1.0 / width
            self._IGORgraph.node[x + "_2"]["cov"] = self._covs[x]["len"] * 1.0 / width


        # stub
        for node in self._IGORgraph.nodes():
            if "cov" not in self._IGORgraph.node[node]:
                self._IGORgraph.node[node]["cov"] = 1


	meancov = 0
        counter = 0
        std   = 0
        covs = []
        
        for x in self._covs.keys():
            meancov += self._IGORgraph.node[x + "_1"]["cov"]
            counter += 1
            covs.append(self._IGORgraph.node[x + "_1"]["cov"])

        meancov *= 1.0
        meancov /= counter

        dispcov = sqrt(sum([(x - meancov) * (x - meancov) for x in covs]) / counter)

        self._settings.set("mean_cov", meancov)
        self._settings.set("disp_cov", dispcov)
        # now the graph should be ready for the next stage
        #self._contigs_coverage2()
        for node in self._IGORgraph.nodes():
            self._IGORgraph.node[node]["mean_cov"] = meancov
            self._IGORgraph.node[node]["disp_cov"] = dispcov
        #nx.write_gpickle(self._IGORgraph, "igor_vasea_graph.cpickle")
        return self._IGORgraph


    def _format_sam_line(self, line1, line2):
        """Take what we need from two lines and format"""
        return "%s %s %s %s %s " % (line1[0], line1[1], line1[2], line1[3], line1[9]) + \
                "%s %s %s %s\n" % (line2[1], line2[2], line2[3], line2[9])




    def _paired_read_to_graph(self, line1, line2, lib_id):
        """
        @param line1: line corresponding to first fragment
        @param line2: line corresponding to second fragment
        @param lib_id: id of the library from which pair comes
        """
        # insert size and standard deviation
        library = self._settings.get("libraries").get(lib_id)
        ins_size, std_dev, pair_mode = library["ins"], library["std"], library["pm"]
 
        # first leg of the read
        oflag1, rname1, lpos1 = line1[1], line1[2], int(line1[3])

        # to remove if something
        #if oflag1 == "0" and lpos1 > ins_size + 3 * std_dev:
        #    return
        #if oflag1 != "0" and lpos1 < len(line1[9]) - 3 * std_dev:
        #    return
        width1 = self._IGORgraph.node[rname1 + "_1"]["width"]
        rpos1 = lpos1 + len(line1[9])
        
        # second leg of the read
        oflag2, rname2, lpos2 = line2[1], line2[2], int(line2[3])
        #if oflag2 == "0" and lpos2 > ins_size + 3 * std_dev:
        #    return
        #if oflag2 != "0" and lpos2 < len(line2[9]) - 3 * std_dev:
        #    return
        width2 = self._IGORgraph.node[rname2 + "_1"]["width"]
        rpos2 = lpos2 + len(line2[9])


        op, oq = self._get_orientation(oflag1, oflag2, int(pair_mode))
        orients = (op, oq)
               

 
        if orients == (0, 0):
            distance = ins_size - (width1 - lpos1) - rpos2
            edge = (rname1 + "_1", rname2 + "_2")
        elif orients == (1, 1):
            distance = ins_size - (width2 - lpos2) - rpos1
            edge = (rname1 + "_2", rname2 + "_1")
        elif orients == (0, 1):
            distance = ins_size - (width1 - lpos1) - (width2 - lpos2)
            edge = (rname1 + "_1", rname2 + "_1")
        elif orients == (1, 0):
            distance = ins_size - rpos1 - rpos2
            edge = (rname1 + "_2", rname2 + "_2")


        #if (rname1 == "CP000144:5:42563:42783" and rname2 == "CP000144:6:42807:43402") or \
        #    (rname1 == "CP000144:6:42807:43402" and rname2 == "CP000144:5:42563:42783"):
        #    print rname1, rname2, lpos1, lpos2, oflag1, oflag2, orients, distance, edge






        #print edge, orients, distance
        if -std_dev * 3 <= distance <= ins_size + 3 * std_dev:
            pair = tuple(sorted(edge))
            if rname1 < rname2:
                pos1, pos2 = lpos1, lpos2
            else:
                pos1, pos2 = lpos2, lpos1
                line1, line2 = line2, line1 # used to calculate entropy
   
	    # removee laterrrrrrrrrrrrr
            para = self._dist.get(pair, [])
            link_dict = {"pos1": pos1, "pos2": pos2, "dist": distance,
                            "ins": ins_size, "std": std_dev}
            link = Link(**link_dict)
            para.append(link)
	    self._dist[pair] = para
                

    def _entropy(self, s):
        """Compute entropy of a read"""
        arr = []
        for i in range(len(s) - 1):
            arr.append(s[i: (i + 2)])
        s = arr
        p, lns = Counter(s), float(len(s))
        return -sum( count/lns * log(count/lns, 2) for count in p.values())   


 
    def _get_orientation(self, oflag1, oflag2, pair_mode):
        """Gets orientation from SAM object for pairs"""
        if oflag1 == "0":
            o1 = 0
        else:
            o1 = 1
        if oflag2 == "0":
            o2 = 0
        else:
            o2 = 1
        if pair_mode == 1:
            o2 = 1 - o2
        if pair_mode == 2:
            o1 = 1 - o1
        return o1, o2


    
    def _contigs_coverage(self):
        """Calculate the coverage and the standard deviation.
        Initialize IGORgraph with nodes"""
        seqs = io.load_fasta(self._settings.get("ctg_fasta"))
        for name, seq in seqs.items():
            tmp = name.split(" ")
            name = tmp[0]
            l = len(seq)
            self._IGORgraph.add_node(name + "_1", {'seq': seq, 'width': l}) # first strand
            self._IGORgraph.add_node(name + "_2", {'seq': seq, 'width': l}) # second strand


   
    def _contigs_coverage2(self):
        """Calculate the coverage and the standard deviation.
        Initialize IGORgraph with nodes"""
        seqs = io.load_fasta(self._settings.get("ctg_fasta"))
        mean_cov = 0
        counter = 0
        covs = []
        for name, seq in seqs.items():
            # skip split names.
            tmp = name.split(" ")
            name = tmp[0]
            l = len(seq)
            cov = self._IGORgraph.degree(name + "_1") * 1.0 / l
            mean_cov += cov
            counter += 1
            covs.append(cov)
            self._IGORgraph.node[name + "_1"]['cov'] = cov
            self._IGORgraph.node[name + "_2"]['cov'] = cov
        mean_cov *= 1.0
        mean_cov /= counter
        
        disp_cov = sqrt(sum([(x - mean_cov) * (x - mean_cov) for x in covs]) / counter)
        
        self._settings.set("mean_cov", mean_cov)
        self._settings.set("disp_cov", disp_cov)
    




    def _build_graph(self):
        """build the graph on the data obtained
        from parsing the libraries"""
        """
        @note: perform bundling step as described in the paper
        "The Greedy Path-Merging Algorithm for Contig Scaffolding"
        by Huson, Reinert and Myers, 2002
        """
        nodedict = {}
        for node1, node2 in self._dist:
            ar = self._dist.get((node1, node2), [])
            ar1 = nodedict.get(node1, {})
            ar2 = nodedict.get(node2, {})
          
            for link in ar:
                pos1 = link.pos1
                aaa = ar1.get(pos1, [])
                aaa.append(node2)
                ar1[pos1] = aaa

                pos2 = link.pos2
                bbb = ar2.get(pos2, [])
                bbb.append(node1)
                ar2[pos2] = bbb

            nodedict[node1] = ar1
            nodedict[node2] = ar2
                 

        for node1, node2 in self._dist:

            links = self._dist[(node1, node2)]
            # ------------------------
            linkdists = {}
            for link in links:
                linkdists[link] = link.dist

            distslinks = {} # identify the link by dist
            for link in links:
                distslinks[link.dist] = link


            sorted_x = sorted(linkdists.items(), key=operator.itemgetter(1))

        
            median = len(sorted_x) / 2
          
            p = 0
            q = 0
            newmean = 0
            newstd = 0
            size = 0
            
            #std_dev = int(self._settings.get("std_dev"))

            lowerbound = sorted_x[median][1] - 3 * distslinks[sorted_x[median][1]].std  #std_dev
            upperbound = sorted_x[median][1] + 3 * distslinks[sorted_x[median][1]].std
            
            iterator = 0
            while iterator < len(linkdists) and sorted_x[iterator][1] < upperbound:
                std_dev = distslinks[sorted_x[iterator][1]].std
                if lowerbound < sorted_x[iterator][1] < upperbound:
                    size += 1
                    p += sorted_x[iterator][1] * 1.0 / pow(std_dev, 2)
                    q += 1.0 / pow(std_dev, 2)
                    iterator += 1
                else:
                    iterator += 1

            newmean = p / q
            self._IGORgraph.add_edge(node1, node2, weight=size, dist=newmean)

            # ------------------------


            # comment it out for a moment while testing bundling step
        """with open("dist_stat.txt", "w") as f:
            for x, y in self._IGORgraph.edges():
                f.write("%s     %s        %s\n" % (x, y, self._IGORgraph.edge[x][y]["weight"]))"""




