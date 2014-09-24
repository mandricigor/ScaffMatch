

from collections import deque
import networkx as nx
import fasta.io as io
from math import sqrt



class Link(object):
    def __init__(self, **args):
        for x, y in args.iteritems():
            self.__dict__[x] = y

    def __str__(self):
        return "distance: %s" % self.dist
    

class GraphConstructor(object):
    
    _settings = None
    
    def __init__(self):
        self._cov = {}
        self._sizes = {}
        self._dist = {}
        self._IGORgraph = nx.Graph() # the main scaffolding graph
        
        
    def set_settings(self, settings):
        self._settings = settings
        
        
    def _read_in_chunks(self, file_object, chunk_size=2000000):
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
        mapped = self._settings.get("mapped_file")
        unmapped = self._settings.get("unmapped_file")
    
        curlen1 = curlen2 = 0
    
        for lib_id in libraries.keys()[0:1]: # stub for multiple libraries
            lib = libraries[lib_id]
            sam1, sam2 = lib["sam1"], lib["sam2"]
            with open(mapped, 'w') as mapped, open(unmapped, 'w') as unmapped, \
                    open(sam1) as xx, open(sam2) as yy:
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
                                base1.popleft()
                            for i in range(c2):
                                base2.popleft()                            
                        else:
                            line1, line2 = base1.popleft().split(), base2.popleft().split()
                            if line1[2] == "*" or line2[2] == "*":
                                unmapped.write(self._format_sam_line(line1, line2))
                            else:
                                if line1[2] != line2[2]: # discard those who has mapped to the same contig
                                    mapped.write(self._format_sam_line(line1, line2))
                                    self._paired_read_to_graph(line1, line2, lib_id)
                    leftover1 = "\n".join(base1) + "\n" + leftover1
                    leftover2 = "\n".join(base2) + "\n" + leftover2 
        self._build_graph()           
        # now the graph should be ready for the next stage
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
        ins_size, std_dev = library["ins"], library["std"]
        
        # first leg of the read
        oflag1, rname1, lpos1 = line1[1], line1[2], int(line1[3])
        width1 = len(line1[9])
        rpos1 = lpos1 + width1
        
        # second leg of the read
        oflag2, rname2, lpos2 = line2[1], line2[2], int(line2[3])
        width2 = len(line2[9])
        rpos2 = lpos2 + width2
        op, oq = self._get_orientation(oflag1, oflag2, int(self._settings.get("pair_mode")))
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

        #print edge, orients, distance, lpos1, lpos2
        if -std_dev * 2 <= distance <= ins_size + 2 * std_dev:
            pair = tuple(sorted(edge))
            # change the order of positions if needed
            if rname1 < rname2:
                pos1, pos2 = lpos1, lpos2
            else:
                pos1, pos2 = lpos2, lpos1
    
            para = self._dist.setdefault(pair, [])
            link_dict = {"pos1": pos1, "pos2": pos2, "dist": distance,
                            "ins": ins_size, "std": std_dev}
            link = Link(**link_dict)
            para.append(link)
                
    
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
        mean_cov = 0
        counter = 0
        covs = []
        for name, seq in seqs.items():
            # skip split names.
            tmp = name.split(" ")
            name = tmp[0]
            l = len(seq)
            try:
                cov = self._settings.get("cov")[name] * 1.0 / l
            except Exception:
                cov = 0
            mean_cov += cov
            counter += 1
            covs.append(cov)
            self._IGORgraph.add_node(name + "_1", {'seq': seq, 'width': l, 'cov': cov}) # first strand
            self._IGORgraph.add_node(name + "_2", {'seq': seq, 'width': l, 'cov': cov}) # second strand
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
        for node1, node2 in self._dist:
            links = self._dist[(node1, node2)]
            dists = [x.dist for x in links]
            if dists:
                linkw = len(dists)
                avgdist = float(sum(dists)) / linkw                
            else:
                avgdist = None
            if avgdist:
                self._IGORgraph.add_edge(node1, node2, weight=linkw, dist=avgdist)