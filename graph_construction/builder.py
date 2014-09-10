
import sys
import os
import logging
import networkx as nx  # @UnresolvedImport
import pickle
import mmap
from math import sqrt  # @UnresolvedImport

import helpers.misc as misc
import helpers.io as io




class SamToken(object):
    QNAME = ""
    OFLAG = ""
    RNAME = ""
    LPOS = 0
    RPOS = 0

    def __str__(self):
        return "|| %s - %s - %s - %s - %s ||" % (self.QNAME, self.OFLAG, self.RNAME, self.LPOS, self.RPOS)


class GraphBuilder(object):
    
    _settings = None
    
    def __init__(self):
        pass
    
        
    def set_settings(self, settings):
        self._settings = settings
        
        
    def build_graph(self):
        
        seqs = io.load_fasta(self._settings.get("ctg_fasta"))
        EG = nx.MultiGraph()
        
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
            print name, cov
            # add node.
            EG.add_node(name, {'seq':seq, 'width':l, 'cov':cov})
            
        mean_cov *= 1.0
        mean_cov /= counter
        
        disp_cov = sqrt(sum([(x - mean_cov) * (x - mean_cov) for x in covs]) / counter)
        
        self._settings.set("mean_cov", mean_cov)
        self._settings.set("disp_cov", disp_cov)
        
        print mean_cov, disp_cov
        
        
        IGORgraph = nx.Graph()
        for node in EG.nodes():
            IGORgraph.add_node(node + "_1", cov=EG.node[node]["cov"])
            IGORgraph.add_node(node + "_2", cov=EG.node[node]["cov"])
            
            
        ins_size = int(self._settings.get("ins_size"))
        std_dev = int(self._settings.get("std_dev"))
    
    
   

        data = {} # data = {pair: {left: {}, right: {}}}
        dist = {} # data = {pair: {(left, right): distance}}
        mapping1 = self._settings.get("mapping1")
        mapping2 = self._settings.get("mapping2")
        
        counter = 0
        counter2 = 0
        
        for sam1, sam2, pos1, pos2 in self.pair_gen(mapping1, mapping2):
            
            counter += 1
            
            p = sam1.RNAME
            q = sam2.RNAME
    
    
            if p == q:
                counter2 += 1
                continue
    
            order = (p, q)
            width1 = EG.node[p]['width']
            width2 = EG.node[q]['width']
            # increment coverage.
            EG.node[p]['cov'] += sam1.RPOS - sam1.LPOS
            EG.node[q]['cov'] += sam2.RPOS - sam2.LPOS
    
    
            # stateA = 0, stateB = 1, stateC = 2, stateD = 3
            op, oq = misc.get_orien(sam1, sam2, int(self._settings.get("pair_mode")))
    
            orients = (op, oq)
            
    
    
            if orients == (0, 0):
                distance = ins_size - (width1 - sam1.LPOS) - (sam2.RPOS)
                edge = (p + "_1", q + "_2")
            elif orients == (1, 1):
                distance2 = ins_size - (width2 - sam2.LPOS) - sam1.RPOS
                edge = (p + "_2", q + "_1")
            elif orients == (0, 1):
                distance = ins_size - (width1 - sam1.LPOS) - (width2 - sam2.LPOS)
                edge = (p + "_1", q + "_1")
            elif orients == (1, 0):
                distance2 = ins_size - sam1.RPOS - sam2.RPOS
                edge = (p + "_2", q + "_2")
            if orients == (1, 1):
                #distdict2 = CORRECT_DIST.get(tuple(sorted((q + "_1", p + "_2"))), []) # this is correct
                if -std_dev * 2 <= distance2 <= ins_size + 2 * std_dev:
                    distanta = distance2
                else:
                    continue
                    #distdict2.append(distance2)
                    #CORRECT_DIST[tuple(sorted((q + "_1", p + "_2")))] = distdict2
            elif orients == (0, 0): # DONE WITH THAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                #distdict = CORRECT_DIST.get(tuple(sorted((p + "_1", q + "_2"))), []) # this is correct
                if -std_dev * 2 <= distance <= ins_size + 2 * std_dev:
                    distanta = distance
                else:
                    continue
                    #distdict.append(distance)
                    #CORRECT_DIST[tuple(sorted((p + "_1", q + "_2")))] = distdict
            elif orients == (0, 1):
                #distdict = CORRECT_DIST.get(tuple(sorted((p + "_1", q + "_1"))), []) # this is correct
                if -std_dev * 2 <= distance <= ins_size + 2 * std_dev:
                    distanta = distance
                else:
                    continue
                    #distdict.append(distance)
                    #CORRECT_DIST[tuple(sorted((p + "_1", q + "_1")))] = distdict
            elif orients == (1, 0):
                #distdict = CORRECT_DIST.get(tuple(sorted((p + "_1", q + "_1"))), [])
                #distdict2 = CORRECT_DIST.get(tuple(sorted((p + "_2", q + "_2"))), []) # this is correct
                if -std_dev * 2 <= distance2 <= ins_size + 2 * std_dev:
                    distanta = distance2
                else:
                    continue
                #distdict2.append(distance2)
                #CORRECT_DIST[tuple(sorted((p + "_2", q + "_2")))] = distdict2

            pair = tuple(sorted(edge))

        if p < q:
            pos1, pos2 = sam1.LPOS, sam2.LPOS
        else:
            pos1, pos2 = sam2.LPOS, sam1.LPOS

        para = dist.get(pair) 


        if para is None:
            data[pair] = {"left": {pos1: pos2}, "right": {pos2: pos1}}
            dist[pair] = {(pos1, pos2): distanta}
        else:
            if 1 != 1: #(pos1, pos2) in para:
                pass
            else:
                """left = data[pair].get("left")
		        right = data[pair].get("right")
		        if left:
		            if pos1 in left:
			        p2 = left[pos1]
			        del data[pair]["left"][pos1]
			        #del data[pair]["right"][pos2]
			        if (pos1, p2) in dist[pair]:
			            del dist[pair][(pos1, p2)]
			        continue	
		        if right:
			        if pos2 in right:
			        p1 = right[pos2]
			        #del data[pair]["left"][pos1]
			        del data[pair]["right"][pos2]
			        if (p1, pos2) in dist[pair]:
			            del dist[pair][(p1, pos2)]
			        continue"""
                dist[pair][(pos1, pos2)] = distanta
                left = data[pair].get("left", {})
                left[pos1] = pos2
                right = data[pair].get("right", {})
                right[pos2] = pos1
                data[pair]["left"], data[pair]["right"] = left, right


        nx.write_gpickle(data, "data.cpickle")
        nx.write_gpickle(dist, "dist.cpickle")


        for x, y in dist:
            distances = dist[(x, y)].values()
            if distances:
                d = sum(distances) * 1.0 / len(distances)
                disp = sqrt(sum([(q - d) * (q - d) for q in distances]) / len(distances))
                IGORgraph.add_edge(x, y, dist=d, weight=len(distances), dispersion=disp)

            #if not IGORgraph.has_edge(*edge):
            #    IGORgraph.add_edge(edge[0], edge[1], weight=1)
            #else:
            #    IGORgraph.edge[edge[0]][edge[1]]["weight"] += 1


        # now add the width of each contig into igor's graph
        for node in IGORgraph.nodes():
            IGORgraph.node[node]["width"] = EG.node[node[:-2]]["width"]
        """for x, y in distances:
            dist, counter = distances[(x, y)]
            dist = dist * 1.0 / counter
            IGORgraph.edge[x][y]["dist"] = dist
    
        CORRECT_DIST2 = {}
        for key in CORRECT_DIST:
            if CORRECT_DIST[key]:
                x1 = CORRECT_DIST[key]
                dist = sum(x1) / len(x1)
	        CORRECT_DIST2[key] = {"dist": sum(x1) / len(x1), "weight": len(x1)}"""
        wdir = self._settings.get("scaff_dir")
        #nx.write_gpickle(CORRECT_DIST2, wdir + "/correct_distances.cpickle")
        nx.write_gpickle(IGORgraph, wdir + "/igor_graph.cpickle")

        return EG

        
        
        
    def pop_sam(self, token, sam):
        ''' populates object '''
    
        sam.QNAME = token[0]
        sam.OFLAG = token[1]
        sam.RNAME = token[2]
        sam.LPOS = int(token[3])
        sam.RPOS = sam.LPOS + len(token[9])
    
    
    def sam_gen(self, file_path):
        ''' generator for sam files '''
    
        # create the SAM object.
        sam = SamToken()
    
        # start the token generator.
        for token, pos in self.token_gen(file_path, "\t"):
    
            # fill the object.
            try:
                self.pop_sam(token, sam)
            except:
                print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                continue
    
            # yield it.
            yield sam, pos
    
    
    def pair_gen(self, file_path1, file_path2):
        ''' generator for sam files '''
    
        # create the SAM object.
        sama = SamToken()
        samb = SamToken()
    
        # start the token generator.
        gena = self.token_gen(file_path1, "\t")
        genb = self.token_gen(file_path2, "\t")
    
        # loop over first iterator.
        for tokena, posa in gena:
            tokenb, posb = genb.next()
    
            # fill the object.
            self.pop_sam(tokena, sama)
            self.pop_sam(tokenb, samb)
    
            # yield it.
            yield sama, samb, posa, posb
    
    
    def openmm(self, file_path):
        fin = open(file_path, "r")
        mmin = mmap.mmap(fin.fileno(), 0, access=mmap.ACCESS_COPY)
        return fin, mmin
    
    
    
    def closemm(self, fin, mmin):
        mmin.close()
        fin.close()
    
    
    def token_gen(self, file_path, delim):
        ''' generates tokens by delim '''
    
        # open the file and memory map.
        with open(file_path, "rb") as fin:
    
            # begin yielding tokens.
            pos = 0
            for line in fin:
    
                # yield the line.
                yield line.strip().split(delim), pos
    
                # update info.
                pos += len(line)
