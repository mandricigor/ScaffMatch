
import time
from os import remove
from fasta.io import write_agp, write_fasta
from alignment.graph import GraphConstructor
from matching.matching import Matcher
import networkx as nx


class Pipeline(object):
    
    def __init__(self):
        pass
        
    def set_settings(self, settings):
        self._settings = settings
        
    def get_settings(self):
        return self._settings
            
    def _graph_contruction(self):
        gconstructor = GraphConstructor()
        gconstructor.set_settings(self._settings)
        return gconstructor.scaffolding_graph()
        
    def _matching(self, graph):
        matcher = Matcher(graph)
        matcher.set_settings(self._settings)
        matcher.match()
        
    def _print_fasta(self):
        wdir = self._settings.get("scaff_dir")
        graph_pickle = wdir + "/final_graph.cpickle"
        agp_file = wdir + "/final_agp.agp"
        ctg_file = self._settings.get("ctg_fasta")
        write_agp(graph_pickle, agp_file)
        scaff_file = wdir + "/scaffolds.fa"
        write_fasta(ctg_file, agp_file, scaff_file)
        try:
            remove(agp_file) # cleaning up the agp_file
        except Exception:
            pass
        
    def scaffold(self):
        print "Construct scaffolding graph"
        start_time = time.time()
        scaffgraph = self._graph_contruction()
        #self._settings.set("mean_cov", scaffgraph.node[scaffgraph.nodes()[0]]["mean_cov"])
        #self._settings.set("disp_cov", scaffgraph.node[scaffgraph.nodes()[0]]["disp_cov"])
        print "Scaffolding graph constructed: %s" % (time.time() - start_time)
        start_time2 = time.time()
        self._matching(scaffgraph)
        print "Scaffolding done: %s" % (time.time() - start_time2)
        start_time3 = time.time()
        self._print_fasta()
        print "Results written to the disk: %s" % (time.time() - start_time3)
        print "Overall time: %s" % (time.time() - start_time)
