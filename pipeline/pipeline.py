
import time
import logging
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
        logger = self._settings.get("logger")
        logger.info("Stage I: Construct scaffolding graph")
        start_time = time.time()
        scaffgraph = self._graph_contruction()
        logger.info("Scaffolding graph constructed: %s" % (time.time() - start_time))
        start_time2 = time.time()
        logger.info("Stage II: Performing matching heuristics (%s)" % self._settings.get("matching"))
        self._matching(scaffgraph)
        logger.info("Scaffolding done: %s" % (time.time() - start_time2))
        start_time3 = time.time()
        self._print_fasta()
        logger.info("Results written to the disk: %s" % (time.time() - start_time3))
        logger.info("Overall time: %s" % (time.time() - start_time))
        logger.info("--**--**--**--**--**--**--**--**--**--**--**--**--**--**--**--**--")
