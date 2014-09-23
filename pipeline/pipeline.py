
from helpers.io import write_agp, write_fasta

from alignment.graph import GraphConstructor
from matching.matching import Matcher



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
        gconstructor.scaffolding_graph()
        
    def _matching(self):
        matcher = Matcher()
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
        
    def scaffold(self):
        self._graph_contruction()
        #import time
        #start_time = time.time()
        #self._matching()
        #self._print_fasta()
        #print("--- %s seconds ---" % (time.time() - start_time))
