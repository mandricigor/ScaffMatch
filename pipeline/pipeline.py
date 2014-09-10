
from helpers.io import write_agp, write_fasta

from alignment.pairer import Pairer
from graph_construction.builder import GraphBuilder
from matching.matching import Matcher



class Pipeline(object):
    
    def __init__(self):
        pass
        
    def set_settings(self, settings):
        self._settings = settings
        
    def get_settings(self):
        return self._settings
        
    def _pair(self):
        pairer = Pairer()
        pairer.set_settings(self._settings)
        pairer.pair()
        
    
    def _scaffold_graph(self):
        graphBuilder = GraphBuilder()
        graphBuilder.set_settings(self._settings)
        graphBuilder.build_graph()
        
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
        self._pair()
        self._scaffold_graph()
        import time
        start_time = time.time()
        self._matching()
        self._print_fasta()
        print("--- %s seconds ---" % (time.time() - start_time))
