#!/usr/bin/env python

'''
Created on Jul 25, 2014

@author: igor
'''

import argparse
import sys
import logging
from pipeline.pipeline import Pipeline
from copy import copy



class MyParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)




class Settings(object):
    _instance = None
    
    _settings = {
                "unmapped_file": "unmapped.txt",
                "mapped_file": "mapped.txt",
                "matching": "max_weight",
                }
    _defaults = copy(_settings)
    
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Settings, cls).__new__(
                                cls, *args, **kwargs)
        return cls._instance
    
    def update(self, new_settings):
        self._settings.update(new_settings)
        
    def get(self, key):
        return self._settings.get(key)
    
    def set(self, key, value):
        self._settings[key] = value
           
    def __str__(self):
        return "\n".join(["%s: %s" % x for x in self._settings.iteritems()])


    def iteritems(self):
        it = []
        for s, v in self._settings.iteritems():
            if s not in self._defaults and s != "logger":
                it.append((s, v))
        return it


def parse_args():
    main_p = MyParser()
    main_p.add_argument('-w', dest='scaff_dir', required=True, help='scaffolding directory (working directory where ScaffMatch files will be stored. The resulting scaffolds are placed in $scaff_dir/scaffolds.fa file')
    main_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta file')
    main_p.add_argument('-m1', dest='mappings1', required=True, help='comma separated list of .sam files (first read in the read pair)')
    main_p.add_argument('-m2', dest='mappings2', required=True, help='comma separated list of .sam files (second read in the read pair)')  
    main_p.add_argument('-i', dest='ins_size', required=True, help='insert sizes (comma separated values)')
    main_p.add_argument('-p', dest='pair_mode', required=True, help='pair modes (0 - SOLiD style -> ->, 1 - innie style -> <-, 2 - outtie style <- ->) (comma separated values)')
    main_p.add_argument('-s', dest='std_dev', required=True, help='libraries standard deviations (comma separated values)')
    main_p.add_argument('-t', dest='bundle_threshold', required=False, default=5, help='bundle threshold: links between contigs with # of reads < bundle_threshold are discarded')
    main_p.add_argument('-g', dest='matching', required=False, default="max_weight", help='matching heuristic: max_weight - Maximum Weight Matching heuristics, greedy - Greedy heuristics, backbone - Maximum Weight Matching heuristics without Insertion step')
    main_p.add_argument('-l', dest='log_file', required=False, default="", help='log file (optional argument)')
    return vars(main_p.parse_args())  


def prepare_libraries(args_dict):
    """If there are more then one single
    library, pair mappings as necessary"""
    mappings1 = args_dict["mappings1"]
    args_dict.pop("mappings1", None)
    mappings2 = args_dict["mappings2"]
    args_dict.pop("mappings2", None)
    insert_sizes = args_dict["ins_size"]
    std_devs = args_dict["std_dev"]
    pair_modes = args_dict["pair_mode"]
    try:
        mappings1 = mappings1.split(",")
        mappings2 = mappings2.split(",")
        insert_sizes = insert_sizes.split(",")
        std_devs = std_devs.split(",")
        pair_modes = pair_modes.split(",")
    except:
        print "You have issues in defining mapping files"
        sys.exit()
    if not (len(mappings1) == len(mappings2) == len(insert_sizes) == len(std_devs) == len(pair_modes)):
        print "You missed some arguments, please check them carefully"
        sys.exit()
    libraries = {}
    libcount = 1
    for sam1, sam2, ins, std, pm in \
                    zip(mappings1, mappings2, insert_sizes, std_devs, pair_modes):
        libraries[libcount] = {"ins": int(ins), "std": int(std), "sam1": sam1, "sam2": sam2, "pm": pm}
        libcount += 1
    args_dict["libraries"] = libraries


if __name__ == '__main__':
    # Parse arguments
    args = parse_args()
    prepare_libraries(args)
    settings = Settings()
    settings.update(args)
    # Set up logging
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    log_file = settings.get("log_file")
    if not log_file or log_file == "NONE":
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(settings.get("log_file"))
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    settings.set("logger", logger)
    # Print out the settings
    logger.info("**--**--**--**--**--**--**--**--**--**--**--**--**--**--**--**--**")
    logger.info("Settings used for this run of ScaffMatch are:")
    for s, v in settings.iteritems():
        if s in ["std_dev", "ins_size", "pair_mode"]:
            continue
        logger.info("    %s  -- %s" % (s, v)) 
    # Feed the settings to the scaffolder pipeline
    scaffolder = Pipeline()
    scaffolder.set_settings(settings)
    # Go!
    scaffolder.scaffold()
    logger.info("Done!")
