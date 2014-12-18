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




class Settings(object):
    _instance = None
    
    _settings = {
                "unmapped_file": "unmapped.txt",
                "mapped_file": "mapped.txt",
                "matching": "max_weight",
                "bundle_threshold": 5
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
    main_p = argparse.ArgumentParser()
    main_p.add_argument('-w', dest='scaff_dir', required=True, help='scaffolding directory')
    main_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    main_p.add_argument('-m1', dest='mappings1', required=True, help='read first file')
    main_p.add_argument('-m2', dest='mappings2', required=True, help='read second file')  
    main_p.add_argument('-i', dest='ins_size', required=True, help='insert size')
    main_p.add_argument('-p', dest='pair_mode', required=True, help='pair mode')
    main_p.add_argument('-s', dest='std_dev', required=True, help='standard deviation')
    main_p.add_argument('-t', dest='bundle_threshold', required=False, default=5, help='bundle threshold')
    main_p.add_argument('-g', dest='matching', required=False, default="max_weight", help='matching heuristic')
    main_p.add_argument('-l', dest='log_file', required=False, help='log file')
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
    try:
        mappings1 = mappings1.split(",")
        mappings2 = mappings2.split(",")
        insert_sizes = insert_sizes.split(",")
        std_devs = std_devs.split(",")
    except:
        print "You have issues in defining mapping files"
        sys.exit()
    if len(mappings1) != len(mappings2):
        print "Lengths of mapping files are not the same"
        sys.exit()
    libraries = {}
    libcount = 1
    for sam1, sam2, ins, std in \
                    zip(mappings1, mappings2, insert_sizes, std_devs):
        libraries[libcount] = {"ins": int(ins), "std": int(std), "sam1": sam1, "sam2": sam2}
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
    if not log_file:
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
        logger.info("    %s  -- %s" % (s, v)) 
    # Feed the settings to the scaffolder pipeline
    scaffolder = Pipeline()
    scaffolder.set_settings(settings)
    # Go!
    scaffolder.scaffold()
