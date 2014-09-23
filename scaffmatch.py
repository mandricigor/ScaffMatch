'''
Created on Jul 25, 2014

@author: igor
'''

import argparse
import sys
from settings import Settings
from pipeline.pipeline import Pipeline


def parse_args():
    main_p = argparse.ArgumentParser()
    main_p.add_argument('-w', dest='scaff_dir', required=True, help='scaffolding directory')
    main_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    main_p.add_argument('-m1', dest='mappings1', required=True, help='read first file')
    main_p.add_argument('-m2', dest='mappings2', required=True, help='read second file')  
    main_p.add_argument('-i', dest='ins_size', required=True, help='insert size')
    main_p.add_argument('-p', dest='pair_mode', required=True, help='pair mode')
    main_p.add_argument('-s', dest='std_dev', required=True, help='standard deviation')
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
    args = parse_args()
    prepare_libraries(args)
    settings = Settings()
    settings.update(args)
    scaffolder = Pipeline()
    scaffolder.set_settings(settings)
    scaffolder.scaffold()
    
    
    
    