'''
Created on Jul 25, 2014

@author: igor
'''

import argparse
from settings import Settings
from pipeline.pipeline import Pipeline


def parse_args():
    main_p = argparse.ArgumentParser()
    main_p.add_argument('-w', dest='scaff_dir', required=True, help='scaffolding directory')
    main_p.add_argument('-c', dest='ctg_fasta', required=True, help='contig fasta')
    main_p.add_argument('-m1', dest='mapping1', required=True, help='read first file')
    main_p.add_argument('-m2', dest='mapping2', required=True, help='read second file')  
    main_p.add_argument('-i', dest='ins_size', required=True, help='insert size')
    main_p.add_argument('-p', dest='pair_mode', required=True, help='pair mode')
    main_p.add_argument('-s', dest='std_dev', required=True, help='standard deviation')
    return vars(main_p.parse_args())  



if __name__ == '__main__':
    args = parse_args()
    settings = Settings()
    settings.update(args)
    scaffolder = Pipeline()
    scaffolder.set_settings(settings)
    scaffolder.scaffold()
    
    
    
    