

import os
import subprocess
import logging
import numpy as np
import sys



class Pairer(object):
    
    _settings = None
    
    def __init__(self):
        self._cov = {}
        self._sizes = {}
    
    def pair(self):
        """ creates the alignment """
    
        # key size.
        key_size = self._settings.get("key_size")
    
        # relavent files.
        base_dir = self._settings.get("scaff_dir")
        ctg_fasta = self._settings.get("ctg_fasta")
        in1_sam = self._settings.get("mapping1")
        in2_sam = self._settings.get("mapping2")
    
        read1_sam = os.path.abspath('%s/read1.sam' % base_dir)
        read2_sam = os.path.abspath('%s/read2.sam' % base_dir)
    
        names1_npy = os.path.abspath('%s/name1.npy' % base_dir)
        names2_npy = os.path.abspath('%s/name2.npy' % base_dir)
        sort1_npy = os.path.abspath('%s/sort1.npy' % base_dir)
        sort2_npy = os.path.abspath('%s/sort2.npy' % base_dir)
    
        # compute name sizes.
        names_size1 = self._name_size(in1_sam)
        names_size2 = self._name_size(in2_sam)
    
    
        # check if sorted is present.
        if os.path.isfile(sort1_npy) == False:
 
            # create / load name array.
            if os.path.isfile(names1_npy) == False:
                logging.info("creating name array 1")
                names1 = self._extract_names(in1_sam, names_size1, key_size)
                self._names(file_name=names1_npy, data=names1)

            else:
                logging.info("loading name array 1")
                names1 = self._names(file_name=names1_npy)
 
     
            # sort it.
            logging.info("sorting name array 1")
            names1.sort(order=['name'])
            self._names(file_name=sort1_npy, data=names1)
            del names1
            subprocess.call(["rm", "-f", names1_npy])
     
        # check if sorted is present.
        if os.path.isfile(sort2_npy) == False:
     
            # create / load name array.
            if os.path.isfile(names2_npy) == False:
                logging.info("creating name array 2")
                names2 = self._extract_names(in2_sam, names_size2, key_size)
                self._names(file_name=names2_npy, data=names2)
            else:
                logging.info("loading name array 2")
                names2 = self._names(file_name=names2_npy)
     
            # sort it.
            logging.info("sorting name array 2")
            names2.sort(order=['name'])
            self._names(file_name=sort2_npy, data=names2)
            del names2
            subprocess.call(["rm", "-f", names2_npy])
    
        # create sizes.
        self._sizes = self.create_sizes(ctg_fasta)
    
    
        # do work.
        self._dual_loop(sort1_npy, sort2_npy, in1_sam, in2_sam, read1_sam, read2_sam)
        
        self._settings.set("cov", self._cov)
        
    
    

        
    def set_settings(self, settings):
        self._settings = settings
        
       
    def create_sizes(self, ctg_fasta):
        sizes = {}
        with open(ctg_fasta) as f:
            lines = f.readlines()
        if lines:
            lines = ''.join(lines)
            lines = lines.split('>')[1:]
            for line in lines:
                line = line.split("\n")
                contig = line[0]
                lens = sum(map(len, line[1:]))
                sizes[contig] = lens
        return sizes
        
            
    def _dual_loop(self, sort1_npy, sort2_npy, in1_sam, in2_sam, out1_sam, out2_sam):
        """ extract unique alignments, pairs them and annotate repeats"""
    
        # open SAM files.
        sam1 = open(in1_sam, "rb")
        sam2 = open(in2_sam, "rb")
        out1 = open(out1_sam, "wb")
        out2 = open(out2_sam, "wb")
    
    
        # create iterators.
        itr1 = self._uniq_gen(sort1_npy, sam1)
        itr2 = self._uniq_gen(sort2_npy, sam2)
    
    
        # first git.
        u1 = itr1.next()
        u2 = itr2.next()
 
        print u1, u2
        
        cnt = 0
        while u1 != None and u2 != None:
            # peek for a match.
            if u1['name'] == u2['name']:
    
                # seek to it.
                sam1.seek(u1['row'])
                sam2.seek(u2['row'])
                out1.write(sam1.readline())
                out2.write(sam2.readline())
    
                # change both.
                u1 = itr1.next()
                u2 = itr2.next()
    
            else:
                # pop smaller.
                if u1['name'] > u2['name']:
                    u2 = itr2.next()
                else:
                    u1 = itr1.next()
    
            # die after 5
            cnt += 1
            #if cnt > 5: break
    
        # close them.
        sam1.close()
        sam2.close()
        out1.close()
        out2.close()
    
    
    
    
    
    def _names(self, file_name=None, data=None, size=None, name_size=None):
        """ returns pointer to mapped file """
    
        if size != None and name_size != None:
            return np.zeros(size,  dtype=np.dtype([('name','S%d' % name_size),('row',np.int)]))
        elif file_name != None and data == None:
            return np.load(file_name)
        elif file_name != None and data != None:
            np.save(file_name, data)
        else:
            logging.error("bad usage")
            sys.exit(1)
    
    def _name_size(self, file_path):
        """ guess string size """
        # determine name size.
        print file_path
        with open(file_path, "rb") as fin:
            for line1 in fin:
                if line1[0] == '@': continue
                name_size = 100 #len(line1.split("\t")[0]) + 10
                break
    
        return name_size
    
    def _extract_names(self, file_name, name_size, key_size):
        """ builds numpy array of name hits"""
    
        # count lines.
        logging.info("reading lines")
        with open(file_name, "rb") as fin:
            size = 0
            for line in fin:
                if line[0] == '@': continue
                size += 1
                
        print size, "THIS IS THE SIZE"
                #if size > 10000000: break
    
        # allocate array.
        names = self._names(size=size, name_size=name_size)
    
        # copy data into array.
        logging.info("copying data")
        with open(file_name, "rb") as fin:
    
            offset = 0
            idx = 0
            for line1 in fin:
    
                # skip header.
                if line1[0] == '@':
                    offset += len(line1)
                    continue
    
                # tokenize.
                tokens = line1.split("\t")
                
    
                # skip no map.
                if tokens[2] == "*":
                    offset += len(line1)
                    print line1, "DDDDDDDDDDDDDDDDDDDDDDDDDDD"
                    continue
    
    
    
    
                # operate.
                if key_size == 0:
                    names[idx]['name'] = tokens[0]
                else:
                    names[idx]['name'] = tokens[0][0:-key_size]
                names[idx]['row'] = offset
    
                # reset.
                idx += 1
                offset += len(line1)
    
        # resize.
        names.resize(idx)
        # return the size.
        return names
    
    
    
    def _uniq_gen(self, names_npy, sam):
        """ generator for unique reads in list """
    
        # create mmap object.
        mmap = np.load(names_npy, mmap_mode='r')
    
        # setup buffered loop.
        buffstep = 10000000
        buffo = 0
        buffs = buffstep
        if buffo + buffstep > mmap.shape[0]:
            buffs = mmap.shape[0] - buffo
    
        # buffer loop.
        while buffo < mmap.shape[0]:
            # make buffer.
            logging.info("unique: buffering: %d %d" % (buffo, buffs))
            names = mmap[buffo:buffs]
            # iterate over non-boundry cases.
            for i in range(1, names.shape[0]-1):
                # must not match its neighbors.
                if names[i-1]['name'] != names[i]['name'] and names[i+1]['name'] != names[i]['name']:
                    tokens = sam.readline().split("\t")
                    if len(tokens) > 9:
                        cov = self._cov.get(tokens[2], 0)
                        self._cov[tokens[2]] = cov + len(tokens[9])
                    yield names[i]
#                 else:
#                     # annotate repeat.
#                     sam.seek(names[i]['row'])
#                     tokens = sam.readline().split("\t")
#                     ctg = tokens[2]
#                     start = int(tokens[3])
#                     stop = start + len(tokens[9])
                buffo += 1
            # check the first one.
            if names[0]['name'] != names[1]['name']:
                tokens = sam.readline().split("\t")
                if len(tokens) > 9:
                    cov = self._cov.get(tokens[2], 0)
                    self._cov[tokens[2]] = cov + len(tokens[9])
                yield names[i]
#             else:
#                 sam.seek(names[i]['row'])
#                 tokens = sam.readline().split("\t")
#                 ctg = tokens[2]
#                 start = int(tokens[3])
#                 stop = start + len(tokens[9])
            buffo += 1
    
            # check the last one.
            if names[-1]['name'] != names[-2]['name']:
                tokens = sam.readline().split("\t")
                if len(tokens) > 9:
                    cov = self._cov.get(tokens[2], 0)
                    self._cov[tokens[2]] = cov + len(tokens[9])
                yield names[i]
#             else:
#                 sam.seek(names[i]['row'])
#                 tokens = sam.readline().split("\t")
#                 ctg = tokens[2]
#                 start = int(tokens[3])
#                 stop = start + len(tokens[9])
            buffo += 1
    
            # update for buffer.
            if buffo + buffstep > mmap.shape[0]:
                buffs = buffo + (mmap.shape[0] - buffo)
            else:
                buffs = buffo + buffstep
                # yield poison pill
        yield None
