#!/usr/bin/env python

import pysam
from time import sleep
from operator import itemgetter
import glob

SAM_FILE = pysam.Samfile("alignments/sorted/EMB1.rm.bam_sorted.bam.bam", "rb")



def lookup_base(scaf, base):
    p_ups = SAM_FILE.pileup(reference=scaf,start=base-150,end=base+150)
    for p_up in p_ups:
        if p_up.pos+1 == base:
            print "Coverage: %s" %(p_up.n)
            #The correct base
            for pileupread in p_up.pileups:
                #print dir(pileupread.alignment)
                #print pileupread.alignment.query[pileupread.qpos]
                #print pileupread.alignment.aligned_pairs
                #print dir(pileupread)
                base = pileupread.alignment.seq[pileupread.qpos]
                if not pileupread.is_del:
                    print pileupread.indel
                    #base = "*"
                    print "Base in read %s = %s" % (pileupread.alignment.qname,\
                        base)



class AlignmentColumn:
    def __init__(self,line):
        try:
            spl_line = line.split()
            self.seq = spl_line[0]
            self.pos = spl_line[1]
            self.ref = spl_line[2]
            self.n = spl_line[3]
            self.c = 0
            self.a = 0
            self.t = 0
            self.g = 0
            self.dot = 0
            self.ins = [] 
            self.dl = []
            self.alt = "N"
            if int(self.n) != 0:
                self.parse_line(spl_line[4])
                self.find_alt()
        except:
            print line
            exit()
            
    def __str__(self):
        x = "%s:%s\t%s\t%s" % (self.seq,self.pos,self.ref,self.alt)
        return x

    def find_alt(self):
        options = [('c',self.c), ('a',self.a), ('t',self.t), ('g',self.g)]
        best = max(options, key = itemgetter(1))
        if float(best[1])/float((self.dot or 1)) >= .25:
            self.alt = best[0]
            

    def parse_line(self, ln):
        i = 0
        while i< len(ln):
            if ln[i] == "+":
                num_char = int(ln[i+1])
                ins = ""
                for offs in range(int(num_char)):
                    ins += ln[i+2+offs]
                i+= num_char
                self.ins.append(ins)
            elif ln[i] == "-":
                num_char = int(ln[i+1])
                dl = ""
                for offs in range(int(num_char)):
                    dl += ln[i+2+offs]
                i+= num_char
                self.dl.append(dl)
            elif ln[i] in ['a','A']:
                self.a += 1
            elif ln[i] in ['c','C']:
                self.c += 1
            elif ln[i] in ['t','T']:
                self.t += 1
            elif ln[i] in ['g','G']:
                self.g += 1
            elif ln[i] in ['.',',']:
                self.dot += 1
            elif ln[i] in ['^','$']:
                self.dot += 1
                i += 1
            i += 1

class PileupReader:
    def __init__(self,f_name):
        self.file_name = f_name
        self.turtle = self.file_name.split('.')[0]
        self.fp = open(self.file_name, 'r')


    def __iter__(self):
        return self

    def next(self):
        ln = self.fp.readline()
        while ln:
            if ln[0] == 'S':
                #GOOD
                align_col = AlignmentColumn(ln)
                return align_col
            self.fp.readline()
        raise StopIteration()


def get_pileup_list():
    fn_list = glob.glob('pileups/*_pileup')
    return fn_list

def main():
    turtles = {}
    for pileup_file in get_pileup_list():
        bar(75)
        print "Begging work on file %s" %(pileup_file)
        cur_arr = [] 
        cur_reader = PileupReader(pileup_file)
        for a_col in cur_reader:
            cur_arr.append(a_col)
        turtles[cur_reader.turtle] = cur_arr

    print turtles['EMB2'][3]

def bar(num):
    print "~" * num

if __name__=="__main__":
    main()



