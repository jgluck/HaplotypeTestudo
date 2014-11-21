#!/usr/bin/env python

import pysam
from time import sleep
from operator import itemgetter
import glob
import pickle
import psycopg2
import os.path
import csv

SAM_FILE = pysam.Samfile("alignments/sorted/EMB1.rm.bam_sorted.bam.bam", "rb")


##Constants
C_C = 'c'
A_C = 'a'
T_C = 't'
G_C = 'g'


A_CONST = ['a', 'A']
C_CONST = ['c', 'C']
G_CONST = ['g', 'G']
T_CONST = ['t', 'T']
DOT_CONST = ['.', ',']
DEL_CONST = ['^', '$']

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
        options = [(C_C,self.c), (A_C,self.a), (T_C,self.t), (G_C,self.g)]
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
            elif ln[i] in A_CONST:
                self.a += 1
            elif ln[i] in C_CONST:
                self.c += 1
            elif ln[i] in T_CONST:
                self.t += 1
            elif ln[i] in G_CONST:
                self.g += 1
            elif ln[i] in DOT_CONST:
                self.dot += 1
            elif ln[i] in DEL_CONST:
                self.dot += 1
                i += 1
            i += 1

class PileupReader:
    def __init__(self,f_name):
        self.file_name = f_name
        self.turtle = os.path.basename(self.file_name).split('.')[0] 
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

def handle_csv_line(csv_line, turtles, cursor):
    turt = csv_line[0]
    seq = csv_line[1]
    pos = csv_line[2]
    ref = csv_line[3]

    a_col = turtles[turt][seq][int(pos)]
    assert int(pos) == int(a_col.pos)
    assert seq == a_col.seq
    query = "INSERT INTO variants (turtle,scaffold,pos,ref,alt,qual,dp,af1,backfilled) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);"
    data = (turt,seq,pos,a_col.ref,a_col.alt,"0","0","0","t")
    cursor.execute(query,data)


def main():
    turtles = {}
    conn=psycopg2.connect("host=biords.cjt20mcdfhfm.us-east-1.rds.amazonaws.com port=5432 user=master password=tomswift dbname=haplotypes")
    cursor = conn.cursor()

    for pileup_file in get_pileup_list():
        bar(75)
        print "Begging work on file %s" %(pileup_file)
        cur_dict = {} 
        cur_reader = PileupReader(pileup_file)
        for a_col in cur_reader:
            if a_col.seq in cur_dict:
                cur_dict[a_col.seq][int(a_col.pos)] = a_col
            else:
                cur_dict[a_col.seq] = [None] * 145000
                cur_dict[a_col.seq][int(a_col.pos)] = a_col
        turtles[cur_reader.turtle] = cur_dict
    with open("to_backfill.csv", "rb") as b_csv:
        b_reader = csv.reader(b_csv, delimiter=',')
        for line in b_reader:
            try:
                handle_csv_line(line,turtles,cursor)
            except Exception, e:
                print "Line was: \n\t%s" %(" ".join(line))
                print turtles.keys()
                import traceback
                print traceback.format_exc()

    conn.commit()

def bar(num):
    print "~" * num

if __name__=="__main__":
    main()



