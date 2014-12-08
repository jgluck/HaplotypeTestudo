#!/usr/bin/env python

import sys
import csv
import os.path
import glob
import pysam
from operator import itemgetter 

variants = {}

haplotype_file_name = ["vote_fix_one_case.csv", "candidates_pg-1.csv"]

haplotypes = {}
haplotypes['1'] = {}
haplotypes['2'] = {}

turtles = {"EMB1":"2","EMB2":"1","EMB3":"1","EMB4":"1","EMB5":"1","EMB7":"1","EMB8":"2","EMB9":"2","EMB10":"2","EMB11":"2","EMB12":"2","EMB13":"2","EMB14":"1","EMB15":"1"}

t_pile = {}

SAM_FILE = pysam.Samfile("alignments/sorted/EMB1.rm.bam_sorted.bam.bam", "rb")

##Constants
C_C = 'C'
A_C = 'A'
T_C = 'T'
G_C = 'G'


A_CONST = ['a', 'A']
C_CONST = ['c', 'C']
G_CONST = ['g', 'G']
T_CONST = ['t', 'T']
DOT_CONST = ['.', ',']
DEL_CONST = ['^', '$']


def main():
    read_variants()
    load_pileup()
    for fn in haplotype_file_name:
        read_haplotypes(fn)
        for turtle in turtles.keys():
            handle_turtle(fn,turtle)


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
        except Exception as e:
            print str(e)
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
    
    def get_alts(self):
        alts = []
        for x in [(C_C,self.c),(A_C,self.a),(T_C,self.t),(G_C,self.g)]:
            cat,quant = x
            if quant > 0:
                alts.append((cat,quant))
        return alts


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
            ln = self.fp.readline()
        raise StopIteration()




def get_pileup_list():
    fn_list = glob.glob('pileups/*_pileup')
    return fn_list



def load_pileup():
    for pileup_file in get_pileup_list():
        bar(75)
        print "Beginning work on file %s" %(pileup_file)
        cur_dict = {} 
        cur_reader = PileupReader(pileup_file)
        for a_col in cur_reader:
            if a_col.seq in cur_dict:
                cur_dict[a_col.seq][int(a_col.pos)] = a_col
            else:
                cur_dict[a_col.seq] = [None] * 145000
                cur_dict[a_col.seq][int(a_col.pos)] = a_col
        t_pile[cur_reader.turtle] = cur_dict


def handle_turtle(fn,turtle):
    dirname = os.path.basename(fn)
    dirname = dirname.split(".")[0]
    t_name = "%s/%s_haplotypes.tsv" % (dirname,turtle)
    t_group = turtles[turtle]
    with open(t_name,'w') as t_file:
        for scaffold in haplotypes[t_group]:
            for pos in haplotypes[t_group][scaffold]:
                h_answer = haplotypes[t_group][scaffold][pos]
                v_answer = "-"
                p_group = []
                try:
                    p_col = t_pile[turtle][scaffold][int(pos)]
                    p_group = p_col.get_alts()
                except Exception, e:
                    import traceback
                    print traceback.format_exc()
                    print "Failed pileup things"
                
                p_group = [ s[0] for s in p_group]
                if h_answer not in p_group:
                    print "Uh oh!: %s not in %s" % (h_answer, str(p_group))
                    v_answer = ",".join(p_group)
                else:
                    p_group.remove(h_answer)
                    v_answer = ",".join(p_group)
                if not v_answer:
                    v_answer = h_answer
                write_line = "%s\t%s\t%s\n" % (scaffold+","+pos,h_answer,v_answer)
                t_file.write(write_line)

                    

def read_variants():
    with open('variants.csv', 'rb') as variant_file:
        variant_reader = csv.reader(variant_file, delimiter=',')
        for row in variant_reader:
            if row[1] not in variants:
                variants[row[1]] = {}
            if row[2] not in variants[row[1]]:
                variants[row[1]][row[2]] = {}
            if row[3] not in variants[row[1]][row[2]]:
                variants[row[1]][row[2]][row[3]] = []
            variants[row[1]][row[2]][row[3]].append(row[5])


def read_haplotypes(name):
    with open(name, 'rb') as haplotype_file:
        haplotype_reader = csv.reader(haplotype_file, delimiter=",")
        for row in haplotype_reader:
            if row[1] not in haplotypes[row[0]]:
                haplotypes[row[0]][row[1]] = {}
            haplotypes[row[0]][row[1]][row[2]] = row[4]
            
def bar(n):
    print "=" * n

if __name__ == "__main__":
    main()
