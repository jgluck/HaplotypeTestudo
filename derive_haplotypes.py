#!/usr/bin/env python

import csv

variants = {}

haplotype_file_name = "vote_thing.csv"

haplotypes = {}
haplotypes['1'] = {}
haplotypes['2'] = {}

turtles = {"EMB1":"2","EMB2":"1","EMB3":"1","EMB4":"1","EMB5":"1","EMB7":"1","EMB8":"2","EMB9":"2","EMB10":"2","EMB11":"2","EMB12":"2","EMB13":"2","EMB14":"1","EMB15":"1"}

def main():
    read_variants()
    read_haplotypes()
    for turtle in turtles.keys():
        handle_turtle(turtle)


def handle_turtle(turtle):
    t_name = "%s_haplotypes.tsv" % (turtle)
    t_group = turtles[turtle]
    with open(t_name,'w') as t_file:
        for scaffold in haplotypes[t_group]:
            for pos in haplotypes[t_group][scaffold]:
                h_answer = haplotypes[t_group][scaffold][pos]
                v_answer = "-"
                try:
                    variant_arr = variants[turtle][scaffold][pos]
                    if h_answer not in variant_arr:
                        print "Uh oh!: %s not in %s" % (h_answer, str(variant_arr))
                    else:
                        variant_arr.remove(h_answer)
                        v_answer = ",".join(variant_arr)
                except:
                    print "Failed something"
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


def read_haplotypes():
    with open(haplotype_file_name, 'rb') as haplotype_file:
        haplotype_reader = csv.reader(haplotype_file, delimiter=",")
        for row in haplotype_reader:
            if row[1] not in haplotypes[row[0]]:
                haplotypes[row[0]][row[1]] = {}
            haplotypes[row[0]][row[1]][row[2]] = row[4]
            

if __name__ == "__main__":
    main()
