#!/usr/bin/env python

import csv

mother = {}
vote_file = open("vote_thing.csv","wb")
vote_writer = csv.writer(vote_file, delimiter = ',')

class vote:

    def __init__(self,spl_ln):
        self.bucket = spl_ln[0]
        self.scaffold = spl_ln[1]
        self.pos = spl_ln[2]
        self.ref = spl_ln[3]
        self.bp = spl_ln[4]
        self.count = int(spl_ln[5])

    def __str__(self):
        s = "\tBucket: %s\n\
                \tScaffold: %s\tPos: %s\n\
                \tRef:%s\tBP:%s\tCount:%d\n" % (self.bucket,self.scaffold,\
                self.pos, self.ref, self.bp, self.count)



def main():
    print "We assume that if we have a vote for one half, we can infer the other"
    read_mother()
    process_child()

def read_mother():
    with open('mother_variants.csv', 'rb') as mothercsv:
        motherreader = csv.reader(mothercsv, delimiter=',')
        for row in motherreader:
            if row[2] not in mother:
                mother[row[2]] = {}
            mother[row[2]][row[3]] = (row[4],row[5])


def read_vote(fp):
    while True:
        v1 = fp.readline()
        v2 = fp.readline()
        print v1
        print v2
        if not v1 or not v2:
            break
        assert v1[0] == '1' and v2[0] == '2'
        v1_split = v1.split(',')
        v2_split = v2.split(',')
        for pos in [1,2]:
            assert v1_split[pos] == v2_split[pos]
        yield vote(v1_split),vote(v2_split)


def output_votes(v1,v2):
    vote_writer.writerow([v1.bucket,v1.scaffold,v1.pos,v1.ref,v1.bp,v1.count])
    vote_writer.writerow([v2.bucket,v2.scaffold,v2.pos,v2.ref,v2.bp,v2.count])


def process_child():
    skipped = 0
    fixed = 0
    fine = 0
    with open('bp_chosen.csv','r') as votes:
        for v1,v2 in read_vote(votes):
            should_output = True
            if v1.bp == '-':
                if v2.bp != '-':
                    if v1.pos in mother[v1.scaffold]:
                        mother_vote1, mother_vote2 = mother[v1.scaffold][v1.pos]
                        if mother_vote1 == v2.bp:
                            v1.bp = mother_vote2
                            v1.count = 0
                            fixed += 1
                        elif mother_vote2 == v2.bp:
                            v1.bp = mother_vote1
                            v1.count = 0
                            fixed += 1
                        else:
                            print "Uh oh at:" + str(v1) + str(v2)
                            should_output = False
                            skipped += 1
                    #get mother's other for v1
                else:
                    print "Skipping, both votes bad"
                    should_output = False
                    skipped += 1
                    #skip both are bad
            elif v2.bp == '-':
                if v1.pos in mother[v1.scaffold]:
                        mother_vote1, mother_vote2 = mother[v1.scaffold][v1.pos]
                        if mother_vote1 == v1.bp:
                            v2.bp = mother_vote2
                            v2.count = 0
                            fixed += 1
                        elif mother_vote2 == v1.bp:
                            v2.bp = mother_vote1
                            v2.count = 0
                            fixed += 1
                        else:
                            print "Uh oh at:" + str(v1) + str(v2)
                            skipped += 1
                            should_output = False
            else:
                fine += 1
                #just output
            if should_output:
                output_votes(v1,v2)
        print "Skipped: %d" % (skipped)
        print "Fixed: %d" % (fixed)
        print "Fine: %d" % (fine)

if __name__=='__main__':
    main()
