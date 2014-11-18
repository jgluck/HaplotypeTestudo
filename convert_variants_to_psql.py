#!/usr/bin/env python
import glob
import os.path as path

def main():
    fn_list = glob.glob('variants/*.vcf')
    for fn in fn_list:
        parse_file(fn)

def parse_file(fn):
    with open(fn, 'r') as v_file:
        child = path.split(fn)[1].split('.')[0] 
        for line in v_file:
            if line[0] == "#":
                continue
            else:
                #process line
                sp_line = line.split("\t")
                info_slot = sp_line[7]
                sp_info = info_slot.split(';')
                scaff = sp_line[0]
                pos = sp_line[1]
                ref = sp_line[3]
                alts = sp_line[4]
                alts = alts.split(",")
                qual = sp_line[5]
                dp = ""
                af1 = ""
                for field in sp_info:
                    sp_field = field.split("=")
                    if sp_field[0] == "DP":
                        dp = sp_field[1]
                    elif sp_field[0] == "AF1":
                        af1 = sp_field[1]

            dlm = '~'
            for alt in alts:
                print (child + dlm +
                    scaff + dlm +
                    pos + dlm +
                    ref + dlm +
                    alt + dlm +
                    qual + dlm +
                    dp + dlm +
                    af1)



if __name__ == "__main__":
    main()
