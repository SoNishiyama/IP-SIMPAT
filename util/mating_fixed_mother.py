#!/usr/bin/env python3
#  2021/01/22 SN

# Simulated crossing between a maternal individual and paternal individuals.
# Output genotype matrix for the specified markers

# 22/04/22 add -m option -m: frequency of masking operation. float value between 0-1, default 0. Converting genotype
# data into missing genotype in a defined frequency. This would be used to test the robustness for missing genotype.

# usage: python this.py genmat.tsv mother-id father-id-list.txt marker-list num-progeny output_prefix [-m float]

import Cervus_prep_util as Cp
import argparse
import random

def output(geno,m):
    out = []
    for gen in geno:
        if gen == 0:
            if m < random.random():
                out.append("1")
                out.append("1")
            else:
                out.append("0")
                out.append("0")
        elif gen == 1:
            if m < random.random():
                out.append("1")
                out.append("2")
            else:
                out.append("0")
                out.append("0")
        elif gen == 2:
            if m < random.random():
                out.append("2")
                out.append("2")
            else:
                out.append("0")
                out.append("0")
        elif gen == 9:
            out.append("0")
            out.append("0")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gen", help="input genotype matrix", type=str)
    parser.add_argument("mother", help="mother ID", type=str)
    parser.add_argument("f", help="list of father ID", type=str)
    parser.add_argument("mar", help="input marker list", type=str)
    parser.add_argument("num_ofs", help="number of offspring per cross", type=int)
    parser.add_argument("out", help="output prefix", type=str)
    parser.add_argument("-m", "--missing", dest="m", type=float, default=0,
                        help="Frequency of masking operation")
    args = parser.parse_args()

    genmat = open(args.gen)
    marlist = open(args.mar)
    father = open(args.f)
    out_gen = open(args.out + ".cervus.gen.txt", 'w')
    out_prg = open(args.out + ".offspring.txt", 'w')
    out_tru = open(args.out + ".truecross.txt", 'w')

    indlist = set()
    indlist.add(args.mother)
    patlist = []
    for i in father:
        indlist.add(i.rstrip())
        patlist.append(i.rstrip())

    prep = Cp.Cervus_prep(genmat, indlist, marlist)
    indv, mar, d = prep.parent_genmat_prep()

    out_head = ["individual"]
    for i in mar:
        out_head.append(str(i) + "a")
        out_head.append(str(i) + "b")
    out_gen.write("\t".join(out_head) + "\n")

    out_prg.write("\t".join(["Offspring ID", "Known mother ID", "Candidate father IDs"+"\n"]))
    out_tru.write("\t".join(["Offspring ID", "true mother ID", "true father IDs"+"\n"]))

    # outputting parent data
    for i in range(len(indv)):
        out = []
        for j in range(len(mar)):
            out.append(d[i][mar[j]])
        out_gen.write(indv[i]+"\t"+"\t".join(output(out,0))+"\n")

    mating = Cp.Mating(d, indv)

    # outputting progeny data
    c = 1
    for i in range(len(patlist)):
        for j in range(args.num_ofs):
            prog = mating.mating(args.mother, patlist[i], mar)
            out = output(prog,args.m)
            prog_name = "OFS"+str(c)
            out_gen.write(prog_name+"\t"+"\t".join(out)+"\n")
            out_prg.write(prog_name+"\t"+args.mother+"\t"+"\t".join(patlist)+"\n")
            out_tru.write(prog_name+"\t"+args.mother+"\t"+patlist[i]+"\n")
            c += 1


if __name__ == '__main__':
    main()
