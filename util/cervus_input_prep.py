#!/usr/bin/env python3
#  2021/01/22 SN

# Convert genotype matrix to the cervus-compatible format.

# 22/04/22 add -m option -m: frequency of masking operation. float value between 0-1, default 0. Converting genotype
# data into missing genotype in a defined frequency. This would be used to test the robustness for missing genotype.

# usage: python this.py input.genmat ind-name-list marker-list output_prefix [-m float]

import Cervus_prep_util as Cp
import argparse
import random


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gen", help="input genotype matrix", type=str)
    parser.add_argument("ind", help="input individual list", type=str)
    parser.add_argument("mar", help="input marker list", type=str)
    parser.add_argument("out", help="output prefix", type=str)
    parser.add_argument("-m", "--missing", dest="m", type=float, default=0,
                        help="Frequency of masking operation")
    args = parser.parse_args()

    genmat = open(args.gen)
    indlist = open(args.ind)
    marlist = open(args.mar)
    o = open(args.out + ".cervus.gen.txt", 'w')

    prep = Cp.Cervus_prep(genmat, indlist, marlist)
    indv, mar, d = prep.parent_genmat_prep()

    out_head = ["individual"]
    for i in mar:
        out_head.append(str(i) + "a")
        out_head.append(str(i) + "b")
    o.write("\t".join(out_head)+"\n")

    for i in range(len(d)):
        out = [indv[i]]
        for j in mar:
            gen = d[i][j]
            if gen == 0:
                if args.m < random.random():
                    out.append("1")
                    out.append("1")
                else:
                    out.append("0")
                    out.append("0")
            elif gen == 1:
                if args.m < random.random():
                    out.append("1")
                    out.append("2")
                else:
                    out.append("0")
                    out.append("0")
            elif gen == 2:
                if args.m < random.random():
                    out.append("2")
                    out.append("2")
                else:
                    out.append("0")
                    out.append("0")
            elif gen == 9:
                out.append("0")
                out.append("0")
        o.write("\t".join(out)+"\n")


if __name__ == '__main__':
    main()
