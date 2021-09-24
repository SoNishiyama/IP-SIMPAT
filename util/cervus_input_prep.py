#!/usr/bin/env python3
#  2021/01/22 SN

# usage: python this.py input.genmat ind-name-list marker-list output_prefix

import sys
import Cervus_prep_util as Cp


def main():
    genmat = open(sys.argv[1])
    indlist = open(sys.argv[2])
    marlist = open(sys.argv[3])
    o = open(sys.argv[4]+".cervus.gen.txt", 'w')

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
                out.append("1")
                out.append("1")
            elif gen == 1:
                out.append("1")
                out.append("2")
            elif gen == 2:
                out.append("2")
                out.append("2")
            elif gen == 9:
                out.append("0")
                out.append("0")
        o.write("\t".join(out)+"\n")


if __name__ == '__main__':
    main()
