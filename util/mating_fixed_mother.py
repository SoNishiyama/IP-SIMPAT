#!/usr/bin/env python3
#  2021/01/22 SN

# usage: python this.py genmat.tsv mother-id father-id-list.txt marker-list num-progeny output_prefix

import sys
import Cervus_prep_util as Cp


def output(geno):
    out = []
    for gen in geno:
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
    return out


def main():
    genmat = open(sys.argv[1])
    mother = sys.argv[2]
    father = open(sys.argv[3])
    marlist = open(sys.argv[4])
    num_prog = int(sys.argv[5])
    out_gen = open(sys.argv[6] + ".cervus.gen.txt", 'w')
    out_prg = open(sys.argv[6] + ".offspring.txt", 'w')
    out_tru = open(sys.argv[6] + ".truecross.txt", 'w')

    indlist = set()
    indlist.add(mother)
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
        out_gen.write(indv[i]+"\t"+"\t".join(output(out))+"\n")

    mating = Cp.Mating(d, indv)

    # outputting progeny data
    c = 1
    for i in range(len(patlist)):
        for j in range(num_prog):
            prog = mating.mating(mother, patlist[i], mar)
            out = output(prog)
            prog_name = "OFS"+str(c)
            out_gen.write(prog_name+"\t"+"\t".join(out)+"\n")
            out_prg.write(prog_name+"\t"+mother+"\t"+"\t".join(patlist)+"\n")
            out_tru.write(prog_name+"\t"+mother+"\t"+patlist[i]+"\n")
            c += 1


if __name__ == '__main__':
    main()
