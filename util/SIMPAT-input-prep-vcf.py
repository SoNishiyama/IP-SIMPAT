#!/usr/bin/env python3

# 2021/01/26 SN

# usage : python this.py input.vcf output-prefix maternal-parent-id

# This script goes through genotypes of input.vcf
# then convert each genotype to numeric representation.
# Output will be formatted for the "maternal-parent-id",
# where all the allele possessed by the mother will be formatted as 0.
# Loci in which the mother possess as heterozygous will be removed from the output.

# The output will be used for the paternity inference with the IP-SIMPAT Matlab codes.

# 0/0 -> 0
# 0/1 -> 1
# 1/1 -> 2
# ./. -> nan

# The other genotype codes won't be recognized and the script will be terminated.


import sys
from statistics import mean


def format_by_mother(gen, mother_id):
    mother_gen = gen[mother_id]
    if mother_gen == 1 or mother_gen == "nan":
        out = 9
    else:
        gen_m = [gen.pop(mother_id)]
        gen_m.extend(gen)
        if mother_gen == 2:
            gen_m = ["t" if i == 2 else i for i in gen_m]
            gen_m = [2 if i == 0 else i for i in gen_m]
            out = [0 if i == "t" else i for i in gen_m]
        else:
            out = gen_m
    return out


def main():
    v = open(sys.argv[1])
    out_mat = open(sys.argv[2]+".num.txt", 'w')
    out_pos = open(sys.argv[2]+".pos.txt", 'w')
    mo = sys.argv[3]

    cnt = 0
    while True:
        x1 = v.readline()
        cnt += 1
        x2 = x1.rstrip().split("\t")
        if x2[0] == "#CHROM":
            libs = x2[9:]
            try:
                mo_id = libs.index(mo)
            except ValueError:
                print(mo+" is not in the input.")
                sys.exit()
            libs_mo = [libs.pop(mo_id)]
            libs_mo.extend(libs)
            out_mat.write("variant\t"+"\t".join(libs_mo)+"\n")
            out_pos.write("variant\tchr\tpos\n")
            break

    prev = ""
    pn = 1
    while True:
        x1 = v.readline()
        if x1 == "":
            break
        cnt += 1
        if cnt % 10000 == 0:
            sys.stderr.write('Checking records: {}\r'.format(cnt))
        x2 = x1.rstrip().split("\t")
        pos = ":".join(x2[0:2])
        if prev != pos:
            pn = 1
        else:
            pn += 1
        var = pos+":"+str(pn)
        out = []
        prev = ":".join(x2[0:2])
        for x3 in x2[9:]:
            if x3[0:3] == "0/0":
                out.append(0)
            elif x3[0:3] == "0/1":
                out.append(1)
            elif x3[0:3] == "1/1":
                out.append(2)
            elif x3[0:3] == "./.":
                out.append("nan")
            else:
                print("can not recognize genotype in Line "+str(cnt)+": "+x3)
                sys.exit()
        gen = format_by_mother(out, mo_id)
        if gen != 9:
            out_mat.write(var+"\t"+"\t".join([str(i) for i in gen])+"\n")
            out_pos.write("\t".join([var, x2[0], x2[1]]) + "\n")

    v.close()
    out_mat.close()
    out_pos.close()


if __name__ == '__main__':
    main()
