#!/usr/bin/env python3
#  2021/02/25 SN

# This will compare the inference of cervus and the true cross combination which is from "cervus_input_prep.py".
# Input files will be specified as a list with two columns.
# The first column is assumed to be the cervus output.csv
# The secound is assumed to be the ".truecross.txt" file from cervus_input_prep.py

# usage: python this.py input.list.txt output.tsv

import sys
import pandas as pd


def validation(cervus, truecross):
    inp_cv = pd.read_csv(cervus, index_col=False)
    inp_tc = pd.read_table(truecross)

    d_tc_pat = dict()
    total = inp_tc.shape[0]
    for i in range(total):
        d_tc_pat[inp_tc.iloc[i, 0]] = str(inp_tc.iloc[i, 2])

    st_TP = 0
    st_FP = 0
    rl_TP = 0
    rl_FP = 0
    ml_TP = 0
    ml_FP = 0

    head = inp_cv.columns.to_list()
    id_ofs = head.index("Offspring ID")
    id_pat = head.index("Candidate father ID")
    id_conf = head.index("Trio confidence")

    for i in range(inp_cv.shape[0]):
        ofs = inp_cv.iloc[i, id_ofs]
        cand_pat = str(inp_cv.iloc[i, id_pat])
        conf = inp_cv.iloc[i, id_conf]
        if conf == "*":
            if d_tc_pat[ofs] == cand_pat:
                st_TP += 1
                rl_TP += 1
                ml_TP += 1
            else:
                st_FP += 1
                rl_FP += 1
                ml_FP += 1
        elif conf == "+":
            if d_tc_pat[ofs] == cand_pat:
                rl_TP += 1
                ml_TP += 1
            else:
                rl_FP += 1
                ml_FP += 1
        elif conf == "-":
            if d_tc_pat[ofs] == cand_pat:
                ml_TP += 1
            else:
                ml_FP += 1
    return [total, st_TP, st_FP, rl_TP, rl_FP, ml_TP, ml_FP]


def main():
    inp_list = open(sys.argv[1])
    out = open(sys.argv[2], 'w')
    head = ["cervus_out", "truecross", "total", "strict-TP", "strict-FP", "relaxed-TP", "relaxed-FP", "most_likely-TP",
            "most_likely-FP"]
    out.write("\t".join(head) + "\n")
    while True:
        line = inp_list.readline()
        if line == "":
            break
        files = line.rstrip().split()
        res = validation(files[0], files[1])
        out.write(files[0]+"\t"+files[1]+"\t")
        out.write("\t".join([str(i) for i in res])+"\n")


if __name__ == '__main__':
    main()