#!/usr/bin/env python3
#  2021/01/22 SN

import numpy as np


class Cervus_prep:
    def __init__(self, genmat, indlist, marlist):
        self.genmat = genmat
        self.indlist = indlist
        self.marlist = marlist

    def ind_prep(self):
        indv = []
        for i in self.indlist:
            indv.append(i.rstrip())
        return indv

    def mar_prep(self):
        mar = []
        for i in self.marlist:
            mar.append(i.rstrip())
        return mar

    def parent_genmat_prep(self):
        indv = self.ind_prep()

        mar = self.mar_prep()

        head = self.genmat.readline().rstrip().split()
        indv_in = []
        indv_out = []
        for i in head[1:]:
            if i in indv:
                indv_in.append(1)
                indv_out.append(i)
            else:
                indv_in.append(0)

        d = [dict() for i in range(len(indv))]

        while True:
            x1 = self.genmat.readline()
            if x1 == "":
                break
            x2 = x1.rstrip().split()
            curr_mar = x2[0]
            curr_gen = x2[1:]
            if curr_mar in mar:
                c = 0
                for i in range(len(curr_gen)):
                    x3 = curr_gen[i]
                    if indv_in[i] == 1:
                        try:
                            d[c][curr_mar] = int(x3)
                        except ValueError:
                            d[c][curr_mar] = 9
                        c += 1
        self.genmat.close()
        return indv_out, mar, d


class Mating:
    def __init__(self, gen_d, indlist):
        self.gen_d = gen_d
        self.indlist = indlist

    def meiosis(self, gen):
        gamete = []
        for i in gen:
            if i == 0:
                gamete.append(0)
            elif i == 1:
                gamete.append(np.random.binomial(1, 0.5))
            elif i == 2:
                gamete.append(1)
            elif i == 9:
                gamete.append(9)
        return gamete

    def mating(self, parent1, parent2, mar):
        gen_d1 = self.gen_d[self.indlist.index(parent1)]
        gen_d2 = self.gen_d[self.indlist.index(parent2)]

        p1 = []
        p2 = []
        for i in mar:
            p1.append(gen_d1[i])
            p2.append(gen_d2[i])
        m1 = self.meiosis(p1)
        m2 = self.meiosis(p2)
        progeny = np.array(m1) + np.array(m2)
        progeny[progeny > 3] = 9
        return progeny




