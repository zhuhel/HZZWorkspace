#!/usr/bin/env python3

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.locator_params(axis='y', nticks=10)

from helper import get_limit_dic

class CompareLimits:
    def __init__(self, file_names, tag_names):
        self.limit_list = []
        for file_name in file_names:
            self.limit_list.append(get_limit_dic(file_name))

        self.tag_names = tag_names

    @staticmethod
    def plot_cmp(a_list, b_list, out_name, tag_names):
        a_x, a_y, a_y_e = a_list
        b_x, b_y, b_y_e = b_list # only toy has error
        print(len(a_x),len(a_y))
        diff_y = [(x-y)/x for x,y in zip(a_y, b_y)]
        diff_y_e = [y/x for x, y in zip(a_y, b_y_e)]
        for x, y in zip(a_x, diff_y): 
            print("{} {:.4f}".format(x,y))

        print(type(a_x),type(a_y))
        line_a, = plt.plot(a_x, a_y, 'bo', lw=2)
        line_b, = plt.plot(b_x, b_y, 'ro', lw=2)
        plt.xlabel("m_{H} [GeV]")
        plt.ylabel("Upper limit")
        plt.legend([line_a, line_b], tag_names)
        plt.savefig(out_name)
        plt.close()

        plt.errorbar(a_x, diff_y, yerr=diff_y_e, fmt='o', lw=2)
        plt.xlabel("$m_H$ [GeV]")
        plt.ylabel("(asym - toys)/asym")
        plt.savefig("diff_"+out_name)
        plt.close()

    def xs_input(self, name):
        try:
            self.xs_dic = {}
            with open(name,'r') as f:
                for line in f:
                    items = line.split()
                    self.xs_dic[int(items[0])] = float(items[1])
        except IOError:
            print(name,"is not there")
            self.xs_dic = None

    def flat_list(self, limit_dic):
        mass_list = list(limit_dic.keys())
        values_list = list(limit_dic.values())
        self.xs_input('xs_input_22.txt')
        if self.xs_dic:
            l_2n = [x[0]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            l_1n = [x[1]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            l_m = [x[2]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            l_1p = [x[3]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            l_2p = [x[4]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            l_ob = [x[5]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            if len(values_list[0]) > 6:
                l_ob_e = [x[6]*self.xs_dic[mass_list[ii]] for ii,x in enumerate(values_list)]
            else:
                l_ob_e = [0]*len(l_ob)
        else:
            l_2n = [x[0] for x in values_list]
            l_1n = [x[1] for x in values_list]
            l_m = [x[2] for x in values_list]
            l_1p = [x[3] for x in values_list]
            l_2p = [x[4] for x in values_list]
            l_ob = [x[5] for x in values_list]
            if len(values_list[0]) > 6:
                l_ob_e = [x[6] for x in values_list]
            else:
                l_ob_e = [0]*len(l_ob)
        return mass_list, l_2n, l_1n, l_m, l_1p, l_2p, l_ob, l_ob_e

    def process(self):
        mass_a, dd, dd, m_a, dd, dd, obs_a, obs_a_e = self.flat_list(self.limit_list[0])
        mass_b, dd, dd, m_b, dd, dd, obs_b, obs_b_e = self.flat_list(self.limit_list[1])
        # self.plot_cmp( (mass_a, m_a), (mass_b, m_b), "cmp_median.pdf", self.tag_names)
        self.plot_cmp( (mass_a, obs_a, obs_a_e), (mass_b, obs_b, obs_b_e), "cmp_obs.pdf", self.tag_names)

if __name__ == "__main__":
    file_names = ["limit_asym.txt", "limit_toys.txt"]
    tag_names = ['asymptotic', 'toys']
    cmp_limit = CompareLimits(file_names, tag_names)
    cmp_limit.process()
