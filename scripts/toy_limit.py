#!/usr/bin/env python
import os

import ROOT

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.locator_params(axis='y', nticks=10)

from scipy import interpolate
import numpy as np

class ToyReader:
    def __init__(self, mass, toy_dir, data_input, poi_list):
        self.toy_dir = toy_dir

        self.data_input = data_input
        self.poi_name = "mu_ggF"
        self.poi_vals = poi_list
        self.poi_list = ["{:.2f}".format(x) for x in poi_list]
        self.pattern = "toys_mH{}_{}".format(mass, self.poi_name)

        self.mass = mass

        self.read_data()

    def toy_output_pattern(self, _str):
        self.pattern = _str

    def get_p0(self, nll_list, nll):
        total = len(nll_list)
        larger = len([x for x in nll_list if x > nll])
        return larger*1.0/total

    def get_cls(self, nll_sb, nll_bonly, nll):
        clsb = self.get_p0(nll_sb, nll)
        clb = self.get_p0(nll_bonly,nll)
        return clsb/clb

    def read_data(self):
        f_data = ROOT.TFile.Open(self.data_input)
        tree = f_data.Get("physics")
        self.obs_nll_dic = {}
        entries = tree.GetEntries()
        if entries != len(self.poi_list):
            print "entries in observed does not match the one in poi list"
            print entries, len(self.poi_list)
            print "--------------------------"
            if entries < len(self.poi_list):
                return

        for ientry in range(entries):
            tree.GetEntry(ientry)
            nll = 2 * (tree.nll_SB_fixed - tree.nll_SB_free)
            self.obs_nll_dic["{:.2f}".format(tree.mu)] = nll

        f_data.Close()
        self.obs_nll_dic

    def read_toys(self, poi):
        chain = ROOT.TChain("physics", "physics")
        chain.Add(os.path.join(self.toy_dir, self.pattern+poi+"_seed*root"))
        nentries = chain.GetEntries()
        print "total entries for: ",poi,nentries
        nll_sb = []
        nll_bonly = []
        for ientry in range(nentries):
            chain.GetEntry(ientry)
            nll_sb.append(2*(chain.nll_SB_fixed - chain.nll_SB_free))
            nll_bonly.append(2*(chain.nll_B_fixed - chain.nll_SB_free))

        return nll_sb, nll_bonly

    def process(self):
        cls_list = []
        for poi in self.poi_list:
            nll_sb, nll_bonly = self.read_toys(poi)
            obs_nll = self.obs_nll_dic[poi]
            cls = self.get_cls(nll_sb, nll_bonly, obs_nll)
            cls_list.append(cls)

            #bins, edges = np.histogram(nll_sb, 20)
            #plt.plot(edges, bins, 'o')
            nbins = 30
            plt.hist(nll_sb, bins=nbins, range=(-5, 35))
            plt.hist(nll_bonly, bins=nbins, histtype='step', range=(-5, 35))
            plt.plot([obs_nll, obs_nll], plt.ylim(), 'r-', lw=2)
            plt.savefig('mH_{}_nll_{}.png'.format(self.mass, poi))
            plt.close()

        for x,y in zip(self.poi_vals, cls_list):
            print x,y

        f = interpolate.interp1d(self.poi_vals, cls_list, 'quadratic', bounds_error=False)
        new_poi = np.arange(self.poi_vals[0]-0.1, self.poi_vals[-1]+0.1, 0.005)
        new_cls = f(new_poi)
        # find the poi that yields the cls that is close to 0.05
        limit = 0.
        min_dis = 999
        for x, y in zip(new_poi, new_cls):
            if abs(y - 0.05) < min_dis:
                min_dis = abs(y - 0.05)
                limit = x
        print "limits:", limit
        plt.plot(self.poi_vals, cls_list, 'o', new_poi, new_cls, '-')
        plt.xlabel('mu')
        plt.ylabel('CLs')
        plt.savefig('mH_{}_cls_test.png'.format(self.mass))
        plt.close()

if __name__ == "__main__":
    mass_dic = {
        # 1600: [6.5, 7.5, 8.5],
        # 1800: [15.0, 17.0, 19.0],
        2000: [38, 42, 46]
    }
    for mass, poi_list in mass_dic.iteritems():
        reader = ToyReader(mass, "data", 'observed_q_{}.root'.format(mass), poi_list)
        reader.process()
