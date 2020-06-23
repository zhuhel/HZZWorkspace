#!/usr/bin/env python3
import os
import sys
import math

import ROOT
from optparse import OptionParser

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.locator_params(axis='y', nticks=10)

from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy.stats import chi2

import numpy as np
from numpy import vectorize


def get_mass(input_str):
    items = input_str.split(',')
    mass_str = 0
    for item in items:
        if "mH" in item:
            mass_str = item
            break

    if ":" in mass_str:
        return mass_str.split(':')[1]
    elif "=" in mass_str:
        return mass_str.split('=')[1]
    else:
        return mass_str


def get_limit_dic(file_name):
    limit_dic = {}
    with open(file_name, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            items = line.split()
            mass_ = int(get_mass(items[0]))
            start_i = 3
            one_entry = [float(x) for x in items[start_i:]]
            limit_dic[mass_] = one_entry
    return limit_dic


def exp_fun(x, a, b):
    return a * np.exp(-b * x)


def clsb_func(x, mu, sigma):
    ratio = mu/sigma
    th = ratio**2
    if x <= 0:
        res = 0.
    elif x <= th:
        res = 0.5* 1/np.sqrt(2*np.pi)/np.sqrt(x)*np.exp(-0.5*x)
    else:
        res = 1/np.sqrt(2*np.pi)/(2*ratio)*np.exp(-0.5*(x + ratio**2)**2/(2*ratio)/(2*ratio))

    return res
vclsb_func = vectorize(clsb_func)


def clb_func(x, mu, sigma):
    ratio = mu/sigma
    th = ratio**2
    if x <= 0:
        res = 0.
    elif abs(x-0.0001)<1E-7:
        res = norm.cdf(-mu/sigma)
    elif x <= th:
        res = 0.5* 1/np.sqrt(2*np.pi)/np.sqrt(x)*np.exp(-0.5*(np.sqrt(x) - ratio)**2)
    else:
        res = 1/np.sqrt(2*np.pi)/(2*ratio)*np.exp(-0.5*(x - ratio**2)**2/(2*ratio)**2)

    return res
vclb_func = vectorize(clb_func)


def get_sigma(qmu, mu_hat, mu):
    num = 0
    if mu_hat < 0:
        num = (mu**2 - 2*mu*mu_hat)
    elif mu_hat <= mu:
        num = (mu - mu_hat)**2
    else:
        print("best fitted should not greater than mu!")
        return 1
    return math.sqrt(num/qmu)


def get_Z(qmu, sigma, mu, mu_prime):
    th = mu**2/sigma**2
    if qmu <= th:
        sig = math.sqrt(qmu) - abs((mu - mu_prime)/sigma)
    else:
        sig = (qmu - (mu**2 - 2*mu*mu_prime)/sigma**2)/(2*abs(mu/sigma))
    return sig


def get_pmu(qmu, sigma, mu):
    return 1 - norm.cdf(get_Z(qmu, sigma, mu, mu))


def get_pb(qmu, sigma, mu):
    return norm.cdf(get_Z(qmu, sigma, mu, 0))


class ToyReader:
    def __init__(self, mass, toy_dir, data_input, exp_input, poi_name, poi_list):
        self.mass = mass
        self.toy_dir = toy_dir

        self.data_input = os.path.join(toy_dir, data_input)
        self.exp_input = os.path.join(toy_dir, exp_input)
        self.poi_name = poi_name

        self.poi_list = ["{:.4f}".format(x) for x in sorted(poi_list)]
        self.pattern = "toys_mH{}_{}".format(mass, self.poi_name)

        self.obs_nll_dic = {}
        self.read_data(self.data_input)
        self.read_exp(self.exp_input)

        self.do_exp = False

    def toy_output_pattern(self, _str):
        self.pattern = _str

    @staticmethod
    def get_p0(nll_list, nll):
        total = len(nll_list)
        larger = len([x for x in nll_list if x > nll])
        f = larger*1.0/total
        e = math.sqrt(larger*(1-f))/total # binominal error
        return f, e

    def get_cls(self, nll_sb, nll_bonly, nll):
        clsb, clsb_e = [x/2.0 for x in self.get_p0(nll_sb, nll)]
        clb, clb_e = self.get_p0(nll_bonly,nll)

        if clb == 0:
            print("Clb is zero!")
            cls = clsb
            cls_e = clsb_e
        else:
            cls = clsb/clb
            if clsb == 0:
                cls_e = cls
            else:
                cls_e = cls*math.sqrt((clsb_e/clsb)**2 + (clb_e/clb)**2)

        return clsb, clb, cls, cls_e

    def read_data(self, data_input):
        """
        reading the NLL of observed data
        """
        if not os.path.exists(data_input):
            print(data_input,"is not there.")
            return

        f_data = ROOT.TFile.Open(data_input)
        if not f_data or f_data.IsZombie():
            return

        tree = f_data.Get("physics")
        entries = tree.GetEntries()
        if entries != len(self.poi_list):
            print("entries in observed does not match the one in poi list")
            print( entries, len(self.poi_list))
            print("--------------------------")
            if entries < len(self.poi_list):
                return

        for ientry in range(entries):
            tree.GetEntry(ientry)
            if tree.nll_SB_free > tree.nll_SB_fixed:
                nll = 0
            else:
                nll = 2 * (tree.nll_SB_fixed - tree.nll_SB_free)
                self.obs_nll_dic["{:.4f}".format(tree.mu)] = (nll, tree.poi_SB_free, tree.qzero, tree.bonly_muhat)

        f_data.Close()

    def read_exp(self, exp_input):
        """ read expected toys
        """
        if not os.path.exists(exp_input):
            print( exp_input,"is not there." )
            return

        f_exp = ROOT.TFile.Open(exp_input)
        if not f_exp or f_exp.IsZombie():
            return

        tree = f_exp.Get("physics")
        entries = tree.GetEntries()
        self.exp_nll_dic = {}

        for ientry in range(entries):
            tree.GetEntry(ientry)
            if tree.nll_SB_free > tree.nll_SB_fixed:
                nll = 0
            else:
                nll = 2 * (tree.nll_SB_fixed - tree.nll_SB_free)
                if nll > 100:
                    continue
            key = "{:.4f}".format(tree.mu)

            if key not in self.exp_nll_dic:
                self.exp_nll_dic[key] = []

            self.exp_nll_dic[key].append(nll)

        print("-----------------------")
        print("summary of expected toys")
        for x,y in self.exp_nll_dic.items():
            print(x,len(y))
            print("-----------------------")
            f_exp.Close()

    def read_toys(self, poi):
        chain = ROOT.TChain("physics", "physics")
        chain.Add(os.path.join(self.toy_dir, self.pattern+poi+"_seed*root"))
        nentries = chain.GetEntries()
        nll_sb = []
        for ientry in range(nentries):
            chain.GetEntry(ientry)
            if ROOT.TMath.IsNaN(chain.nll_SB_fixed) or\
               ROOT.TMath.IsNaN(chain.nll_SB_free):
                continue

            if chain.nll_SB_free > chain.nll_SB_fixed:
                nll = 0
            else:
                nll = 2*(chain.nll_SB_fixed - chain.nll_SB_free)

            nll_sb.append(nll)

        print("total entries for S+B toys: ",len(nll_sb))
        return nll_sb

    def get_limit(self, cls_list, out_name, cls_e_list=None):
        kind = 'quadratic' # quadratic, linear
        if cls_e_list is not None:
            for x,y,z in zip(self.poi_vals, cls_list, cls_e_list):
                print("{:.4f} {:.4f} {:.4f}".format(x, y, z))

        steps = (self.poi_vals[-1] - self.poi_vals[0])/1000.
        new_poi = np.arange(self.poi_vals[0]-0.1, self.poi_vals[-1]+0.1, steps)
        use_extrapolate = False
        popt = None
        perr = None
        if use_extrapolate:
            f = interpolate.interp1d(self.poi_vals, cls_list, kind, bounds_error=False)
            new_cls = f(new_poi)
        else:
            try:
                if cls_e_list is None:
                    popt, pcov = curve_fit(exp_fun, np.array(self.poi_vals), np.array(cls_list), p0=(1.67, 0.54))
                else:
                    popt, pcov = curve_fit(exp_fun, np.array(self.poi_vals), np.array(cls_list),
                                     sigma=np.array(cls_e_list), absolute_sigma=True, p0=[1.67, 0.54])
                # print(popt)
                # print("cov:", pcov)
                perr = np.sqrt(np.diag(pcov))
                # print (perr)
                new_cls = exp_fun(new_poi, *popt)
            except RuntimeError as e:
                print("cannot perform the fit")
                return 0., 0

        # find the poi that yields the cls that is close to 0.05
        limit = 0.
        min_dis = 999
        for x, y in zip(new_poi, new_cls):
            if abs(y - 0.05) < min_dis:
                min_dis = abs(y - 0.05)
                limit = x

        limit = self.find_x(new_poi, new_cls, 0.05)
        limit_e = 0. # error of the limit
        if use_extrapolate:
            plt.plot(self.poi_vals, cls_list, 'o', new_poi, new_cls, '-')
        else:
            # plt.plot(self.poi_vals, cls_list, 'o', new_poi, new_cls, '-')
            plt.semilogy(self.poi_vals, cls_list, 'o', new_poi, new_cls, '-')
            up_cls = exp_fun(new_poi, *(popt + perr))
            down_cls = exp_fun(new_poi, *(popt - perr))
            limit_up = self.find_x(new_poi, up_cls, 0.05)
            limit_down = self.find_x(new_poi, down_cls, 0.05)
            limit_e = abs(limit_up - limit_down)/2.

        plt.xlabel('mu')
        plt.ylabel('CLs')
        plt.plot(plt.xlim(), [0.05, 0.05], 'r-', lw=2)
        plt.plot([limit, limit], plt.ylim(), 'r-', lw=2)

        plt.annotate('({:.4f},{:.2f})'.format(limit, 0.05),
                     xy=(limit, 0.05),
                     xytext=(limit*1.1, 0.1),
                     arrowprops=dict(facecolor='black', shrink=0.05),
                    )
        if use_extrapolate:
            plt.savefig('mH_{}_{}_expo{}.pdf'.format(self.mass, out_name, kind))
        else:
            plt.savefig('mH_{}_{}_expFun.pdf'.format(self.mass, out_name))

        plt.close()
        print("Limit: {:.4f} {:.4f}", limit, limit_e)
        return limit, min_dis, limit_e

    @staticmethod
    def find_x(x_list, y_list, th):
        res = 0;
        min_dis = 9999.
        for x, y in zip(x_list, y_list):
            if abs(y - th) < min_dis:
                min_dis = abs(y - th)
                res = x
        return res

    def get_expected_nll(self, poi):
        try:
            nll_list = self.exp_nll_dic[poi]
        except (AttributeError, KeyError) as e :
            return None

        th_2n = ROOT.Math.gaussian_cdf(-2)
        th_1n = ROOT.Math.gaussian_cdf(-1)
        th_median = ROOT.Math.gaussian_cdf(0)
        th_1p= ROOT.Math.gaussian_cdf(1)
        th_2p= ROOT.Math.gaussian_cdf(2)

        bins = 80
        hist, X1 = np.histogram(nll_list, bins=bins, density=True)
        dx = X1[1] - X1[0]
        cdf_list = np.cumsum(hist)*dx
        kind = "cubic" # 'quadratic'
        f = interpolate.interp1d(X1[:-1], cdf_list, kind, bounds_error=False)

        #new_nll = np.arange(X1[0], X1[-1], 1)
        new_nll = np.arange(0, X1[-1], 0.01)
        new_cdf = f(new_nll)

        # find the NLL gives +/-1 1/2 sigma
        val_2n = self.find_x(new_nll, new_cdf, th_2n)
        val_1n = self.find_x(new_nll, new_cdf, th_1n)
        val_m = self.find_x(new_nll, new_cdf, th_median)
        val_1p = self.find_x(new_nll, new_cdf, th_1p)
        val_2p = self.find_x(new_nll, new_cdf, th_2p)
        print(val_1p,val_2p,val_m,val_1n,val_2n)

        # make plots for debugging
        fig, ax1 = plt.subplots()
        ax1.hist(nll_list, bins=bins)
        l_2n, = ax1.plot([val_2n, val_2n], plt.ylim(), 'y-', lw=2)
        l_1n, = ax1.plot([val_1n, val_1n], plt.ylim(), 'y-', lw=2)
        l_m,  = ax1.plot([val_m, val_m], plt.ylim(), 'g-', lw=2)
        l_1p, = ax1.plot([val_1p, val_1p], plt.ylim(), 'r-', lw=2)
        l_2p, = ax1.plot([val_2p, val_2p], plt.ylim(), 'r-', lw=2)

        ax2 = ax1.twinx()
        scaled_new_cdf = new_cdf * len(nll_list)
        l_cdf, = ax2.plot(new_nll, scaled_new_cdf, 'm--', lw=2)
        ax2.set_ylim(0, len(nll_list)*1.1)

        leg = ax1.legend([l_2n, l_1n, l_m, l_1p, l_2p, l_cdf], 
                         ["-2$\sigma$", '-1$\sigma$', "median", '+1$\sigma$', '+$2\sigma$',"CDF"],
                        framealpha=0.0)
        ax1.set_xlabel('-2 $\ln\Lambda$')
        ax1.set_ylabel('Events')
        ax1.set_title('background only toys')
        fig.tight_layout()
        #plt.xlabel('-2 $\ln\Lambda$')
        #plt.ylabel('Events')
        #plt.title('background only toys')
        #plt.tight_layout()

        #plt.savefig('mH_{}_nll_bonlyToys_{}.pdf'.format(self.mass, poi))
        plt.close()
        return (val_2n, val_1n, val_m, val_1p, val_2p)

    def process(self):
        cls_obs = []
        cls_obs_e = []
        cls_2n = []
        cls_1n = []
        cls_m = []
        cls_1p = []
        cls_2p = []
        self.poi_vals = []
        for poi in self.poi_list:
            mu_val = float(poi)
            nll_sb = self.read_toys(poi)
            nll_bonly = self.exp_nll_dic[poi]
            if not hasattr(self, "obs_nll_dic"):
                print("a DUMMY observed NLL used!!")
                obs_nll = 1.
                muhat = mu_val
                qzero = bonly_muhat = 0
            else:
                obs_nll = self.obs_nll_dic[poi][0]
                muhat = self.obs_nll_dic[poi][1]
                qzero = self.obs_nll_dic[poi][2]
                bonly_muhat = self.obs_nll_dic[poi][3]

            clsb, clb, cls, cls_e = self.get_cls(nll_sb, nll_bonly, obs_nll)
            if  cls < 0.0001:
                continue
            self.poi_vals.append(mu_val)
            cls_obs.append(cls)
            cls_obs_e.append(cls_e)

            nbins = 80
            is_log = True
            bin_values, bin_bound, l_sb = plt.hist(nll_sb, bins=nbins, range=(-5, 35), log=is_log)
            dd, bin_bonly, l_bonly = plt.hist(nll_bonly, bins=nbins, histtype='step',
                                       range=(-5, 35), log=is_log, color='m')
            if is_log:
                plt.ylim(0.1, 1E4)
            else:
                plt.ylim(0, 5000)

            l_obs_nll, = plt.plot([obs_nll, obs_nll], plt.ylim(), 'r-', lw=2)

            # put chi-squared distribution on top of p_mu.
            center_list = [0.5*(bin_bound[i+1] + bin_bound[i]) for i in range(len(bin_bound)-1)]
            obs_sigma = get_sigma(obs_nll, muhat, mu_val)
            chi_val = vclsb_func(np.array(center_list), mu_val, obs_sigma)
            chi_val *= len(nll_sb)
            l_chi2, = plt.plot(center_list, chi_val, 'y-', lw=2)

            # put clb function on top of p_b.
            sigma_bkg = get_sigma(qzero, bonly_muhat, mu_val)
            bonly_curve = vclb_func(np.array(center_list), mu_val, sigma_bkg)
            bonly_curve *= 0.5*len(nll_bonly)
            l_nochi2, = plt.plot(center_list, bonly_curve, 'g-', lw=2)

            # add legend..
            leg = plt.legend([l_sb[0], l_bonly[0], l_obs_nll, l_chi2, l_nochi2],
                             ["S+B", "B-only", "observed", "$p_\mu$ Asymp.", "$p_b$ Asymp."],
                             framealpha=0.)

            # add asymptotic info
            print("observed fitted info:", obs_nll, obs_sigma, sigma_bkg, mu_val)
            asym_clsb = get_pmu(obs_nll, sigma_bkg, mu_val)
            asym_clb = 1 - get_pb(obs_nll, sigma_bkg, mu_val)
            asym_cls = asym_clsb/asym_clb
            height = 0.95
            x_loc = 0.4
            ax = plt.gca()
            asymp_summary = "Asymptotic Summary\n$p_\mu$\t = {:.4f}\n$p_b$\t = {:.4f}\nCL$_s$\t = {:.4f}".format(
                        asym_clsb, asym_clb, asym_cls)
            plt.text(x_loc, height, asymp_summary,
                     verticalalignment='top', horizontalalignment='left',
                     transform=ax.transAxes,
                     bbox=dict(facecolor='none', alpha=0.5),
                     multialignment='center'
                    )
            #plt.text(x_loc, height-, r'CLb:{:.4f}'.format(asym_clb))
            #plt.text(x_loc, height, r'CLs:{:.4f}'.format(asym_cls))

            #plt.xlabel("-2 $\ln\Lambda$")
            plt.xlabel(r"""$\tilde{q}_{\mu}$""")
            plt.ylabel("Events")
            plt.title('$m_H$ = {:.0f}, $\mu$ = {:.4f}, $p_\mu$ = {:.4f}, $p_b$ = {:.4f}, CL$_s$ = {:.4f}'.format(self.mass, float(poi), clsb, clb, cls))

            plt.savefig('mH_{}_nll_{}.pdf'.format(self.mass, poi))
            plt.close()

            # find the median and +/- 1/2 sigma NLL from b-only fit
            try:
                nll_2n,nll_1n,nll_m,nll_1p,nll_2p = self.get_expected_nll(poi)
            except (TypeError, ValueError):
                continue
            #print "\tpoi:{}, observed_NLL:{:.4f}, expected:{:.4f}".format(poi, obs_nll, nll_m)
            if self.do_exp:
                cls_2n.append(self.get_cls(nll_sb, nll_bonly, nll_2n)[2])
                cls_1n.append(self.get_cls(nll_sb, nll_bonly, nll_1n)[2])
                cls_m.append(self.get_cls(nll_sb, nll_bonly, nll_m)[2])
                cls_1p.append(self.get_cls(nll_sb, nll_bonly, nll_1p)[2])
                cls_2p.append(self.get_cls(nll_sb, nll_bonly, nll_2p)[2])
            else:
                print("not doing expected limits")

        obs, dis, obs_e = self.get_limit(cls_obs, "obs", cls_obs_e)
        if self.do_exp and len(cls_2n) > 0:
            exp_m, dis, exp_m_e = self.get_limit(cls_m, "exp_median")
            exp_2n, dis, exp_2n_e = self.get_limit(cls_2n, "exp_2n")
            exp_1n, dis, exp_1n_e = self.get_limit(cls_1n, "exp_1n")
            exp_1p, dis, exp_1p_e = self.get_limit(cls_1p, "exp_1p")
            exp_2p, dis, exp_2p_e = self.get_limit(cls_2p, "exp_2p")
        else:
            exp_2n = exp_1n = exp_m = exp_1p = exp_2p = 0.
        return exp_2p, exp_1p, exp_m, exp_1n, exp_2n, obs, obs_e

if __name__ == "__main__":
    usage = "%prog [options]"
    version = "%prog 1.0"
    parser = OptionParser(usage=usage, description="process toy data", version=version)
    parser.add_option('-p', '--poi_name', default='XS_ggF')
    parser.add_option('-t', '--test', default=None)
    parser.add_option('-e', '--exp', default=False, action='store_true')

    options, args = parser.parse_args()
    poi_name = options.poi_name
    print("poi name", poi_name)
    mass_dic = {}
    if options.test is None:
        org_mass_dic = get_limit_dic("limit_asym.txt")
        for x,y in org_mass_dic.items():
            new_y = list(map(lambda x: round(x, 4), y))
            #new_y = new_y[:-1]
            mass_dic[x] = new_y
    else:
        items = options.test.split(':')
        mass_val = int(items[0])
        poi_list = [float(x) for x in items[1].split(',')]
        mass_dic = {mass_val: poi_list}

    print(mass_dic)
    results = ""
    for mass, poi_list in mass_dic.items():
        reader = ToyReader(mass, "data",
                           'observed_mH{}.root'.format(mass),
                           'expected_mH{}.root'.format(mass),
                           poi_name,
                           poi_list)
        reader.do_exp = options.exp
        limits = reader.process()
        results += "mH="+str(mass)+" -1 -1 "+" ".join(map(lambda x: "{:.4f}".format(x), limits)) + "\n"
        del reader
    print(results)
