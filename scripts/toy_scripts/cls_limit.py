#!/usr/bin/env python

import math
from scipy.stats import norm
from scipy.optimize import curve_fit

import ROOT
from ROOT import RooFit

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.locator_params(axis='y', nticks=10)

def exp_fun(x, a, b):
    return a * np.exp(-b * x)

class AsympLimit:
    def __init__(self):
        print "cls limit"
        self.do_tilde = True

        # (mu_hat, NLL)
        self.best_fit = (0, 0)

        # for each mu, calculate qmu
        self.mu_dic = {}

    def read_workspace(self, file_name):
        self.fin = ROOT.TFile.Open(file_name)
        self.ws = self.fin.Get("combined")
        self.mc = self.ws.obj("ModelConfig")
        self.data = self.ws.obj("obsData")
        self.poi = self.ws.obj("mu_ggF")
        self.simPdf = self.mc.GetPdf()

    def test(self, mu):
        if mu not in self.mu_dic:
            #qmu, sigma, pmu, pb, cls
            self.mu = mu
            self.qmu = self.cal_qmu()
            self.sigma = self.get_sigma()
            self.get_cls_IN()
            self.mu_dic[mu] = (self.qmu, self.sigma, self.pmu, self.pb, self.cls)

        return self.mu_dic[mu]

    def cal_qmu(self):
        """perform Fit and get NLL"""
        # signal + background fit
        self.poi.setConstant(0)
        nll_free = self.fit()
        self.best_fit = (self.poi.getVal(), nll_free)

        if self.poi.getVal() < 0:
            self.poi.setVal(0)
            self.poi.setConstant(1)
            nll_free = self.fit()

        self.poi.setVal(self.mu)
        self.poi.setConstant(1)
        nll_fixed = self.fit()
        qmu = 2*(nll_fixed - nll_free)
        return qmu


    def fit(self):
        nuisance = self.mc.GetNuisanceParameters()
        nll = self.simPdf.createNLL(
            self.data,
            RooFit.Constrain(nuisance),
            RooFit.GlobalObservables(self.mc.GetGlobalObservables())
        )
        nll.enableOffsetting(True)
        minim = ROOT.RooMinimizer(nll)
        minim.optimizeConst(2)
        minim.setStrategy(1)
        status = minim.minimize(
            "Minuit2",
            ROOT.Math.MinimizerOptions.DefaultMinimizerAlgo()
        )
        if status != 0:
            print "status is not zero!"

        return nll.getVal()

    def get_sigma(self):
        num = 0.
        mu_hat = self.best_fit[0]
        mu = self.mu
        if mu_hat < 0 and self.do_tilde:
            num = (mu**2 - 2*mu*mu_hat)
        elif mu_hat < mu:
            num = (mu - mu_hat)**2
        else:
            print "best fitted should not greater than mu!"
            return 1

        print num, type(num), self.qmu, type(self.qmu)
        return math.sqrt(num/self.qmu)

    def get_Z(self, qmu, sigma, mu, mu_prime):
        th = mu**2/sigma**2
        if qmu  < th or (not self.do_tilde):
            sig = math.sqrt(qmu) - abs((mu - mu_prime)/sigma)
        else:
            sig = (qmu - (mu**2 - 2*mu*mu_prime)/sigma**2)/(2*abs(mu/sigma))
        return sig

    def get_pmu(self):
        return 1 - norm.cdf(self.get_Z(
            self.qmu, self.sigma, self.mu, self.mu
        ))

    def get_pb(self):
        return norm.cdf(self.get_Z(
            self.qmu, self.sigma, self.mu, 0
        ))

    def get_cls(self, mu_hat, qmu, mu):
        self.best_fit = (mu_hat, 0)
        self.mu = mu
        self.qmu = qmu
        self.sigma = self.get_sigma()
        self.get_cls_IN()
        self.print_output()

    def print_output(self):
        return "mu: {:.4f}\nqmu: {:.4f}\nsigma: {:.4f}\
        \nCLs+b: {:.4f}\nClb: {:.4f}\nCls: {:.4f}\n".format(
            self.mu, self.qmu, self.sigma, self.pmu, self.pb, self.cls
        )

    def get_cls_IN(self):
        self.pmu = self.get_pmu()
        self.pb = self.get_pb()
        self.cls = self.pmu/(1-self.pb)

    def process(self):
        self.read_workspace("HZZllvv_Graviton800_combined_HZZllvv_Graviton800_model_forLimits.root")
        mu_vals = [0.1031, 0.1921, 0.1383]
        cls_vals = []
        out_text = ""
        for mu in mu_vals:
            cls_vals.append(self.test(mu)[4])
            out_text += self.print_output()

        print out_text
        steps = (mu_vals[-1] - mu_vals[0])/1000.
        new_poi = np.arange(mu_vals[0]-0.1, mu_vals[-1]+0.1, steps)
        popt = curve_fit(exp_fun, np.array(mu_vals), np.array(cls_vals), p0=(1.67, 0.54))[0]
        new_cls = exp_fun(new_poi, *popt)

        # find the poi that yields the cls that is close to 0.05
        limit = 0.
        min_dis = 999
        for x, y in zip(new_poi, new_cls):
            if abs(y - 0.05) < min_dis:
                min_dis = abs(y - 0.05)
                limit = x

        print "limits:", limit
        plt.plot(mu_vals, cls_vals, 'o', new_poi, new_cls, '-')
        plt.xlabel('mu')
        plt.ylabel('CLs')
        plt.plot(plt.xlim(), [0.05, 0.05], 'r-', lw=2)
        plt.plot([limit, limit], plt.ylim(), 'r-', lw=2)
        plt.savefig('TEST.pdf')

        self.fin.Close()

if __name__ == "__main__":
    limit = AsympLimit()
    #limit.get_cls(5.4739, 3.4366, 38.6059)
    #limit.get_cls(4.22642e-06, 4.78328, 0.1383)

    limit.process()
