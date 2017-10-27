#!/usr/bin/env python
"""
obtain global pvalue according to
- number of upper crosssing, w.r.t 0.7sigma
- local pvalue
ref: http://arxiv.org/abs/1005.1891
"""
import ROOT
import math
import sys

class AnalyzeGlobal(object):
    def __init__(self, file_name):
        import AtlasStyle
        print "Input file:", file_name
        self.f1 = ROOT.TFile.Open(file_name)
        self.tree = self.f1.Get("physics")

    def process(self):
        self.draw_sig("sigma_s1", "most significant excess", 0.70)
        self.draw_sig("sigma_s2", "second-to-most significant excess", 0.65)

    def draw_sig(self, br_name, title, x_off):
        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        ROOT.gStyle.SetOptStat(1111)
        h_sigma_s1 = ROOT.TH1F("h_"+br_name, "s1", 100, 1, 5)
        self.tree.Draw(br_name+">>h_"+br_name, "")

        h_sigma_s1.GetXaxis().SetTitle("local significance")
        ROOT.gPad.Update()
        h_sigma_s1.SetName(title)
        st = h_sigma_s1.FindObject("stats")
        st.SetX1NDC(x_off)
        st.SetY1NDC(x_off)
        st.SetX2NDC(0.9)
        st.SetY2NDC(0.9)
        st.SetShadowColor(0)
        canvas.SaveAs(br_name+".eps")
        canvas.SaveAs(br_name+".pdf")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print sys.argv[0]," file_name"
        sys.exit(1)

    file_ = sys.argv[1]
    ag = AnalyzeGlobal(file_)
    ag.process()
