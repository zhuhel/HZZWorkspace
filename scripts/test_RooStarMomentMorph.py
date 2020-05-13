#!/usr/bin/env python
import ROOT
from ROOT import RooRealVar
from ROOT import RooArgList
from ROOT import std
from ROOT import RooGaussian

ROOT.gROOT.SetBatch()

def createMorph(obs, nuis_var, norm, up, down):
    pdf_list = RooArgList()
    pdf_list.add(down)
    pdf_list.add(up)
    pdf_list.add(norm) # nominal pdf

    nref_points = std.vector('double')()
    nref_points.push_back(-1)
    nref_points.push_back(1)

    nnuis_points = std.vector('int')()
    nnuis_points.push_back(2)

    par_list = RooArgList()
    par_list.add(nuis_var)

    morph = ROOT.RooStarMomentMorph("morph", "morph",
                                    par_list, obs, pdf_list,
                                    nnuis_points, nref_points,
                                    ROOT.RooStarMomentMorph.Linear)
    return morph

def test_1d():

    m1 = RooRealVar("m1", "m1", 125.0, 110, 140)
    s1 = RooRealVar("s1", "s1", 2.5, 0, 10)

    m2 = RooRealVar("m2", "m2", 124.0, 110, 140)
    s2 = RooRealVar("s2", "s2", 2.0, 0, 10)
    m3 = RooRealVar("m3", "m3", 126.0, 110, 140)
    s3 = RooRealVar("s3", "s3", 3.0, 0, 10)

    m4l = RooRealVar("m4l", "m4l", 110, 140)
    g1 = RooGaussian("g1", "g1", m4l, m1, s1)
    g2 = RooGaussian("g2", "g2", m4l, m2, s2)
    g3 = RooGaussian("g3", "g3", m4l, m3, s3)

    nuis_var = RooRealVar('nuis_var', "nuis_var", 0., -5., 5.)
    morph = createMorph(m4l, nuis_var, g1, g3, g2)

    canvas = ROOT.TCanvas("canvas", "canvas", 450, 450)
    frame = m4l.frame()
    morph.plotOn(frame, ROOT.RooFit.LineColor(2))
    nuis_var.setVal(1.0)
    morph.plotOn(frame, ROOT.RooFit.LineColor(4))

    nuis_var.setVal(0.5)
    morph.plotOn(frame, ROOT.RooFit.LineColor(8),
                ROOT.RooFit.LineStyle(2))
    nuis_var.setVal(-0.5)
    morph.plotOn(frame, ROOT.RooFit.LineColor(8),
                ROOT.RooFit.LineStyle(2))
    nuis_var.setVal(-1.0)
    morph.plotOn(frame, ROOT.RooFit.LineColor(4))

    frame.Draw()
    canvas.SaveAs("gaussian_morph_1d.pdf")

def test_2d():
    obs_x = RooRealVar("obs_x", "obs_x", 0, 10)
    obs_y = RooRealVar("obs_y", "obs_y", 0, 10)
    obs_list = RooArgList(obs_x, obs_y) 
    
    formular = "obs_x + obs_y*0.3"
    pdf_sum = ROOT.RooGenericPdf("pdf_sum", "sum pdf", formular,
                                 obs_list)
    f_up = "obs_x*1.8 + obs_y*0.3"
    f_down = "obs_x*0.2 + obs_y*0.3"
    pdf_up = ROOT.RooGenericPdf("pdf_up", "sum pdf", f_up,
                                 obs_list)
    pdf_down = ROOT.RooGenericPdf("pdf_down", "sum pdf", f_down,
                                 obs_list)
    
    nuis_var = RooRealVar('nuis_var', "nuis_var", 0., -5., 5.)
    morph = createMorph(obs_list, nuis_var, pdf_sum, pdf_up, pdf_down)
    morph.Print()

    th2f = pdf_sum.createHistogram("th2f", obs_x, 
                                   ROOT.RooFit.YVar(obs_y))
    
    th2f_morph = morph.createHistogram("th2f_morph", obs_x, 
                                   ROOT.RooFit.YVar(obs_y))

    canvas = ROOT.TCanvas("canvas", "canvas", 450, 450)

    th2f.Draw('colz')
    canvas.SaveAs("sum_hist.pdf")
    ##project to X
    th1_x = th2f.ProjectionX("th1_x", 1)
    th1_x.Draw()
    th1_x_morph = th2f_morph.ProjectionX("th1_x_morph", 1)
    th1_x_morph.Draw("same")
    canvas.SaveAs("sum_hist_x.pdf")

if __name__ == "__main__":
    print "hello"
    test_2d()
