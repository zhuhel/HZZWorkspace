#!/usr/bin/env python
import ROOT
import os
from optparse import OptionParser



def plot(flist,
         poi_name="mu", x_title="#mu",
         tree_name="physics", nll_name="NLL",
        ):
    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    ROOT.SetAtlasStyleCanvas(canvas, True)
    mg = ROOT.TMultiGraph()
    color = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kOrange]
    ic = 0

    legend = ROOT.AtlasUtils.myLegend(0.6, 0.7, 0.8, 0.9)

    ## get each graph one by one
    try:inputlist=open(flist, "r")
    except IOError:
        print "Error: file could not be opened! ", inputlist
        sys.exit(1)
    for f1_name in inputlist.readlines():
        if f1_name[-1]=="\n": f1_name=f1_name[:-1]
        minMu = ROOT.Double(-999)
        mulow = ROOT.Double(-999)
        muHi =  ROOT.Double(-999)
        file_name, l1_name = f1_name.split(':')

        graph = ROOT.getGraphFromFile(
            file_name,
            tree_name, nll_name, poi_name,
            minMu, mulow, muHi)
        graph.SetLineWidth(2)
        graph.SetLineColor(color[ic])
        ic += 1

        mg.Add(graph, "L")
        legend.AddEntry(graph, l1_name, "l")


    mg.Draw("A")
    mg.GetXaxis().SetTitle(x_title)
    mg.GetYaxis().SetTitle("-2ln#Lambda")
    xmin = 0.9
    xmax = 1.1
    mg.GetXaxis().SetRangeUser(xmin,xmax)
    mg.GetYaxis().SetRangeUser(0,10)
    #xmin = mg.GetXaxis().GetXmin()
    #xmax = mg.GetXaxis().GetXmax()

    onesigma = ROOT.TLine()
    onesigma.SetLineStyle(2)
    onesigma.SetLineWidth(2)
    y_var = 1.0
    x_off_sigma = xmax*1.05
    y_off_sigma = 0.02
    onesigma.DrawLine(xmin, y_var, xmax, y_var)
    ROOT.AtlasUtils.myText(x_off_sigma, y_var-y_off_sigma, 1, "1#sigma", False)
    y_var = 4.0
    onesigma.DrawLine(xmin, y_var, xmax, y_var)
    ROOT.AtlasUtils.myText(x_off_sigma, y_var-y_off_sigma, 1, "2#sigma", False)
    #y_var = 9.0
    #onesigma.DrawLine(xmin, y_var, xmax, y_var)
    #ROOT.myText(x_off_sigma, y_var-y_off_sigma, 1, "3#sigma", False)

    x_off_title = 0.185
    lumi_weight = 139.0
    ROOT.AtlasUtils.myText(x_off_title, 0.85, 1, "#bf{#it{ATLAS}} Internal")
    ROOT.AtlasUtils.myText(x_off_title, 0.80, 1, "mH = 300GeV")
    ROOT.AtlasUtils.myText(x_off_title, 0.75, 1, "13 TeV, {:.1f} fb^{{-1}}".format(lumi_weight))
    legend.Draw()

    canvas.SaveAs(flist.replace("txt",'pdf'))
    #canvas.SaveAs('test.pdf')

if __name__ == "__main__":
    # plot("test_mass_new.root:observed m_{H}", "mH", "m_{H} [GeV]")
    usage = "%prog filelist poi_name x_title"
    version="%prog 1.0"
    parser = OptionParser(usage=usage, description="plot 1D NLL scan", version=version)
    (options,args) = parser.parse_args()

    if len(args) < 3:
        print parser.print_help()
        exit(1)
    f1_name = args[0]
    poi_name = args[1]
    x_title = args[2]

    if not hasattr(ROOT, "getGraphFromFile"):
        ROOT.gROOT.LoadMacro("draw1DNLL.cxx")

    if not hasattr(ROOT, "loader"):
        ROOT.gROOT.LoadMacro("loader.c")

    from HZZWorkspace import AtlasStyle
    ROOT.gROOT.SetBatch()

    # call real program
    plot(f1_name, poi_name, x_title)
