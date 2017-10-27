#!/usr/bin/env python
__author__ = "Xiangyang Ju"
__version__ = "0.1"
import ROOT
from array import array
import AtlasStyle

import os
import errno

from adder import adder

class Ploter:
    def __init__(self, status="Internal", lumi=36.1):
        self.status = status
        self.lumi = lumi

        # few options
        self.add_ratio = True

        self.has_data = False

        # predefined colors
        self.COLORS = [206, 64, 95, 28, 29, 209, 5, 432, 433, 434, 435, 436, 8, 6]
        self.LINE_STYLE = [1]*10
        self.VerticalCanvasSplit = 0.4

        # parameters for label/title size
        self.t_size = 0.04
        self.x_title_size = 0.05
        self.text_size = 0.05

        # my canvas
        self.can = None
        self.pad1 = None
        self.pad2 = None

        # legend
        self.legend = None

        # atlas and legend offset
        self.x_offset = 0.20
        self.x_off_atlas = 0.65
        self.y_offset = 0.85

        # show sum of background for x-check
        self.show_sum_bkg = True

        self.totalObj = []


    def prepare_2pad_canvas(self, cname, width=600, height=600):

        self.can = ROOT.TCanvas(cname,cname,width,height)
        self.pad1 = ROOT.TPad("p1_"+cname,cname,0.0, self.VerticalCanvasSplit,1.0,1.0)
        self.pad2 = ROOT.TPad("p2_"+cname,cname,0.0,0.0, 1.0, self.VerticalCanvasSplit)
        self.pad1.SetBottomMargin(0)
        self.pad1.SetTopMargin(0.09)
        self.pad1.SetLeftMargin(0.17)
        self.pad2.Draw()
        self.pad2.SetTopMargin(0)
        self.pad2.SetBottomMargin(0.4)
        self.pad2.SetLeftMargin(0.17)
        self.pad2.SetGridy()
        self.pad2.SetGridx()
        self.can.cd()
        self.pad1.Draw()
        self.pad2.Draw()

    def add_ratio_panel(self, hist_list, y_title,
                        y_min, y_max, reverse=False):
        """
        hist_list = [Data, MC1, MC2]
        plot Data/MC1, Data/MC2
        @para=reverse, plot MC1/Data
        """
        self.add_ratio = True
        #y_title = "Data/MC"
        #if reverse:
        #    y_title = "MC/Data"

        if len(hist_list) < 2:
            print "less than 2 histograms, kidding?"
            return None

        hist_list_cp = [x.Clone(x.GetName()+"_cloneRatio") for x in hist_list]
        self.totalObj.append(hist_list_cp)

        h_refer = hist_list_cp[0].Clone("Histreference")
        h_refer.Sumw2()
        self.totalObj.append(h_refer)
        # print "REFER:", h_refer.Integral()
        for i, hist in enumerate(hist_list_cp):
            if i==0:

                hist.Sumw2()
                hist.Divide(h_refer)
                hist.SetFillColor(1)
                hist.SetFillStyle(3010)
                hist.SetMarkerSize(0.001)

                labelscalefact = 1. / (1. - self.VerticalCanvasSplit)
                hist.GetYaxis().SetTitle(y_title)
                hist.GetYaxis().SetTitleSize(self.t_size*labelscalefact)
                hist.GetYaxis().SetLabelSize(self.t_size*labelscalefact)
                hist.GetYaxis().SetRangeUser(y_min, y_max)
                hist.GetYaxis().SetTitleOffset(0.8)

                hist.GetXaxis().SetLabelSize(self.t_size*labelscalefact)
                hist.GetXaxis().SetTitleSize(self.x_title_size*labelscalefact)
                hist.GetXaxis().SetTitleOffset(1.4)

                #hist.Draw("AXIS")
                hist.Draw("E2")
            else:
                # start to calculate the ratio
                if reverse: # MC/Data
                    this_hist = hist.Clone(hist.GetName()+"_cp")
                    this_hist.Divide(h_refer)
                else: # Data/MC
                    this_hist = h_refer.Clone(hist.GetName()+"_cpDI")
                    this_hist.SetLineColor(hist.GetLineColor())
                    this_hist.Divide(hist)
                    print "Yields:",hist.Integral(), h_refer.Integral()

                self.totalObj.append(this_hist)
                this_hist.Draw("EP SAME")


    def stack_hists(self,
        hist_list, tag_list, out_name,
        x_title, y_title,
        is_log=False, has_data=True,
    ):
        ## not saving the plots..

        # In hist_list, the data should be first element, if has_data
        if len(hist_list) > len(self.COLORS):
            print "I don't have enough colors", len(hist_list), len(self.COLORS)
            return

        # clone current histograms, so that inputs are untouched.
        hist_list_cp = [] # a list of non-data histograms
        h_data = None  # the first element is assumed to be data

        hist_sum = None
        for i, hist in enumerate(hist_list):
            new_hist = hist.Clone(hist.GetName()+"_clone")
            if has_data and i!=0:
                color = self.COLORS[i-1]
            elif has_data:
                color = 1
            else:
                color = self.COLORS[i]

            new_hist.SetLineColor(color)

            if i==0 and has_data:
                # decorate data points
                new_hist.SetMarkerStyle(20)
                new_hist.SetMarkerSize(1.2)
                h_data = new_hist
                self.get_offset(h_data)
                continue
            elif i==0 or hist_sum is None:
                hist_sum = hist
            else:
                hist_sum.Add(hist)

            new_hist.SetFillColor(color)
            hist_list_cp.append(new_hist)

        # always plot the smallest component in the bottom
        hist_sorted_list = sorted(hist_list_cp, key=lambda k:k.Integral())
        hs = ROOT.THStack("hs", "")
        for hist in hist_sorted_list:
            hs.Add(hist)

        # start to plot them
        if self.add_ratio and has_data:
            self.prepare_2pad_canvas("canvas", 600, 600)
            self.pad2.cd()
            self.pad2.SetGridy()
            hist_sum.SetLineColor(8)
            new_data_copy = h_data.Clone("data_copy")
            #self.add_ratio_panel([new_data_copy, hist_sum], y_title, 0.5, 1.52)
            self.add_ratio_panel([new_data_copy, hist_sum], "Data/MC", 0.55, 1.42)
            self.pad1.cd()
        else:
            self.can = ROOT.TCanvas("canvas", "canvas", 600, 600)

        y_max = hs.GetMaximum()
        y_min = hs.GetMinimum()
        if has_data and h_data is not None:
            if y_max < h_data.GetMaximum():
                y_max = h_data.GetMaximum()
            if y_min > h_data.GetMinimum():
                y_min = h_data.GetMinimum()

        if has_data:
            this_hist = h_data
        else:
            this_hist = hs

        if is_log:
            if self.add_ratio:
                self.pad1.SetLogy()
            else:
                self.can.SetLogy()
            this_hist.GetYaxis().SetRangeUser(4E-3, y_max*1e3)
        else:
            this_hist.GetYaxis().SetRangeUser(1E-3, y_max*1.1)

        #this_hist.SetNdivisions(8, "X")
        this_hist.SetXTitle(x_title)
        this_hist.SetYTitle(y_title)

        if has_data:
            h_data.Draw("EP")
            hs.Draw("HISTsame")
            h_data.Draw("AXISsame")
            h_data.Draw("EPsame")
        else:
            hs.Draw("HIST")

        # add legend
        hist_id = 0
        if has_data:
            hist_all = [h_data] + hist_list_cp
        else:
            hist_all = hist_list_cp

        if self.show_sum_bkg:
            legend = self.get_legend(len(hist_all)+1)
        else:
            legend = self.get_legend(len(hist_all))

        for hist, tag in zip(hist_all, tag_list):
            if has_data and hist_id == 0:
                legend.AddEntry(hist, tag+" {:.0f}".format(hist.Integral()), "LP")
                # add sum of background...
                if self.show_sum_bkg:
                    hist_sum.SetFillColor(0)
                    hist_sum.SetLineColor(0)
                    legend.AddEntry(hist_sum, "Total Bkg {:.0f}".format(hist_sum.Integral()), "F")
            else:
                legend.AddEntry(hist, tag+" {:.1f}".format(hist.Integral()), "F")
            hist_id += 1

        legend.Draw("same")
        self.add_atlas()
        self.add_lumi()

        if out_name is not None:
            self.can.SaveAs(out_name)


    def get_legend(self, nentries):
        x_min = self.x_offset
        x_max = x_min + 0.3
        # y_max = self.y_offset-self.text_size*2-0.001
        # y_min = y_max - self.t_size*nentries
        y_max = self.y_offset
        y_min = y_max - self.t_size*nentries

        legend = ROOT.TLegend(x_min, y_min, x_max, y_max)
        legend.SetFillColor(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        # legend.SetTextSize(0.035)
        legend.SetTextSize(self.t_size)

        return legend


    def add_atlas(self):
        adder.add_text(self.x_off_atlas, self.y_offset-0.03,
                      1, "#bf{#it{ATLAS}} "+self.status,
                      self.text_size)

    def add_lumi(self):
        adder.add_text(self.x_off_atlas,
                      self.y_offset - self.text_size - 0.007-0.03,
                      1, "13 TeV, "+str(self.lumi)+" fb^{-1}",
                      self.text_size)

    def get_offset(self, hist):
        max_bin = hist.GetMaximumBin()

        last_bin = hist.GetXaxis().GetLast()
        first_bin = hist.GetXaxis().GetFirst()
        if max_bin < first_bin + (last_bin - first_bin)/2.:
            self.x_offset = 0.60
            self.x_off_atlas = 0.2


    def stack(self, hist_list):
        hist_list_cp = [] # a list of non-data histograms
        for hist in hist_list:
            hist_list_cp.append( hist.Clone(hist.GetName()+"stackClone"))

        self.totalObj.append(hist_list_cp)

        hist_sorted_list = sorted(hist_list_cp, key=lambda k:k.Integral())
        hs = ROOT.THStack("hs", "")
        hist_sum = None
        for hist in hist_sorted_list:
            hs.Add(hist)
            if hist_sum is None:
                hist_sum = hist.Clone(hist.GetName()+"sumClone")
            else:
                hist_sum.Add(hist)

        self.totalObj.append(hist_sum)
        return hist_sum, hs

    def color(self, hist_list, no_fill=False):
        for i, hist in enumerate(hist_list):
            color = self.COLORS[i]
            hist.SetLineColor(color)
            if not no_fill:
                hist.SetFillColor(color)
            hist.SetMarkerColor(color)
            hist.SetLineStyle(self.LINE_STYLE[i])

    def set_y_range(self, hist_list, is_logY):
        hist = hist_list[0]
        y_max = max(hist_list, key=lambda x: x.GetMaximum()).GetMaximum()
        y_min = min(hist_list, key=lambda x: x.GetMinimum()).GetMinimum()

        if is_logY:
            if self.add_ratio:
                self.pad1.SetLogy()
            else:
                self.can.SetLogy()

            hist.GetYaxis().SetRangeUser(4E-3, y_max*1e2)
        else:
            if y_min < 0:
                y_min *= 1.1
            else:
                y_min *= 0.9
            hist.GetYaxis().SetRangeUser(y_min, y_max*1.1)
            pass

    def compare_hists(self, hist_list, tag_list, **kwargs):
        """
        a list of histograms,
        Key words include:
            ratio_title, ratio_range, logY, out_name
            no_fill, x_offset, draw_option,
            add_yields, add_ratio,
            out_folder, label
        """
        self.del_obj()

        if len(hist_list) < 2:
            print "not enough hitograms for comparison"
            return
        try:
            no_fill = kwargs["no_fill"]
        except KeyError:
            no_fill = False

        try:
            add_ratio = self.add_ratio = kwargs["add_ratio"]
        except KeyError:
            add_ratio = self.add_ratio


        self.color(hist_list, no_fill)

        if add_ratio:
            self.prepare_2pad_canvas('canvas', 600, 600)
            self.pad2.cd()
            try:
                ratio_title = kwargs["ratio_title"]
            except KeyError:
                ratio_title = "MC/Data"

            try:
                ratio_x, ratio_y = kwargs["ratio_range"]
            except KeyError:
                ratio_x, ratio_y = 0.55, 1.42

            self.add_ratio_panel(hist_list, ratio_title, ratio_x, ratio_y, True)
            self.pad1.cd()
        else:
            self.text_size = 0.04
            self.can = ROOT.TCanvas("canvas", "canvas", 600, 600)

        self.set_y_offset()
        try:
            self.x_offset = kwargs["x_offset"]
        except KeyError:
            self.get_offset(hist_list[0])

        try:
            is_logy = kwargs["logY"]
        except KeyError:
            is_logy = False

        legend = self.get_legend(len(hist_list))

        self.set_y_range(hist_list, is_logy)
        try:
            draw_option = kwargs["draw_option"]
        except KeyError:
            draw_option = "HIST"

        try:
            add_yield = kwargs['add_yields']
        except KeyError:
            add_yield = False

        for i, hist in enumerate(hist_list):
            if add_yield:
                legend.AddEntry(hist, "{}: {:.3f}".format(tag_list[i], hist.Integral()))
            else:
                legend.AddEntry(hist, tag_list[i])

            if i==0:
                hist.Draw(draw_option)
            else:
                hist.Draw(draw_option+" SAME")

        legend.Draw("same")
        self.add_atlas()
        self.add_lumi()
        if 'label' in kwargs.keys():
            adder.add_text(self.x_off_atlas, 
                           self.y_offset - self.text_size*2 - 0.007-0.03,
                           1, kwargs['label'], self.text_size)

        try:
            out_name = kwargs["out_name"]
        except KeyError:
            out_name = "TEST"

        try:
            out_folder = kwargs["out_folder"]
        except KeyError:
            out_folder = "./"

        self.mkdir_p(out_folder)

        if is_logy:
            self.can.SaveAs(out_folder+"/"+out_name+"_Log.eps")
            self.can.SaveAs(out_folder+"/"+out_name+"_Log.pdf")
        else:
            self.can.SaveAs(out_folder+"/"+out_name+".eps")
            self.can.SaveAs(out_folder+"/"+out_name+".pdf")

    def set_y_offset(self):
        if not self.add_ratio:
            self.y_offset = 0.88

    def del_obj(self):
        for obj in self.totalObj:
            del obj

        if self.can:
            del self.can

        if self.pad1:
            del self.pad1

        if self.pad2:
            del self.pad2

    def set_palette(self):
        red_list = [1., 1., 0.]
        white_list = [0., 1.0, 0.0]
        blue_list = [0., 1., 1.]
        length_list = [0.0, 0.5, 1.0]
        nb = 50
        ROOT.TColor.CreateGradientColorTable(
            3, 
            array('d', length_list),
            array('d', red_list),
            array('d', white_list),
            array('d', blue_list),
            nb)

    def plot_correlation(self, corr_hist, out_name, thre=0.05):
        # plot the 2D,
        if self.can is None:
            self.can = ROOT.TCanvas("canvas", "canvas", 1200, 800)
            self.can.SetRightMargin(0.08)
            self.can.SetLeftMargin(0.18)

        h_new = self.remove_weak_correlation_bins(corr_hist, thre)

        #ROOT.gStyle.SetPalette(55)
        # palette = array('d', [15, 20, 23, 30, 32])
        # ROOT.gStyle.SetPalette(5, palette)
        self.set_palette()
        h_new.Draw("colz text")
        h_new.SetMarkerSize(0.7)
        ROOT.gStyle.SetPaintTextFormat("3.2f")
        h_new.GetXaxis().SetLabelSize(0.02)
        h_new.GetXaxis().LabelsOption('u')
        h_new.GetYaxis().SetLabelSize(0.02)
        h_new.GetZaxis().SetLabelSize(0.02)
        h_new.GetZaxis().SetRangeUser(-1, 1)
        h_new.GetXaxis().SetTickSize(0.)
        h_new.GetYaxis().SetTickSize(0.)

        self.can.SaveAs(out_name+".eps")
        self.can.SaveAs(out_name+".pdf")

    def remove_weak_correlation_bins(self, h2d, thre=0.05):
        empty_xbins = []
        for xbin in range(h2d.GetNbinsX()):
            is_empty = True
            for ybin in range(h2d.GetNbinsY()):
                value = h2d.GetBinContent(xbin+1, ybin+1)
                if abs(value) > thre and abs(value - 1) > 1E-3:
                    is_empty = False
                    break
            if is_empty:
                empty_xbins.append(xbin+1)

        final_bins = h2d.GetNbinsX() - len(empty_xbins)
        org_bins = h2d.GetNbinsX()
        print "original bins", org_bins
        print "final bins: ", final_bins
        #print empty_xbins
        if True:
            h2d_new = ROOT.TH2D("reduced_correlation", "reduced_correlation",
                                final_bins, 0.5, final_bins+0.5,
                                final_bins, 0.5, final_bins+0.5)
            new_x = 0
            for xbin in range(h2d.GetNbinsX()):
                if (xbin+1) in empty_xbins:
                    continue
                new_x += 1 
                new_y = 0
                for ybin in range(h2d.GetNbinsY()):
                    if (org_bins - ybin) in empty_xbins:
                        continue
                    new_y += 1
                    value = h2d.GetBinContent(xbin+1, ybin+1)
                    # print xbin,ybin,new_x,new_y, value

                    h2d_new.SetBinContent(new_x, new_y, value)
                    h2d_new.GetXaxis().SetBinLabel(new_x, h2d.GetXaxis().GetBinLabel(xbin+1))
                    h2d_new.GetYaxis().SetBinLabel(new_y, h2d.GetYaxis().GetBinLabel(ybin+1))

            return h2d_new

        return None

    def create_graph_pulls(self, input_list):
        if input_list is None or len(input_list) < 0:
            return None
        
        info = input_list[0]

        if len(info) != 4:
            print "input should be a list of tuple of size of 4"
            return None
    
        label_list = []
        x_list = []
        y_nom_list = []
        y_up_list = []
        y_down_list = []
        x_val = 0.5
        for item in input_list:
            name, nom, up, down = item
            label_list.append(name)
            x_list.append(x_val)
            y_nom_list.append(nom)
            y_up_list.append(up)
            y_down_list.append(abs(down))
            x_val += 1

        nbins = len(x_list)
        gr = ROOT.TGraphAsymmErrors(
            nbins,
            array('d', x_list),
            array('d', y_nom_list),
            array('d', [0.]*nbins), # error of x (low)
            array('d', [0.]*nbins), # error of x (high)
            array('d', y_down_list), # error of y (low)
            array('d', y_up_list) # error of y (high)
        )
        nbins += 1
        gr_one = ROOT.TGraphErrors(
            nbins,
            array('d', [x-0.5 for x in x_list]+[nbins]),
            array('d', [0.]*nbins),
            array('d', [0.]*nbins), # error of x
            array('d', [1.]*nbins) # error of y
        )
        gr_two = ROOT.TGraphErrors(
            nbins,
            array('d', [x-0.5 for x in x_list]+[nbins]),
            array('d', [0.]*nbins),
            array('d', [0.]*nbins), # error of x
            array('d', [2.]*nbins) # error of y
        )
        return gr,gr_one,gr_two

    @staticmethod
    def mkdir_p(path):
        try:
            os.makedirs(path)
            print path,"is created"
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
