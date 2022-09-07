#! /usr/bin/env python
from optparse import OptionParser
import sys
import os
import datetime
from array import *
from ROOT import *
import numpy
import math
import glob
from plotUtil import *

gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)

TickSize = 0.03
AxisTitleSize = 0.05
AxisLabelSize = 0.04
LeftMargin = 0.15
RightMargin = 0.08
TopMargin = 0.08
BottomMargin = 0.11

def plot_Stack(totalhist, list_hist, list_legtext, logy, ymaxscale, XaxisName, Ytitle_unit, plotname):
    color = ['#073b4c','#118ab2','#06d6a0','#ffd166','#ef476f']
    binwidth = totalhist.GetXaxis().GetBinWidth(1)
    hs = THStack('hs','hs');
    for i, hist in enumerate(list_hist):
        hist.SetLineColor(1)
        hist.SetLineWidth(1)
        hist.SetFillColor(TColor.GetColor(color[i]))
        hs.Add(hist)

    c = TCanvas('c', 'c', 800, 700)
    if logy:
        c.SetLogy()
    c.cd()
    hs.Draw()
    hs.GetXaxis().SetTitle(XaxisName)
    if Ytitle_unit == '':
        hist.GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth))
    else:
        hist.GetYaxis().SetTitle('Entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))
    hs.GetXaxis().SetTitleSize(AxisTitleSize)
    hs.GetYaxis().SetTitleSize(AxisTitleSize)
    hs.GetXaxis().SetTickSize(TickSize)
    hs.GetYaxis().SetTickSize(TickSize)
    hs.GetXaxis().SetLabelSize(AxisLabelSize)
    hs.GetYaxis().SetLabelSize(AxisLabelSize)
    hs.GetYaxis().SetTitleOffset(1.3)
    hs.SetMaximum(totalhist.GetMaximum() * ymaxscale)
    hs.SetMinimum(0.001)
    totalhist.SetLineWidth(3)
    totalhist.Draw('histsame')
    c.Update()
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.15,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    leg.AddEntry('', 'Au+Au #sqrt{s_{NN}}=200 GeV', '')
    leg.Draw()
    leg1 = TLegend(LeftMargin+0.05, (1-TopMargin)-0.22,
                   LeftMargin+0.3, (1-TopMargin)-0.01)
    # leg1.SetNColumns(3)
    leg1.SetTextSize(0.03)
    leg1.SetFillStyle(0)
    for i, text in enumerate(list_legtext):
        leg1.AddEntry(list_hist[i], text, 'f')
    leg1.Draw()
    c.Update()
    c.SaveAs(plotname+'.png')
    c.SaveAs(plotname+'.pdf')


if __name__ == '__main__':
    plotpath = './TrackletAna/'
    os.makedirs(plotpath, exist_ok=True)

    hM_Eta_proto = GetHist('hM_Eta_proto', 0)
    hM_Eta_proto_PVzRange1 = GetHist('hM_Eta_proto_PVzRange1', 0)
    hM_Eta_proto_PVzRange2 = GetHist('hM_Eta_proto_PVzRange2', 0)
    hM_Eta_proto_PVzRange3 = GetHist('hM_Eta_proto_PVzRange3', 0)
    hM_Eta_proto_PVzRange4 = GetHist('hM_Eta_proto_PVzRange4', 0)
    hM_Eta_proto_PVzRange5 = GetHist('hM_Eta_proto_PVzRange5', 0)

    list_hM_Eta_proto_layer12 = [hM_Eta_proto_PVzRange1[0], hM_Eta_proto_PVzRange2[0], hM_Eta_proto_PVzRange3[0], hM_Eta_proto_PVzRange4[0], hM_Eta_proto_PVzRange5[0]]
    list_hM_Eta_proto_layer23 = [hM_Eta_proto_PVzRange1[1], hM_Eta_proto_PVzRange2[1], hM_Eta_proto_PVzRange3[1], hM_Eta_proto_PVzRange4[1], hM_Eta_proto_PVzRange5[1]]
    list_hM_Eta_proto_layer13 = [hM_Eta_proto_PVzRange1[2], hM_Eta_proto_PVzRange2[2], hM_Eta_proto_PVzRange3[2], hM_Eta_proto_PVzRange4[2], hM_Eta_proto_PVzRange5[2]]

    l_legtext = ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm']

    plot_Stack(hM_Eta_proto[0], list_hM_Eta_proto_layer12, l_legtext, False, 1.5, 'Proto-tracklet #eta', '', plotpath+'TrackletAna_MisAlignNum0/ProtoTklEta_layer12_Stack')
    plot_Stack(hM_Eta_proto[1], list_hM_Eta_proto_layer23, l_legtext, False, 1.5, 'Proto-tracklet #eta', '', plotpath+'TrackletAna_MisAlignNum0/ProtoTklEta_layer23_Stack')
    plot_Stack(hM_Eta_proto[2], list_hM_Eta_proto_layer13, l_legtext, False, 1.5, 'Proto-tracklet #eta', '', plotpath+'TrackletAna_MisAlignNum0/ProtoTklEta_layer13_Stack')

    