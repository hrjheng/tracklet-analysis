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
from plotUtil import GetHistogram

gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)

TickSize = 0.03
AxisTitleSize = 0.05
AxisLabelSize = 0.04
LeftMargin = 0.15
RightMargin = 0.08
TopMargin = 0.08
BottomMargin = 0.13

def Draw_1DhistsComp(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, outname):
    # norm1 is to normalize case0 to 1 
    # others are normalized w.r.t case0
    color = ['#1B1A17', '#086972', '#035397', '#9B0000']
    legtext = ['Nominal', '+1% splitting', '+2% splitting', '+5% splitting']
    ymax = 0;
    norm_case0 = lhist[0].Integral(-1, -1)
    for i, h in enumerate(lhist):
        h.Sumw2()
        if norm1:
            h.Scale(1. / norm_case0)

        if h.GetMaximum() > ymax:
            ymax = h.GetMaximum()

    binwidth_bin1 = lhist[0].GetXaxis().GetBinWidth(1)
    binwidth_bin2 = lhist[0].GetXaxis().GetBinWidth(2)
    printbinwidth = True
    if binwidth_bin1 != binwidth_bin2:
        printbinwidth = False

    c = TCanvas('c', 'c', 800, 700)
    if logx:
        c.SetLogx()
    if logy:
        c.SetLogy()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    if printbinwidth:
        if norm1:
            if Ytitle_unit == '':
                lhist[0].GetYaxis().SetTitle('Normalized entries / ({:g})'.format(binwidth_bin1))
            else:
                lhist[0].GetYaxis().SetTitle('Normalized entries / ({:g} {unit})'.format(binwidth_bin1, unit=Ytitle_unit))
        else:
            if Ytitle_unit == '':
                lhist[0].GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth_bin1))
            else:
                lhist[0].GetYaxis().SetTitle('Entries / ({:g} {unit})'.format(binwidth_bin1, unit=Ytitle_unit))
    else:
        if norm1:
            if Ytitle_unit == '':
                lhist[0].GetYaxis().SetTitle('Normalized entries')
            else:
                lhist[0].GetYaxis().SetTitle('Normalized entries {unit})'.format(unit=Ytitle_unit))
        else:
            if Ytitle_unit == '':
                lhist[0].GetYaxis().SetTitle('Entries')
            else:
                lhist[0].GetYaxis().SetTitle('Entries {unit}'.format(unit=Ytitle_unit))
        # lhist[0].GetXaxis().SetRangeUser(lhist[0].GetBinLowEdge(1)-binwidth, lhist[0].GetBinLowEdge(lhist[0].GetNbinsX())+2*binwidth)
    if logy:
        lhist[0].GetYaxis().SetRangeUser(lhist[0].GetMinimum(0)* 0.5, ymax * 30)
    else:
        lhist[0].GetYaxis().SetRangeUser(
            0., ymax * ymaxscale)
    lhist[0].GetXaxis().SetTitle(XaxisName)
    lhist[0].GetXaxis().SetTickSize(TickSize)
    lhist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    lhist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    lhist[0].GetYaxis().SetTickSize(TickSize)
    lhist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    lhist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    lhist[0].GetXaxis().SetTitleOffset(1.1)
    lhist[0].GetYaxis().SetTitleOffset(1.4)
    for i, h in enumerate(lhist):
        if i == 0:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('hist')
        else:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('histsame')

    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.15,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    leg.AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "")
    leg.Draw()

    leg1 = TLegend(LeftMargin+0.04, (1-TopMargin)-0.21,
                   LeftMargin+0.34, (1-TopMargin)-0.03)
    leg1.SetTextSize(0.035)
    leg1.SetFillStyle(0)
    for i, h in enumerate(lhist):
        leg1.AddEntry(h, legtext[i], "l")
    leg1.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')

if __name__ == '__main__':
    plotpath = './ClusterSplit/'
    os.makedirs(plotpath, exist_ok=True)

    histfilepath = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/'
    filelist = [histfilepath + 'Hists_Tracklets_RandhitCase0_ClusSplitCase0_MisAlignNum0_dRcut0p5.root',
                histfilepath + 'Hists_Tracklets_RandhitCase0_ClusSplitCase1_MisAlignNum0_dRcut0p5.root',
                histfilepath + 'Hists_Tracklets_RandhitCase0_ClusSplitCase2_MisAlignNum0_dRcut0p5.root',
                histfilepath + 'Hists_Tracklets_RandhitCase0_ClusSplitCase3_MisAlignNum0_dRcut0p5.root']
    
    l_hM_CluspairDr_l1 = []
    l_hM_CluspairDr_l2 = []
    l_hM_CluspairDr_l3 = []

    for i, fname in enumerate(filelist):
        hM_CluspairDr_l1 = GetHistogram(fname, 'hM_ClusPairDr_layer1')
        hM_CluspairDr_l2 = GetHistogram(fname, 'hM_ClusPairDr_layer2')
        hM_CluspairDr_l3 = GetHistogram(fname, 'hM_ClusPairDr_layer3')
        l_hM_CluspairDr_l1.append(hM_CluspairDr_l1)
        l_hM_CluspairDr_l2.append(hM_CluspairDr_l2)
        l_hM_CluspairDr_l3.append(hM_CluspairDr_l3)

    # print(len(l_hM_CluspairDr), l_hM_CluspairDr[1].Integral(-1,-1))
    # Draw_1DhistsComp(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, outname)
    Draw_1DhistsComp(l_hM_CluspairDr_l1, False, False, False, 1.5, 'Cluster-Cluster #DeltaR (layer 1)', '', plotpath+'ClusterPairDr_l1_comp')
    Draw_1DhistsComp(l_hM_CluspairDr_l2, False, False, False, 1.5, 'Cluster-Cluster #DeltaR (layer 2)', '', plotpath+'ClusterPairDr_l2_comp')
    Draw_1DhistsComp(l_hM_CluspairDr_l3, False, False, False, 1.5, 'Cluster-Cluster #DeltaR (layer 3)', '', plotpath+'ClusterPairDr_l3_comp')
