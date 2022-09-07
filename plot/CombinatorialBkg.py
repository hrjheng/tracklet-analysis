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
# from plotUtil import *

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
    legtext = ['Nominal', '+1% random hits', '+5% random hits', '+10% random hits']
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
    plotpath = './CombinatorialBkg/'
    os.makedirs(plotpath, exist_ok=True)

    minitreepath = '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal/'
    chain_Case1 = TChain('minitree_12')
    for f in glob.glob(minitreepath + 'TrackletAna_minitree_layer12_Evt0to4000_RandhitCase0_MisAlignNum0.root'):
        chain_Case1.Add(f)

    chain_Case2 = TChain('minitree_12')
    for f in glob.glob(minitreepath + 'TrackletAna_minitree_layer12_Evt0to4000_RandhitCase1_MisAlignNum0.root'):
        chain_Case2.Add(f)

    chain_Case3 = TChain('minitree_12')
    for f in glob.glob(minitreepath +  'TrackletAna_minitree_layer12_Evt0to4000_RandhitCase2_MisAlignNum0.root'):
        chain_Case3.Add(f)
    
    chain_Case4 = TChain('minitree_12')
    for f in glob.glob(minitreepath +  'TrackletAna_minitree_layer12_Evt0to4000_RandhitCase3_MisAlignNum0.root'):
        chain_Case4.Add(f)

    list_chain = [chain_Case1, chain_Case2, chain_Case3, chain_Case4]

    NBins = 100
    xbin = numpy.logspace(-4, math.log10(0.5), NBins+1)

    l_hM_dR_reco_logx = []
    l_hM_dR_reco = []
    for ichain, chain in enumerate(list_chain):
        hM_dR_reco_logx = TH1F('hM_dR_reco_logx', 'hM_dR_reco_logx', NBins, xbin)
        hM_dR_reco = TH1F('hM_dR_reco', 'hM_dR_reco', 200, 0, 0.5)
        for index, ev in enumerate(chain):
            for elem in ev.recotkl_dR:
                # print(elem)
                hM_dR_reco_logx.Fill(elem)
                hM_dR_reco.Fill(elem)

        l_hM_dR_reco_logx.append(hM_dR_reco_logx)
        l_hM_dR_reco.append(hM_dR_reco)

    Draw_1DhistsComp(l_hM_dR_reco_logx, True, True, True, 1.3,
                     'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_randhit_dR_logX_comp')
    Draw_1DhistsComp(l_hM_dR_reco, True, False, True, 1.3,
                     'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_randhit_dR_comp')