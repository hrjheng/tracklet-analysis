#! /usr/bin/env python
from optparse import OptionParser
import sys
import os
import datetime
from array import *
from ROOT import *
import numpy as np
from plotUtil import *

gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)
gStyle.SetPalette(kLightTemperature)
# TColor.InvertPalette()
gStyle.SetPaintTextFormat(".3g");

TickSize = 0.03
AxisTitleSize = 0.05
AxisLabelSize = 0.04
LeftMargin = 0.15
RightMargin = 0.08
TopMargin = 0.08
BottomMargin = 0.13

def Draw_1Dhist_cumulative(hist, logy, ymaxscale, XaxisName, YaxisName, Ytitle_unit, outname):
    hist.Sumw2()
    binwidth = hist.GetXaxis().GetBinWidth(1)
    c = TCanvas('c', 'c', 800, 700)
    if logy:
        c.SetLogy()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    if Ytitle_unit == '':
        hist.GetYaxis().SetTitle('{} / ({:g})'.format(YaxisName, binwidth))
    else:
        hist.GetYaxis().SetTitle('{} / ({:g} {unit})'.format(YaxisName, binwidth, unit=Ytitle_unit))

    # hist.GetXaxis().SetRangeUser(hist.GetBinLowEdge(1)-binwidth, hist.GetBinLowEdge(hist.GetNbinsX())+2*binwidth)
    if logy:
        hist.GetYaxis().SetRangeUser(hist.GetMinimum(0)*0.5, (hist.GetMaximum()) * ymaxscale)
    else:
        hist.GetYaxis().SetRangeUser(0., (hist.GetMaximum()) * ymaxscale)
    hist.GetXaxis().SetTitle(XaxisName)
    hist.GetXaxis().SetTickSize(TickSize)
    hist.GetXaxis().SetTitleSize(AxisTitleSize)
    hist.GetXaxis().SetLabelSize(AxisLabelSize)
    hist.GetYaxis().SetTickSize(TickSize)
    hist.GetYaxis().SetTitleSize(AxisTitleSize)
    hist.GetYaxis().SetLabelSize(AxisLabelSize)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.35)
    hist.SetLineColor(1)
    hist.SetLineWidth(2)
    hist.Draw('hist')
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.1,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    # leg.AddEntry("","Au+Au #sqrt{s_{NN}}=200 GeV","")
    leg.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_1DhistsComp2_wRatio(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname):
    # ratio = TGraphAsymmErrors()
    ratio = lhist[1].Clone();
    ratio.Divide(lhist[0])
    xrange_low = lhist[0].GetBinLowEdge(1)
    xrange_high = lhist[0].GetBinLowEdge(lhist[0].GetNbinsX())+lhist[0].GetXaxis().GetBinWidth(lhist[0].GetNbinsX())

    color = ['#035397', '#9B0000']
    ymax = -1
    ymin = 10e10
    for h in lhist:
        h.Sumw2()
        if h.GetMaximum() > ymax:
            ymax = h.GetMaximum()
        if h.GetMinimum(0) < ymin:
            ymin = h.GetMinimum(0)
        if norm1:
            h.Scale(1. / h.Integral(-1, -1))
            ymax = h.GetMaximum()
            ymin = h.GetMinimum(0)

    binwidth_bin1 = lhist[0].GetXaxis().GetBinWidth(1)
    binwidth_bin2 = lhist[0].GetXaxis().GetBinWidth(2)
    printbinwidth = True
    if binwidth_bin1 != binwidth_bin2:
        printbinwidth = False

    c = TCanvas('c', 'c', 800, 700)
    pad1 = TPad( 'pad1', ' ', 0, 0.3, 1, 1)
    pad2 = TPad( 'pad2', ' ', 0, 0, 1, 0.3)
    if logx:
        pad1.SetLogx()
        pad2.SetLogy()
    if logy:
        pad1.SetLogy()

    pad1.SetRightMargin(RightMargin)
    pad1.SetTopMargin(TopMargin)
    pad1.SetLeftMargin(LeftMargin)
    pad1.SetBottomMargin(0.01)
    pad1.Draw()
    pad2.SetGridy()
    pad2.SetRightMargin(RightMargin)
    pad2.SetLeftMargin(LeftMargin)
    pad2.SetTopMargin(0.01)
    pad2.SetBottomMargin(0.35)
    pad2.Draw() # Draw the TPad on the TCanvas before plotting something on TPad  
    # pad2.SetFrameFillColor(0)
    # pad2.SetFrameBorderMode(0)
    # pad2.SetBorderMode(0)
    # pad2.SetBorderSize(0)
    pad1.cd()
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
        lhist[0].GetYaxis().SetRangeUser(ymin * 0.1, ymax * 100)
    else:
        lhist[0].GetYaxis().SetRangeUser(0., ymax * ymaxscale)
    textscale = 1.2
    lhist[0].GetXaxis().SetTickSize(TickSize)
    lhist[0].GetYaxis().SetTickSize(TickSize)
    lhist[0].GetYaxis().SetTitleSize(AxisTitleSize*textscale)
    lhist[0].GetYaxis().SetLabelSize(AxisLabelSize*textscale)
    lhist[0].GetYaxis().SetTitleOffset(1.1)
    for i, h in enumerate(lhist):
        if i == 0:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('hist')
        else:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('histsame')

    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.19,
                  (1-RightMargin)-0.1/textscale, (1-TopMargin)-0.05)
    leg.SetTextSize(0.045*textscale)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    leg.AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "")
    leg.Draw()

    leg1 = TLegend(LeftMargin+0.04, (1-TopMargin)-0.17,
                   LeftMargin+0.34, (1-TopMargin)-0.05)
    leg1.SetTextSize(0.035*textscale)
    leg1.SetFillStyle(0)
    for i, h in enumerate(lhist):
        leg1.AddEntry(h, legtext[i], "l")
    leg1.Draw()
    c.Update()

    pad2.cd()
    ratio.GetYaxis().SetTitle('Ratio')
    ratio.GetYaxis().SetTitleSize(AxisTitleSize*2.7)
    ratio.GetYaxis().SetTitleOffset(0.5)
    ratio.GetYaxis().SetLabelSize(AxisLabelSize*2.7)
    ratio.GetYaxis().SetNdivisions(205)
    ratio.GetYaxis().SetRangeUser(0, ratio.GetMaximum()*1.1)
    ratio.GetXaxis().SetTitle(XaxisName)
    ratio.GetXaxis().SetTitleSize(AxisTitleSize*2.7)
    ratio.GetXaxis().SetTitleOffset(1)
    ratio.GetXaxis().SetLabelSize(AxisLabelSize*2.7)
    ratio.GetXaxis().SetTickSize(TickSize)
    ratio.GetXaxis().SetRangeUser(xrange_low, xrange_high)
    # ratio.SetMarkerColor(1)
    # ratio.SetMarkerSize(0.6)
    # ratio.SetMarkerStyle(20)
    ratio.SetLineColor(1)
    ratio.SetLineWidth(2)
    ratio.Draw("hist")
    c.Update()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

if __name__ == '__main__':
    plotpath = './Bkg_Beta/'
    os.makedirs(plotpath, exist_ok=True)

    # df = ROOT.RDataFrame('minitree_12', '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer12_Evt0to2000_RandhitCase0_MisAlignNum0.root')
    # hM_recotkl_dR = df.Histo1D(("hM_recotkl_dR", " ", 50, 0, 0.05), "recotkl_dR")
    # hM_recotkl_genmatched_dR = df.Histo1D(("hM_recotkl_genmatched_dR", " ", 50, 0, 0.05), "recotkl_genmatched_dR")
    # hM_recotkl_dR_fr = df.Histo1D(("hM_recotkl_dR_fr", " ", 100, 0, 0.5), "recotkl_dR")
    # hM_recotkl_genmatched_dR_fr = df.Histo1D(("hM_recotkl_genmatched_dR_fr", " ", 100, 0, 0.5), "recotkl_genmatched_dR")
    # hM_recotkl_deta = df.Histo1D(("hM_recotkl_deta", " ", 100, -0.05, 0.05), "recotkl_deta")
    # hM_recotkl_genmatched_deta = df.Histo1D(("hM_recotkl_genmatched_deta", " ", 100, -0.05, 0.05), "recotkl_genmatched_deta")
    # hM_recotkl_dphi = df.Histo1D(("hM_recotkl_dphi", " ", 100, -0.05, 0.05), "recotkl_dphi")
    # hM_recotkl_genmatched_dphi = df.Histo1D(("hM_recotkl_genmatched_dphi", " ", 100, -0.05, 0.05), "recotkl_genmatched_dphi")
    hM_recotkl_dR = GetHist('hM_recotkl_dR', 0)
    hM_recotkl_dR_fr = GetHist('hM_recotkl_dR_fr', 0)
    hM_recotkl_genmatched_dR_ghadPtIncl = GetHist('hM_recotkl_genmatched_dR_ghadPtIncl', 0)
    hM_recotkl_genmatched_dR_fr_ghadPtIncl = GetHist('hM_recotkl_genmatched_dR_fr_ghadPtIncl', 0)
    hM_recotkl_deta = GetHist('hM_recotkl_deta', 0)
    hM_recotkl_genmatched_deta_ghadPtIncl = GetHist('hM_recotkl_genmatched_deta_ghadPtIncl', 0)
    hM_recotkl_dphi = GetHist('hM_recotkl_dphi', 0)
    hM_recotkl_genmatched_dphi_ghadPtIncl = GetHist('hM_recotkl_genmatched_dphi_ghadPtIncl', 0)

    hM_Eta_vtxZ_reco = []
    hM_Eta_vtxZ_reco_genmatched = []
    hM_Eta_vtxZ_anabin_reco = []
    hM_Eta_vtxZ_anabin_reco_genmatched = []
    hM_anabin_1minusBeta = []
    layer = [12,23,13]
    
    for i in range(20):
        hM_Eta_vtxZ_reco = GetHist('hM_Eta_vtxZ_reco_centbin{}'.format(i), 0)
        hM_Eta_vtxZ_anabin_reco = GetHist('hM_Eta_vtxZ_anabin_reco_centbin{}'.format(i), 0)
        hM_Eta_vtxZ_reco_genmatched = GetHist('hM_Eta_vtxZ_reco_genmatched_centbin{}'.format(i), 0)
        hM_Eta_vtxZ_anabin_reco_genmatched = GetHist('hM_Eta_vtxZ_anabin_reco_genmatched_centbin{}'.format(i), 0)

        hM_Eta_vtxZ_anabin_ratio = GetHist('hM_Eta_vtxZ_anabin_reco_genmatched_centbin{}'.format(i), 0)
        for j in range(len(hM_Eta_vtxZ_anabin_ratio)):
            hM_Eta_vtxZ_anabin_ratio[j].Divide(hM_Eta_vtxZ_anabin_reco[j])
            hM_Eta_vtxZ_anabin_ratio[j].SetNameTitle('hM_anabin_1minusBeta_centbin{}_Layer{}'.format(i,layer[j]),'hM_anabin_1minusBeta_centbin{}_Layer{}'.format(i,layer[j]))

        Draw_2Dhist(hM_Eta_vtxZ_reco[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_centbin{}_layer12'.format(i))
        Draw_2Dhist(hM_Eta_vtxZ_reco_genmatched[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_genmatched_centbin{}_layer12'.format(i))
        Draw_2Dhist(hM_Eta_vtxZ_anabin_reco[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colztext45', plotpath+'RecoTracklet_Eta_vtxZ_anabin_centbin{}_layer12'.format(i))
        Draw_2Dhist(hM_Eta_vtxZ_anabin_reco_genmatched[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colztext45', plotpath+'RecoTracklet_Eta_vtxZ_anabin_genmatched_centbin{}_layer12'.format(i))
        Draw_2Dhist(hM_Eta_vtxZ_anabin_ratio[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colztext45', plotpath+'RecoTracklet_Eta_vtxZ_anabin_beta_centbin{}_layer12'.format(i))

        hM_anabin_1minusBeta.append(hM_Eta_vtxZ_anabin_ratio)

    # Draw_1DhistsComp2(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname)
    # Draw_1DhistsComp2_wRatio([hM_recotkl_dR_fr.GetValue(),hM_recotkl_genmatched_dR_fr.GetValue()], False, False, True, 200, 'Reco tracklet #DeltaR', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_dR_fr_GenmatchedComp')
    # Draw_1DhistsComp2_wRatio([hM_recotkl_dR.GetValue(),hM_recotkl_genmatched_dR.GetValue()], False, False, True, 200, 'Reco tracklet #DeltaR', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_dR_GenmatchedComp')
    # Draw_1DhistsComp2_wRatio([hM_recotkl_deta.GetValue(),hM_recotkl_genmatched_deta.GetValue()], False, False, True, 200, 'Reco tracklet #Delta#eta', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_deta_GenmatchedComp')
    # Draw_1DhistsComp2_wRatio([hM_recotkl_dphi.GetValue(),hM_recotkl_genmatched_dphi.GetValue()], False, False, True, 200, 'Reco tracklet #Delta#phi', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_dphi_GenmatchedComp')
    Draw_1DhistsComp2_wRatio([hM_recotkl_dR_fr[0], hM_recotkl_genmatched_dR_fr_ghadPtIncl[0]], False, False, True, 200, 'Reco tracklet #DeltaR', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_dR_fr_GenmatchedComp')
    Draw_1DhistsComp2_wRatio([hM_recotkl_dR[0], hM_recotkl_genmatched_dR_ghadPtIncl[0]], False, False, True, 200, 'Reco tracklet #DeltaR', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_dR_GenmatchedComp')
    Draw_1DhistsComp2_wRatio([hM_recotkl_deta[0], hM_recotkl_genmatched_deta_ghadPtIncl[0]], False, False, True, 200, 'Reco tracklet #Delta#eta', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_deta_GenmatchedComp')
    Draw_1DhistsComp2_wRatio([hM_recotkl_dphi[0], hM_recotkl_genmatched_dphi_ghadPtIncl[0]], False, False, True, 200, 'Reco tracklet #Delta#phi', '', ['Before gen-matched','Gen-matched'], plotpath+'RecoTracklet_dphi_GenmatchedComp')

    # hM_recotkl_dR_diff_fr = hM_recotkl_dR_fr.GetValue().Clone()
    # hM_recotkl_dR_diff_fr.Add(hM_recotkl_genmatched_dR_fr.GetValue(), -1)
    # hM_recotkl_dR_diff_fr_cumu = hM_recotkl_dR_diff_fr.GetCumulative();
    # hM_recotkl_dR_diff_fr_cumu.Scale(1./hM_recotkl_dR_fr.GetValue().Integral(-1,-1))
    # hM_recotkl_dR_diff = hM_recotkl_dR.GetValue().Clone()
    # hM_recotkl_dR_diff.Add(hM_recotkl_genmatched_dR.GetValue(), -1)
    # hM_recotkl_dR_diff_cumu = hM_recotkl_dR_diff.GetCumulative();
    # hM_recotkl_dR_diff_cumu.Scale(1./hM_recotkl_dR.GetValue().Integral(-1,-1))
    hM_recotkl_dR_diff_fr = hM_recotkl_dR_fr[0].Clone()
    hM_recotkl_dR_diff_fr.Add(hM_recotkl_genmatched_dR_fr_ghadPtIncl[0], -1)
    hM_recotkl_dR_diff_fr_cumu = hM_recotkl_dR_diff_fr.GetCumulative();
    hM_recotkl_dR_diff_fr_cumu.Scale(1./hM_recotkl_dR_fr[0].Integral(-1,-1))
    hM_recotkl_dR_diff = hM_recotkl_dR[0].Clone()
    hM_recotkl_dR_diff.Add(hM_recotkl_genmatched_dR_ghadPtIncl[0], -1)
    hM_recotkl_dR_diff_cumu = hM_recotkl_dR_diff.GetCumulative();
    hM_recotkl_dR_diff_cumu.Scale(1./hM_recotkl_dR[0].Integral(-1,-1))
    # Draw_1Dhist_cumulative(hist, logy, ymaxscale, XaxisName, YaxisName, Ytitle_unit, outname)
    Draw_1Dhist_cumulative(hM_recotkl_dR_diff_fr_cumu, False, 1.2, 'Reco tracklet #DeltaR', 'Cumulative fraction', '', plotpath+'RecoTracklet_dR_FracNotMatched_Cumulative_FullRange')
    Draw_1Dhist_cumulative(hM_recotkl_dR_diff_cumu, False, 1.2, 'Reco tracklet #DeltaR', 'Cumulative fraction', '', plotpath+'RecoTracklet_dR_FracNotMatched_Cumulative')

    fout = TFile(plotpath+'BkgFraction.root','recreate')
    fout.cd()
    for i in range(20): # Centrality bin
        for j in range(3): # Layer combination
            hM_anabin_1minusBeta[i][j].Write()

    fout.Close()
