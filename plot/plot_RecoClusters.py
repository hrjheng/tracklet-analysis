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

def Draw_1Dhist_fitPossion(hist, norm1, logy, ymaxscale, XaxisName, Ytitle_unit, fitparams, outname):
    hist.Sumw2()
    binwidth = hist.GetXaxis().GetBinWidth(1)
    c = TCanvas('c', 'c', 800, 700)
    if norm1:
        hist.Scale(1. / hist.Integral(-1, -1))
    if logy:
        c.SetLogy()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    if norm1:
        if Ytitle_unit == '':
            hist.GetYaxis().SetTitle(
                'Normalized entries / ({:g})'.format(binwidth))
        else:
            hist.GetYaxis().SetTitle(
                'Normalized entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))
    else:
        if Ytitle_unit == '':
            hist.GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth))
        else:
            hist.GetYaxis().SetTitle(
                'Entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))

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
    f1 = TF1('f1', '[0] * TMath::Poisson(x, [1])', hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    f1.SetParameters(fitparams[0], fitparams[1])
    f1.SetLineColorAlpha(TColor.GetColor('#F54748'), 0.9)
    hist.Fit(f1)
    f1.Draw('same')
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.1,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    # leg.AddEntry("","Au+Au #sqrt{s_{NN}}=200 GeV","")
    leg.Draw()
    leg2 = TLegend(0.54, 0.7, 0.88, 0.8)
    leg2.SetTextSize(0.034)
    leg2.SetFillStyle(0)
    leg2.AddEntry(f1,'Poisson fit','l')
    # leg2.AddEntry('', '#mu={0:.2e}#pm{1:.2e}'.format(f1.GetParameter(0), f1.GetParError(0)), '')
    leg2.AddEntry('', '#mu={0:.2f}#pm{1:.3g}'.format(f1.GetParameter(1), f1.GetParError(1)), '')
    leg2.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_1Dhist_PerEv(hist, truthPV_z, NClusters, norm1, logy, ymaxscale, XaxisName, Ytitle_unit, outname):
    hist.Sumw2()
    binwidth = hist.GetXaxis().GetBinWidth(1)
    c = TCanvas('c', 'c', 800, 700)
    if norm1:
        hist.Scale(1. / hist.Integral(-1, -1))
    if logy:
        c.SetLogy()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    if norm1:
        if Ytitle_unit == '':
            hist.GetYaxis().SetTitle(
                'Normalized entries / ({:g})'.format(binwidth))
        else:
            hist.GetYaxis().SetTitle(
                'Normalized entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))
    else:
        if Ytitle_unit == '':
            hist.GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth))
        else:
            hist.GetYaxis().SetTitle(
                'Entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))

    # hist.GetXaxis().SetRangeUser(hist.GetBinLowEdge(1)-binwidth, hist.GetBinLowEdge(hist.GetNbinsX())+2*binwidth)
    if logy:
        hist.GetYaxis().SetRangeUser(hist.GetMinimum(0)*0.5, (hist.GetMaximum()) * 150)
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
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.1, (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    # leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
    leg2 = TLegend(RightMargin+0.08, (1-TopMargin)-0.12, RightMargin+0.18, (1-TopMargin)-0.05)
    leg2.SetTextSize(0.035)
    leg2.SetFillStyle(0)
    leg2.AddEntry('','Truth PV Z pos = {:.4g} cm'.format(truthPV_z),'')
    leg2.AddEntry('','# of clusters = {}'.format(NClusters),'')
    leg2.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_2Dhist_alt(hist, logz, norm1, rmargin, XaxisName, YaxisName, outname):
    c = TCanvas('c', 'c', 1400, 700)
    if logz:
        c.SetLogz()
    c.cd()
    gPad.SetRightMargin(rmargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    if norm1:
        hist.Scale(1. / hist.Integral(-1, -1, -1, -1))
    hist.GetXaxis().SetTitle(XaxisName)
    hist.GetYaxis().SetTitle(YaxisName)
    hist.GetXaxis().SetTickSize(TickSize)
    hist.GetYaxis().SetTickSize(TickSize)
    hist.GetXaxis().SetTitleSize(AxisTitleSize)
    hist.GetYaxis().SetTitleSize(AxisTitleSize)
    hist.GetXaxis().SetLabelSize(AxisLabelSize)
    hist.GetYaxis().SetLabelSize(AxisLabelSize)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.15)
    hist.GetZaxis().SetLabelSize(AxisLabelSize)
    hist.SetContour(1000)
    hist.Draw('COLZ')

    profile = hist.ProfileX()
    profile.SetMarkerColor(2)
    profile.SetLineColor(2)
    profile.Draw('sameE1')

    leg = TLegend(LeftMargin, 1-TopMargin*1.1, LeftMargin+0.01, 0.98)
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

if __name__ == '__main__':
    plotpath = './MVTXRecoClusters_Basic/'
    os.makedirs(plotpath, exist_ok=True)
    os.makedirs(plotpath+'PerEvent', exist_ok=True)

    hM_NClusLayer1 = TH1F('hM_NClusLayer1', 'hM_NClusLayer1', 130, 0, 13000)
    hM_NClusLayer1_NTruthVtx1 = TH1F('hM_NClusLayer1_NTruthVtx1', 'hM_NClusLayer1_NTruthVtx1', 100, 0, 10000)
    hM_NClusLayer1_NTruthVtxg1 = TH1F('hM_NClusLayer1_NTruthVtxg1', 'hM_NClusLayer1_NTruthVtxg1', 130, 0, 13000)
    # hM_NtruthVtx = TH1F('hM_NtruthVtx', 'hM_NtruthVtx', 11, -0.5, 10.5)
    hM_NtruthVtx = TH1F('hM_NtruthVtx', 'hM_NtruthVtx', 11, 0, 11)
    hM_NTruthVtx_MVTXwindow = TH1F('hM_NTruthVtx_MVTXwindow','hM_NTruthVtx_MVTXwindow', 11, 0, 11)
    hM_NtruthVtx_NClusL1l10k = TH1F('hM_NtruthVtx_NClusL1l10k', 'hM_NtruthVtx_NClusL1l10k', 11, 0, 11)
    hM_NtruthVtx_NClusL1g10k = TH1F('hM_NtruthVtx_NClusL1g10k', 'hM_NtruthVtx_NClusL1g10k', 11, 0, 11)
    hM_NtruthVtx_NClusL1g15k = TH1F('hM_NtruthVtx_NClusL1g15k', 'hM_NtruthVtx_NClusL1g15k', 11, 0, 11)
    hM_NClusLayer1_NTruthVtx = TH2F('hM_NClusLayer1_NTruthVtx','hM_NClusLayer1_NTruthVtx',130, 0, 13000, 11, 0, 11)

    hM_TruthVtxTime = TH1F('hM_TruthVtxTime','hM_TruthVtxTime',60,-15000,15000)
    hM_TruthVtxTime_Embedneq0 = TH1F('hM_TruthVtxTime_Embedneq0', 'hM_TruthVtxTime_Embedneq0', 60,-15000,15000)
    hM_TruthVtxTime_Embedneq0_binw2 = TH1F('hM_TruthVtxTime_Embedneq0_binw2', 'hM_TruthVtxTime_Embedneq0_binw2', 30,-15000,15000)
    hM_TruthVtxTime_altrange1 = TH1F('hM_TruthVtxTime_altrange1','hM_TruthVtxTime_altrange1',60,-0.003,0.003)
    hM_TruthVtxTime_NClusL1g10k = TH1F('hM_TruthVtxTime_NClusL1g10k','hM_TruthVtxTime_NClusL1g10k', 60,-15000,15000)
    hM_TruthVtxTime_NClusL1g10k_altrange1 = TH1F('hM_TruthVtxTime_NClusL1g10k_altrange1','hM_TruthVtxTime_NClusL1g10k_altrange1',60,-0.003,0.003)

    hM_NMvtxClusVtx = TH1F('hM_NMvtxClusVtx', 'hM_NMvtxClusVtx', 100,0,10000)
    hM_NMvtxClusVtx_VtxTime = TH2F('hM_NMvtxClusVtx_VtxTime','hM_NMvtxClusVtx_VtxTime', 150, 0, 15000, 150, -15000, 15000)
    hM_NMvtxClusVtx_VtxTime_altrange = TH2F('hM_NMvtxClusVtx_VtxTime_altrange','hM_NMvtxClusVtx_VtxTime_altrange', 100, 0, 1000, 100, -10000, 10000)
    hM_NMvtxClusVtx_VtxTime_altrange2 = TH2F('hM_NMvtxClusVtx_VtxTime_altrange2','hM_NMvtxClusVtx_VtxTime_altrange2', 100, 0, 100, 100, -10000, 10000)
    hM_NMvtxClusVtx_VtxTime_Embedneq0 = TH2F('hM_NMvtxClusVtx_VtxTime_Embedneq0','hM_NMvtxClusVtx_VtxTime_Embedneq0',100, 0, 1000, 100, -10000, 10000)
    hM_NMvtxClusVtx_NMvtxHitVtx = TH2F('hM_NMvtxClusVtx_NMvtxHitVtx','hM_NMvtxClusVtx_NMvtxHitVtx', 120, 0, 12000, 120, 0, 12000)
    hM_NMvtxClusVtx_NMvtxHitVtx_Embedneq0 = TH2F('hM_NMvtxClusVtx_NMvtxHitVtx_Embedneq0','hM_NMvtxClusVtx_NMvtxHitVtx_Embedneq0', 120, 0, 12000, 120, 0, 12000)
    hM_NMvtxClusVtx_NMvtxHitVtx_Embedeq0 = TH2F('hM_NMvtxClusVtx_NMvtxHitVtx_Embedeq0','hM_NMvtxClusVtx_NMvtxHitVtx_Embedeq0', 120, 0, 12000, 120, 0, 12000)
    hM_NMvtxClusVtx_NG4ParticleVtx = TH2F('hM_NMvtxClusVtx_NG4ParticleVtx','hM_NMvtxClusVtx_NG4ParticleVtx',120, 0, 12000, 120, 0, 12000)
    hM_NMvtxClusVtx_NG4ParticleVtx_Embedneq0 = TH2F('hM_NMvtxClusVtx_NG4ParticleVtx_Embedneq0','hM_NMvtxClusVtx_NG4ParticleVtx_Embedneq0',120, 0, 12000, 120, 0, 12000)
    hM_NMvtxClusVtx_NG4ParticleVtx_Embedeq0 = TH2F('hM_NMvtxClusVtx_NG4ParticleVtx_Embedeq0','hM_NMvtxClusVtx_NG4ParticleVtx_Embedeq0',120, 0, 12000, 120, 0, 12000)

    f = TFile('/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/MVTXRecoClusters_Nevt4000.root', 'r')
    tree = f.Get('EventTree')
    for idx in range(tree.GetEntries()):
    # for idx in range(5):
        tree.GetEntry(idx)
        # hM_ClusZ_all_ev.Reset('ICESM')

        NLayer1 = 0
        for i in range(len(tree.ClusLayer)):
            if tree.ClusLayer[i] == 0:
                NLayer1 += 1
        
        hM_NtruthVtx.Fill(tree.NTruthVtx)
        if NLayer1 < 10000:
            hM_NtruthVtx_NClusL1l10k.Fill(tree.NTruthVtx)
        if NLayer1 > 10000:
            # print('ev={}, total number of clusters={}, number of 1st layer clusters={}'.format(idx,len(tree.ClusLayer),NLayer1))
            hM_NtruthVtx_NClusL1g10k.Fill(tree.NTruthVtx)
        if NLayer1 > 15000:
            hM_NtruthVtx_NClusL1g15k.Fill(tree.NTruthVtx)

        NTruthVtx_MVTXwindow = 0
        for t in tree.TruthPV_t:
            if abs(t) < 5000.:
                NTruthVtx_MVTXwindow += 1

        hM_NTruthVtx_MVTXwindow.Fill(NTruthVtx_MVTXwindow)

        hM_NClusLayer1.Fill(NLayer1)
        hM_NClusLayer1_NTruthVtx.Fill(NLayer1, tree.NTruthVtx)
        if tree.NTruthVtx == 1:
            hM_NClusLayer1_NTruthVtx1.Fill(NLayer1)
        else:
            hM_NClusLayer1_NTruthVtxg1.Fill(NLayer1)

        for i in range(len(tree.TruthPV_t)):
            hM_NMvtxClusVtx.Fill(tree.TruthPV_NClus[i])
            hM_NMvtxClusVtx_VtxTime.Fill(tree.TruthPV_NClus[i], tree.TruthPV_t[i])
            hM_NMvtxClusVtx_VtxTime_altrange.Fill(tree.TruthPV_NClus[i], tree.TruthPV_t[i])
            hM_NMvtxClusVtx_VtxTime_altrange2.Fill(tree.TruthPV_NClus[i], tree.TruthPV_t[i])
            hM_NMvtxClusVtx_NMvtxHitVtx.Fill(tree.TruthPV_NClus[i], tree.TruthPV_Nhits[i])
            hM_NMvtxClusVtx_NG4ParticleVtx.Fill(tree.TruthPV_NClus[i], tree.TruthPV_Npart[i])
            hM_TruthVtxTime.Fill(tree.TruthPV_t[i])
            hM_TruthVtxTime_altrange1.Fill(tree.TruthPV_t[i])
            if tree.TruthPV_embed[i] != 0:
                # print(idx, i, tree.TruthPV_t[i], tree.TruthPV_Npart[i], tree.TruthPV_Nhits[i], tree.TruthPV_NClus[i])
                hM_TruthVtxTime_Embedneq0.Fill(tree.TruthPV_t[i])
                hM_TruthVtxTime_Embedneq0_binw2.Fill(tree.TruthPV_t[i])
                hM_NMvtxClusVtx_VtxTime_Embedneq0.Fill(tree.TruthPV_NClus[i], tree.TruthPV_t[i])
                hM_NMvtxClusVtx_NMvtxHitVtx_Embedneq0.Fill(tree.TruthPV_NClus[i], tree.TruthPV_Nhits[i])
                hM_NMvtxClusVtx_NG4ParticleVtx_Embedneq0.Fill(tree.TruthPV_NClus[i], tree.TruthPV_Npart[i])
            else:
                hM_NMvtxClusVtx_NMvtxHitVtx_Embedeq0.Fill(tree.TruthPV_NClus[i], tree.TruthPV_Nhits[i])
                hM_NMvtxClusVtx_NG4ParticleVtx_Embedeq0.Fill(tree.TruthPV_NClus[i], tree.TruthPV_Npart[i])

            if NLayer1 > 10000:
                hM_TruthVtxTime_NClusL1g10k.Fill(tree.TruthPV_t[i])
                hM_TruthVtxTime_NClusL1g10k_altrange1.Fill(tree.TruthPV_t[i])
            
            

    hM_NClusLayer1.GetXaxis().SetMaxDigits(3)
    hM_NClusLayer1_NTruthVtx1.GetXaxis().SetMaxDigits(3)
    hM_NClusLayer1_NTruthVtxg1.GetXaxis().SetMaxDigits(3)
    hM_TruthVtxTime_altrange1.GetXaxis().SetMaxDigits(3)
    hM_TruthVtxTime_NClusL1g10k.GetXaxis().SetMaxDigits(3)
    hM_TruthVtxTime_NClusL1g10k_altrange1.GetXaxis().SetMaxDigits(3)
    # hM_NClusLayer1_NTruthVtx.GetXaxis().SetMaxDigits(3)

    Draw_1Dhist(hM_NClusLayer1, False, True, 150, 'Number of clusters (1st layer)', '', plotpath+'NClusLayer1')
    Draw_1Dhist(hM_NClusLayer1_NTruthVtx1, False, True, 150, 'Number of clusters (1st layer)', '', plotpath+'NClusLayer1_NTruthVtx1')
    Draw_1Dhist(hM_NClusLayer1_NTruthVtxg1, False, True, 150, 'Number of clusters (1st layer)', '', plotpath+'NClusLayer1_NTruthVtxg1')

    Draw_1Dhist_fitPossion(hM_NtruthVtx, False, False, 1.3, 'Number of truth vertices', '', [1200,1], plotpath+'NTruthVtx_Inclusive')
    Draw_1Dhist_fitPossion(hM_NTruthVtx_MVTXwindow, False, False, 1.5, 'Number of truth vertices (|t|<5#mus)', '', [2500,1], plotpath+'NTruthVtx_MVTXwindow')

    Draw_1Dhist(hM_NtruthVtx_NClusL1l10k, False, False, 1.3, 'Number of truth vertices', '', plotpath+'NTruthVtx_NClusL1l10k')
    Draw_1Dhist(hM_NtruthVtx_NClusL1g10k, False, False, 1.3, 'Number of truth vertices', '', plotpath+'NTruthVtx_NClusL1g10k')
    Draw_1Dhist(hM_NtruthVtx_NClusL1g15k, False, False, 1.3, 'Number of truth vertices', '', plotpath+'NTruthVtx_NClusL1g15k')

    Draw_1Dhist(hM_TruthVtxTime, False, True, 100, 'Truth vertex time (ns)', 'ns', plotpath+'TruthVtxTime')
    Draw_1Dhist(hM_TruthVtxTime_altrange1, False, True, 100, 'Truth vertex time (ns)', 'ns', plotpath+'TruthVtxTime_altrange1')
    Draw_1Dhist(hM_TruthVtxTime_Embedneq0, False, False, 1.3, 'Truth vertex time (ns)', 'ns', plotpath+'TruthVtxTime_Embedneq0')
    Draw_1Dhist(hM_TruthVtxTime_Embedneq0_binw2, False, False, 1.3, 'Truth vertex time (ns)', 'ns', plotpath+'TruthVtxTime_Embedneq0_binw2')
    Draw_1Dhist(hM_TruthVtxTime_NClusL1g10k, False, True, 100, 'Truth vertex time (ns)', 'ns', plotpath+'TruthVtxTime_NClusL1g10k')
    Draw_1Dhist(hM_TruthVtxTime_NClusL1g10k_altrange1, False, True, 100, 'Truth vertex time (ns)', 'ns', plotpath+'TruthVtxTime_NClusL1g10k_altrange1')

    Draw_2Dhist_alt(hM_NClusLayer1_NTruthVtx, False, False, 0.13, 'Number of clusters (1st layer)', 'Number of truth vertices', plotpath+'NClusLayer1_NTruthVtx')

    # Draw_2Dhist(hist, logz, norm1, rmargin, XaxisName, YaxisName, outname)
    Draw_1Dhist(hM_NMvtxClusVtx, False, False, 1.3, 'Number of MVTX clusters from the vtx', '', plotpath+'NMvtxClusVtx')
    # hM_NMvtxClusVtx_VtxTime.GetYaxis().SetMaxDigits(3)
    Draw_2Dhist(hM_NMvtxClusVtx_VtxTime, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Truth vtx time (ns)', plotpath+'NMvtxClusVtx_VtxTime')
    Draw_2Dhist(hM_NMvtxClusVtx_VtxTime_altrange, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Truth vtx time (ns)', plotpath+'NMvtxClusVtx_VtxTime_altrange')
    Draw_2Dhist(hM_NMvtxClusVtx_VtxTime_altrange2, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Truth vtx time (ns)', plotpath+'NMvtxClusVtx_VtxTime_altrange2')
    Draw_2Dhist(hM_NMvtxClusVtx_VtxTime_Embedneq0, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Truth vtx time (ns)', plotpath+'NMvtxClusVtx_VtxTime_Embedneq0')

    Draw_2Dhist(hM_NMvtxClusVtx_NMvtxHitVtx, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Number of MVTX truth hits from the vtx', plotpath+'NMvtxClusVtx_NMvtxHitVtx')
    Draw_2Dhist(hM_NMvtxClusVtx_NG4ParticleVtx, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Number of MVTX truth particles from the vtx', plotpath+'NMvtxClusVtx_NG4ParticleVtx')
    Draw_2Dhist(hM_NMvtxClusVtx_NG4ParticleVtx_Embedneq0, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Number of MVTX truth particles from the vtx', plotpath+'NMvtxClusVtx_NG4ParticleVtx_Embedneq0')
    Draw_2Dhist(hM_NMvtxClusVtx_NG4ParticleVtx_Embedeq0, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Number of MVTX truth particles from the vtx', plotpath+'NMvtxClusVtx_NG4ParticleVtx_Embedeq0')
    Draw_2Dhist(hM_NMvtxClusVtx_NMvtxHitVtx_Embedneq0, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Number of MVTX truth hits from the vtx', plotpath+'NMvtxClusVtx_NMvtxHitVtx_Embedneq0')
    Draw_2Dhist(hM_NMvtxClusVtx_NMvtxHitVtx_Embedeq0, False, False, 0.13, 'Number of MVTX clusters from the vtx', 'Number of MVTX truth hits from the vtx', plotpath+'NMvtxClusVtx_NMvtxHitVtx_Embedeq0')
