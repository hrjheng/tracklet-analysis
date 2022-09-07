#! /usr/bin/env python
from optparse import OptionParser
import sys
import os
import datetime
from array import *
from ROOT import *
import numpy as np
import math
import glob
from plotUtil import Draw_1Dhist

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

def Draw_1Dhist_wTF1(hist, norm1, logy, ymaxscale, XaxisName, Ytitle_unit, outname):
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
    f1 = TF1('f1', 'gaus', hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    f1.SetLineColorAlpha(TColor.GetColor('#F54748'), 0.9)
    hist.Fit(f1)
    f1.Draw('same')
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.1,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    # leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
    leg2 = TLegend(0.54, 0.67, 0.88, 0.8)
    leg2.SetTextSize(0.033)
    leg2.SetFillStyle(0)
    leg2.AddEntry(f1,'Gaussian fit','l')
    leg2.AddEntry('', '#mu={0:.2e}#pm{1:.2e}'.format(f1.GetParameter(1), f1.GetParError(1)), '')
    leg2.AddEntry('', '#sigma={0:.2e}#pm{1:.2e}'.format(f1.GetParameter(2), f1.GetParError(2)), '')
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
    return f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(2), f1.GetParError(2)

def Draw_2Dhist(hist, logz, norm1, rmargin, XaxisName, YaxisName, outname):
    c = TCanvas('c', 'c', 800, 700)
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
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetZaxis().SetLabelSize(AxisLabelSize)
    hist.SetContour(1000)
    hist.Draw('colz')

    # l=TLine()
    # l.SetLineStyle(2)
    # l.DrawLine(-25, -24, 25, 26)
    # l.DrawLine(-25, -26, 25, 24)

    leg = TLegend(LeftMargin, 1-TopMargin*1.1, LeftMargin+0.01, 0.98)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    # leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
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

def Draw2Dhisttable(hist, XaxisName, YaxisName, DrawOpt, outname):
    c = TCanvas('c', 'c', 3300, 3200)
    c.cd()
    gPad.SetRightMargin(0.15)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    # gStyle.SetPaintTextFormat('1.3f')
    hist.GetXaxis().SetTitle(XaxisName)
    hist.GetYaxis().SetTitle(YaxisName)
    hist.GetXaxis().SetTickSize(TickSize)
    hist.GetYaxis().SetTickSize(TickSize)
    hist.GetXaxis().SetTitleSize(AxisTitleSize)
    hist.GetYaxis().SetTitleSize(AxisTitleSize)
    hist.GetXaxis().SetLabelSize(AxisLabelSize)
    hist.GetYaxis().SetLabelSize(AxisLabelSize)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetZaxis().SetLabelSize(AxisLabelSize)
    hist.LabelsOption("v")
    hist.SetContour(10000)
    hist.Draw(DrawOpt)
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def main(dphicutbin, dzcutbin):
    plotpath = './PV_TruthReco/optimization/dPhiCutBin{}_dZCutBin{}/'.format(dphicutbin,dzcutbin)
    os.makedirs(plotpath, exist_ok=True)

    hM_DiffVtxZ_Nvtx1 = TH1F('hM_DiffVtxZ_Nvtx1', 'hM_DiffVtxZ_Nvtx1_altrange', 40, -0.1, 0.1)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 100, -0.1, 0.1, 100, 0, 10000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2 = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 50, -0.05, 0.05, 50, 0, 2000)

    df = ROOT.RDataFrame('minitree', '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin{}_dZCutbin{}_Nevt4000.root'.format(dphicutbin,dzcutbin))
    data = df.AsNumpy(columns=['NTruthVtx','PV_z','TruthPV_trig_z','NhitsLayer1'])

    for i, ntruthvtx in enumerate(data['NTruthVtx']):
        if ntruthvtx == 1:
            hM_DiffVtxZ_Nvtx1.Fill((data['PV_z'][i]-data['TruthPV_trig_z'][i]))
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange.Fill((data['PV_z'][i] - data['TruthPV_trig_z'][i]), data['NhitsLayer1'][i])
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2.Fill((data['PV_z'][i] - data['TruthPV_trig_z'][i]), data['NhitsLayer1'][i])
    

    out = TFile(plotpath+'fitresult_dPhiCutbin{}_dZCutbin{}.root'.format(dphicutbin, dzcutbin), 'recreate') 
    tree = TTree('tree', 'tree')
    arrmean = array('f', [0])
    arrmeanerr = array('f', [0])
    arrsigma = array('f', [0])
    arrsigmaerr = array('f', [0])
    # cmd = 'struct tree_str {Float_t mean, mean_err, sigma, sigma_err;};' 
    # gROOT.ProcessLine(cmd)
    # tree_out = tree_str()
    tree.Branch('mean', arrmean, 'mean/F') 
    tree.Branch('mean_err', arrmeanerr, 'mean_err/F')
    tree.Branch('sigma', arrsigma, 'sigma/F') 
    tree.Branch('sigma_err', arrsigmaerr, 'sigma_err/F')
    arrmean[0], arrmeanerr[0], arrsigma[0], arrsigmaerr[0] = Draw_1Dhist_wTF1(hM_DiffVtxZ_Nvtx1, False, True, 50, '#Deltaz(PV_{Truth}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1')
    tree.Fill()
    tree.Write('', TObject.kOverwrite) 
    out.Close()

    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2')

    del out, df, data
    

if __name__ == '__main__':
    plotpath = './PV_TruthReco/optimization/'
    dofit = True

    hM_mean_dPhidZCut = TH2F('hM_mean_dPhidZCut','hM_mean_dPhidZCut',40,0,40,40,0,40)
    hM_sigma_dPhidZCut = TH2F('hM_sigma_dPhidZCut','hM_sigma_dPhidZCut',40,0,40,40,0,40)
    
    for i in range(0, 40):
        for j in range(0, 40):
            if dofit == True:
                main(i+1,j+1)
            
            fitrs = ROOT.RDataFrame('tree', './PV_TruthReco/optimization/dPhiCutBin{}_dZCutBin{}/fitresult_dPhiCutbin{}_dZCutbin{}.root'.format(i+1,j+1,i+1,j+1))
            res = fitrs.AsNumpy(columns=['mean','sigma'])

            hM_mean_dPhidZCut.SetBinContent(i + 1, j + 1, res['mean'][0])
            hM_sigma_dPhidZCut.SetBinContent(i + 1, j + 1, res['sigma'][0])
            hM_mean_dPhidZCut.GetXaxis().SetBinLabel(i+1, '{:.2f}'.format((i+1)*0.01))
            hM_mean_dPhidZCut.GetYaxis().SetBinLabel(j+1, '{:.2f}'.format((j+1)*0.01))
            hM_sigma_dPhidZCut.GetXaxis().SetBinLabel(i+1, '{:.2f}'.format((i+1)*0.01))
            hM_sigma_dPhidZCut.GetYaxis().SetBinLabel(j+1, '{:.2f}'.format((j+1)*0.01))

            del fitrs, res
    
    # hM_mean_dPhidZCut.SetMinimum(1e-10)
    # hM_sigma_dPhidZCut.SetMinimum(1e-6)
    Draw2Dhisttable(hM_mean_dPhidZCut, '#Delta#phi cut', '|#Delta z| cut', 'colz0', plotpath+'Optimization_GaussianMean_Table')
    Draw2Dhisttable(hM_sigma_dPhidZCut, '#Delta#phi cut', '|#Delta z| cut', 'colz0', plotpath+'Optimization_GaussianSigma_Table')


    