#! /usr/bin/env python
from optparse import OptionParser
import sys
import os
import datetime
from array import *
import ROOT
import numpy as np
import pandas as pd
import math
import glob
import ctypes
from plotUtil import Draw_1Dhist

ROOT.gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
ROOT.gROOT.ProcessLine('SetsPhenixStyle()')
ROOT.gROOT.SetBatch(True)

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
    c = ROOT.TCanvas('c', 'c', 800, 700)
    if norm1:
        hist.Scale(1. / hist.Integral(-1, -1))
    if logy:
        c.SetLogy()
    c.cd()
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
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
    # Gaussian fit
    # f1 = ROOT.TF1('f1', 'gaus', hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    # f1.SetParameter(1,0.)
    # f1.SetParLimits(1,-0.1,0.1)
    # f1.SetParameter(2,0.0025)
    # f1.SetParLimits(2,-0.2,0.2)
    # f1.SetLineColorAlpha(ROOT.TColor.GetColor('#F54748'), 0.9)
    # hist.Fit(f1,'B')
    # f1.Draw('same')
    leg = ROOT.TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.1,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    # leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
    leg2 = ROOT.TLegend(0.54, 0.67, 0.88, 0.8)
    leg2.SetTextSize(0.033)
    leg2.SetFillStyle(0)
    # leg2.AddEntry(f1,'Gaussian fit','l')
    # leg2.AddEntry('', '#mu={0:.2e}#pm{1:.2e}'.format(f1.GetParameter(1), f1.GetParError(1)), '')
    # leg2.AddEntry('', '#sigma={0:.2e}#pm{1:.2e}'.format(f1.GetParameter(2), f1.GetParError(2)), '')
    leg2.AddEntry('', '#mu={0:.2e}#pm{1:.2e}'.format(hist.GetMean(), hist.GetMeanError()), '')
    leg2.AddEntry('', '#sigma={0:.2e}#pm{1:.2e}'.format(hist.GetStdDev(), hist.GetStdDevError()), '')
    leg2.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        ROOT.gSystem.ProcessEvents()
        del c
        c = 0
    # return f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(2), f1.GetParError(2)
    return hist.GetMean(), hist.GetMeanError(), hist.GetStdDev(), hist.GetStdDevError()

def Draw_2Dhist(hist, logz, norm1, rmargin, XaxisName, YaxisName, outname):
    c = ROOT.TCanvas('c', 'c', 800, 700)
    if logz:
        c.SetLogz()
    c.cd()
    ROOT.gPad.SetRightMargin(rmargin)
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
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

    leg = ROOT.TLegend(LeftMargin, 1-TopMargin*1.1, LeftMargin+0.01, 0.98)
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
        ROOT.gSystem.ProcessEvents()
        del c
        c = 0

def Draw2Dhisttable(hist, XaxisName, YaxisName, ZaxisName, DrawOpt, outname):
    c = ROOT.TCanvas('c', 'c', 4000, 3700)
    c.cd()
    ROOT.gPad.SetRightMargin(0.18)
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gStyle.SetPaintTextFormat('1.5f')
    hist.SetMarkerSize(0.4)
    hist.GetXaxis().SetTitle(XaxisName)
    hist.GetYaxis().SetTitle(YaxisName)
    hist.GetZaxis().SetTitle(ZaxisName)
    hist.GetXaxis().SetTickSize(TickSize)
    hist.GetYaxis().SetTickSize(TickSize)
    hist.GetXaxis().SetTitleSize(AxisTitleSize)
    hist.GetYaxis().SetTitleSize(AxisTitleSize)
    hist.GetZaxis().SetTitleSize(AxisTitleSize)
    hist.GetXaxis().SetLabelSize(AxisLabelSize)
    hist.GetYaxis().SetLabelSize(AxisLabelSize)
    hist.GetXaxis().SetTitleOffset(1.1)
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetZaxis().SetTitleOffset(1.2)
    hist.GetZaxis().SetLabelSize(AxisLabelSize)
    hist.GetZaxis().SetRangeUser(hist.GetMinimum(), hist.GetMaximum())
    hist.LabelsOption("v")
    hist.SetContour(10000)
    hist.Draw(DrawOpt)
    # bx,by,bz = (ctypes.c_int(),ctypes.c_int(),ctypes.c_int()) # Ref: https://github.com/svenkreiss/decouple/blob/master/Decouple/plot_utils.py#L44-L52
    # hist.GetBinXYZ(hist.GetMinimumBin(),bx,by,bz)
    # binxCenter = hist.GetXaxis().GetBinCenter(bx.value)
    # binyCenter = hist.GetYaxis().GetBinCenter(by.value)
    # el = ROOT.TEllipse(binxCenter,binyCenter,.1,.1)
    # el.SetFillColor(2)
    # el.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        ROOT.gSystem.ProcessEvents()
        del c
        c = 0

def doFitSaveFitresult(reshist, fitresname, histwTF1name, dphicutbin, dzcutbin):
    print ('In doFitSaveFitresult()')
    out = ROOT.TFile(fitresname, 'recreate') 
    tree = ROOT.TTree('tree', 'tree')
    arrmean = array('f', [0])
    arrmeanerr = array('f', [0])
    arrsigma = array('f', [0])
    arrsigmaerr = array('f', [0])
    tree.Branch('mean', arrmean, 'mean/F') 
    tree.Branch('mean_err', arrmeanerr, 'mean_err/F')
    tree.Branch('sigma', arrsigma, 'sigma/F') 
    tree.Branch('sigma_err', arrsigmaerr, 'sigma_err/F')
    arrmean[0], arrmeanerr[0], arrsigma[0], arrsigmaerr[0] = Draw_1Dhist_wTF1(reshist, False, True, 50, '#Deltaz(PV_{Truth}, PV_{Reco}) (cm)', 'cm', histwTF1name)
    tree.Fill()
    tree.Write('', ROOT.TObject.kOverwrite) 
    out.Close()
    del reshist, out, tree, arrmean, arrmeanerr, arrsigma, arrsigmaerr

def main(dphicutbin, dzcutbin):
    print ('In main()')
    # plotpath = './PV_TruthReco/optimization/dPhiCutBin{}_dZCutBin{}/'.format(dphicutbin,dzcutbin)
    plotpath = './RecoPV_optimization/dPhiCutBin{}_dZCutBin{}/'.format(dphicutbin,dzcutbin)
    os.makedirs(plotpath, exist_ok=True)

    hM_DiffVtxZ_Nvtx1 = ROOT.TH1F('hM_DiffVtxZ_Nvtx1', 'hM_DiffVtxZ_Nvtx1_altrange', 80, -0.1, 0.1)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange = ROOT.TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 100, -0.1, 0.1, 100, 0, 5000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2 = ROOT.TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 100, -0.05, 0.05, 100, 0, 2000)
    hM_DiffVtxZ_Nvtx1_NClusL1less100 = ROOT.TH1F('hM_DiffVtxZ_Nvtx1_NClusL1less100', 'hM_DiffVtxZ_Nvtx1_NClusL1less100', 80, -0.2, 0.2)

    df = ROOT.RDataFrame('minitree', '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin{}_dZCutbin{}.root'.format(dphicutbin,dzcutbin))
    data = df.AsNumpy(columns=['NTruthVtx','PV_z','TruthPV_trig_z','NhitsLayer1'])

    for i, ntruthvtx in enumerate(data['NTruthVtx']):
        # if ntruthvtx == 1 and abs(data['TruthPV_trig_z'][i]) < 10:
        if ntruthvtx == 1:
            hM_DiffVtxZ_Nvtx1.Fill((data['PV_z'][i]-data['TruthPV_trig_z'][i]))
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange.Fill((data['PV_z'][i] - data['TruthPV_trig_z'][i]), data['NhitsLayer1'][i])
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2.Fill((data['PV_z'][i] - data['TruthPV_trig_z'][i]), data['NhitsLayer1'][i])
            if data['NhitsLayer1'][i] <= 100:
                hM_DiffVtxZ_Nvtx1_NClusL1less100.Fill((data['PV_z'][i]-data['TruthPV_trig_z'][i]))

    doFitSaveFitresult(hM_DiffVtxZ_Nvtx1, plotpath+'fitresult_dPhiCutbin{}_dZCutbin{}.root'.format(dphicutbin, dzcutbin), plotpath+'RecoClusters_DiffVtxZ_Nvtx1', dphicutbin, dzcutbin)
    doFitSaveFitresult(hM_DiffVtxZ_Nvtx1_NClusL1less100, plotpath+'fitresult_NClusL1less100_dPhiCutbin{}_dZCutbin{}.root'.format(dphicutbin, dzcutbin), plotpath+'RecoClusters_DiffVtxZ_Nvtx1_NClusL1less100', dphicutbin, dzcutbin)
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2')

    del df, data
    
def get2DTable(fsaddname, histaddname):
    print('In get2DTable')
    hM_mean_dPhidZCut = ROOT.TH2F('hM_mean_dPhidZCut{}'.format(histaddname), 'hM_mean_dPhidZCut{}'.format(histaddname), 30, 0, 30, 30, 0, 30)
    hM_sigma_dPhidZCut = ROOT.TH2F('hM_sigma_dPhidZCut{}'.format(histaddname), 'hM_sigma_dPhidZCut{}'.format(histaddname), 30, 0, 30, 30, 0, 30)

    l_sigma = []
    l_mean = []
    l_dPhiBin = []
    l_dZBin = []
    for i in range(0, 30):
        for j in range(0, 30):
            # print (i,j)
            try:
                fitrs = ROOT.RDataFrame('tree', './RecoPV_optimization/dPhiCutBin{}_dZCutBin{}/fitresult{}_dPhiCutbin{}_dZCutbin{}.root'.format(i+1,j+1,fsaddname,i+1,j+1))
                res = fitrs.AsNumpy(columns=['mean','sigma'])

                l_dPhiBin.append(i+1)
                l_dZBin.append(j+1)
                l_mean.append(res['mean'][0])
                l_sigma.append(res['sigma'][0])

                hM_mean_dPhidZCut.SetBinContent(i + 1, j + 1, res['mean'][0])
                hM_sigma_dPhidZCut.SetBinContent(i + 1, j + 1, res['sigma'][0])
                hM_mean_dPhidZCut.GetXaxis().SetBinLabel(i+1, '{:.2f}'.format((i+1)*0.01))
                hM_mean_dPhidZCut.GetYaxis().SetBinLabel(j+1, '{:.2f}'.format((j+1)*0.01))
                hM_sigma_dPhidZCut.GetXaxis().SetBinLabel(i+1, '{:.2f}'.format((i+1)*0.01))
                hM_sigma_dPhidZCut.GetYaxis().SetBinLabel(j+1, '{:.2f}'.format((j+1)*0.01))
                del fitrs, res

                # plot_check = (i+1 == 8 and j+1 == 2)
                plot_check = False
                if plot_check:
                    for evt in range(2000):
                        os.makedirs('./RecoPV_optimization/dPhiCutBin{}_dZCutBin{}/event/rmsAll'.format(i+1,j+1), exist_ok=True)
                        os.makedirs('./RecoPV_optimization/dPhiCutBin{}_dZCutBin{}/event/vtxclusterZ'.format(i+1,j+1), exist_ok=True)
                        f = ROOT.TFile('/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/hists/dPhiCutbin{}_dZCutbin{}/VtxCluster_hist_dPhiCutbin{}_dZCutbin{}_event{}.root'.format(i+1, j+1, i+1, j+1, evt), 'r')
                        hM_rmsAll = f.Get('hM_rmsAll')
                        hM_vtxclusterZ = f.Get('hM_vtxclusterZ')
                        Draw_1Dhist(hM_rmsAll, False, True, 50, 'vertex cluster std (cm)', 'cm', './RecoPV_optimization/dPhiCutBin{}_dZCutBin{}/event/rmsAll/hM_rmsAll_evt{}'.format(i+1,j+1,evt))
                        Draw_1Dhist(hM_vtxclusterZ, False, True, 50, 'vertex candidate Z position (cm)', 'cm', './RecoPV_optimization/dPhiCutBin{}_dZCutBin{}/event/vtxclusterZ/hM_vtxclusterZ_evt{}'.format(i+1,j+1,evt))
                        del f, hM_rmsAll, hM_vtxclusterZ

            except Exception as e:
                print(e)

    df_opt = pd.DataFrame(list(zip(l_dPhiBin, l_dZBin, l_mean, l_sigma)), columns =['dPhiBin', 'dZBin', 'mean', 'sigma'])
    df_opt = df_opt.sort_values(by='sigma', ascending=True, ignore_index=True)        

    del l_dPhiBin, l_dZBin, l_mean, l_sigma
    return hM_mean_dPhidZCut, hM_sigma_dPhidZCut, df_opt

if __name__ == '__main__':
    plotpath = './RecoPV_optimization/'
    dofit = False

    for i in range(0, 30):
        for j in range(0, 30):
            try:
                if dofit == True:
                    main(i+1,j+1)
            except Exception as e:
                print(e)

    hM_mean_dPhidZCut, hM_sigma_dPhidZCut, df_opt = get2DTable('', '')
    # hM_mean_dPhidZCut_NClusL1less100, hM_sigma_dPhidZCut_NClusL1less100, df_opt_NClusL1less100 = get2DTable('_NClusL1less100', '_NClusL1less100')

    # hM_mean_dPhidZCut.SetMinimum(1e-10)
    # hM_sigma_dPhidZCut.SetMinimum(0)
    # hM_sigma_dPhidZCut_NClusL1less100.SetMinimum(0)
    ROOT.gStyle.SetPalette(ROOT.kThermometer)
    # ROOT.TColor.InvertPalette()
    ROOT.TGaxis.SetMaxDigits(2)

    Draw2Dhisttable(hM_mean_dPhidZCut, '#Delta#phi cut [radian]', '|#Delta z| cut [cm]', '#mu_{Gaussian fit} [cm]', 'colztext45', plotpath+'Optimization_GaussianMean_Table')
    Draw2Dhisttable(hM_sigma_dPhidZCut, '#Delta#phi cut [radian]', '|#Delta z| cut [cm]', '#sigma_{Gaussian fit} [cm]', 'colztext45', plotpath+'Optimization_GaussianSigma_Table')
    # Draw2Dhisttable(hM_mean_dPhidZCut_NClusL1less100, '#Delta#phi cut [radian]', '|#Delta z| cut [cm]', 'Gaussian width', 'colztext45', plotpath+'Optimization_NClusL1less100_GaussianMean_Table')
    # Draw2Dhisttable(hM_sigma_dPhidZCut_NClusL1less100, '#Delta#phi cut [radian]', '|#Delta z| cut [cm]', 'Gaussian width' , 'colztext45', plotpath+'Optimization_NClusL1less100_GaussianSigma_Table')

    print ('Inclusive')
    for i in range(5):
        sigmainfo = df_opt.iloc[i]
        print(sigmainfo)

    # print ('Peripheral 1: NClusL1 < 100')
    # for i in range(5):
    #     sigmainfo = df_opt_NClusL1less100.iloc[i]
    #     print(sigmainfo)

    