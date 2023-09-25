import os
import datetime
from array import *
from ROOT import *
import numpy as np
import math
import glob
from plotUtil import *
import itertools

gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)

NCentBins = 20

def centbin(fname='/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root', NpercentileDiv=20):
    df = ROOT.RDataFrame('minitree', fname)
    np_NhitsL1 = df.AsNumpy(columns=['NhitsLayer1'])
    Binedge_NhitsL1_percentile = [0]
    for i in range(NpercentileDiv-1):
        Binedge_NhitsL1_percentile.append(np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5))

    Binedge_NhitsL1_percentile.append(10000)
    return Binedge_NhitsL1_percentile


def hists_PV(vtxfname):
    hM_RecoPV = TH1F('hM_RecoPV', 'hM_RecoPV', 100, -25, 25)
    hM_PVdiff = TH1F('hM_PVdiff', 'hM_PVdiff', 100, -0.1, 0.1)
    CentBin = centbin()
    hM_NClusLayer1_NTurhtVtx1 = TH1F('hM_NClusLayer1_NTurhtVtx1','hM_NClusLayer1_NTurhtVtx1', NCentBins, np.asarray(CentBin))
    hM_NClusLayer1_NTruthVtx1_HasRecoPV = TH1F('hM_NClusLayer1_NTruthVtx1_HasRecoPV','hM_NClusLayer1_NTruthVtx1_HasRecoPV', NCentBins, np.asarray(CentBin))

    f = TFile(vtxfname, 'READ')
    tree = f.Get('minitree')
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)
        hM_RecoPV.Fill(tree.PV_z)
        hM_PVdiff.Fill(tree.PV_z - tree.TruthPV_trig_z)

        if tree.NTruthVtx == 1:
            hM_NClusLayer1_NTurhtVtx1.Fill(tree.NhitsLayer1)
            if tree.PV_z > -99.:
                hM_NClusLayer1_NTruthVtx1_HasRecoPV.Fill(tree.NhitsLayer1)

    return hM_RecoPV, hM_PVdiff, hM_NClusLayer1_NTurhtVtx1, hM_NClusLayer1_NTruthVtx1_HasRecoPV


def hists_PV_twohalves(vtxfname):
    hM_RecoPV = TH1F('hM_RecoPV', 'hM_RecoPV', 100, -25, 25)
    hM_PVdiff = TH1F('hM_PVdiff', 'hM_PVdiff', 100, -0.1, 0.1)
    CentBin = centbin()
    hM_NClusLayer1_NTurhtVtx1 = TH1F('hM_NClusLayer1_NTurhtVtx1','hM_NClusLayer1_NTurhtVtx1', NCentBins, np.asarray(CentBin))
    hM_NClusLayer1_NTruthVtx1_HasRecoPV = TH1F('hM_NClusLayer1_NTruthVtx1_HasRecoPV','hM_NClusLayer1_NTruthVtx1_HasRecoPV', NCentBins, np.asarray(CentBin))

    f = TFile(vtxfname, 'READ')
    tree = f.Get('minitree')
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)
        finalPVz = (tree.PV_z_halves1+tree.PV_z_halves2) / 2.
        hM_RecoPV.Fill(finalPVz)
        hM_PVdiff.Fill(finalPVz - tree.TruthPV_trig_z)

        if tree.NTruthVtx == 1:
            hM_NClusLayer1_NTurhtVtx1.Fill(tree.NhitsLayer1)
            if finalPVz > -99.:
                hM_NClusLayer1_NTruthVtx1_HasRecoPV.Fill(tree.NhitsLayer1)

    return hM_RecoPV, hM_PVdiff, hM_NClusLayer1_NTurhtVtx1, hM_NClusLayer1_NTruthVtx1_HasRecoPV


def Draw_1Dhists_MisalignComp(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, doGaussFit, legtext, outname):
    color = []
    legtextsize = 0.035
    legylowedge = 0.21
    if len(lhist) == 2:
        legtextsize = 0.035
        color = ['#1B1A17','#9B0000']
    elif len(lhist) == 3:
        legtextsize = 0.035
        legylowedge = 0.26
        color = ['#1B1A17','#9B0000','#035397']
    elif len(lhist) == 4:
        legtextsize = 0.035
        legylowedge = 0.3
        color = ['#1B1A17','#9B0000','#035397','#377D71']
    elif len(lhist) == 5:
        legtextsize = 0.03
        legylowedge = 0.35
        color = ['#1B1A17','#9B0000','#035397','#377D71','#FF9F29']
    elif len(lhist) == 6:
        legtextsize = 0.03
        legylowedge = 0.4
        color = ['#1B1A17','#9B0000','#035397','#377D71','#FF9F29','#AA96DA']
    else:
        print('[WARNING] There are {} histograms, but there are not so many colors defined yet'.format(len(lhist)))

    # legtext = ['1st+2nd layers', '2nd+3rd layers', '1st+3rd layers']
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
        lhist[0].GetYaxis().SetRangeUser(ymin * 0.05, ymax * 100)
    else:
        lhist[0].GetYaxis().SetRangeUser(0., ymax * ymaxscale)
    lhist[0].GetXaxis().SetTitle(XaxisName)
    lhist[0].GetXaxis().SetTickSize(TickSize)
    lhist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    lhist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    lhist[0].GetYaxis().SetTickSize(TickSize)
    lhist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    lhist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    lhist[0].GetXaxis().SetTitleOffset(1.1)
    lhist[0].GetYaxis().SetTitleOffset(1.4)

    fitresult = []
    for i, h in enumerate(lhist):
        if i == 0:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('hist')
        else:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('histsame')
        
        if (doGaussFit):
            f = TF1('f_{}'.format(i), 'gaus', h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
            f.SetParameter(1,0.)
            f.SetParLimits(1,-0.1,0.1)
            f.SetParameter(2,0.005)
            f.SetParLimits(2,0,0.2)
            f.SetLineColorAlpha(TColor.GetColor(color[i]), 0.9)
            h.Fit(f, 'B')
            f.Draw('same')
            fitresult.append(f)

    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.15,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    leg.AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "")
    leg.Draw()

    leg1 = TLegend(LeftMargin+0.03, (1-TopMargin)-legylowedge,
                   LeftMargin+0.33, (1-TopMargin)-0.03)
    leg1.SetTextSize()
    leg1.SetFillStyle(0)
    for i, h in enumerate(lhist):
        leg1.AddEntry(h, legtext[i], "l")
        if (doGaussFit):
            leg1.AddEntry(0, "#sigma_{Gauss}" +" = {:.2e} cm".format(fitresult[i].GetParameter(2)), "")
    leg1.SetTextSize(legtextsize)
    leg1.SetFillStyle(0)
    leg1.Draw()
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_1DEffComp(leff, lcolor, logx, XaxisName, legtext, axesrange, outname):
    c = TCanvas('c', 'c', 800, 700)
    if logx:
        c.SetLogx()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(0.14)
    for i,eff in enumerate(leff):
        if i == 0:
            leff[i].GetXaxis().SetTitle(XaxisName)
            leff[i].GetYaxis().SetTitle('Efficiency')
            leff[i].GetXaxis().SetRangeUser(axesrange[0], axesrange[1])
            leff[i].GetYaxis().SetRangeUser(axesrange[2], axesrange[3])
            leff[i].SetLineColor(TColor.GetColor(lcolor[i]))
            leff[i].SetMarkerColor(TColor.GetColor(lcolor[i]))
            leff[i].SetMarkerSize(1.3)
            leff[i].SetMarkerStyle(20)
            leff[i].Draw('AP')
        else:
            leff[i].SetLineColor(TColor.GetColor(lcolor[i]))
            leff[i].SetMarkerColor(TColor.GetColor(lcolor[i]))
            leff[i].SetMarkerSize(1.3)
            leff[i].SetMarkerStyle(20)
            leff[i].Draw('Psame')

    leg = TLegend(LeftMargin, 1-TopMargin*1.1, LeftMargin+0.01, 0.98)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    # leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
    leg1 = TLegend(0.65, 0.23, 0.85, 0.43)
    leg1.SetTextSize(0.035)
    leg1.SetFillStyle(0)
    for i, eff in enumerate(leff):
        leg1.AddEntry(eff, legtext[i], "pe")
    leg1.Draw()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0


if __name__ == '__main__':
    plotpath = './TrackletAna_MisAlignComp/'
    os.makedirs(plotpath, exist_ok=True)

    hM_RecoPV_nominal, hM_PVdiff_nominal, hM_NClusLayer1_NTurhtVtx1_nominal, hM_NClusLayer1_NTruthVtx1_HasRecoPV_nominal = hists_PV('/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root')
    hM_RecoPV_misalign, hM_PVdiff_misalign, hM_NClusLayer1_NTurhtVtx1_misalign, hM_NClusLayer1_NTruthVtx1_HasRecoPV_misalign = hists_PV_twohalves('/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth2p9_3DVertex_twohalves.root')
    hM_RecoPV_misalign_extreme, hM_PVdiff_misalign_extreme, hM_NClusLayer1_NTurhtVtx1_misalign_extreme, hM_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme = hists_PV_twohalves('/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth3p5_3DVertex_twohalves.root')
    hM_RecoPV_misalign_extreme2, hM_PVdiff_misalign_extreme2, hM_NClusLayer1_NTurhtVtx1_misalign_extreme2, hM_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme2 = hists_PV_twohalves('/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth7_3DVertex_twohalves.root')

    Draw_1Dhists_MisalignComp(lhist=[hM_RecoPV_nominal,hM_RecoPV_misalign, hM_RecoPV_misalign_extreme], norm1=False, logx=False, logy=False, ymaxscale=1.5, XaxisName='Primary vertex V_{z} (cm)', Ytitle_unit='cm', doGaussFit=False, legtext=['Nominal', 'Misaligned', 'Misaligned (extreme)'], outname=plotpath+'MisalignComp_twohalfangle_RecoPV')
    Draw_1Dhists_MisalignComp(lhist=[hM_PVdiff_nominal,hM_PVdiff_misalign, hM_PVdiff_misalign_extreme], norm1=False, logx=False, logy=True, ymaxscale=100, XaxisName='#Deltaz(PV_{Truth}, PV_{Reco}) (cm)', Ytitle_unit='cm', doGaussFit=True, legtext=['Nominal', 'Misaligned', 'Misaligned (extreme)'], outname=plotpath+'MisalignComp_twohalfangle_PVdiff')

    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_nominal = TGraphAsymmErrors()
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign = TGraphAsymmErrors()
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme = TGraphAsymmErrors()
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme2 = TGraphAsymmErrors()
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_nominal.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_nominal, hM_NClusLayer1_NTurhtVtx1_nominal)
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_misalign, hM_NClusLayer1_NTurhtVtx1_misalign)
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme, hM_NClusLayer1_NTurhtVtx1_misalign_extreme)
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme2.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme2, hM_NClusLayer1_NTurhtVtx1_misalign_extreme2)

    print('Nominal:', err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_nominal.Print())
    print('Misalign:', err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign.Print())
    print('Misalign (extreme):', err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme.Print())
    print('Misalign (extreme2):', err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme2.Print())

    Draw_1DEffComp(leff=[err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_nominal, err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign, err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme, err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_misalign_extreme2], 
                   lcolor=['#03001C', '#035397', '#9B0000', '#377D71'], 
                   logx=True, XaxisName='Number of clusters (1st layer)', 
                   legtext=['Nominal', 'Gap=2.9mm', 'Gap=3.5mm', 'Gap=7mm'], 
                   axesrange=[0,10100,0.65,1.1], 
                   outname=plotpath+'RecoPVEff_NClusLayer1_EffMisalignComp')
