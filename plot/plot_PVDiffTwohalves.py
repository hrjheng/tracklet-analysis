import os
import datetime
from array import *
from ROOT import *
import numpy as np
import math
import glob
from plotUtil import *
import itertools
import statistics

gROOT.LoadMacro("./sPHENIXStyle/sPhenixStyle.C")
gROOT.ProcessLine("SetsPhenixStyle()")
gROOT.SetBatch(True)

def sigmaEff(v, threshold, xmin, xmax):
    v = np.sort(v)

    total = np.size(v)
    max = threshold * total

    start = []
    stop = []
    width = []

    i = 0
    while i != np.size(v) - 1:
        count = 0
        j = i
        while j != np.size(v) - 1 and count < max:
            count += 1
            j += 1

        if j != np.size(v) - 1:
            start.append(v[i])
            stop.append(v[j])
            width.append(v[j] - v[i])

        i += 1

    npminwidth = np.array(width)

    minwidth = np.amin(npminwidth)
    pos = np.argmin(npminwidth)

    xmin += [start[pos]]
    xmax += [stop[pos]]

    return minwidth


def getPVzDiffSigma(fname):
    list_diffPVz = []

    f = TFile(fname, "READ")
    tree = f.Get("minitree")
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)
        list_diffPVz.append(tree.PV_z_halves1 - tree.PV_z_halves2)
    
    f.Close()
    del f

    data = np.array(list_diffPVz)
    datamin = [-999.]
    datamax = [999.]
    sigmaeff = sigmaEff(data, 0.685, datamin, datamax)

    med = statistics.median(list_diffPVz)

    return datamin, datamax, med, sigmaeff

def getPVDiffHists(fname, DiffPVx_binparam, DiffPVy_binparam, DiffPVz_binparam):
    hM_RecoPVz_DiffPVx = TH2F("hM_RecoPVz_DiffPVx", "hM_RecoPVz_DiffPVx", 100, -25, 25, DiffPVx_binparam[0], DiffPVx_binparam[1], DiffPVx_binparam[2])
    hM_RecoPVz_DiffPVy = TH2F("hM_RecoPVz_DiffPVy", "hM_RecoPVz_DiffPVy", 100, -25, 25, DiffPVy_binparam[0], DiffPVy_binparam[1], DiffPVy_binparam[2])
    hM_RecoPVz_DiffPVz = TH2F("hM_RecoPVz_DiffPVz", "hM_RecoPVz_DiffPVz", 100, -25, 25, DiffPVz_binparam[0], DiffPVz_binparam[1], DiffPVz_binparam[2])
    hM_RecoPVx_Res = TH1F("hM_RecoPVx_Res", "hM_RecoPVx_Res", 100, -0.5, 0.5)
    hM_RecoPVy_Res = TH1F("hM_RecoPVy_Res", "hM_RecoPVy_Res", 100, -0.5, 0.5)
    hM_RecoPVz_Res = TH1F("hM_RecoPVz_Res", "hM_RecoPVz_Res", 100, -0.5, 0.5)
    # hM_RecoPVx_Res = TH1F("hM_RecoPVx_Res", "hM_RecoPVx_Res", 100, -5, 5)
    # hM_RecoPVy_Res = TH1F("hM_RecoPVy_Res", "hM_RecoPVy_Res", 100, -5, 5)
    # hM_RecoPVz_Res = TH1F("hM_RecoPVz_Res", "hM_RecoPVz_Res", 100, -5, 5)
    
    f = TFile(fname, "READ")
    tree = f.Get("minitree")
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)
        if tree.NTruthVtx != 1:
            continue

        finalPVx = (tree.PV_x_halves1 + tree.PV_x_halves2) / 2.
        finalPVy = (tree.PV_y_halves1 + tree.PV_y_halves2) / 2.
        finalPVz = (tree.PV_z_halves1 + tree.PV_z_halves2) / 2.

        hM_RecoPVz_DiffPVx.Fill(finalPVz, tree.PV_x_halves1 - tree.PV_x_halves2)
        hM_RecoPVz_DiffPVy.Fill(finalPVz, tree.PV_y_halves1 - tree.PV_y_halves2)
        hM_RecoPVz_DiffPVz.Fill(finalPVz, tree.PV_z_halves1 - tree.PV_z_halves2)
        # hM_RecoPVx_Res.Fill(finalPVx - tree.TruthPV_trig_x)
        # hM_RecoPVy_Res.Fill(finalPVy - tree.TruthPV_trig_y)
        # hM_RecoPVz_Res.Fill(finalPVz - tree.TruthPV_trig_z)
        hM_RecoPVx_Res.Fill(finalPVx - tree.TruthPV_mostNpart_x)
        hM_RecoPVy_Res.Fill(finalPVy - tree.TruthPV_mostNpart_y)
        hM_RecoPVz_Res.Fill(finalPVz - tree.TruthPV_mostNpart_z)


    f.Close()
    del f

    return hM_RecoPVz_DiffPVx, hM_RecoPVz_DiffPVy, hM_RecoPVz_DiffPVz, hM_RecoPVx_Res, hM_RecoPVy_Res, hM_RecoPVz_Res


def Draw_Graph_wFit(xarray, yarray, lmargin, xaxistitle, yaxistitle, ytitleoffset, xlim, ylim, dofit, outname):
    c = TCanvas('c1', 'c1', 800, 700)
    c.cd()
    gPad.SetLeftMargin(lmargin)
    gr = TGraph(len(xarray), xarray, yarray)
    # gr.SetMarkerStyle(21)
    gr.GetXaxis().SetTitle(xaxistitle)
    gr.GetYaxis().SetTitle(yaxistitle)
    gr.GetXaxis().SetLimits(xlim[0], xlim[1]);
    gr.SetMinimum(ylim[0])
    gr.SetMaximum(ylim[1])
    gr.Draw('AP')
    gr.GetYaxis().SetTitleOffset(ytitleoffset)
    if dofit:
        fit = TF1('fit', 'pol5', xarray[0], xarray[-1])
        fit.SetLineColorAlpha(TColor.GetColor("#F54748"), 0.8)
        gr.Fit(fit,'R')
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_2Dhist_wProfile(hist, logz, norm1, rmargin, XaxisName, YaxisName, drawopt, outname):
    gStyle.SetPalette(kBird)
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
    hist.Draw(drawopt)
    profile = hist.ProfileX()
    profile.SetLineColor(TColor.GetColor("#F54748"))
    profile.SetMarkerColor(TColor.GetColor("#F54748"))
    profile.SetMarkerStyle(20)
    profile.SetMarkerSize(0.5)
    profile.SetLineWidth(2)
    profile.Draw('same')

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


def Draw2Dhisttable(hist, XaxisName, YaxisName, DrawOpt, outname):
    gStyle.SetPalette(kBird)
    c = TCanvas('c', 'c', 4000, 3700)
    c.cd()
    gPad.SetRightMargin(0.18)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    # ROOT.gStyle.SetPaintTextFormat('1.5f')
    # hist.SetMarkerSize(0.4)
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
    hist.GetZaxis().SetRangeUser(hist.GetMinimum(), hist.GetMaximum())
    # hist.LabelsOption("v")
    hist.SetContour(10000)
    hist.Draw(DrawOpt)
    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        ROOT.gSystem.ProcessEvents()
        del c
        c = 0

# draw a list of tprofiles 
def DrawTProfileList(tprofilelist, xaxistitle, yaxistitle, ylim, legpos, legtitle, legtext, outname):
# def DrawTProfileList(tprofilelist, xaxistitle, yaxistitle, ylim, legtitle, legtext, outname):
    # gStyle.SetPalette(kThermometer)
    gStyle.SetPalette(kLightTemperature)
    c = TCanvas('c', 'c', 800, 700)
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    tprofilelist[0].GetXaxis().SetTitle(xaxistitle)
    tprofilelist[0].GetYaxis().SetTitle(yaxistitle)
    tprofilelist[0].GetXaxis().SetTickSize(TickSize)
    tprofilelist[0].GetYaxis().SetTickSize(TickSize)
    tprofilelist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    tprofilelist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    tprofilelist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    tprofilelist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    tprofilelist[0].GetXaxis().SetTitleOffset(1.1)
    tprofilelist[0].GetYaxis().SetRangeUser(ylim[0], ylim[1])
    for i in range(len(tprofilelist)):
        tprofilelist[i].SetMarkerSize(0.8)
        if i == 0:
            tprofilelist[i].Draw('PLCPMC')
        else:
            tprofilelist[i].Draw('SAME PLCPMC')

    leg = TLegend(legpos[0], legpos[1], legpos[2], legpos[3])
    leg.SetHeader(legtitle, 'l')
    leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    for i in range(len(tprofilelist)):
        leg.AddEntry(tprofilelist[i], legtext[i], 'lep')
    leg.Draw()

    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')

    # canv_leg = TCanvas('canv_leg', 'canv_leg', 500, 250)
    # canv_leg.cd()
    # leg = TLegend(0.05, 0.05, 0.95, 0.95)
    # leg.SetHeader(legtitle,'l')
    # leg.SetNColumns(5)
    # leg.SetFillStyle(0)
    # leg.SetLineColor(1)
    # leg.SetLineWidth(2)
    # leg.SetTextSize(0.04)
    # for i in range(len(tprofilelist)):
    #     leg.AddEntry(tprofilelist[i], legtext[i], 'lep')
    # leg.Draw()
    # canv_leg.SaveAs(outname+'_leg.pdf')
    # canv_leg.SaveAs(outname+'_leg.png')
    

# draw a list of TH1F with color palette
def DrawTH1FList(th1flist, logy, xaxistitle, ytitle_unit, legpos, legtitle, legtext, outname):
    gStyle.SetPalette(kLightTemperature)
    ymax = -1
    ymin = 10e10
    for h in th1flist:
        h.Sumw2()
        if h.GetMaximum() > ymax:
            ymax = h.GetMaximum()
        if h.GetMinimum(0) < ymin:
            ymin = h.GetMinimum(0)

    binwidth = th1flist[0].GetXaxis().GetBinWidth(1)

    c = TCanvas('c', 'c', 800, 700)
    c.cd()
    if logy:
        gPad.SetLogy()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)

    th1flist[0].GetXaxis().SetTitle(xaxistitle)
    if ytitle_unit == '':
        th1flist[0].GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth))
    else:
        th1flist[0].GetYaxis().SetTitle('Entries / ({:g} {unit})'.format(binwidth, unit=ytitle_unit))
    th1flist[0].GetXaxis().SetTickSize(TickSize)
    th1flist[0].GetYaxis().SetTickSize(TickSize)
    th1flist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    th1flist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    th1flist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    th1flist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    th1flist[0].GetXaxis().SetTitleOffset(1.1)
    th1flist[0].GetYaxis().SetTitleOffset(1.35)
    if logy:
        th1flist[0].GetYaxis().SetRangeUser(ymin*0.5, ymax*10)
    else:
        th1flist[0].GetYaxis().SetRangeUser(0, ymax*1.3)

    for i in range(len(th1flist)):
        th1flist[i].SetLineWidth(2)
        if i == 0:
            th1flist[i].Draw('HISTPLC')
        else:
            th1flist[i].Draw('SAME HISTPLC')

    leg = TLegend(legpos[0], legpos[1], legpos[2], legpos[3])
    leg.SetHeader(legtitle, 'l')
    leg.SetNColumns(2)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    for i in range(len(th1flist)):
        leg.AddEntry(th1flist[i], legtext[i], 'l')
    leg.Draw()

    c.RedrawAxis()
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')


if __name__ == "__main__":
    plotpath = "./PV_TruthReco/"
    os.makedirs(plotpath, exist_ok=True)

    # listfname1 = ['/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth2p0_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth2p5_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth2p9_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth3p0_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth3p5_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth4p0_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth5_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth6_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth7_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth8_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth9_3DVertex_twohalves.root',
    #              '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth10_3DVertex_twohalves.root']
    
    # list_gap1 = [2, 2.5, 2.9, 3, 3.5, 4, 5, 6, 7, 8, 9, 10]
    # list_sigmaeff1 = []
    # list_median1 = []
    
    # for i, fname in enumerate(listfname1):
    #     diffmin, diffmax, mdn, se = getPVzDiffSigma(fname)
    #     list_sigmaeff1.append(se)
    #     list_median1.append(mdn)
    #     print("GapNorth = {}: median = {}, sigmaeff = {}".format(list_gap1[i], mdn, se))

    # Draw_Graph_wFit(array('d', list_gap1), array('d', list_sigmaeff1), 0.15, 'Gap at the north end (mm)', '#sigma_{eff} (mm)', 1.3, [min(list_gap1)-1, max(list_gap1)+1], [0, max(list_sigmaeff1)*1.2], True, plotpath+'PVzDiff_sigmaeff_asFuncOfGapNorth')
    # Draw_Graph_wFit(array('d', list_gap1), array('d', list_median1), 0.15, 'Gap at the north end (mm)', 'Median (mm)', 1.3, [min(list_gap1)-1, max(list_gap1)+1], [-0.05, 0.05], False, plotpath+'PVzDiff_median_asFuncOfGapNorth')

    # list_gap2 = ['{:.1f}'.format(1.6+0.1*i) for i in range(25)]
    # listfname2 = ['/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth{}_3DVertex_twohalves.root'.format(str_pttop(i)) for i in list_gap2]
    # list_sigmaeff2 = []
    # list_median2 = []

    # for i, fname in enumerate(listfname2):
    #     diffmin, diffmax, mdn, se = getPVzDiffSigma(fname)
    #     list_sigmaeff2.append(se)
    #     list_median2.append(mdn)
    #     print("GapNorth = {}: median = {}, sigmaeff = {}".format(list_gap2[i], mdn, se))

    # Draw_Graph_wFit(array('d', [1.6+0.1*i for i in range(25)]), array('d', list_sigmaeff2), 0.15, 'Gap at the north end (mm)', '#sigma_{eff} (mm)', 1.3, [min([1.6+0.1*i for i in range(25)])-0.1, max([1.6+0.1*i for i in range(25)])+0.1], [0, max(list_sigmaeff2)*1.2], True, plotpath+'PVzDiff_sigmaeff_asFuncOfGapNorth_range2')
    # Draw_Graph_wFit(array('d', [1.6+0.1*i for i in range(25)]), array('d', list_median2), 0.15, 'Gap at the north end (mm)', 'Median (mm)', 1.3, [min([1.6+0.1*i for i in range(25)])-0.1, max([1.6+0.1*i for i in range(25)])+0.1], [-0.05, 0.05], False, plotpath+'PVzDiff_median_asFuncOfGapNorth_range2')

    plotpath_scan = "./PV_TruthReco/AsymScan_FieldOff/"
    os.makedirs(plotpath_scan, exist_ok=True)
    RecoPVfilepath = '/sphenix/user/hjheng/TrackletAna/minitree/HijingAuAuMB_NoPileup_0T_RecoVtx_Optimization/'
    list_gapupper = ['{:.2f}'.format(0.75+0.1*i).replace('.', 'p') for i in range(21)]
    list_profile_DiffPVx = []
    list_profile_DiffPVy = []
    list_profile_DiffPVz = []
    list_hist_RecoPVxRes = []
    list_hist_RecoPVyRes = []
    list_hist_RecoPVzRes = []
    for i, gapupper_str in enumerate(list_gapupper):
        hM_RecoPVz_DiffPVx, hM_RecoPVz_DiffPVy, hM_RecoPVz_DiffPVz, hM_RecoPVx_Res, hM_RecoPVy_Res, hM_RecoPVz_Res = getPVDiffHists(fname=RecoPVfilepath+'TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth3p5_GapUpper{}_CentShift0p00_3DVertex_twohalves.root'.format(gapupper_str), DiffPVx_binparam=[100,-0.5,0.5], DiffPVy_binparam=[100,-0.5,0.5], DiffPVz_binparam=[100,-0.5,0.5])
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVx, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{x}^{1} - PV_{x}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVx_GapNorth3p5_GapUpper{}_CentShift0p00'.format(gapupper_str))
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVy, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{y}^{1} - PV_{y}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVy_GapNorth3p5_GapUpper{}_CentShift0p00'.format(gapupper_str))
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVz, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{z}^{1} - PV_{z}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVz_GapNorth3p5_GapUpper{}_CentShift0p00'.format(gapupper_str))
        profile_DiffPVx = hM_RecoPVz_DiffPVx.ProfileX('DiffPVx_ProfileX_GapUpper{}'.format(gapupper_str))
        profile_DiffPVy = hM_RecoPVz_DiffPVy.ProfileX('DiffPVy_ProfileX_GapUpper{}'.format(gapupper_str))
        profile_DiffPVz = hM_RecoPVz_DiffPVz.ProfileX('DiffPVz_ProfileX_GapUpper{}'.format(gapupper_str))
        list_profile_DiffPVx.append(profile_DiffPVx)
        list_profile_DiffPVy.append(profile_DiffPVy)
        list_profile_DiffPVz.append(profile_DiffPVz)
        list_hist_RecoPVxRes.append(hM_RecoPVx_Res)
        list_hist_RecoPVyRes.append(hM_RecoPVy_Res)
        list_hist_RecoPVzRes.append(hM_RecoPVz_Res)
    
    DrawTProfileList(tprofilelist=list_profile_DiffPVx, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{x}^{1} - PV_{x}^{2} (cm)', ylim=[-0.25,0.25], legpos=[0.6, 0.16, 0.88, 0.44], legtitle='East-half gap', legtext=['{}mm'.format(list_gapupper[i]).replace('p', '.') for i in range(len(list_gapupper))], outname=plotpath_scan+'RecoPV_DiffPVx_TProfile_GapNorth3p5_ScanGapUpper_CentShift0p00')
    DrawTProfileList(tprofilelist=list_profile_DiffPVy, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{y}^{1} - PV_{y}^{2} (cm)', ylim=[-0.2,0.2], legpos=[0.4, 0.155, 0.65, 0.43], legtitle='East-half gap', legtext=['{}mm'.format(list_gapupper[i]).replace('p', '.') for i in range(len(list_gapupper))], outname=plotpath_scan+'RecoPV_DiffPVy_TProfile_GapNorth3p5_ScanGapUpper_CentShift0p00')
    DrawTProfileList(tprofilelist=list_profile_DiffPVz, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{z}^{1} - PV_{z}^{2} (cm)', ylim=[-0.45,0.45], legpos=[0.4, 0.155, 0.65, 0.43], legtitle='East-half gap', legtext=['{}mm'.format(list_gapupper[i]).replace('p', '.') for i in range(len(list_gapupper))], outname=plotpath_scan+'RecoPV_DiffPVz_TProfile_GapNorth3p5_ScanGapUpper_CentShift0p00')
    DrawTH1FList(th1flist=list_hist_RecoPVxRes, logy=True, xaxistitle='#Delta(PV_{x}^{Reco}, PV_{x}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.58, 0.89, 0.88], legtitle='East-half gap', legtext=['{}mm'.format(list_gapupper[i]).replace('p', '.') for i in range(len(list_gapupper))], outname=plotpath_scan+'RecoPVxRes_GapNorth3p5_ScanGapUpper_CentShift0p00')
    DrawTH1FList(th1flist=list_hist_RecoPVyRes, logy=True, xaxistitle='#Delta(PV_{y}^{Reco}, PV_{y}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.58, 0.89, 0.88], legtitle='East-half gap', legtext=['{}mm'.format(list_gapupper[i]).replace('p', '.') for i in range(len(list_gapupper))], outname=plotpath_scan+'RecoPVyRes_GapNorth3p5_ScanGapUpper_CentShift0p00')
    DrawTH1FList(th1flist=list_hist_RecoPVzRes, logy=True, xaxistitle='#Delta(PV_{z}^{Reco}, PV_{z}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.58, 0.89, 0.88], legtitle='East-half gap', legtext=['{}mm'.format(list_gapupper[i]).replace('p', '.') for i in range(len(list_gapupper))], outname=plotpath_scan+'RecoPVzRes_GapNorth3p5_ScanGapUpper_CentShift0p00')

    list_profile_DiffPVx.clear()
    list_profile_DiffPVy.clear()
    list_profile_DiffPVz.clear()
    list_hist_RecoPVxRes.clear()
    list_hist_RecoPVyRes.clear()
    list_hist_RecoPVzRes.clear()
    list_centshift = ['{:.2f}'.format(-0.5+0.1*i).replace('.', 'p') for i in range(11)]
    # binparam_diffz = [[90, -3, 3], [90, -3, 3], [90,-1.5, 1.5], [90, -0.9, 0.9], [90, -0.6, 0.6], [90, -0.3, 0.3], [90, -0.6, 0.6], [90, -0.9, 0.9], [90,-1.5, 1.5], [90, -3, 3], [90, -3, 3]]
    # binparam_diffx = [[90, -3, 3], [90, -3, 3], [90,-1.5, 1.5], [90, -0.9, 0.9], [90, -0.6, 0.6], [90, -0.3, 0.3], [90, -0.6, 0.6], [90, -0.9, 0.9], [90,-1.5, 1.5], [90, -3, 3], [90, -3, 3]]
    for i, centshift_str in enumerate(list_centshift):
        hM_RecoPVz_DiffPVx, hM_RecoPVz_DiffPVy, hM_RecoPVz_DiffPVz, hM_RecoPVx_Res, hM_RecoPVy_Res, hM_RecoPVz_Res = getPVDiffHists(fname=RecoPVfilepath+'TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth3p5_GapUpper1p75_CentShift{}_3DVertex_twohalves.root'.format(centshift_str), DiffPVx_binparam=[100,-0.5,0.5], DiffPVy_binparam=[100,-0.5,0.5], DiffPVz_binparam=[100,-5,5])
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVx, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{x}^{1} - PV_{x}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVx_GapNorth3p5_GapUpper1p75_CentShift{}'.format(centshift_str))
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVy, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{y}^{1} - PV_{y}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVy_GapNorth3p5_GapUpper1p75_CentShift{}'.format(centshift_str))
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVz, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{z}^{1} - PV_{z}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVz_GapNorth3p5_GapUpper1p75_CentShift{}'.format(centshift_str))
        profile_DiffPVx = hM_RecoPVz_DiffPVx.ProfileX('DiffPVx_ProfileX_CentShift{}'.format(centshift_str))
        profile_DiffPVy = hM_RecoPVz_DiffPVy.ProfileX('DiffPVy_ProfileX_CentShift{}'.format(centshift_str))
        profile_DiffPVz = hM_RecoPVz_DiffPVz.ProfileX('DiffPVz_ProfileX_CentShift{}'.format(centshift_str))
        list_profile_DiffPVx.append(profile_DiffPVx)
        list_profile_DiffPVy.append(profile_DiffPVy)
        list_profile_DiffPVz.append(profile_DiffPVz)
        list_hist_RecoPVxRes.append(hM_RecoPVx_Res)
        list_hist_RecoPVyRes.append(hM_RecoPVy_Res)
        list_hist_RecoPVzRes.append(hM_RecoPVz_Res)

    DrawTProfileList(tprofilelist=list_profile_DiffPVx, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{x}^{1} - PV_{x}^{2} (cm)', ylim=[-0.6,0.6], legpos=[0.4, 0.17, 0.65, 0.4], legtitle='Shift in MVTX Center', legtext=['{}mm'.format(list_centshift[i]).replace('p', '.') for i in range(len(list_centshift))], outname=plotpath_scan+'RecoPV_DiffPVx_TProfile_GapNorth3p5_GapUpper1p75_ScanCentShift')
    DrawTProfileList(tprofilelist=list_profile_DiffPVy, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{y}^{1} - PV_{y}^{2} (cm)', ylim=[-0.6,0.6], legpos=[0.4, 0.17, 0.65, 0.4], legtitle='Shift in MVTX Center', legtext=['{}mm'.format(list_centshift[i]).replace('p', '.') for i in range(len(list_centshift))], outname=plotpath_scan+'RecoPV_DiffPVy_TProfile_GapNorth3p5_GapUpper1p75_ScanCentShift')
    DrawTProfileList(tprofilelist=list_profile_DiffPVz, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{z}^{1} - PV_{z}^{2} (cm)', ylim=[-3,3], legpos=[0.4, 0.17, 0.65, 0.4], legtitle='Shift in MVTX Center', legtext=['{}mm'.format(list_centshift[i]).replace('p', '.') for i in range(len(list_centshift))], outname=plotpath_scan+'RecoPV_DiffPVz_TProfile_GapNorth3p5_GapUpper1p75_ScanCentShift')
    DrawTH1FList(th1flist=list_hist_RecoPVxRes, logy=True, xaxistitle='#Delta(PV_{x}^{Reco}, PV_{x}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.68, 0.89, 0.88], legtitle='Shift in MVTX Center', legtext=['{}mm'.format(list_centshift[i]).replace('p', '.') for i in range(len(list_centshift))], outname=plotpath_scan+'RecoPVxRes_GapNorth3p5_GapUpper1p75_ScanCentShift')
    DrawTH1FList(th1flist=list_hist_RecoPVyRes, logy=True, xaxistitle='#Delta(PV_{y}^{Reco}, PV_{y}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.68, 0.89, 0.88], legtitle='Shift in MVTX Center', legtext=['{}mm'.format(list_centshift[i]).replace('p', '.') for i in range(len(list_centshift))], outname=plotpath_scan+'RecoPVyRes_GapNorth3p5_GapUpper1p75_ScanCentShift')
    DrawTH1FList(th1flist=list_hist_RecoPVzRes, logy=True, xaxistitle='#Delta(PV_{z}^{Reco}, PV_{z}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.68, 0.89, 0.88], legtitle='Shift in MVTX Center', legtext=['{}mm'.format(list_centshift[i]).replace('p', '.') for i in range(len(list_centshift))], outname=plotpath_scan+'RecoPVzRes_GapNorth3p5_GapUpper1p75_ScanCentShift')

    list_profile_DiffPVx.clear()
    list_profile_DiffPVy.clear()
    list_profile_DiffPVz.clear()
    list_hist_RecoPVxRes.clear()
    list_hist_RecoPVyRes.clear()
    list_hist_RecoPVzRes.clear()
    list_gapnorth = [1.5 + i * 0.1 for i in range(21)]
    for i, gapnorth in enumerate(list_gapnorth):
        gapupper = gapnorth / 2.
        gapnorth_str = '{:.1f}'.format(gapnorth).replace('.', 'p')
        gapupper_str = '{:.2f}'.format(gapupper).replace('.', 'p')
        hM_RecoPVz_DiffPVx, hM_RecoPVz_DiffPVy, hM_RecoPVz_DiffPVz, hM_RecoPVx_Res, hM_RecoPVy_Res, hM_RecoPVz_Res = getPVDiffHists(fname=RecoPVfilepath+'TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth{}_GapUpper{}_CentShift0p0_3DVertex_twohalves.root'.format(gapnorth_str, gapupper_str), DiffPVx_binparam=[100,-0.5,0.5], DiffPVy_binparam=[100,-0.5,0.5], DiffPVz_binparam=[100,-0.5,0.5])
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVx, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{x}^{1} - PV_{x}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVx_GapNorth{}_GapUpper{}_CentShift0p0'.format(gapnorth_str, gapupper_str))
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVy, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{y}^{1} - PV_{y}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVy_GapNorth{}_GapUpper{}_CentShift0p0'.format(gapnorth_str, gapupper_str))
        Draw_2Dhist_wProfile(hist=hM_RecoPVz_DiffPVz, logz=False, norm1=False, rmargin=0.15, XaxisName='Reco PV_{z} (cm)', YaxisName='PV_{z}^{1} - PV_{z}^{2} (cm)', drawopt='colz', outname=plotpath_scan+'RecoPV_DiffPVz_GapNorth{}_GapUpper{}_CentShift0p0'.format(gapnorth_str, gapupper_str))
        profile_DiffPVx = hM_RecoPVz_DiffPVx.ProfileX('DiffPVx_ProfileX_GapNorth{}_GapUpper{}_CentShift0p0'.format(gapnorth_str, gapupper_str))
        profile_DiffPVy = hM_RecoPVz_DiffPVy.ProfileX('DiffPVy_ProfileX_GapNorth{}_GapUpper{}_CentShift0p0'.format(gapnorth_str, gapupper_str))
        profile_DiffPVz = hM_RecoPVz_DiffPVz.ProfileX('DiffPVz_ProfileX_GapNorth{}_GapUpper{}_CentShift0p0'.format(gapnorth_str, gapupper_str))
        list_profile_DiffPVx.append(profile_DiffPVx)
        list_profile_DiffPVy.append(profile_DiffPVy)
        list_profile_DiffPVz.append(profile_DiffPVz)
        list_hist_RecoPVxRes.append(hM_RecoPVx_Res)
        list_hist_RecoPVyRes.append(hM_RecoPVy_Res)
        list_hist_RecoPVzRes.append(hM_RecoPVz_Res)

    DrawTProfileList(tprofilelist=list_profile_DiffPVx, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{x}^{1} - PV_{x}^{2} (cm)', ylim=[-0.3, 0.3], legpos=[0.4, 0.16, 0.65, 0.43], legtitle='Gap North', legtext=['{:.2f}mm'.format(list_gapnorth[i]).replace('p', '.') for i in range(len(list_gapnorth))], outname=plotpath_scan+'RecoPV_DiffPVx_TProfile_GapNorthScan_GapUpperHalfGapNorth_CentShift0p0')
    DrawTProfileList(tprofilelist=list_profile_DiffPVy, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{y}^{1} - PV_{y}^{2} (cm)', ylim=[-0.3, 0.3], legpos=[0.4, 0.16, 0.65, 0.43], legtitle='Gap North', legtext=['{:.2f}mm'.format(list_gapnorth[i]).replace('p', '.') for i in range(len(list_gapnorth))], outname=plotpath_scan+'RecoPV_DiffPVy_TProfile_GapNorthScan_GapUpperHalfGapNorth_CentShift0p0')
    DrawTProfileList(tprofilelist=list_profile_DiffPVz, xaxistitle='Reco PV_{z} (cm)', yaxistitle='PV_{z}^{1} - PV_{z}^{2} (cm)', ylim=[-0.3, 0.3], legpos=[0.4, 0.16, 0.65, 0.43], legtitle='Gap North', legtext=['{:.2f}mm'.format(list_gapnorth[i]).replace('p', '.') for i in range(len(list_gapnorth))], outname=plotpath_scan+'RecoPV_DiffPVz_TProfile_GapNorthScan_GapUpperHalfGapNorth_CentShift0p0')
    DrawTH1FList(th1flist=list_hist_RecoPVxRes, logy=True, xaxistitle='#Delta(PV_{x}^{Reco}, PV_{x}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.58, 0.89, 0.88], legtitle='Gap North', legtext=['{:.2f}mm'.format(list_gapnorth[i]).replace('p', '.') for i in range(len(list_gapnorth))], outname=plotpath_scan+'RecoPVxRes_GapNorthScan_GapUpperHalfGapNorth_CentShift0p0')
    DrawTH1FList(th1flist=list_hist_RecoPVyRes, logy=True, xaxistitle='#Delta(PV_{y}^{Reco}, PV_{y}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.58, 0.89, 0.88], legtitle='Gap North', legtext=['{:.2f}mm'.format(list_gapnorth[i]).replace('p', '.') for i in range(len(list_gapnorth))], outname=plotpath_scan+'RecoPVyRes_GapNorthScan_GapUpperHalfGapNorth_CentShift0p0')
    DrawTH1FList(th1flist=list_hist_RecoPVzRes, logy=True, xaxistitle='#Delta(PV_{z}^{Reco}, PV_{z}^{Truth}) (cm)', ytitle_unit='cm', legpos=[0.65, 0.58, 0.89, 0.88], legtitle='Gap North', legtext=['{:.2f}mm'.format(list_gapnorth[i]).replace('p', '.') for i in range(len(list_gapnorth))], outname=plotpath_scan+'RecoPVzRes_GapNorthScan_GapUpperHalfGapNorth_CentShift0p0')
