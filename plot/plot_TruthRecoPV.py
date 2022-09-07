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

def Draw_Graph(xarray, yarray, lmargin, xaxistitle, yaxistitle, ytitleoffset, xlim, ylim, outname):
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
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_GraphError(xarray, xerrarray, yarray, yerrarray, lmargin, xaxistitle, yaxistitle, ytitleoffset, xlim, ylim, outname):
    c = TCanvas('c1', 'c1', 800, 700)
    c.cd()
    gPad.SetLeftMargin(lmargin)
    grerr = TGraphErrors(len(xarray), xarray, yarray, xerrarray, yerrarray)
    # gr.SetMarkerStyle(21)
    grerr.GetXaxis().SetTitle(xaxistitle)
    grerr.GetYaxis().SetTitle(yaxistitle)
    grerr.GetXaxis().SetLimits(xlim[0], xlim[1]);
    grerr.SetMinimum(ylim[0])
    grerr.SetMaximum(ylim[1])
    grerr.SetLineWidth(2)
    grerr.SetMarkerSize(2)
    grerr.GetXaxis().SetTitleOffset(1.2)
    grerr.GetYaxis().SetTitleOffset(ytitleoffset)
    grerr.GetXaxis().SetNdivisions(10, 0, 0, kTRUE);
    grerr.Draw('AP')
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_1Dhist_wPercentile(hist, l_percentile, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, outname):
    hist.Sumw2()
    binwidth = hist.GetXaxis().GetBinWidth(1)
    binwidth2 = hist.GetXaxis().GetBinWidth(2)
    printbinwidth = True
    if binwidth != binwidth2:
        printbinwidth = False
    
    c = TCanvas('c', 'c', 800, 700)
    if norm1:
        hist.Scale(1. / hist.Integral(-1, -1))
    if logy:
        c.SetLogy()
    if logx:
        c.SetLogx()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    if printbinwidth:
        if norm1:
            if Ytitle_unit == '':
                hist.GetYaxis().SetTitle('Normalized entries / ({:g})'.format(binwidth))
            else:
                hist.GetYaxis().SetTitle('Normalized entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))
        else:
            if Ytitle_unit == '':
                hist.GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth))
            else:
                hist.GetYaxis().SetTitle('Entries / ({:g} {unit})'.format(binwidth, unit=Ytitle_unit))
    else:
        if norm1:
            if Ytitle_unit == '':
                hist.GetYaxis().SetTitle('Normalized entries')
            else:
                hist.GetYaxis().SetTitle('Normalized entries {unit})'.format(unit=Ytitle_unit))
        else:
            if Ytitle_unit == '':
                hist.GetYaxis().SetTitle('Entries')
            else:
                hist.GetYaxis().SetTitle('Entries {unit}'.format(unit=Ytitle_unit))

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
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.13,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
    linestorage = []
    textstorage = []
    for i,p in enumerate(l_percentile):
        if (logx == False and hist.GetXaxis().GetXmax() <= 1000 and i > 9):
            continue
        pline = TLine(p, 0, p, (hist.GetMaximum()/ymaxscale)*0.8)
        pline.SetLineWidth(1)
        pline.SetLineStyle(kDashed)
        pline.SetLineColor(2)
        if len(linestorage)==0:
            pline.Draw()
        else:
            pline.Draw('same')
        gPad.Update()
        linestorage.append(pline)

        if (logx == False and hist.GetXaxis().GetXmax() > 1000 and i > 9) or (logx == False and hist.GetXaxis().GetXmax() <= 1000 and i <= 9) or logx == True:
            ptext = TText(p,(hist.GetMaximum()/ymaxscale)*0.1,'{}-{}%'.format(5*((len(l_percentile)+1)-(i+2)),5*((len(l_percentile)+1)-(i+1))))
            ptext.SetTextAlign(13)
            ptext.SetTextSize(0.02)
            ptext.SetTextColor(2)
            ptext.SetTextAngle(90)
            if len(textstorage)==0:
                ptext.Draw()
            else:
                ptext.Draw('same')
            gPad.Update()
            textstorage.append(ptext)

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
    plotpath = './PV_TruthReco/'
    os.makedirs(plotpath, exist_ok=True)

    df = ROOT.RDataFrame('minitree', '/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_RecoClusters_RecoVtx_TklCluster_Nevt4000.root')
    np_NhitsL1 = df.AsNumpy(columns=['NhitsLayer1'])
    NhitsL1_percentile = []
    NhitsL1_percentile_cut = [0]
    NpercentileDiv = 20
    for i in range(NpercentileDiv-1):
        # print('percentile={}%, Nhits={}'.format((i+1)*5, np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5)))
        NhitsL1_percentile.append(np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5))
        if i % 2 == 1:
            NhitsL1_percentile_cut.append(np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5))
    NhitsL1_percentile_cut.append(10000)

    Ncuts = 50
    incre = 0.01
    list_absdiff = array('d', [i*incre for i in range(Ncuts+1)])
    list_count_absdiff_trig = array('d', [0 for i in range(Ncuts+1)])
    list_count_absdiff_MostNpart = array('d', [0 for i in range(Ncuts+1)])
    list_count_absdiff_Nvtx1 = array('d', [0 for i in range(Ncuts+1)])
    # print (list_absdiff)

    hM_PVz = TH1F('hM_PVz', 'hM_PVz', 60, -15, 15);
    hM_NClusLayer1_NTruthVtx1 = TH1F('hM_NClusLayer1_NTruthVtx1', 'hM_NClusLayer1_NTruthVtx1', 100, 0, 10000)
    hM_NClusLayer1_NTruthVtx1_altrange1 = TH1F('hM_NClusLayer1_NTruthVtx1_altrange1', 'hM_NClusLayer1_NTruthVtx1_altrange1', 50, 0, 1000)
    NBins = 50
    xbin = np.logspace(0, 4, NBins+1)
    hM_NClusLayer1_NTruthVtx1_logX = TH1F('hM_NClusLayer1_NTruthVtx1_logX', 'hM_NClusLayer1_NTruthVtx1_logX', NBins, xbin)

    hM_DiffVtxZ_trig_altrange = TH1F('hM_DiffVtxZ_trig_altrange', 'hM_DiffVtxZ_trig_altrange', 100, -0.5, 0.5)
    hM_DiffVtxZ_mostNpart_altrange = TH1F('hM_DiffVtxZ_mostNpart_altrange', 'hM_DiffVtxZ_mostNpart_altrange', 100, -0.5, 0.5)
    hM_DiffVtxZ_Nvtx1_altrange = TH1F('hM_DiffVtxZ_Nvtx1_altrange', 'hM_DiffVtxZ_Nvtx1_altrange', 40, -0.1, 0.1)
    list_hM_DiffVtxZ_Nvtx1 = []
    for i in range(int(NpercentileDiv/2)):
        list_hM_DiffVtxZ_Nvtx1.append(TH1F('hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 40, -0.1, 0.1))

    hM_DiffVtxZTrig_NhitsLayer1_altrange = TH2F('hM_DiffVtxZTrig_NhitsLayer1_altrange', 'hM_DiffVtxZTrig_NhitsLayer1_altrange', 100, -0.1, 0.1, 200, 0, 20000)
    hM_DiffVtxZMostNpart_NhitsLayer1_altrange = TH2F('hM_DiffVtxZMostNpart_NhitsLayer1_altrange', 'hM_DiffVtxZMostNpart_NhitsLayer1_altrange', 100, -0.1, 0.1, 200, 0, 20000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 100, -0.1, 0.1, 100, 0, 10000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2 = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 50, -0.05, 0.05, 50, 0, 2000)

    Nevt_Nvtx1 = 0
    f = TFile('/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_RecoClusters_RecoVtx_TklCluster_Nevt4000.root', 'r')
    tree = f.Get('minitree')
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)

        hM_PVz.Fill(tree.PV_z)
        hM_DiffVtxZ_trig_altrange.Fill(tree.PV_z - tree.TruthPV_trig_z)
        hM_DiffVtxZ_mostNpart_altrange.Fill(tree.PV_z - tree.TruthPV_mostNpart_z)
        hM_DiffVtxZTrig_NhitsLayer1_altrange.Fill((tree.PV_z - tree.TruthPV_trig_z), tree.NhitsLayer1)
        hM_DiffVtxZMostNpart_NhitsLayer1_altrange.Fill((tree.PV_z - tree.TruthPV_mostNpart_z), tree.NhitsLayer1)
        if tree.NTruthVtx == 1:
            Nevt_Nvtx1 += 1
            hM_NClusLayer1_NTruthVtx1.Fill(tree.NhitsLayer1)
            hM_NClusLayer1_NTruthVtx1_altrange1.Fill(tree.NhitsLayer1)
            hM_NClusLayer1_NTruthVtx1_logX.Fill(tree.NhitsLayer1)
            for i in range(len(NhitsL1_percentile_cut)):
                if tree.NhitsLayer1 > NhitsL1_percentile_cut[i] and tree.NhitsLayer1 <= NhitsL1_percentile_cut[i+1]:
                    list_hM_DiffVtxZ_Nvtx1[i].Fill(tree.PV_z - tree.TruthPV_trig_z)

            hM_DiffVtxZ_Nvtx1_altrange.Fill(tree.PV_z - tree.TruthPV_trig_z)
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange.Fill((tree.PV_z - tree.TruthPV_trig_z), tree.NhitsLayer1)
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2.Fill((tree.PV_z - tree.TruthPV_trig_z), tree.NhitsLayer1)

            if abs(tree.TruthPV_trig_z - tree.PV_z) > 1 and tree.PV_z > -999:
                print (tree.TruthPV_trig_z, tree.PV_z, tree.NhitsLayer1)
        
        for i, diff in enumerate(list_absdiff):
            if abs(tree.TruthPV_trig_z - tree.PV_z) > diff:
                list_count_absdiff_trig[i] += 1
            
            if abs(tree.TruthPV_mostNpart_z - tree.PV_z) > diff:
                list_count_absdiff_MostNpart[i] += 1

            if tree.NTruthVtx == 1 and abs(tree.TruthPV_trig_z - tree.PV_z) > diff:
                list_count_absdiff_Nvtx1[i] += 1

    list_count_absdiff_trig = array('d', [i/tree.GetEntries() for i in list_count_absdiff_trig])
    list_count_absdiff_MostNpart = array('d', [i/tree.GetEntries() for i in list_count_absdiff_MostNpart])
    list_count_absdiff_Nvtx1 = array('d', [i/Nevt_Nvtx1 for i in list_count_absdiff_MostNpart])
    # print (list_count_absdiff_trig)

    hM_NTruthVtx = f.Get('hM_NTruthVtx')
    hM_DiffVtxZ_trig = f.Get('hM_DiffVtxZ_trig')
    hM_DiffVtxZ_mostNpart = f.Get('hM_DiffVtxZ_mostNpart')
    hM_DiffVtxZ_Nvtx1 = f.Get('hM_DiffVtxZ_Nvtx1')
    hM_TruthPVzTrig_RecoPVz = f.Get('hM_TruthPVzTrig_RecoPVz')
    hM_TruthPVzMostNpart_RecoPVz = f.Get('hM_TruthPVzMostNpart_RecoPVz')
    hM_TruthPVz_RecoPVz_Nvtx1 = f.Get('hM_TruthPVz_RecoPVz_Nvtx1')
    hM_DiffVtxZTrig_NhitsLayer1 = f.Get('hM_DiffVtxZTrig_NhitsLayer1')
    hM_DiffVtxZMostNpart_NhitsLayer1 = f.Get('hM_DiffVtxZMostNpart_NhitsLayer1')
    hM_DiffVtxZ_NhitsLayer1_Nvtx1 = f.Get('hM_DiffVtxZ_NhitsLayer1_Nvtx1')

    Draw_1Dhist(hM_PVz, False, False, 1.3, 'Primary vertex V_{z} (cm)', 'cm', plotpath+'PVz_woCuts')
    Draw_1Dhist(hM_NTruthVtx, False, False, 1.2, 'Number of truth vertices', '', plotpath+'RecoClusters_NTruthVtx')
    Draw_1Dhist(hM_DiffVtxZ_trig, False, True, 150, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_trig')
    Draw_1Dhist(hM_DiffVtxZ_mostNpart, False, True, 150, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_mostNpart')
    Draw_1Dhist(hM_DiffVtxZ_Nvtx1, False, True, 150, '#Deltaz(PV_{Truth}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1')
    Draw_1Dhist_wTF1(hM_DiffVtxZ_trig_altrange, False, True, 150, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_trig_altrange')
    Draw_1Dhist_wTF1(hM_DiffVtxZ_mostNpart_altrange, False, True, 150, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_mostNpart_altrange')
    Draw_1Dhist_wTF1(hM_DiffVtxZ_Nvtx1_altrange, False, True, 50, '#Deltaz(PV_{Truth}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1_altrange')
    l_res = []
    l_reserr = []
    for i, hist in enumerate(list_hM_DiffVtxZ_Nvtx1):
        mean, meanerr, sigma, sigmaerr = Draw_1Dhist_wTF1(hist, False, True, 50, '#Deltaz(PV_{Truth}, PV_{Reco}) (cm)', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1))
        l_res.append(sigma)
        l_reserr.append(sigmaerr)

    list_res = array('d', l_res)
    list_reserr = array('d', l_reserr)
    Draw_GraphError(array('d', [i+1 for i in range(len(list_res))]), array('d', [0 for i in range(len(list_res))]), list_res, list_reserr, 0.15, 'Bin', 'Vertex resolution (cm)', 1.6, [0.5, len(list_res)+0.5], [0, 0.02], plotpath+'VtxReolution_Bin')

    Draw_1Dhist_wPercentile(hM_NClusLayer1_NTruthVtx1, NhitsL1_percentile, False, False, True, 5, 'Number of clusters (1st layer)', '', plotpath+'NClusLayer1_NTruthVtx1_wPercentile')
    Draw_1Dhist_wPercentile(hM_NClusLayer1_NTruthVtx1_altrange1, NhitsL1_percentile, False, False, True, 5, 'Number of clusters (1st layer)', '', plotpath+'NClusLayer1_NTruthVtx1_wPercentile_altrange1')
    Draw_1Dhist_wPercentile(hM_NClusLayer1_NTruthVtx1_logX, NhitsL1_percentile, False, True, True, 5, 'Number of clusters (1st layer)', '', plotpath+'NClusLayer1_NTruthVtx1_wPercentile_logX')

    Draw_2Dhist(hM_TruthPVzTrig_RecoPVz, False, False, 0.11, 'Truth PV (trig) Z pos (cm)', 'Reco PV Z pos (cm)', plotpath+'RecoClusters_TruthTrigRecoPVz')
    Draw_2Dhist(hM_TruthPVzMostNpart_RecoPVz, False, False, 0.11, 'Truth PV Z (most N_{particles}) pos (cm)', 'Reco PV Z pos (cm)', plotpath+'RecoClusters_TruthMostNpartRecoPVz')
    Draw_2Dhist(hM_TruthPVz_RecoPVz_Nvtx1, False, False, 0.11, 'Truth PV Z pos (cm)', 'Reco PV Z pos (cm)', plotpath+'RecoClusters_TruthRecoPVz_Nvtx1')
    Draw_2Dhist(hM_DiffVtxZTrig_NhitsLayer1, False, False, 0.11, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZTrig_NhitsLayer1')
    Draw_2Dhist(hM_DiffVtxZMostNpart_NhitsLayer1, False, False, 0.11, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZMostNpart_NhitsLayer1')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1, False, False, 0.11, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1')
    Draw_2Dhist(hM_DiffVtxZTrig_NhitsLayer1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZTrig_NhitsLayer1_altrange')
    Draw_2Dhist(hM_DiffVtxZMostNpart_NhitsLayer1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZMostNpart_NhitsLayer1_altrange')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2')

    Draw_Graph(list_absdiff, list_count_absdiff_trig, 0.15, '#Deltaz_{Threshold} (cm)', 'Fraction (|#Deltaz(PV_{Truth}^{trig}, PV_{Reco})| > #Deltaz_{Threshold})', 1.3, [0, float(Ncuts*incre)], [0, 1.1], plotpath+'RecoClusters_TruthTrigRecoPVz_AbsDiffFrac')
    Draw_Graph(list_absdiff, list_count_absdiff_MostNpart, 0.15, '#Deltaz_{Threshold} (cm)', 'Fraction (|#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco})| > #Deltaz_{Threshold})', 1.3, [0, float(Ncuts*incre)], [0, 1.1], plotpath+'RecoClusters_TruthMostNpartRecoPVz_AbsDiffFrac')
    Draw_Graph(list_absdiff, list_count_absdiff_Nvtx1, 0.2, '#Deltaz_{Threshold} (cm)', 'Fraction (|#Deltaz(PV_{Truth}, PV_{Reco})| > #Deltaz_{Threshold})', 2, [0, float(Ncuts*incre)], [0, 0.005], plotpath+'RecoClusters_TruthRecoPVz_Nvtx1_AbsDiffFrac')