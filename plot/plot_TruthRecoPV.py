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
    f1.SetParameter(1,0.)
    f1.SetParLimits(1,-0.1,0.1)
    f1.SetParameter(2,0.005)
    f1.SetParLimits(2,0,0.2)
    f1.SetLineColorAlpha(TColor.GetColor('#F54748'), 0.9)
    hist.Fit(f1, 'B')
    f1.Draw('same')
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.1,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.040)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
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

    # leg = TLegend(LeftMargin, 1-TopMargin*1.1, LeftMargin+0.01, 0.98)
    leg = TLegend(LeftMargin+0.03, 1-TopMargin-0.15, LeftMargin+0.1, 1-TopMargin-0.035)
    leg.SetTextSize(0.040)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
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
    grerr.SetMarkerSize(1.5)
    grerr.GetXaxis().SetTitleOffset(1.2)
    grerr.GetYaxis().SetTitleOffset(ytitleoffset)
    grerr.GetXaxis().SetNdivisions(10, 0, 0, kTRUE);
    grerr.Draw('AP')
    leg = TLegend(1 - RightMargin - 0.4, 1-TopMargin-0.15, 1 - RightMargin - 0.08, 1-TopMargin-0.035)
    leg.SetTextSize(0.040)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
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
    leg = TLegend(1 - RightMargin - 0.45, 1-TopMargin-0.15, 1 - RightMargin - 0.08, 1-TopMargin-0.035)
    leg.SetTextSize(0.040)
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

def Draw_1DhistsComp2(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname):
    color = ['#035397', '#9B0000']
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
            ymax = h.GetMaximum()
            ymin = h.GetMinimum(0)

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
    leg.SetTextSize(0.040)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    leg.AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "")
    leg.Draw()

    leg1 = TLegend(LeftMargin+0.04, (1-TopMargin)-0.15,
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
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_1DEff(EffErr, logx, XaxisName, axesrange, outname):
    c = TCanvas('c', 'c', 800, 700)
    if logx:
        c.SetLogx()
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
    EffErr.GetXaxis().SetTitle(XaxisName)
    EffErr.GetYaxis().SetTitle('Efficiency')
    EffErr.GetXaxis().SetRangeUser(axesrange[0], axesrange[1])
    # EffErr.GetXaxis().SetMoreLogLabels()
    EffErr.GetYaxis().SetRangeUser(axesrange[2], axesrange[3])
    EffErr.SetMarkerColor(1)
    EffErr.SetMarkerSize(1.5)
    EffErr.SetMarkerStyle(20)
    EffErr.GetXaxis().SetTitleOffset(1.2)
    # EffErr.GetYaxis().SetTitleOffset(ytitleoffset)
    # EffErr.GetXaxis().SetNdivisions(10, 0, 0, kTRUE)
    EffErr.Draw('AP')
    leg = TLegend(LeftMargin+0.03, 1-TopMargin-0.15, LeftMargin+0.1, 1-TopMargin-0.035)
    leg.SetTextSize(0.040)
    leg.SetFillStyle(0)
    leg.AddEntry('', '#it{#bf{sPHENIX}} Simulation', '')
    leg.AddEntry('','Au+Au #sqrt{s_{NN}}=200 GeV','')
    leg.Draw()
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
    gPad.SetBottomMargin(BottomMargin)
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
    leg1 = TLegend(LeftMargin+0.04, (1-TopMargin)-0.15,
                   LeftMargin+0.34, (1-TopMargin)-0.03)
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
    plotpath = './PV_TruthReco/FieldOn/'
    os.makedirs(plotpath, exist_ok=True)
    # Field off
    # fname = '/sphenix/user/hjheng/TrackletAna/minitree/HijingAuAuMB_NoPileup_0T_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth1p5_GapUpper0p75_CentShift0p0_3DVertex_twohalves.root'
    # Field on
    fname = '/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth1p5_GapUpper0p75_CentShift0p0_3DVertex_twohalves.root'

    df = ROOT.RDataFrame('minitree', fname)
    np_NhitsL1 = df.AsNumpy(columns=['NhitsLayer1'])
    NhitsL1_percentile = []
    Binedge_NhitsL1_percentile = [0]
    NhitsL1_percentile_cut = [0]
    NpercentileDiv = 20
    for i in range(NpercentileDiv-1):
        # print('percentile={}%, Nhits={}'.format((i+1)*5, np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5)))
        NhitsL1_percentile.append(np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5))
        Binedge_NhitsL1_percentile.append(np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5))
        if i % 2 == 1:
            NhitsL1_percentile_cut.append(np.percentile(np_NhitsL1['NhitsLayer1'], (i+1)*5))
    NhitsL1_percentile_cut.append(10000)
    Binedge_NhitsL1_percentile.append(10000)
    
    with open("Centrality_bin.txt", 'w') as f:
        for i in Binedge_NhitsL1_percentile:
            print('{:3g}'.format(i), file=f)

    # print(Binedge_NhitsL1_percentile, file=f)

    Ncuts = 50
    incre = 0.01
    list_absdiff = array('d', [i*incre for i in range(Ncuts+1)])
    list_count_absdiff_trig = array('d', [0 for i in range(Ncuts+1)])
    list_count_absdiff_MostNpart = array('d', [0 for i in range(Ncuts+1)])
    list_count_absdiff_Nvtx1 = array('d', [0 for i in range(Ncuts+1)])
    # print (list_absdiff)

    hM_NTruthVtx = TH1F('hM_NTruthVtx', 'hM_NTruthVtx', 11, 0, 11)
    hM_PVz = TH1F('hM_PVz', 'hM_PVz', 50, -25, 25)
    hM_TruthTrigPVz = TH1F('hM_TruthTrigPVz', 'hM_TruthTrigPVz', 50, -25, 25)
    hM_TrithTrigPVz_Nvtx1 = TH1F('hM_TrithTrigPVz_Nvtx1', 'hM_TrithTrigPVz_Nvtx1', 50, -25, 25)
    hM_NClusLayer1_NTruthVtx1 = TH1F('hM_NClusLayer1_NTruthVtx1', 'hM_NClusLayer1_NTruthVtx1', 60, 0, 6000)
    hM_NClusLayer1_NTruthVtx1_altrange1 = TH1F('hM_NClusLayer1_NTruthVtx1_altrange1', 'hM_NClusLayer1_NTruthVtx1_altrange1', 50, 0, 1000)
    NBins = 50
    xbin = np.logspace(0, 4, NBins+1)
    hM_NClusLayer1_NTruthVtx1_logX = TH1F('hM_NClusLayer1_NTruthVtx1_logX', 'hM_NClusLayer1_NTruthVtx1_logX', NBins, xbin)

    hM_NClusLayer1_all_PercBin = TH1F('hM_NClusLayer1_all_PercBin','hM_NClusLayer1_all_PercBin', NpercentileDiv, np.asarray(Binedge_NhitsL1_percentile))
    hM_NClusLayer1_NTurhtVtx1_PervBin = TH1F('hM_NClusLayer1_NTurhtVtx1_PervBin','hM_NClusLayer1_NTurhtVtx1_PervBin', NpercentileDiv, np.asarray(Binedge_NhitsL1_percentile))
    hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin = TH1F('hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin','hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin', NpercentileDiv, np.asarray(Binedge_NhitsL1_percentile))
    hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin = TH1F('hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin','hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin', NpercentileDiv, np.asarray(Binedge_NhitsL1_percentile))
    RecoPVEff_NClusLayer1_NTruthVtx1 = TH1F('RecoPVEff_NClusLayer1_NTruthVtx1', 'RecoPVEff_NClusLayer1_NTruthVtx1',NpercentileDiv, np.asarray(Binedge_NhitsL1_percentile))
    RecoPVEff_passPVSel_NClusLayer1_NTruthVtx1 = TH1F('RecoPVEff_passPVSel_NClusLayer1_NTruthVtx1', 'RecoPVEff_passPVSel_NClusLayer1_NTruthVtx1',NpercentileDiv, np.asarray(Binedge_NhitsL1_percentile))

    hM_DiffVtxZ_trig_altrange = TH1F('hM_DiffVtxZ_trig_altrange', 'hM_DiffVtxZ_trig_altrange', 100, -0.5, 0.5)
    hM_DiffVtxZ_mostNpart_altrange = TH1F('hM_DiffVtxZ_mostNpart_altrange', 'hM_DiffVtxZ_mostNpart_altrange', 100, -0.5, 0.5)
    hM_DiffVtxZ_Nvtx1_altrange = TH1F('hM_DiffVtxZ_Nvtx1_altrange', 'hM_DiffVtxZ_Nvtx1_altrange', 40, -0.1, 0.1)
    list_hM_DiffVtxX_Nvtx1 = []
    list_hM_DiffVtxY_Nvtx1 = []
    list_hM_DiffVtxZ_Nvtx1 = []
    for i in range(int(NpercentileDiv/2)):
        if i < 2:
            list_hM_DiffVtxX_Nvtx1.append(TH1F('hM_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.1, 0.1))
            list_hM_DiffVtxY_Nvtx1.append(TH1F('hM_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.1, 0.1))
            list_hM_DiffVtxZ_Nvtx1.append(TH1F('hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.05, 0.05))
        elif i >= 2 and i < 5:
            list_hM_DiffVtxX_Nvtx1.append(TH1F('hM_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.05, 0.05))
            list_hM_DiffVtxY_Nvtx1.append(TH1F('hM_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.05, 0.05))
            list_hM_DiffVtxZ_Nvtx1.append(TH1F('hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.03, 0.03))
        elif i >= 5: 
            list_hM_DiffVtxX_Nvtx1.append(TH1F('hM_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.03, 0.03))
            list_hM_DiffVtxY_Nvtx1.append(TH1F('hM_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.03, 0.03))
            list_hM_DiffVtxZ_Nvtx1.append(TH1F('hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 'hM_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1), 50, -0.01, 0.01))
        else:
            print('Weird! Check!!')

    hM_TruthPVzTrig_RecoPVz = TH2F('hM_TruthPVzTrig_RecoPVz', 'hM_TruthPVzTrig_RecoPVz', 120, -15, 15, 120, -15, 15)
    hM_TruthPVzMostNpart_RecoPVz = TH2F('hM_TruthPVzMostNpart_RecoPVz', 'hM_TruthPVzMostNpart_RecoPVz', 120, -15, 15, 120, -15, 15)
    hM_TruthPVz_RecoPVz_Nvtx1 = TH2F('hM_TruthPVz_RecoPVz_Nvtx1', 'hM_TruthPVz_RecoPVz_Nvtx1', 120, -15, 15, 120, -15, 15)
    hM_TruthPVx_RecoPVx_Nvtx1 = TH2F('hM_TruthPVx_RecoPVx_Nvtx1', 'hM_TruthPVx_RecoPVx_Nvtx1', 100, -0.1, 0.1, 100, -0.1, 0.1)
    hM_TruthPVy_RecoPVy_Nvtx1 = TH2F('hM_TruthPVy_RecoPVy_Nvtx1', 'hM_TruthPVy_RecoPVy_Nvtx1', 100, -0.1, 0.1, 100, -0.1, 0.1)
    hM_TruthPVz_NhitsLayer1 = TH2F('hM_TruthPVz_NhitsLayer1','hM_TruthPVz_NhitsLayer1', 50, -25, 25, 200, 0, 10000)
    hM_DiffVtxZTrig_NhitsLayer1 = TH2F('hM_DiffVtxZTrig_NhitsLayer1', 'hM_DiffVtxZTrig_NhitsLayer1', 100, -5, 5, 200, 0, 10000)
    hM_DiffVtxZMostNpart_NhitsLayer1 = TH2F('hM_DiffVtxZMostNpart_NhitsLayer1', 'hM_DiffVtxZMostNpart_NhitsLayer1', 100, -5, 5, 200, 0, 10000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1 = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1', 100, -5, 5, 200, 0, 10000)
    hM_DiffVtxZTrig_NhitsLayer1_altrange = TH2F('hM_DiffVtxZTrig_NhitsLayer1_altrange', 'hM_DiffVtxZTrig_NhitsLayer1_altrange', 100, -0.1, 0.1, 250, 0, 5000)
    hM_DiffVtxZMostNpart_NhitsLayer1_altrange = TH2F('hM_DiffVtxZMostNpart_NhitsLayer1_altrange', 'hM_DiffVtxZMostNpart_NhitsLayer1_altrange', 100, -0.1, 0.1, 250, 0, 5000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange', 100, -0.1, 0.1, 250, 0, 5000)
    hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2 = TH2F('hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 'hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2', 50, -0.05, 0.05, 50, 0, 2000)

    Nevt_Nvtx1 = 0
    Nevt_NotRecoPV = 0
    Nevt_NotRecoPV_TruthPVless10 = 0
    Nevt_TruthPVless10 = 0
    levt_NotRecoPV = []
    levt_NotPropRecoPV = []
    lNClusL1_NotRecoPV = []
    lNClusL1_NotPropRecoPV = []
    f = TFile(fname, 'r')
    tree = f.Get('minitree')
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)

        finalPVx = (tree.PV_x_halves1+tree.PV_x_halves2)/2. 
        finalPVy = (tree.PV_y_halves1+tree.PV_y_halves2)/2.
        finalPVz = (tree.PV_z_halves1+tree.PV_z_halves2)/2.

        if finalPVz < -99.:
            Nevt_NotRecoPV += 1

        hM_NTruthVtx.Fill(tree.NTruthVtx)
        hM_PVz.Fill(finalPVz)
        hM_TruthTrigPVz.Fill(tree.TruthPV_mostNpart_z) # change to TruthPV_mostNpart_z 
        hM_NClusLayer1_all_PercBin.Fill(tree.NhitsLayer1) 
        if tree.NTruthVtx == 1:
            hM_NClusLayer1_NTurhtVtx1_PervBin.Fill(tree.NhitsLayer1)
            hM_TrithTrigPVz_Nvtx1.Fill(tree.TruthPV_trig_z)
            hM_TruthPVz_NhitsLayer1.Fill(tree.TruthPV_mostNpart_z, tree.NhitsLayer1) # change to TruthPV_mostNpart_z 
            if finalPVz > -99:
                hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin.Fill(tree.NhitsLayer1)
                if finalPVz >= -10. and finalPVz <= 10.:
                    hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin.Fill(tree.NhitsLayer1)

        if abs(tree.TruthPV_mostNpart_z) > 10.: # change to TruthPV_mostNpart_z 
            continue

        Nevt_TruthPVless10 += 1

        if finalPVz < -99:
            print ('-->Reco PV not reconstructed: event={},NTruthVtx={},TruthPV_trig_z={:.3f},NhitsLayer1={}'.format(tree.event, tree.NTruthVtx, tree.TruthPV_mostNpart_z, tree.NhitsLayer1)) # change to TruthPV_mostNpart_z 
            levt_NotRecoPV.append(tree.event)
            lNClusL1_NotRecoPV.append(tree.NhitsLayer1)
            Nevt_NotRecoPV_TruthPVless10 += 1
        #     continue

        hM_DiffVtxZ_trig_altrange.Fill(finalPVz - tree.TruthPV_mostNpart_z) # change to TruthPV_mostNpart_z 
        hM_DiffVtxZ_mostNpart_altrange.Fill(finalPVz - tree.TruthPV_mostNpart_z)
        hM_DiffVtxZTrig_NhitsLayer1.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1) # change to TruthPV_mostNpart_z 
        hM_DiffVtxZTrig_NhitsLayer1_altrange.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1) # change to TruthPV_mostNpart_z 
        hM_DiffVtxZMostNpart_NhitsLayer1.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1)
        hM_DiffVtxZMostNpart_NhitsLayer1_altrange.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1)
        hM_TruthPVzTrig_RecoPVz.Fill(tree.TruthPV_mostNpart_z, finalPVz) # change to TruthPV_mostNpart_z 
        hM_TruthPVzMostNpart_RecoPVz.Fill(tree.TruthPV_mostNpart_z, finalPVz)
        if tree.NTruthVtx == 1:
            Nevt_Nvtx1 += 1
            hM_NClusLayer1_NTruthVtx1.Fill(tree.NhitsLayer1)
            hM_NClusLayer1_NTruthVtx1_altrange1.Fill(tree.NhitsLayer1)
            hM_NClusLayer1_NTruthVtx1_logX.Fill(tree.NhitsLayer1)
            # hM_NClusLayer1_all_PercBin.Fill(tree.NhitsLayer1)
            for i in range(len(NhitsL1_percentile_cut)):
                if tree.NhitsLayer1 > NhitsL1_percentile_cut[i] and tree.NhitsLayer1 <= NhitsL1_percentile_cut[i+1]:
                    list_hM_DiffVtxX_Nvtx1[i].Fill(finalPVx - tree.TruthPV_mostNpart_x) # change to TruthPV_mostNpart_x
                    list_hM_DiffVtxY_Nvtx1[i].Fill(finalPVy - tree.TruthPV_mostNpart_y) # change to TruthPV_mostNpart_y
                    list_hM_DiffVtxZ_Nvtx1[i].Fill(finalPVz - tree.TruthPV_mostNpart_z) # change to TruthPV_mostNpart_z 

            hM_DiffVtxZ_Nvtx1_altrange.Fill(finalPVz - tree.TruthPV_mostNpart_z) # change to TruthPV_mostNpart_z 
            hM_DiffVtxZ_NhitsLayer1_Nvtx1.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1) # change to TruthPV_mostNpart_z 
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1) # change to TruthPV_mostNpart_z 
            hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2.Fill((finalPVz - tree.TruthPV_mostNpart_z), tree.NhitsLayer1) # change to TruthPV_mostNpart_z 
            hM_TruthPVz_RecoPVz_Nvtx1.Fill(tree.TruthPV_mostNpart_z, finalPVz) # change to TruthPV_mostNpart_z 
            hM_TruthPVx_RecoPVx_Nvtx1.Fill(tree.TruthPV_mostNpart_x, finalPVx) # change to TruthPV_mostNpart_x
            hM_TruthPVy_RecoPVy_Nvtx1.Fill(tree.TruthPV_mostNpart_y, finalPVy) # change to TruthPV_mostNpart_y

            # if tree.PV_z > -999:
            #     hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin.Fill(tree.NhitsLayer1)
            #     if tree.PV_z >= -10. and tree.PV_z <= 10.:
            #         hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin.Fill(tree.NhitsLayer1)

            if abs(tree.TruthPV_mostNpart_z - finalPVz) > 1 and finalPVz > -99: # change to TruthPV_mostNpart_z 
                levt_NotPropRecoPV.append(tree.event)
                lNClusL1_NotPropRecoPV.append(tree.NhitsLayer1)
                print ('--->Reco PV not properly reconstructed: event={},TruthPV_trig_z(TruthPV_mostNpart_z)={:.3f},PV_z={:.3f},NhitsLayer1={}'.format(tree.event, tree.TruthPV_mostNpart_z, finalPVz, tree.NhitsLayer1)) # change to TruthPV_mostNpart_z 
        
        for i, diff in enumerate(list_absdiff):
            if abs(tree.TruthPV_trig_z - finalPVz) > diff:
                list_count_absdiff_trig[i] += 1
            
            if abs(tree.TruthPV_mostNpart_z - finalPVz) > diff:
                list_count_absdiff_MostNpart[i] += 1

            if tree.NTruthVtx == 1 and abs(tree.TruthPV_trig_z - finalPVz) > diff:
                list_count_absdiff_Nvtx1[i] += 1

    list_count_absdiff_trig = array('d', [i/tree.GetEntries() for i in list_count_absdiff_trig])
    list_count_absdiff_MostNpart = array('d', [i/tree.GetEntries() for i in list_count_absdiff_MostNpart])
    list_count_absdiff_Nvtx1 = array('d', [i/Nevt_Nvtx1 for i in list_count_absdiff_MostNpart])
    # print (list_count_absdiff_trig)

    err_RecoPVEff_NClusLayer1_NTruthVtx1 = TGraphAsymmErrors()
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV = TGraphAsymmErrors()
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm = TGraphAsymmErrors()
    err_EventSelEff_NClusLayer1_absPVle10cm = TGraphAsymmErrors() # numerator: events with |PV_z| < 10cm; denominator: events with N_TruthVtx = 1
    # err_RecoPVEff_NClusLayer1_NTruthVtx1.BayesDivide(hM_NClusLayer1_NTurhtVtx1_PervBin, hM_NClusLayer1_all_PercBin)
    err_RecoPVEff_NClusLayer1_NTruthVtx1.BayesDivide(hM_NClusLayer1_NTurhtVtx1_PervBin, hM_NClusLayer1_NTurhtVtx1_PervBin)
    # err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin, hM_NClusLayer1_all_PercBin)
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin, hM_NClusLayer1_NTurhtVtx1_PervBin)
    # err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin, hM_NClusLayer1_all_PercBin)
    err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin, hM_NClusLayer1_NTurhtVtx1_PervBin)
    err_EventSelEff_NClusLayer1_absPVle10cm.BayesDivide(hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin, hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin)
    Draw_1DEff(err_RecoPVEff_NClusLayer1_NTruthVtx1, True, 'Number of clusters (Layer 0)', [0,10100,0,1.35], plotpath+'RecoPVEff_NClusLayer1_NTruthVtx1')
    Draw_1DEff(err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV, True, 'Number of clusters (Layer 0)', [0,10100,0,1.35], plotpath+'RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV')
    Draw_1DEff(err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm, True, 'Number of clusters (Layer 0)', [0,10100,0,1.35], plotpath+'RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm')
    Draw_1DEff(err_EventSelEff_NClusLayer1_absPVle10cm, True, 'Number of clusters (Layer 0)', [0,10100,0,1.35], plotpath+'EventSelEff_NClusLayer1_NTruthVtx1_absPVle10cm')
    list_eff = [err_RecoPVEff_NClusLayer1_NTruthVtx1, err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV, err_RecoPVEff_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm]
    list_color = ['#03001C', '#035397', '#9B0000']
    list_leg = ['N_{truth vtx}=1', 'N_{truth vtx}=1 & Has PV_{Reco}', 'N_{truth vtx}=1 & Has PV_{Reco} & |PV^{Reco}_{z}|<10cm']
    Draw_1DEffComp(list_eff, list_color, True, 'Number of clusters (Layer 0)', list_leg, [0,10100,0,1.4], plotpath+'RecoPVEff_NClusLayer1_EffComp')
    
    fout = TFile(plotpath+'EventSelEff.root','recreate')
    fout.cd()
    err_EventSelEff_NClusLayer1_absPVle10cm.Write()
    fout.Close()

    Draw_1Dhist(hM_PVz, False, False, 1.3, 'Primary vertex V_{z} [cm]', 'cm', plotpath+'PVz_woCuts')
    Draw_1Dhist(hM_TruthTrigPVz, False, False, 1.3, 'Truth vertex V_{z} [cm]', 'cm', plotpath+'TruthTrigPVz_woCuts')
    Draw_1Dhist(hM_TrithTrigPVz_Nvtx1, False, False, 1.3, 'Truth vertex V_{z} [cm]', 'cm', plotpath+'TruthTrigPVz_Nvtx1')
    # Draw_1DhistsComp2(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname)
    Draw_1DhistsComp2([hM_PVz, hM_TruthTrigPVz], False, False, False, 1.2, 'Vertex V_{z} [cm]', 'cm', ['Reco PV', 'Truth PV (triggered)'], plotpath+'RecoTruthPVzComp_woCuts')
    Draw_1Dhist(hM_NTruthVtx, False, False, 1.2, 'Number of truth vertices', '', plotpath+'RecoClusters_NTruthVtx')
    # Draw_1Dhist(hM_DiffVtxZ_trig, False, True, 150, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_trig')
    # Draw_1Dhist(hM_DiffVtxZ_mostNpart, False, True, 150, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_mostNpart')
    # Draw_1Dhist(hM_DiffVtxZ_Nvtx1, False, True, 150, '#Deltaz(PV_{Truth}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1')
    Draw_1Dhist_wTF1(hM_DiffVtxZ_trig_altrange, False, True, 150, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_trig_altrange')
    Draw_1Dhist_wTF1(hM_DiffVtxZ_mostNpart_altrange, False, True, 150, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_mostNpart_altrange')
    Draw_1Dhist_wTF1(hM_DiffVtxZ_Nvtx1_altrange, False, True, 50, '#Deltaz(PV_{Truth}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1_altrange')
    l_res_vtxX = []
    l_reserr_vtxX = []
    l_res_vtxY = []
    l_reserr_vtxY = []
    l_res_vtxZ = []
    l_reserr_vtxZ = []
    for i, hist in enumerate(list_hM_DiffVtxZ_Nvtx1):
        mean, meanerr, sigma, sigmaerr = Draw_1Dhist_wTF1(hist, False, True, 50, '#Deltaz(PV_{Truth}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxZ_Nvtx1_NhitsBin{}'.format(i+1))
        l_res_vtxZ.append(sigma)
        l_reserr_vtxZ.append(sigmaerr)
        mean, meanerr, sigma, sigmaerr = Draw_1Dhist_wTF1(list_hM_DiffVtxX_Nvtx1[i], False, True, 50, '#Deltax(PV_{Truth}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxX_Nvtx1_NhitsBin{}'.format(i+1))
        l_res_vtxX.append(sigma)
        l_reserr_vtxX.append(sigmaerr)
        mean, meanerr, sigma, sigmaerr = Draw_1Dhist_wTF1(list_hM_DiffVtxY_Nvtx1[i], False, True, 50, '#Deltay(PV_{Truth}, PV_{Reco}) [cm]', 'cm', plotpath+'RecoClusters_DiffVtxY_Nvtx1_NhitsBin{}'.format(i+1))
        l_res_vtxY.append(sigma)
        l_reserr_vtxY.append(sigmaerr)


    list_res_vtxX = array('d', l_res_vtxX)
    list_reserr_vtxX = array('d', l_reserr_vtxX)
    Draw_GraphError(array('d', [i+1 for i in range(len(list_res_vtxX))]), array('d', [0 for i in range(len(list_res_vtxX))]), list_res_vtxX, list_reserr_vtxX, 0.17, 'MVTX cluster-multiplicity bin', 'Vertex x-position resolution [cm]', 1.7, [0.5, len(list_res_vtxX)+0.5], [0, 0.042], plotpath+'VtxXReolution_Bin')
    list_res_vtxY = array('d', l_res_vtxY)
    list_reserr_vtxY = array('d', l_reserr_vtxY)
    Draw_GraphError(array('d', [i+1 for i in range(len(list_res_vtxY))]), array('d', [0 for i in range(len(list_res_vtxY))]), list_res_vtxY, list_reserr_vtxY, 0.17, 'MVTX cluster-multiplicity bin', 'Vertex y-position resolution [cm]', 1.7, [0.5, len(list_res_vtxY)+0.5], [0, 0.042], plotpath+'VtxYReolution_Bin')
    list_res_vtxZ = array('d', l_res_vtxZ)
    list_reserr_vtxZ = array('d', l_reserr_vtxZ)
    Draw_GraphError(array('d', [i+1 for i in range(len(list_res_vtxZ))]), array('d', [0 for i in range(len(list_res_vtxZ))]), list_res_vtxZ, list_reserr_vtxZ, 0.17, 'MVTX cluster-multiplicity bin', 'Vertex z-position resolution [cm]', 1.7, [0.5, len(list_res_vtxZ)+0.5], [0, 0.02], plotpath+'VtxZReolution_Bin')

    Draw_1Dhist_wPercentile(hM_NClusLayer1_NTruthVtx1, NhitsL1_percentile, False, False, True, 5, 'Number of clusters (Layer 0)', '', plotpath+'NClusLayer1_NTruthVtx1_wPercentile')
    Draw_1Dhist_wPercentile(hM_NClusLayer1_NTruthVtx1_altrange1, NhitsL1_percentile, False, False, True, 5, 'Number of clusters (Layer 0)', '', plotpath+'NClusLayer1_NTruthVtx1_wPercentile_altrange1')
    Draw_1Dhist_wPercentile(hM_NClusLayer1_NTruthVtx1_logX, NhitsL1_percentile, False, True, True, 5, 'Number of clusters (Layer 0)', '', plotpath+'NClusLayer1_NTruthVtx1_wPercentile_logX')

    Draw_2Dhist(hM_TruthPVz_NhitsLayer1, False, False, 0.11, 'Truth PV (trig) Z pos [cm]', 'Number of clusters (Layer 0)', plotpath+'RecoClusters_TruthPVz_NhitsLayer1')
    Draw_2Dhist(hM_TruthPVzTrig_RecoPVz, False, False, 0.11, 'Truth PV (trig) Z pos [cm]', 'Reco PV Z pos [cm]', plotpath+'RecoClusters_TruthTrigRecoPVz')
    Draw_2Dhist(hM_TruthPVzMostNpart_RecoPVz, False, False, 0.11, 'Truth PV Z (most N_{particles}) pos [cm]', 'Reco PV Z pos [cm]', plotpath+'RecoClusters_TruthMostNpartRecoPVz')
    Draw_2Dhist(hM_TruthPVz_RecoPVz_Nvtx1, False, False, 0.11, 'Truth PV_{Z} [cm]', 'Reconstructed PV_{Z} [cm]', plotpath+'RecoClusters_TruthRecoPVz_Nvtx1')
    Draw_2Dhist(hM_TruthPVx_RecoPVx_Nvtx1, False, False, 0.11, 'Truth PV_{X} [cm]', 'Reconstructed PV_{X} [cm]', plotpath+'RecoClusters_TruthRecoPVx_Nvtx1')
    Draw_2Dhist(hM_TruthPVy_RecoPVy_Nvtx1, False, False, 0.11, 'Truth PV_{Y} [cm]', 'Reconstructed PV_{Y} [cm]', plotpath+'RecoClusters_TruthRecoPVy_Nvtx1')
    Draw_2Dhist(hM_DiffVtxZTrig_NhitsLayer1, False, False, 0.11, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZTrig_NhitsLayer1')
    Draw_2Dhist(hM_DiffVtxZMostNpart_NhitsLayer1, False, False, 0.11, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZMostNpart_NhitsLayer1')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1, False, False, 0.11, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1')
    Draw_2Dhist(hM_DiffVtxZTrig_NhitsLayer1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}^{trig}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZTrig_NhitsLayer1_altrange')
    Draw_2Dhist(hM_DiffVtxZMostNpart_NhitsLayer1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZMostNpart_NhitsLayer1_altrange')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange')
    Draw_2Dhist(hM_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2, False, False, 0.13, '#Deltaz(PV_{Truth}, PV_{Reco})', 'Number of clusters on 1st layer', plotpath+'RecoClusters_DiffVtxZ_NhitsLayer1_Nvtx1_altrange2')

    Draw_Graph(list_absdiff, list_count_absdiff_trig, 0.15, '#Deltaz_{Threshold} [cm]', 'Fraction (|#Deltaz(PV_{Truth}^{trig}, PV_{Reco})| > #Deltaz_{Threshold})', 1.3, [0, float(Ncuts*incre)], [0, 1.1], plotpath+'RecoClusters_TruthTrigRecoPVz_AbsDiffFrac')
    Draw_Graph(list_absdiff, list_count_absdiff_MostNpart, 0.15, '#Deltaz_{Threshold} [cm]', 'Fraction (|#Deltaz(PV_{Truth}^{most Npart}, PV_{Reco})| > #Deltaz_{Threshold})', 1.3, [0, float(Ncuts*incre)], [0, 1.1], plotpath+'RecoClusters_TruthMostNpartRecoPVz_AbsDiffFrac')
    Draw_Graph(list_absdiff, list_count_absdiff_Nvtx1, 0.2, '#Deltaz_{Threshold} [cm]', 'Fraction (|#Deltaz(PV_{Truth}, PV_{Reco})| > #Deltaz_{Threshold})', 2, [0, float(Ncuts*incre)], [0, 0.005], plotpath+'RecoClusters_TruthRecoPVz_Nvtx1_AbsDiffFrac')

    print(hM_NClusLayer1_all_PercBin.Integral(-1,-1), hM_NClusLayer1_NTruthVtx1_HasRecoPV_PercBin.Integral(-1,-1), hM_NClusLayer1_NTruthVtx1_HasRecoPV_absPVle10cm_PercBin.Integral(-1,-1))
    print('Number of events where PV is not reconstructed =', Nevt_NotRecoPV)
    print('Number of events where |truth PV (trig) z| < 10 cm =', Nevt_TruthPVless10)
    print('Number of events where |truth PV (trig) z| < 10 cm but not reconstructed =', Nevt_NotRecoPV_TruthPVless10, '; list of events =', levt_NotRecoPV, '; list of NClusNlayer1 =', lNClusL1_NotRecoPV)
    print('Number of events where PV is reconstructed but abs(tree.TruthPV_trig_z - tree.PV_z) > 1 =', len(levt_NotPropRecoPV), '; list of events =', levt_NotPropRecoPV, '; list of NClusNlayer1 =', lNClusL1_NotPropRecoPV)