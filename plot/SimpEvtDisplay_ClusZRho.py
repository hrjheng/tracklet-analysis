#! /usr/bin/env python
from optparse import OptionParser
import sys
import os
import datetime
from array import *
from ROOT import *
import numpy as np

gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)
gStyle.SetPalette(kLightTemperature)
# TColor.InvertPalette()

TickSize = 0.03
AxisTitleSize = 0.05
AxisLabelSize = 0.04
LeftMargin = 0.15
RightMargin = 0.08
TopMargin = 0.08
BottomMargin = 0.13

# MVTX dimension
MVTXLength = 27.1 # cm 
MVTXL1_min = 2.461 
MVTXL1_mid = 2.523
MVTXL1_max = 2.793
MVTXL2_min = 3.198 
MVTXL2_mid = 3.335
MVTXL2_max = 3.625
MVTXL3_min = 3.993 
MVTXL3_mid = 4.148
MVTXL3_max = 4.426

def Draw_2Dhist_scat(hist, XaxisName, YaxisName, outname):
    c = TCanvas('c', 'c', 800, 700)
    c.cd()
    gPad.SetRightMargin(RightMargin)
    gPad.SetTopMargin(TopMargin)
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetBottomMargin(BottomMargin)
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
    hist.Draw('SCAT')

    MVTXl1_in = TLine(-(MVTXLength/2),MVTXL1_min,MVTXLength/2,MVTXL1_min)
    MVTXl1_out = TLine(-(MVTXLength/2),MVTXL1_max,MVTXLength/2,MVTXL1_max)
    MVTXl2_in = TLine(-(MVTXLength/2),MVTXL2_min,MVTXLength/2,MVTXL2_min)
    MVTXl2_out = TLine(-(MVTXLength/2),MVTXL2_max,MVTXLength/2,MVTXL2_max)
    MVTXl3_in = TLine(-(MVTXLength/2),MVTXL3_min,MVTXLength/2,MVTXL3_min)
    MVTXl3_out = TLine(-(MVTXLength/2),MVTXL3_max,MVTXLength/2,MVTXL3_max)
    MVTXl1_in.SetLineStyle(kDashed)
    MVTXl1_out.SetLineStyle(kDashed)
    MVTXl2_in.SetLineStyle(kDashed)
    MVTXl2_out.SetLineStyle(kDashed)
    MVTXl3_in.SetLineStyle(kDashed)
    MVTXl3_out.SetLineStyle(kDashed)
    MVTXl1_in.Draw('same')
    MVTXl1_out.Draw('same')
    MVTXl2_in.Draw('same')
    MVTXl2_out.Draw('same')
    MVTXl3_in.Draw('same')
    MVTXl3_out.Draw('same')

    t_l1 = TText(MVTXLength/2,MVTXL1_min,'Layer 1')
    t_l2 = TText(MVTXLength/2,MVTXL2_min,'Layer 2')
    t_l3 = TText(MVTXLength/2,MVTXL3_min,'Layer 3')
    t_l1.SetTextSize(0.04)
    t_l2.SetTextSize(0.04)
    t_l3.SetTextSize(0.04)
    t_l1.SetTextAlign(11)
    t_l2.SetTextAlign(11)
    t_l3.SetTextAlign(11)
    t_l1.Draw('same')
    t_l2.Draw('same')
    t_l3.Draw('same')

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
    plotpath = './MVTXRecoClusters_SimpEvtDisplay/'
    os.makedirs(plotpath, exist_ok=True)

    f = TFile('/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters_20221021_Nevt2000.root','r')
    tree = f.Get('EventTree')
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)

        NLayer1 = sum(map(lambda x : x == 0, tree.ClusLayer))

        if NLayer1 < 50:
            hM_ClusZ_ClusR = TH2F('hM_ClusZ_ClusR_ev{}'.format(tree.event), 'hM_ClusZ_ClusR_ev{}'.format(tree.event), 200,-20,20,100,2,5)
            for i in range(len(tree.ClusLayer)):
                hM_ClusZ_ClusR.Fill(tree.ClusZ[i],tree.ClusR[i])

            Draw_2Dhist_scat(hM_ClusZ_ClusR, 'Cluster position Z (cm)', 'Cluster radius (cm)', plotpath+'RecoClusZRho_ev{}'.format(tree.event))
            
