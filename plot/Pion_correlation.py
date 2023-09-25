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

if __name__ == '__main__':
    plotpath = './pion_correlation/'
    os.makedirs(plotpath, exist_ok=True)
    
    hM_recotkldr_genpartpt_gm = TH2F('hM_recotkldr_genpartpt_gm', 'hM_recotkldr_genpartpt_gm', 100, 0, 0.025, 100, 0, 1)
    hM_recotkldr_raw = TH1F('hM_recotkldr_raw', 'hM_recotkldr_raw', 100, 0, 0.025)

    f_pion = TFile('/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer12_Evt0to500_RandhitCase0_MisAlignNum0_dRcut0p5.root','r')
    tree_pion = f_pion.Get('minitree_12')
    for idx in range(tree_pion.GetEntries()):
        tree_pion.GetEntry(idx)
        # print(tree_pion.event, tree_pion.NhitsLayer1, tree_pion.NGenHadron, tree_pion.CentralityBin, len(tree_pion.GenHadron_eta))
        if len(tree_pion.recotklgm_dR) != len(tree_pion.GenHadron_matched_Pt):
            print('recotklgm_dR and GenHadron_matched_Pt not match!')
            continue
        
        for i in range(len(tree_pion.recotklraw_dR)):
            hM_recotkldr_raw.Fill(tree_pion.recotklraw_dR[i])

        for i in range(len(tree_pion.recotklgm_dR)):
            hM_recotkldr_genpartpt_gm.Fill(tree_pion.recotklgm_dR[i], tree_pion.GenHadron_matched_Pt[i])

    f_pion.Close()

    Draw_2Dhist(hM_recotkldr_genpartpt_gm, False, False, 0.13, 'Reco-tracklet #DeltaR', 'Generated pion p_{T} (GeV)', 'colz', plotpath+'recotklgmdr_genpionpt_pion')

    hM_recotkldr_gm = hM_recotkldr_genpartpt_gm.ProjectionX('hM_recotkldr_gm')
    
    Draw_1Dhist(hM_recotkldr_raw, False, False, 1.3, 'Reco-tracklet #DeltaR', '', plotpath+'recotkldr_raw_pion')
    Draw_1Dhist(hM_recotkldr_gm, False, False, 1.3, 'Reco-tracklet #DeltaR', '', plotpath+'recotkldr_gm_pion')

    # Normalize the distributions to have the same value at [normatdR]
    normatdR = 0.001
    hM_recotkldr_raw.Scale(1./hM_recotkldr_raw.GetBinContent(hM_recotkldr_raw.FindBin(normatdR)))
    hM_recotkldr_gm.Scale(1./hM_recotkldr_gm.GetBinContent(hM_recotkldr_gm.FindBin(normatdR)))

    c1 = TCanvas('c1', 'c1', 800, 600)
    pad1 = TPad('pad1', 'pad1', 0, 0.3, 1, 1)
    pad2 = TPad('pad2', 'pad2', 0, 0, 1, 0.3)
    pad1.SetBottomMargin(0)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad1.Draw()
    pad2.Draw()
    pad1.cd()
    pad1.SetLogy()
    hM_recotkldr_raw.GetYaxis().SetTitleOffset(1.1)
    hM_recotkldr_raw.GetXaxis().SetTitleSize(AxisTitleSize*1.2)
    hM_recotkldr_raw.GetYaxis().SetTitleSize(AxisTitleSize*1.2)
    hM_recotkldr_raw.GetXaxis().SetLabelSize(AxisLabelSize*1.2)
    hM_recotkldr_raw.GetYaxis().SetLabelSize(AxisLabelSize*1.2)
    hM_recotkldr_raw.SetLineColor(38)
    hM_recotkldr_raw.SetLineWidth(2)
    hM_recotkldr_raw.Draw('hist')
    hM_recotkldr_gm.SetLineColor(46)
    hM_recotkldr_gm.SetLineWidth(2)
    hM_recotkldr_gm.Draw('histsame')
    leg = TLegend(0.55, 0.75, 0.9, 0.9)
    leg.AddEntry(hM_recotkldr_gm, 'Signal tracklets', 'l')
    leg.AddEntry(hM_recotkldr_raw, 'Reconstructed tracklets', 'l')
    leg.Draw()
    c1.Update()
    pad2.cd()
    hM_recotkldr_signal = hM_recotkldr_gm.Clone('hM_recotkldr_signal')
    hM_recotkldr_signal.Divide(hM_recotkldr_raw)
    hM_recotkldr_signal.GetXaxis().SetTitle('Reco-tracklet #DeltaR')
    hM_recotkldr_signal.GetYaxis().SetTitle('Ratio')
    hM_recotkldr_signal.GetXaxis().SetTitleOffset(1)
    hM_recotkldr_signal.GetYaxis().SetTitleOffset(0.5)
    hM_recotkldr_signal.GetXaxis().SetTitleSize(AxisTitleSize*2.6)
    hM_recotkldr_signal.GetYaxis().SetTitleSize(AxisTitleSize*2.6)
    hM_recotkldr_signal.GetXaxis().SetLabelSize(AxisLabelSize*2.6)
    hM_recotkldr_signal.GetYaxis().SetLabelSize(AxisLabelSize*2.6)
    hM_recotkldr_signal.SetLineColor(kBlack)
    hM_recotkldr_signal.Draw('hist')
    c1.SaveAs(plotpath+'bkgRelContribution_recotkldr.png')
    c1.SaveAs(plotpath+'bkgRelContribution_recotkldr.pdf')