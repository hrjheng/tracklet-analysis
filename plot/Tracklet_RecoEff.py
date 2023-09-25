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

def Draw_simpleEff(eff, xrange, yrange, xtitle, ytitle, outname):
    c = TCanvas('c', 'c', 800, 700)
    c.cd()
    eff.GetXaxis().SetTitle(xtitle)
    eff.GetYaxis().SetTitle(ytitle)
    eff.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    eff.GetYaxis().SetRangeUser(yrange[0], yrange[1])
    eff.SetMarkerColor(1)
    eff.SetMarkerSize(0.9)
    eff.SetMarkerStyle(20)
    eff.Draw('AP')
    c.Draw()
    c.SaveAs(outname+'.pdf')
    c.SaveAs(outname+'.png')
    if(c):
        c.Close()
        gSystem.ProcessEvents()
        del c
        c = 0

def Draw_EffComp(effs, colors, legtext, xrange, yrange, xtitle, ytitle, legpos, outname):
    c = TCanvas('c', 'c', 800, 700)
    c.cd()
    effs[0].GetXaxis().SetTitle(xtitle)
    effs[0].GetYaxis().SetTitle(ytitle)
    effs[0].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    effs[0].GetYaxis().SetRangeUser(yrange[0], yrange[1])
    for i in range(len(effs)):
        effs[i].SetMarkerColor(TColor.GetColor(colors[i]))
        effs[i].SetLineColor(TColor.GetColor(colors[i]))
        effs[i].SetMarkerSize(0.9)
        effs[i].SetMarkerStyle(20)
        if i == 0:
            effs[i].Draw('AP')
        else:
            effs[i].Draw('Psame')

    leg = TLegend(legpos[0], legpos[1], legpos[2], legpos[3])
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    for i, t in enumerate(effs):
        leg.AddEntry(t, legtext[i], 'pl')
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

def GetRecoEff_TGAError(drtext, layer):
    hM_genpartpt_all_AuAuMB_pt0to1 = TH1F('hM_genpartpt_all_AuAuMB_pt0to1_dr{}'.format(drtext), 'hM_genpartpt_all_AuAuMB_pt0to1_dr{}'.format(drtext), 50, 0, 1)
    hM_genpartpt_matched_AuAuMB_pt0to1 = TH1F('hM_genpartpt_matched_AuAuMB_pt0to1_dr{}'.format(drtext), 'hM_genpartpt_matched_AuAuMB_pt0to1_dr{}'.format(drtext), 50, 0, 1)
    hM_genpartpt_all_AuAuMB_pt0to5 = TH1F('hM_genpartpt_all_AuAuMB_pt0to5_dr{}'.format(drtext), 'hM_genpartpt_all_AuAuMB_pt0to5_dr{}'.format(drtext), 50, 0, 5)
    hM_genpartpt_matched_AuAuMB_pt0to5 = TH1F('hM_genpartpt_matched_AuAuMB_pt0to5_dr{}'.format(drtext), 'hM_genpartpt_matched_AuAuMB_pt0to5_dr{}.format(drtext)', 50, 0, 5)
    hM_genpartpt_all_AuAuMB_pt0to0p5 = TH1F('hM_genpartpt_all_AuAuMB_pt0to0p5_dr{}'.format(drtext), 'hM_genpartpt_all_AuAuMB_pt0to0p5_dr{}'.format(drtext), 50, 0, 0.5)
    hM_genpartpt_matched_AuAuMB_pt0to0p5 = TH1F('hM_genpartpt_matched_AuAuMB_pt0to0p5_dr{}'.format(drtext), 'hM_genpartpt_matched_AuAuMB_pt0to0p5_dr{}'.format(drtext), 50, 0, 0.5)
    hM_genpartEta_all_AuAuMB = TH1F('hM_genpartEta_all_AuAuMB_dr{}'.format(drtext), 'hM_genpartEta_all_AuAuMB_dr{}'.format(drtext), 120, -3, 3)
    hM_genpartEta_matched_AuAuMB = TH1F('hM_genpartEta_matched_AuAuMB_dr{}'.format(drtext), 'hM_genpartEta_matched_AuAuMB_dr{}'.format(drtext), 120, -3, 3)
    hM_genpartPhi_all_AuAuMB = TH1F('hM_genpartPhi_all_AuAuMB_dr{}'.format(drtext), 'hM_genpartPhi_all_AuAuMB_dr{}'.format(drtext), 66, -3.3, 3.3)
    hM_genpartPhi_matched_AuAuMB = TH1F('hM_genpartPhi_matched_AuAuMB_dr{}'.format(drtext), 'hM_genpartPhi_matched_AuAuMB_dr{}'.format(drtext), 66, -3.3, 3.3)

    f = TFile('/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer{}_Evt0to500_RandhitCase0_MisAlignNum0_dRcut{}.root'.format(layer,drtext),'r')
    tree = f.Get('minitree_{}'.format(layer))
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)
        if len(tree.recotklgm_dR) != len(tree.GenHadron_matched_Pt):
            print('recotklgm_dR and GenHadron_matched_Pt not match!')
            continue

        for i in range(len(tree.GenHadron_Pt)):
            hM_genpartpt_all_AuAuMB_pt0to1.Fill(tree.GenHadron_Pt[i])
            hM_genpartpt_all_AuAuMB_pt0to5.Fill(tree.GenHadron_Pt[i])
            hM_genpartpt_all_AuAuMB_pt0to0p5.Fill(tree.GenHadron_Pt[i])
            hM_genpartEta_all_AuAuMB.Fill(tree.GenHadron_eta[i])
            hM_genpartPhi_all_AuAuMB.Fill(tree.GenHadron_phi[i])
        
        for i in range(len(tree.GenHadron_matched_Pt)):
            hM_genpartpt_matched_AuAuMB_pt0to1.Fill(tree.GenHadron_matched_Pt[i])
            hM_genpartpt_matched_AuAuMB_pt0to5.Fill(tree.GenHadron_matched_Pt[i])
            hM_genpartpt_matched_AuAuMB_pt0to0p5.Fill(tree.GenHadron_matched_Pt[i])
            hM_genpartEta_matched_AuAuMB.Fill(tree.GenHadron_matched_eta[i])
            hM_genpartPhi_matched_AuAuMB.Fill(tree.GenHadron_matched_phi[i])

    f.Close()

    print(drtext, hM_genpartpt_matched_AuAuMB_pt0to1.Integral(-1,-1), hM_genpartpt_all_AuAuMB_pt0to1.Integral(-1,-1), hM_genpartpt_matched_AuAuMB_pt0to1.Integral(-1,-1)/hM_genpartpt_all_AuAuMB_pt0to1.Integral(-1,-1))

    err_TklRecoEff_pt0to1 = TGraphAsymmErrors()
    err_TklRecoEff_pt0to1.BayesDivide(hM_genpartpt_matched_AuAuMB_pt0to1, hM_genpartpt_all_AuAuMB_pt0to1)
    err_TklRecoEff_pt0to5 = TGraphAsymmErrors()
    err_TklRecoEff_pt0to5.BayesDivide(hM_genpartpt_matched_AuAuMB_pt0to5, hM_genpartpt_all_AuAuMB_pt0to5)
    err_TklRecoEff_pt0to0p5 = TGraphAsymmErrors()
    err_TklRecoEff_pt0to0p5.BayesDivide(hM_genpartpt_matched_AuAuMB_pt0to0p5, hM_genpartpt_all_AuAuMB_pt0to0p5)
    err_TklRecoEff_Eta = TGraphAsymmErrors()
    err_TklRecoEff_Eta.BayesDivide(hM_genpartEta_matched_AuAuMB, hM_genpartEta_all_AuAuMB)
    err_TklRecoEff_Phi = TGraphAsymmErrors()
    err_TklRecoEff_Phi.BayesDivide(hM_genpartPhi_matched_AuAuMB, hM_genpartPhi_all_AuAuMB)

    return err_TklRecoEff_pt0to1, err_TklRecoEff_pt0to5, err_TklRecoEff_pt0to0p5, err_TklRecoEff_Eta, err_TklRecoEff_Phi

if __name__ == '__main__':
    plotpath = './Traklet_RecoEff/'
    os.makedirs(plotpath, exist_ok=True)

    err_TklRecoEff_pt0to1_dr0p4, err_TklRecoEff_pt0to5_dr0p4, err_TklRecoEff_pt0to0p5_dr0p4, err_TklRecoEff_Eta_dr0p4, err_TklRecoEff_Phi_dr0p4 = GetRecoEff_TGAError('0p4', 12)
    err_TklRecoEff_pt0to1_dr0p5, err_TklRecoEff_pt0to5_dr0p5, err_TklRecoEff_pt0to0p5_dr0p5, err_TklRecoEff_Eta_dr0p5, err_TklRecoEff_Phi_dr0p5 = GetRecoEff_TGAError('0p5', 12)
    err_TklRecoEff_pt0to1_dr0p6, err_TklRecoEff_pt0to5_dr0p6, err_TklRecoEff_pt0to0p5_dr0p6, err_TklRecoEff_Eta_dr0p6, err_TklRecoEff_Phi_dr0p6 = GetRecoEff_TGAError('0p6', 12)
    err_TklRecoEff_pt0to1_dr999, err_TklRecoEff_pt0to5_dr999, err_TklRecoEff_pt0to0p5_dr999, err_TklRecoEff_Eta_dr999, err_TklRecoEff_Phi_dr999 = GetRecoEff_TGAError('999', 12)

    Draw_simpleEff(err_TklRecoEff_pt0to1_dr0p4, [0,1.], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to1_dr0p4')
    Draw_simpleEff(err_TklRecoEff_pt0to5_dr0p4, [0,5.], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to5_dr0p4')
    Draw_simpleEff(err_TklRecoEff_pt0to0p5_dr0p4, [0,0.5], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to0p5_dr0p4')
    Draw_simpleEff(err_TklRecoEff_Eta_dr0p4, [-3,3], [0,1.05], '#eta', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartEta_dr0p4')
    Draw_simpleEff(err_TklRecoEff_Phi_dr0p4, [-3.3,3.3], [0,1.05], '#phi', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPhi_dr0p4')

    Draw_simpleEff(err_TklRecoEff_pt0to1_dr0p5, [0,1.], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to1_dr0p5')
    Draw_simpleEff(err_TklRecoEff_pt0to5_dr0p5, [0,5.], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to5_dr0p5')
    Draw_simpleEff(err_TklRecoEff_pt0to0p5_dr0p5, [0,0.5], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to0p5_dr0p5')
    Draw_simpleEff(err_TklRecoEff_Eta_dr0p5, [-3,3], [0,1.05], '#eta', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartEta_dr0p5')
    Draw_simpleEff(err_TklRecoEff_Phi_dr0p5, [-3.3,3.3], [0,1.05], '#phi', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPhi_dr0p5')

    Draw_simpleEff(err_TklRecoEff_pt0to1_dr0p6, [0,1.], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to1_dr0p6')
    Draw_simpleEff(err_TklRecoEff_pt0to5_dr0p6, [0,5.], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to5_dr0p6')
    Draw_simpleEff(err_TklRecoEff_pt0to0p5_dr0p6, [0,0.5], [0,1.05], 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPt_pt0to0p5_dr0p6')
    Draw_simpleEff(err_TklRecoEff_Eta_dr0p6, [-3,3], [0,1.05], '#eta', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartEta_dr0p6')
    Draw_simpleEff(err_TklRecoEff_Phi_dr0p6, [-3.3,3.3], [0,1.05], '#phi', 'Tracklet reconstruction efficiency', plotpath + 'TrackletRecoEff_GenPartPhi_dr0p6')

    Draw_EffComp([err_TklRecoEff_pt0to1_dr0p4, err_TklRecoEff_pt0to1_dr0p5, err_TklRecoEff_pt0to1_dr0p6, err_TklRecoEff_pt0to1_dr999], 
                 ['#9B0000', '#1B1A17', '#035397', '#285430'], 
                 ['#DeltaR=0.4','#DeltaR=0.5','#DeltaR=0.6', 'w/o #DeltaR cut'], [0,1.], [0,1.05], 
                 'p_{T} (GeV)', 'Tracklet reconstruction efficiency', 
                 [0.65,0.25,0.88,0.4],
                 plotpath + 'TrackletRecoEff_GenPartPt_pt0to1_comparison')

    Draw_EffComp([err_TklRecoEff_Eta_dr0p4, err_TklRecoEff_Eta_dr0p5, err_TklRecoEff_Eta_dr0p6, err_TklRecoEff_Eta_dr999], 
                 ['#9B0000', '#1B1A17', '#035397', '#285430'], 
                 ['#DeltaR=0.4','#DeltaR=0.5','#DeltaR=0.6', 'w/o #DeltaR cut'], [-3,3], [0,1.05], 
                 '#eta', 'Tracklet reconstruction efficiency', 
                 [0.45,0.25,0.65,0.4],
                 plotpath + 'TrackletRecoEff_GenPartEta_comparison')
    
    Draw_EffComp([err_TklRecoEff_Phi_dr0p4, err_TklRecoEff_Phi_dr0p5, err_TklRecoEff_Phi_dr0p6, err_TklRecoEff_Phi_dr999],
                 ['#9B0000', '#1B1A17', '#035397', '#285430'],
                 ['#DeltaR=0.4','#DeltaR=0.5','#DeltaR=0.6', 'w/o #DeltaR cut'], [-3.3,3.3], [0,1.05],
                 '#phi', 'Tracklet reconstruction efficiency',
                 [0.45,0.25,0.65,0.4],
                 plotpath + 'TrackletRecoEff_GenPartPhi_comparison')
