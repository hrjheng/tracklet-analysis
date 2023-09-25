import os
import datetime
from array import *
from ROOT import *
import numpy
import math
import glob
from plotUtil import *
import itertools

gROOT.LoadMacro('./sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)

def getHistToComp(lfname, histname):
    hM_nominal = GetHistogram(lfname[0], histname)
    hM_misalign = GetHistogram(lfname[1], histname)
    hM_misalign_extreme = GetHistogram(lfname[2], histname)
    return hM_nominal, hM_misalign, hM_misalign_extreme

def Draw_1Dhists_MisalignComp(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname):
    color = []
    legtextsize = 0.035
    if len(lhist) == 2:
        legtextsize = 0.035
        color = ['#1B1A17','#9B0000']
    elif len(lhist) == 3:
        legtextsize = 0.035
        color = ['#1B1A17','#9B0000','#035397']
    elif len(lhist) == 4:
        legtextsize = 0.035
        color = ['#1B1A17','#9B0000','#035397','#377D71']
    elif len(lhist) == 5:
        legtextsize = 0.03
        color = ['#1B1A17','#9B0000','#035397','#377D71','#FF9F29']
    elif len(lhist) == 6:
        legtextsize = 0.03
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
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    leg.AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "")
    leg.Draw()

    leg1 = TLegend(LeftMargin+0.03, (1-TopMargin)-0.21,
                   LeftMargin+0.33, (1-TopMargin)-0.03)
    leg1.SetTextSize()
    leg1.SetFillStyle(0)
    for i, h in enumerate(lhist):
        leg1.AddEntry(h, legtext[i], "l")
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

def Draw_1Dhists_MisalignComp_tracklet(lhist_nominal, lhist_misalign, lhist_misalign_extreme, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname):
    if len(lhist_nominal) != len(lhist_misalign) or len(lhist_nominal) != len(lhist_misalign_extreme):
        print('Error: len(lhist_nominal) != len(lhist_misalign) or len(lhist_nominal) != len(lhist_misalign_extreme)')
        print (len(lhist_nominal), len(lhist_misalign), len(lhist_misalign_extreme))
        return
    
    color = ['#1B1A17', '#035397', '#9B0000']

    legtextsize = 0.035
    
    # legtext = ['1st+2nd layers', '2nd+3rd layers', '1st+3rd layers']
    ymax = -1
    ymin = 10e10
    for h in itertools.chain(lhist_nominal, lhist_misalign):
        h.Sumw2()
        if h.GetMaximum() > ymax:
            ymax = h.GetMaximum()
        if h.GetMinimum(0) < ymin:
            ymin = h.GetMinimum(0)
        if norm1:
            h.Scale(1. / h.Integral(-1, -1))

    binwidth_bin1 = lhist_nominal[0].GetXaxis().GetBinWidth(1)
    binwidth_bin2 = lhist_nominal[0].GetXaxis().GetBinWidth(2)
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
                lhist_nominal[0].GetYaxis().SetTitle('Normalized entries / ({:g})'.format(binwidth_bin1))
            else:
                lhist_nominal[0].GetYaxis().SetTitle('Normalized entries / ({:g} {unit})'.format(binwidth_bin1, unit=Ytitle_unit))
        else:
            if Ytitle_unit == '':
                lhist_nominal[0].GetYaxis().SetTitle('Entries / ({:g})'.format(binwidth_bin1))
            else:
                lhist_nominal[0].GetYaxis().SetTitle('Entries / ({:g} {unit})'.format(binwidth_bin1, unit=Ytitle_unit))
    else:
        if norm1:
            if Ytitle_unit == '':
                lhist_nominal[0].GetYaxis().SetTitle('Normalized entries')
            else:
                lhist_nominal[0].GetYaxis().SetTitle('Normalized entries {unit})'.format(unit=Ytitle_unit))
        else:
            if Ytitle_unit == '':
                lhist_nominal[0].GetYaxis().SetTitle('Entries')
            else:
                lhist_nominal[0].GetYaxis().SetTitle('Entries {unit}'.format(unit=Ytitle_unit))
    if logy:
        lhist_nominal[0].GetYaxis().SetRangeUser(ymin * 0.05, ymax * 100)
    else:
        lhist_nominal[0].GetYaxis().SetRangeUser(0., ymax * ymaxscale)
    lhist_nominal[0].GetXaxis().SetTitle(XaxisName)
    lhist_nominal[0].GetXaxis().SetTickSize(TickSize)
    lhist_nominal[0].GetXaxis().SetTitleSize(AxisTitleSize)
    lhist_nominal[0].GetXaxis().SetLabelSize(AxisLabelSize)
    lhist_nominal[0].GetYaxis().SetTickSize(TickSize)
    lhist_nominal[0].GetYaxis().SetTitleSize(AxisTitleSize)
    lhist_nominal[0].GetYaxis().SetLabelSize(AxisLabelSize)
    lhist_nominal[0].GetXaxis().SetTitleOffset(1.1)
    lhist_nominal[0].GetYaxis().SetTitleOffset(1.4)
    for i, h in enumerate(lhist_nominal):
        if i == 0:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('hist')
        else:
            h.SetLineColor(TColor.GetColor(color[i]))
            h.SetLineWidth(2)
            h.Draw('histsame')

    for i, h in enumerate(lhist_misalign):
        h.SetLineColor(TColor.GetColor(color[i]))
        h.SetLineWidth(2)
        h.SetLineStyle(2)
        h.Draw('histsame')

    for i, h in enumerate(lhist_misalign_extreme):
        h.SetLineColor(TColor.GetColor(color[i]))
        h.SetLineWidth(2)
        h.SetLineStyle(3)
        h.Draw('histsame')
    
    leg = TLegend((1-RightMargin)-0.45, (1-TopMargin)-0.15,
                  (1-RightMargin)-0.1, (1-TopMargin)-0.03)
    leg.SetTextSize(0.045)
    leg.SetFillStyle(0)
    leg.AddEntry("", "#it{#bf{sPHENIX}} Simulation", "")
    leg.AddEntry("", "Au+Au #sqrt{s_{NN}}=200 GeV", "")
    leg.Draw()

    leg1 = TLegend(LeftMargin+0.03, (1-TopMargin)-0.15,
                   LeftMargin+0.33, (1-TopMargin)-0.05)
    leg1.SetHeader('#splitline{Solid: Nominal}{#splitline{Dashed: Misaligned}{Dotted: Misaligned (extreme)}}')
    leg1.SetTextSize(legtextsize)
    leg1.SetFillStyle(0)
    leg1.Draw()
    leg2 = TLegend(LeftMargin+0.03, (1-TopMargin)-0.25,
                   LeftMargin+0.37, (1-TopMargin)-0.16)
    leg2.SetNColumns(2)
    leg2.SetFillStyle(0)
    for i, h in enumerate(lhist_nominal):
        leg2.AddEntry(h, legtext[i], "l")
    leg2.SetTextSize(legtextsize)
    leg2.SetFillStyle(0)
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

if __name__ == '__main__':
    plotpath = './TrackletAna_MisAlignComp/'
    os.makedirs(plotpath, exist_ok=True)

    flist = [
        './hists/HijingAuAuNoPileup_ana325private/Hists_RecoClusters.root', 
        './hists/HijingAuAuNoPileup_ana325private/Hists_RecoClusters_twohalfangle_GapNorth2p9_updated.root', 
        './hists/HijingAuAuNoPileup_ana325private/Hists_RecoClusters_twohalfangle_GapNorth3p5_updated.root']
    legendentry = ['Nominal', 'Misaligned', 'Misaligned (extreme)']

    hM_ClusZ_all_nominal, hM_ClusZ_all_misalign, hM_ClusZ_all_misalign_extreme = getHistToComp(flist, 'hM_ClusZ_all')
    Draw_1Dhists_MisalignComp(lhist=[hM_ClusZ_all_nominal, hM_ClusZ_all_misalign, hM_ClusZ_all_misalign_extreme], norm1=False, logx=False, logy=False, ymaxscale=1.5, XaxisName='Cluster Z pos (cm)', Ytitle_unit='cm', legtext=legendentry, outname=plotpath+'MisalignComp_twohalfangle_ClusZPos_all')

    hM_ClusPhi_all_nominal, hM_ClusPhi_all_misalign, hM_ClusPhi_all_misalign_extreme = getHistToComp(flist, 'hM_ClusPhi_all')
    Draw_1Dhists_MisalignComp([hM_ClusPhi_all_nominal, hM_ClusPhi_all_misalign, hM_ClusPhi_all_misalign_extreme], False, False, False, 1.5, 'Cluster #phi', '', legendentry, plotpath+'MisalignComp_twohalfangle_ClusPhi_all')

    hM_ClusEtaPV_all_nominal, hM_ClusEtaPV_all_misalign, hM_ClusEtaPV_all_misalign_extreme = getHistToComp(flist, 'hM_ClusEtaPV_all')
    Draw_1Dhists_MisalignComp([hM_ClusEtaPV_all_nominal, hM_ClusEtaPV_all_misalign, hM_ClusEtaPV_all_misalign_extreme], False, False, False, 1.5, 'Cluster #eta', '', legendentry, plotpath+'MisalignComp_twohalfangle_ClusEtaPV_all')


    f_tkl_nominal = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase0_ClusSplitCase0_MisAlignNum0_dRcut0p5.root'
    f_tkl_misalign = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase0_ClusSplitCase0_MisAlignNum100_dRcut0p5_twohalfangle_GapNorth2p9.root'
    f_tkl_misalign_extreme = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase0_ClusSplitCase0_MisAlignNum101_dRcut0p5_twohalfangle_GapNorth3p5.root'
    hM_dEta_reco_nominal = [GetHistogram(f_tkl_nominal, 'hM_dEta_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_dEta_reco_misalign = [GetHistogram(f_tkl_misalign, 'hM_dEta_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_dEta_reco_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_dEta_reco_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_dEta_reco_nominal, lhist_misalign=hM_dEta_reco_misalign, lhist_misalign_extreme=hM_dEta_reco_misalign_extreme, norm1=False, logx=False, logy=True, ymaxscale=150, XaxisName='Reco-tracklet #Delta#eta', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_dEta_reco')

    hM_dEta_reco_altrange_nominal = [GetHistogram(f_tkl_nominal, 'hM_dEta_reco_altrange_layer{}'.format(l)) for l in [12,23,13]]
    hM_dEta_reco_altrange_misalign = [GetHistogram(f_tkl_misalign, 'hM_dEta_reco_altrange_layer{}'.format(l)) for l in [12,23,13]]
    hM_dEta_reco_altrange_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_dEta_reco_altrange_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_dEta_reco_altrange_nominal, lhist_misalign=hM_dEta_reco_altrange_misalign, lhist_misalign_extreme=hM_dEta_reco_altrange_misalign_extreme, norm1=False, logx=False, logy=True, ymaxscale=150, XaxisName='Reco-tracklet #Delta#eta', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_dEta_reco_altrange')

    hM_dPhi_reco_nominal = [GetHistogram(f_tkl_nominal, 'hM_dPhi_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_dPhi_reco_misalign = [GetHistogram(f_tkl_misalign, 'hM_dPhi_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_dPhi_reco_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_dPhi_reco_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_dPhi_reco_nominal, lhist_misalign=hM_dPhi_reco_misalign, lhist_misalign_extreme=hM_dPhi_reco_misalign_extreme, norm1=False, logx=False, logy=True, ymaxscale=150, XaxisName='Reco-tracklet #Delta#phi', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_dPhi_reco')

    hM_dPhi_reco_altrange_nominal = [GetHistogram(f_tkl_nominal, 'hM_dPhi_reco_altrange_layer{}'.format(l)) for l in [12,23,13]]
    hM_dPhi_reco_altrange_misalign = [GetHistogram(f_tkl_misalign, 'hM_dPhi_reco_altrange_layer{}'.format(l)) for l in [12,23,13]]
    hM_dPhi_reco_altrange_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_dPhi_reco_altrange_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_dPhi_reco_altrange_nominal, lhist_misalign=hM_dPhi_reco_altrange_misalign, lhist_misalign_extreme=hM_dPhi_reco_altrange_misalign_extreme, norm1=False, logx=False, logy=True, ymaxscale=150, XaxisName='Reco-tracklet #Delta#phi', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_dPhi_reco_altrange')

    hM_dR_reco_nominal = [GetHistogram(f_tkl_nominal, 'hM_dR_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_dR_reco_misalign = [GetHistogram(f_tkl_misalign, 'hM_dR_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_dR_reco_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_dR_reco_layer{}'.format(l)) for l in [12,23,13]] 
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_dR_reco_nominal, lhist_misalign=hM_dR_reco_misalign, lhist_misalign_extreme=hM_dR_reco_misalign_extreme, norm1=False, logx=False, logy=True, ymaxscale=150, XaxisName='Reco-tracklet #DeltaR', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_dR_reco')

    hM_dR_reco_LogX_nominal = [GetHistogram(f_tkl_nominal, 'hM_dR_reco_LogX_layer{}'.format(l)) for l in [12,23,13]]
    hM_dR_reco_logX_misalign = [GetHistogram(f_tkl_misalign, 'hM_dR_reco_LogX_layer{}'.format(l)) for l in [12,23,13]]
    hM_dR_reco_logX_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_dR_reco_LogX_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_dR_reco_LogX_nominal, lhist_misalign=hM_dR_reco_logX_misalign, lhist_misalign_extreme=hM_dR_reco_logX_misalign_extreme, norm1=False, logx=True, logy=True, ymaxscale=150, XaxisName='Reco-tracklet #DeltaR', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_dR_reco_LogX')

    hM_Eta_reco_nominal = [GetHistogram(f_tkl_nominal, 'hM_Eta_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_Eta_reco_misalign = [GetHistogram(f_tkl_misalign, 'hM_Eta_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_Eta_reco_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_Eta_reco_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_Eta_reco_nominal, lhist_misalign=hM_Eta_reco_misalign, lhist_misalign_extreme=hM_Eta_reco_misalign_extreme, norm1=False, logx=False, logy=False, ymaxscale=1.5, XaxisName='Reco-tracklet #eta', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_Eta_reco')

    hM_Phi_reco_nominal = [GetHistogram(f_tkl_nominal, 'hM_Phi_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_Phi_reco_misalign = [GetHistogram(f_tkl_misalign, 'hM_Phi_reco_layer{}'.format(l)) for l in [12,23,13]]
    hM_Phi_reco_misalign_extreme = [GetHistogram(f_tkl_misalign_extreme, 'hM_Phi_reco_layer{}'.format(l)) for l in [12,23,13]]
    Draw_1Dhists_MisalignComp_tracklet(lhist_nominal=hM_Phi_reco_nominal, lhist_misalign=hM_Phi_reco_misalign, lhist_misalign_extreme=hM_Phi_reco_misalign_extreme, norm1=False, logx=False, logy=False, ymaxscale=1.5, XaxisName='Reco-tracklet #phi', Ytitle_unit='', legtext=['layer1+2','layer2+3','layer1+3'], outname=plotpath+'MisalignComp_twohalfangle_Phi_reco')