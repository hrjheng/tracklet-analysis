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
# gStyle.SetPalette(kSunset)
# TColor.InvertPalette()

def Draw_1Dhists_MisalignComp(lhist, norm1, logx, logy, ymaxscale, XaxisName, Ytitle_unit, legtext, outname):
    color = []
    legtextsize = 0.035
    if len(lhist) == 3:
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

def GetHistList(histname, list_misalignnum, layercomb):
    list_hist = []
    for num in list_misalignnum:
        list_hist.append(GetHist(histname, num)[layercomb])
    
    return list_hist

if __name__ == '__main__':
    plotpath = './TrackletAna_MisAlignComp/'
    os.makedirs(plotpath, exist_ok=True)

    # layernum = 0 -> layer1+2; layernum = 1 -> layer2+3; layernum = 2 -> layer1+3
    list_legend_shiftZ = ['#DeltaZ_{layer 1}=0#mum', 
                          '#DeltaZ_{layer 1}=+10#mum', 
                          '#DeltaZ_{layer 1}=+50#mum', 
                          '#DeltaZ_{layer 1}=+100#mum']
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusZ_layer1', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Cluster Z pos (cm)', 'cm', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ClusZ_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusEta_layer1', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Cluster #eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ClusEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange2', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_altrange', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_LogX', [1, 11, 31, 51], 0), False, True, True, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Eta_proto', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Proto-tracklet #eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_ProtoTracklet_Eta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco', [1, 11, 31, 51], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange', [1, 11, 31, 51], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange2', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco', [1, 11, 31, 51], 0), False, False, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_altrange', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_LogX', [1, 11, 31, 51], 0), False, True, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Eta_reco', [1, 11, 31, 51], 0), False, False, False, 1.5, 'Reco-tracklet #eta', '', list_legend_shiftZ, plotpath+'MisalignComp_shiftZ_RecoTracklet_Eta_layer1')

    list_legend_shiftR_diffangle = ['#Delta r_{layer 1}=0#mum, #psi=0', 
                                    '#Delta r_{layer 1}=50#mum, #psi=0', 
                                    '#Delta r_{layer 1}=50#mum, #psi=#pi/6', 
                                    '#Delta r_{layer 1}=50#mum, #psi=#pi/4', 
                                    '#Delta r_{layer 1}=50#mum, #psi=#pi/3', 
                                    '#Delta r_{layer 1}=50#mum, #psi=#pi/2']
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusX_layer1', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Cluster X pos (cm)', 'cm', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ClusX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusY_layer1', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Cluster Y pos (cm)', 'cm', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ClusY_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusPhi_layer1', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Cluster #phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ClusPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto_altrange', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange2', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_altrange', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_LogX', [1, 151, 161, 171, 181, 191], 0), False, True, True, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_proto', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Proto-tracklet #phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_ProtoTracklet_Phi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco', [1, 151, 161, 171, 181, 191], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco_altrange', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco', [1, 151, 161, 171, 181, 191], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange', [1, 151, 161, 171, 181, 191], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange2', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco', [1, 151, 161, 171, 181, 191], 0), False, False, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_altrange', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_LogX', [1, 151, 161, 171, 181, 191], 0), False, True, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_reco', [1, 151, 161, 171, 181, 191], 0), False, False, False, 1.5, 'Reco-tracklet #phi', '', list_legend_shiftR_diffangle, plotpath+'MisalignComp_shiftR_diffangle_RecoTracklet_Phi_layer1')

    list_legend_shiftR_diffangle_v2 = ['#Delta r_{layer 1}=0#mum, #psi=0', 
                                       '#Delta r_{layer 1}=50#mum, #psi=#pi/4', 
                                       '#Delta r_{layer 1}=50#mum, #psi=3#pi/4',
                                       '#Delta r_{layer 1}=50#mum, #psi=5#pi/4',
                                       '#Delta r_{layer 1}=50#mum, #psi=7#pi/4']
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusX_layer1', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Cluster X pos (cm)', 'cm', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ClusX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusY_layer1', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Cluster Y pos (cm)', 'cm', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ClusY_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusPhi_layer1', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Cluster #phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ClusPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto_altrange', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange2', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_altrange', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_LogX', [1, 171, 201, 211, 221], 0), False, True, True, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_proto', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Proto-tracklet #phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_ProtoTracklet_Phi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco', [1, 171, 201, 211, 221], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco_altrange', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco', [1, 171, 201, 211, 221], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange', [1, 171, 201, 211, 221], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange2', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco', [1, 171, 201, 211, 221], 0), False, False, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_altrange', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_LogX', [1, 171, 201, 211, 221], 0), False, True, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_reco', [1, 171, 201, 211, 221], 0), False, False, False, 1.5, 'Reco-tracklet #phi', '', list_legend_shiftR_diffangle_v2, plotpath+'MisalignComp_shiftR_diffangle_v2_RecoTracklet_Phi_layer1')

    list_legend_shiftR_diffR = ['#Delta r_{layer 1}=0#mum, #psi=0', 
                                '#Delta r_{layer 1}=10#mum, #psi=#pi/4', 
                                '#Delta r_{layer 1}=50#mum, #psi=#pi/4', 
                                '#Delta r_{layer 1}=100#mum, #psi=#pi/4']
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusX_layer1', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Cluster X pos (cm)', 'cm', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ClusX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusY_layer1', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Cluster Y pos (cm)', 'cm', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ClusY_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusPhi_layer1', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Cluster #phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ClusPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto_altrange', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange2', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_altrange', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_LogX', [1, 91, 171, 251], 0), False, True, True, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_proto', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Proto-tracklet #phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_ProtoTracklet_Phi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco', [1, 91, 171, 251], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco_altrange', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco', [1, 91, 171, 251], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange', [1, 91, 171, 251], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange2', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco', [1, 91, 171, 251], 0), False, False, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_altrange', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_LogX', [1, 91, 171, 251], 0), False, True, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_reco', [1, 91, 171, 251], 0), False, False, False, 1.5, 'Reco-tracklet #phi', '', list_legend_shiftR_diffR, plotpath+'MisalignComp_shiftR_diffR_RecoTracklet_Phi_layer1')

    list_legend_shiftZandR = ['#DeltaZ_{l1}=10#mum, #Delta r_{l1}=0cm',
                              '#DeltaZ_{l1}=10#mum, #Delta r_{l1}=10#mum',
                              '#DeltaZ_{l1}=10#mum, #Delta r_{l1}=50#mum',
                              '#DeltaZ_{l1}=10#mum, #Delta r_{l1}=100#mum']
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusX_layer1', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Cluster X pos (cm)', 'cm', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ClusX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusY_layer1', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Cluster Y pos (cm)', 'cm', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ClusY_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_ClusPhi_layer1', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Cluster #phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ClusPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_proto_altrange', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_proto_altrange2', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #Delta#eta', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_altrange', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_proto_LogX', [11, 311, 321, 331], 0), False, True, True, 1.5, 'Proto-tracklet #DeltaR', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_proto', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Proto-tracklet #phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_ProtoTracklet_Phi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco', [11, 311, 321, 331], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dPhi_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dPhi_reco_altrange', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dPhi_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco', [11, 311, 321, 331], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dEta_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange', [11, 311, 321, 331], 0), False, False, True, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dEta_altrange1_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dEta_reco_altrange2', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Reco-tracklet #Delta#eta', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dEta_altrange2_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco', [11, 311, 321, 331], 0), False, False, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dR_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_altrange', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dR_altrange_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_dR_reco_LogX', [11, 311, 321, 331], 0), False, True, True, 1.5, 'Reco-tracklet #DeltaR', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_dR_LogX_layer1')
    Draw_1Dhists_MisalignComp(GetHistList('hM_Phi_reco', [11, 311, 321, 331], 0), False, False, False, 1.5, 'Reco-tracklet #phi', '', list_legend_shiftZandR, plotpath+'MisalignComp_shiftZandR_RecoTracklet_Phi_layer1')