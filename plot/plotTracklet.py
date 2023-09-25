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
    parser = OptionParser(usage="usage: %prog ver [options -n]")
    parser.add_option("-f", "--histfile", dest="histfname", type="string", default='/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase0_ClusSplitCase0_MisAlignNum0_dRcut0p5.root', help="Input file name")
    parser.add_option("-r", "--randomhitcase", dest="randcase", type="int", default=0, help="Random hit case [0: nominal, 1: 1%, 2: 5%, 3: 10%]")
    parser.add_option("-s", "--clustersplit", dest="clustersplit", type="int", default=0, help="Cluster split [0: nominal, 1: 1%, 2: 2%, 3: 5%]")
    parser.add_option("-m", "--misalignnum", dest='misalignnum', type="int", default=0, help="Mis-alignment number (see src/misalignment.h)")
    parser.add_option("-d", "--drcut", dest='drcut', type="string", default='0p5', help="dR cut (default: 0.5)")
    (opt, args) = parser.parse_args()

    histfname = opt.histfname
    randcase = opt.randcase
    clustersplit = opt.clustersplit
    misalignnum = opt.misalignnum
    drcut = opt.drcut

    layer = [12, 23, 13]

    plotpath = './TrackletAna/TrackletAna_RandhitCase{}_ClusSplitCase{}_MisAlignNum{}_DrCut{}/'.format(randcase, clustersplit, misalignnum, drcut)
    os.makedirs(plotpath, exist_ok=True)

    # TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/Hists_Tracklets_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase,  clussplitcase, misalignnum, drcut);
    # histfname = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/Hists_Tracklets_RandhitCase{}_ClusSplitCase{}_MisAlignNum{}_dRcut{}.root'.format(randcase, clustersplit, misalignnum, drcut)

    print('Plotting: [random case, cluster split case, mis-alignment number, Tracklet reco dr cut]=[{}, {}, {}, {}]'.format(randcase, clustersplit, misalignnum, drcut))

    hM_dEta_proto = [GetHistogram(histfname, 'hM_dEta_proto_layer{}'.format(l)) for l in layer]
    hM_dEta_proto_altrange = [GetHistogram(histfname, 'hM_dEta_proto_altrange_layer{}'.format(l)) for l in layer]
    hM_dEta_proto_altrange2 = [GetHistogram(histfname, 'hM_dEta_proto_altrange2_layer{}'.format(l)) for l in layer]
    hM_dPhi_proto = [GetHistogram(histfname, 'hM_dPhi_proto_layer{}'.format(l)) for l in layer]
    hM_dPhi_proto_altrange = [GetHistogram(histfname, 'hM_dPhi_proto_altrange_layer{}'.format(l)) for l in layer]
    hM_dR_proto = [GetHistogram(histfname, 'hM_dR_proto_layer{}'.format(l)) for l in layer]
    hM_dR_proto_altrange = [GetHistogram(histfname, 'hM_dR_proto_altrange_layer{}'.format(l)) for l in layer]
    hM_dR_proto_LogX = [GetHistogram(histfname, 'hM_dR_proto_LogX_layer{}'.format(l)) for l in layer]
    hM_Eta_proto = [GetHistogram(histfname, 'hM_Eta_proto_layer{}'.format(l)) for l in layer]
    hM_Phi_proto = [GetHistogram(histfname, 'hM_Phi_proto_layer{}'.format(l)) for l in layer]
    hM_NTracklet_proto = [GetHistogram(histfname, 'hM_NTracklet_proto_layer{}'.format(l)) for l in layer]
    hM_Eta_proto_PVzRange1 = [GetHistogram(histfname, 'hM_Eta_proto_PVzRange1_layer{}'.format(l)) for l in layer]
    hM_Eta_proto_PVzRange2 = [GetHistogram(histfname, 'hM_Eta_proto_PVzRange2_layer{}'.format(l)) for l in layer]
    hM_Eta_proto_PVzRange3 = [GetHistogram(histfname, 'hM_Eta_proto_PVzRange3_layer{}'.format(l)) for l in layer]
    hM_Eta_proto_PVzRange4 = [GetHistogram(histfname, 'hM_Eta_proto_PVzRange4_layer{}'.format(l)) for l in layer]
    hM_Eta_proto_PVzRange5 = [GetHistogram(histfname, 'hM_Eta_proto_PVzRange5_layer{}'.format(l)) for l in layer]

    hM_dEta_reco = [GetHistogram(histfname, 'hM_dEta_reco_layer{}'.format(l)) for l in layer]
    hM_dEta_reco_altrange = [GetHistogram(histfname, 'hM_dEta_reco_altrange_layer{}'.format(l)) for l in layer]
    hM_dEta_reco_altrange2 = [GetHistogram(histfname, 'hM_dEta_reco_altrange2_layer{}'.format(l)) for l in layer]
    hM_dPhi_reco = [GetHistogram(histfname, 'hM_dPhi_reco_layer{}'.format(l)) for l in layer]
    hM_dPhi_reco_altrange = [GetHistogram(histfname, 'hM_dPhi_reco_altrange_layer{}'.format(l)) for l in layer]
    hM_dR_reco = [GetHistogram(histfname, 'hM_dR_reco_layer{}'.format(l)) for l in layer]
    hM_dR_reco_altrange = [GetHistogram(histfname, 'hM_dR_reco_altrange_layer{}'.format(l)) for l in layer]
    hM_dR_reco_LogX = [GetHistogram(histfname, 'hM_dR_reco_LogX_layer{}'.format(l)) for l in layer]
    hM_Eta_reco = [GetHistogram(histfname, 'hM_Eta_reco_layer{}'.format(l)) for l in layer]
    hM_Phi_reco = [GetHistogram(histfname, 'hM_Phi_reco_layer{}'.format(l)) for l in layer]
    hM_NTracklet_reco = [GetHistogram(histfname, 'hM_NTracklet_reco_layer{}'.format(l)) for l in layer]
    hM_Eta_reco_PVzRange1 = [GetHistogram(histfname, 'hM_Eta_reco_PVzRange1_layer{}'.format(l)) for l in layer]
    hM_Eta_reco_PVzRange2 = [GetHistogram(histfname, 'hM_Eta_reco_PVzRange2_layer{}'.format(l)) for l in layer]
    hM_Eta_reco_PVzRange3 = [GetHistogram(histfname, 'hM_Eta_reco_PVzRange3_layer{}'.format(l)) for l in layer]
    hM_Eta_reco_PVzRange4 = [GetHistogram(histfname, 'hM_Eta_reco_PVzRange4_layer{}'.format(l)) for l in layer]
    hM_Eta_reco_PVzRange5 = [GetHistogram(histfname, 'hM_Eta_reco_PVzRange5_layer{}'.format(l)) for l in layer]

    hM_dEta_recogmatch_ghadPtIncl = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPtIncl_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPt0to0p5 = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPt0to0p5_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPt0p5to1 = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPt0p5to1_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPt1to2 = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPt1to2_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPt2to3 = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPt2to3_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPtg3 = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPtg3_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPtIncl_altrange = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPtIncl_altrange_layer{}'.format(l)) for l in layer]
    hM_dEta_recogmatch_ghadPtIncl_altrange2 = [GetHistogram(histfname, 'hM_dEta_recogmatch_ghadPtIncl_altrange2_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPtIncl = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPtIncl_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPt0to0p5 = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPt0to0p5_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPt0p5to1 = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPt0p5to1_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPt1to2 = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPt1to2_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPt2to3 = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPt2to3_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPtg3 = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPtg3_layer{}'.format(l)) for l in layer]
    hM_dPhi_recogmatch_ghadPtIncl_altrange = [GetHistogram(histfname, 'hM_dPhi_recogmatch_ghadPtIncl_altrange_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPtIncl = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPtIncl_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt0to0p5 = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt0to0p5_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt0p5to1 = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt0p5to1_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt1to2 = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt1to2_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt2to3 = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt2to3_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPtg3 = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPtg3_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPtIncl_LogX = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPtIncl_LogX_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPtIncl_fr = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPtIncl_fr_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt0to0p5_fr = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt0to0p5_fr_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt0p5to1_fr = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt0p5to1_fr_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt1to2_fr = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt1to2_fr_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPt2to3_fr = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPt2to3_fr_layer{}'.format(l)) for l in layer]
    hM_dR_recogmatch_ghadPtg3_fr = [GetHistogram(histfname, 'hM_dR_recogmatch_ghadPtg3_fr_layer{}'.format(l)) for l in layer]
    hM_Eta_recogmatch = [GetHistogram(histfname, 'hM_Eta_recogmatch_layer{}'.format(l)) for l in layer]
    hM_Phi_recogmatch = [GetHistogram(histfname, 'hM_Phi_recogmatch_layer{}'.format(l)) for l in layer]
    hM_NTracklet_recogmatch = [GetHistogram(histfname, 'hM_NTracklet_recogmatch_layer{}'.format(l)) for l in layer]

    hM_hit1Eta_hit2Eta_proto = [GetHistogram(histfname, 'hM_hit1Eta_hit2Eta_proto_layer{}'.format(l)) for l in layer]
    hM_hit1Phi_hit2Phi_proto = [GetHistogram(histfname, 'hM_hit1Phi_hit2Phi_proto_layer{}'.format(l)) for l in layer]
    hM_dPhi_dEta_proto = [GetHistogram(histfname, 'hM_dPhi_dEta_proto_layer{}'.format(l)) for l in layer]
    hM_dPhi_dEta_proto_altrange = [GetHistogram(histfname, 'hM_dPhi_dEta_proto_altrange_layer{}'.format(l)) for l in layer]
    hM_NTracklet_NClusLayer1_proto = [GetHistogram(histfname, 'hM_NTracklet_NClusLayer1_proto_layer{}'.format(l)) for l in layer]

    hM_hit1Eta_hit2Eta_reco = [GetHistogram(histfname, 'hM_hit1Eta_hit2Eta_reco_layer{}'.format(l)) for l in layer]
    hM_hit1Phi_hit2Phi_reco = [GetHistogram(histfname, 'hM_hit1Phi_hit2Phi_reco_layer{}'.format(l)) for l in layer]
    hM_dPhi_dEta_reco = [GetHistogram(histfname, 'hM_dPhi_dEta_reco_layer{}'.format(l)) for l in layer]
    hM_dPhi_dEta_reco_altrange = [GetHistogram(histfname, 'hM_dPhi_dEta_reco_altrange_layer{}'.format(l)) for l in layer]
    hM_NTracklet_NClusLayer1_reco = [GetHistogram(histfname, 'hM_NTracklet_NClusLayer1_reco_layer{}'.format(l)) for l in layer]

    hM_NTracklet_NClusLayer1_recogmatch = [GetHistogram(histfname, 'hM_NTracklet_NClusLayer1_recogmatch_layer{}'.format(l)) for l in layer]

    hM_Eta_vtxZ_proto_incl = [GetHistogram(histfname, 'hM_Eta_vtxZ_proto_incl_layer{}'.format(l)) for l in layer]
    hM_Eta_vtxZ_reco_incl = [GetHistogram(histfname, 'hM_Eta_vtxZ_reco_incl_layer{}'.format(l)) for l in layer]
    hM_Eta_vtxZ_recogmatch_incl = [GetHistogram(histfname, 'hM_Eta_vtxZ_recogmatch_incl_layer{}'.format(l)) for l in layer]

    hM_GenmatchedRecotkldR_GenhadPt = [GetHistogram(histfname, 'hM_GenmatchedRecotkldR_GenhadPt_layer{}'.format(l)) for l in layer]
    hM_GenmatchedRecotkldR_GenhadPt_altrange2 = [GetHistogram(histfname, 'hM_GenmatchedRecotkldR_GenhadPt_altrange2_layer{}'.format(l)) for l in layer]

    # Do the plots
    Draw_1DhistsComp(hM_dEta_proto, False, False, False, 1.3, 'Proto-tracklet #Delta#eta', '', plotpath+'ProtoTracklet_dEta_comp')
    hM_dEta_proto_altrange[0].GetXaxis().SetMaxDigits(2)
    Draw_1DhistsComp(hM_dEta_proto_altrange, False, False, False, 1.3, 'Proto-tracklet #Delta#eta', '', plotpath+'ProtoTracklet_dEta_altrange_comp')
    hM_dEta_proto_altrange2[0].GetXaxis().SetMaxDigits(2)
    Draw_1DhistsComp(hM_dEta_proto_altrange2, False, False, False, 1.8, 'Proto-tracklet #Delta#eta', '', plotpath+'ProtoTracklet_dEta_altrange2_comp')
    Draw_1DhistsComp(hM_dPhi_proto, False, False, False, 1.3, 'Proto-tracklet #Delta#phi', '', plotpath+'ProtoTracklet_dPhi_comp')
    Draw_1DhistsComp(hM_dPhi_proto_altrange, False, False, False, 1.3, 'Proto-tracklet #Delta#phi', '', plotpath+'ProtoTracklet_dPhi_altrange_comp')
    Draw_1DhistsComp(hM_Eta_proto, False, False, False, 1.5, 'Proto-tracklet #eta', '', plotpath+'ProtoTracklet_Eta_comp')
    Draw_1DhistsComp(hM_Phi_proto, False, False, False, 1.5, 'Proto-tracklet #phi', '', plotpath+'ProtoTracklet_Phi_comp')
    Draw_1DhistsComp(hM_NTracklet_proto, False, False, True, 150, 'Number of proto-tracklets', '', plotpath+'ProtoTracklet_NTracklet_comp')
    Draw_1DhistsComp(hM_dR_proto, False, False, True, 1.3, 'Proto-tracklet #DeltaR', '', plotpath+'ProtoTracklet_dR_comp')
    Draw_1DhistsComp(hM_dR_proto_altrange, False, False, True, 1.3, 'Proto-tracklet #DeltaR', '', plotpath+'ProtoTracklet_dR_altrange_comp')
    Draw_1DhistsComp(hM_dR_proto_LogX, False, True, True, 1.3, 'Proto-tracklet #DeltaR', '', plotpath+'ProtoTracklet_dR_logX_comp')

    plot_Stack(hM_Eta_proto[0], [hM_Eta_proto_PVzRange1[0], hM_Eta_proto_PVzRange2[0], hM_Eta_proto_PVzRange3[0], hM_Eta_proto_PVzRange4[0], hM_Eta_proto_PVzRange5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'],
               ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm'], False, 1, 'Proto-tracklet #eta', '', plotpath+'ProtoTracklet_Eta_layer12_Stack')
    plot_Stack(hM_Eta_proto[1], [hM_Eta_proto_PVzRange1[1], hM_Eta_proto_PVzRange2[1], hM_Eta_proto_PVzRange3[1], hM_Eta_proto_PVzRange4[1], hM_Eta_proto_PVzRange5[1]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'],
               ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm'], False, 1.5, 'Proto-tracklet #eta', '', plotpath+'ProtoTracklet_Eta_layer23_Stack')
    plot_Stack(hM_Eta_proto[2], [hM_Eta_proto_PVzRange1[2], hM_Eta_proto_PVzRange2[2], hM_Eta_proto_PVzRange3[2], hM_Eta_proto_PVzRange4[2], hM_Eta_proto_PVzRange5[2]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'],
               ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm'], False, 1.5, 'Proto-tracklet #eta', '', plotpath+'ProtoTracklet_Eta_layer13_Stack')

    Draw_1DhistsComp(hM_dEta_reco, False, False, True, 1.3, 'Reco-tracklet #Delta#eta', '', plotpath+'RecoTracklet_dEta_comp')
    hM_dEta_reco_altrange[0].GetXaxis().SetMaxDigits(2)
    Draw_1DhistsComp(hM_dEta_reco_altrange, False, False, False, 1.3, 'Reco-tracklet #Delta#eta', '', plotpath+'RecoTracklet_dEta_altrange_comp')
    hM_dEta_reco_altrange2[0].GetXaxis().SetMaxDigits(2)
    Draw_1DhistsComp(hM_dEta_reco_altrange2, False, False, False, 1.3, 'Reco-tracklet #Delta#eta', '', plotpath+'RecoTracklet_dEta_altrange2_comp')
    Draw_1DhistsComp(hM_dPhi_reco, False, False, False, 1.3, 'Reco-tracklet #Delta#phi', '', plotpath+'RecoTracklet_dPhi_comp')
    Draw_1DhistsComp(hM_dPhi_reco_altrange, False, False, False, 1.3, 'Reco-tracklet #Delta#phi', '', plotpath+'RecoTracklet_dPhi_altrange_comp')
    Draw_1DhistsComp(hM_Eta_reco, False, False, False, 1.5, 'Reco-tracklet #eta', '', plotpath+'RecoTracklet_Eta_comp')
    Draw_1DhistsComp(hM_Phi_reco, False, False, False, 1.5, 'Reco-tracklet #phi', '', plotpath+'RecoTracklet_Phi_comp')
    Draw_1DhistsComp(hM_NTracklet_reco, False, False, True, 150, 'Number of reco-tracklets', '', plotpath+'RecoTracklet_NTracklet_comp')
    Draw_1DhistsComp(hM_dR_reco, False, False, True, 1.3, 'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_dR_comp')
    Draw_1DhistsComp(hM_dR_reco_altrange, False, False, True, 1.3, 'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_dR_altrange_comp')
    Draw_1DhistsComp(hM_dR_reco_LogX, False, True, True, 1.3, 'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_dR_logX_comp')

    plot_Stack(hM_Eta_reco[0], [hM_Eta_reco_PVzRange1[0], hM_Eta_reco_PVzRange2[0], hM_Eta_reco_PVzRange3[0], hM_Eta_reco_PVzRange4[0], hM_Eta_reco_PVzRange5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'],
               ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm'], False, 1, 'Reco-tracklet #eta', '', plotpath+'RecoTracklet_Eta_layer12_Stack')
    plot_Stack(hM_Eta_reco[1], [hM_Eta_reco_PVzRange1[1], hM_Eta_reco_PVzRange2[1], hM_Eta_reco_PVzRange3[1], hM_Eta_reco_PVzRange4[1], hM_Eta_reco_PVzRange5[1]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'],
               ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm'], False, 1.5, 'Reco-tracklet #eta', '', plotpath+'RecoTracklet_Eta_layer23_Stack')
    plot_Stack(hM_Eta_reco[2], [hM_Eta_reco_PVzRange1[2], hM_Eta_reco_PVzRange2[2], hM_Eta_reco_PVzRange3[2], hM_Eta_reco_PVzRange4[2], hM_Eta_reco_PVzRange5[2]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'],
               ['PV_{z} < -7 cm', '-7 < PV_{z} < -2.5 cm', '-2.5 < PV_{z} < 2.5 cm', '2.5 < PV_{z} < 7 cm', 'PV_{z} > 7 cm'], False, 1.5, 'Reco-tracklet #eta', '', plotpath+'RecoTracklet_Eta_layer13_Stack')

    Draw_1DhistsComp(hM_NTracklet_recogmatch, False, False, True, 150, 'Number of gen-matched reco-tracklets', '', plotpath+'RecoTracklet_genmatched_NTracklet_comp')

    plot_Stack(hM_dR_recogmatch_ghadPtIncl[0], [hM_dR_recogmatch_ghadPtg3[0], hM_dR_recogmatch_ghadPt2to3[0], hM_dR_recogmatch_ghadPt1to2[0], hM_dR_recogmatch_ghadPt0p5to1[0], hM_dR_recogmatch_ghadPt0to0p5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'], [
               'p_{T}^{Gen hadron} > 3 GeV', '2 < p_{T}^{Gen hadron} < 3 GeV', '1 < p_{T}^{Gen hadron} < 2 GeV', '0.5 < p_{T}^{Gen hadron} < 1 GeV', '0 < p_{T}^{Gen hadron} < 0.5 GeV'], True, 2000, 'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_genmatched_dR_layer12_Stack_logy')
    plot_Stack(hM_dR_recogmatch_ghadPtIncl[0], [hM_dR_recogmatch_ghadPtg3[0], hM_dR_recogmatch_ghadPt2to3[0], hM_dR_recogmatch_ghadPt1to2[0], hM_dR_recogmatch_ghadPt0p5to1[0], hM_dR_recogmatch_ghadPt0to0p5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'], [
               'p_{T}^{Gen hadron} > 3 GeV', '2 < p_{T}^{Gen hadron} < 3 GeV', '1 < p_{T}^{Gen hadron} < 2 GeV', '0.5 < p_{T}^{Gen hadron} < 1 GeV', '0 < p_{T}^{Gen hadron} < 0.5 GeV'], False, 1.5, 'Reco-tracklet #DeltaR', '', plotpath+'RecoTracklet_genmatched_dR_layer12_Stack')
    plot_Stack(hM_dEta_recogmatch_ghadPtIncl[0], [hM_dEta_recogmatch_ghadPtg3[0], hM_dEta_recogmatch_ghadPt2to3[0], hM_dEta_recogmatch_ghadPt1to2[0], hM_dEta_recogmatch_ghadPt0p5to1[0], hM_dEta_recogmatch_ghadPt0to0p5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'], [
               'p_{T}^{Gen hadron} > 3 GeV', '2 < p_{T}^{Gen hadron} < 3 GeV', '1 < p_{T}^{Gen hadron} < 2 GeV', '0.5 < p_{T}^{Gen hadron} < 1 GeV', '0 < p_{T}^{Gen hadron} < 0.5 GeV'], True, 2000, 'Reco-tracklet #Delta#eta', '', plotpath+'RecoTracklet_genmatched_deta_layer12_Stack_logy')
    plot_Stack(hM_dEta_recogmatch_ghadPtIncl[0], [hM_dEta_recogmatch_ghadPtg3[0], hM_dEta_recogmatch_ghadPt2to3[0], hM_dEta_recogmatch_ghadPt1to2[0], hM_dEta_recogmatch_ghadPt0p5to1[0], hM_dEta_recogmatch_ghadPt0to0p5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'], [
               'p_{T}^{Gen hadron} > 3 GeV', '2 < p_{T}^{Gen hadron} < 3 GeV', '1 < p_{T}^{Gen hadron} < 2 GeV', '0.5 < p_{T}^{Gen hadron} < 1 GeV', '0 < p_{T}^{Gen hadron} < 0.5 GeV'], False, 1.5, 'Reco-tracklet #Delta#eta', '', plotpath+'RecoTracklet_genmatched_deta_layer12_Stack')
    plot_Stack(hM_dPhi_recogmatch_ghadPtIncl[0], [hM_dPhi_recogmatch_ghadPtg3[0], hM_dPhi_recogmatch_ghadPt2to3[0], hM_dPhi_recogmatch_ghadPt1to2[0], hM_dPhi_recogmatch_ghadPt0p5to1[0], hM_dPhi_recogmatch_ghadPt0to0p5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'], [
               'p_{T}^{Gen hadron} > 3 GeV', '2 < p_{T}^{Gen hadron} < 3 GeV', '1 < p_{T}^{Gen hadron} < 2 GeV', '0.5 < p_{T}^{Gen hadron} < 1 GeV', '0 < p_{T}^{Gen hadron} < 0.5 GeV'], True, 2000, 'Reco-tracklet #Delta#phi', '', plotpath+'RecoTracklet_genmatched_dphi_layer12_Stack_logy')
    plot_Stack(hM_dPhi_recogmatch_ghadPtIncl[0], [hM_dPhi_recogmatch_ghadPtg3[0], hM_dPhi_recogmatch_ghadPt2to3[0], hM_dPhi_recogmatch_ghadPt1to2[0], hM_dPhi_recogmatch_ghadPt0p5to1[0], hM_dPhi_recogmatch_ghadPt0to0p5[0]], ['#073b4c', '#118ab2', '#06d6a0', '#ffd166', '#ef476f'], [
               'p_{T}^{Gen hadron} > 3 GeV', '2 < p_{T}^{Gen hadron} < 3 GeV', '1 < p_{T}^{Gen hadron} < 2 GeV', '0.5 < p_{T}^{Gen hadron} < 1 GeV', '0 < p_{T}^{Gen hadron} < 0.5 GeV'], False, 1.5, 'Reco-tracklet #Delta#phi', '', plotpath+'RecoTracklet_genmatched_dphi_layer12_Stack')

    Draw_2Dhist(hM_hit1Eta_hit2Eta_proto[0], False, False, 0.16, '#eta_{cluster 1}', '#eta_{cluster 2}', 'colz', plotpath+'ProtoTracklet_Clus1EtaClus2Eta_layer12')
    Draw_2Dhist(hM_hit1Phi_hit2Phi_proto[0], False, False, 0.16, '#phi_{cluster 1}', '#phi_{cluster 2}', 'colz', plotpath+'ProtoTracklet_Clus1PhiClus2Phi_layer12')
    Draw_2Dhist(hM_hit1Eta_hit2Eta_reco[0], False, False, 0.16, '#eta_{cluster 1}', '#eta_{cluster 2}', 'colz', plotpath+'RecoTracklet_Clus1EtaClus2Eta_layer12')
    Draw_2Dhist(hM_hit1Phi_hit2Phi_reco[0], False, False, 0.16, '#phi_{cluster 1}', '#phi_{cluster 2}', 'colz', plotpath+'RecoTracklet_Clus1PhiClus2Phi_layer12')
    Draw_2Dhist(hM_dPhi_dEta_proto[0], False, False, 0.15, 'Proto-tracklet #Delta#phi', 'Proto-tracklet #Delta#eta', 'colz', plotpath+'ProtoTracklet_dPhi_dEta_layer12')
    Draw_2Dhist(hM_dPhi_dEta_proto_altrange[0], True, False, 0.15, 'Proto-tracklet #Delta#phi', 'Proto-tracklet #Delta#eta', 'colz', plotpath+'ProtoTracklet_dPhi_dEta_altrange_LogZ_layer12')
    Draw_2Dhist(hM_dPhi_dEta_proto_altrange[0], False, False, 0.15, 'Proto-tracklet #Delta#phi', 'Proto-tracklet #Delta#eta', 'colz', plotpath+'ProtoTracklet_dPhi_dEta_altrange_LinZ_layer12')
    Draw_2Dhist(hM_dPhi_dEta_reco[0], False, False, 0.13, 'Reco-tracklet #Delta#phi', 'Reco-tracklet #Delta#eta', 'colz', plotpath+'RecoTracklet_dPhi_dEta_layer12')
    Draw_2Dhist(hM_dPhi_dEta_reco_altrange[0], True, False, 0.17, 'Reco-tracklet #Delta#phi', 'Reco-tracklet #Delta#eta', 'colz', plotpath+'RecoTracklet_dPhi_dEta_altrange_LogZ_layer12')
    Draw_2Dhist(hM_dPhi_dEta_reco_altrange[0], False, False, 0.17, 'Reco-tracklet #Delta#phi', 'Reco-tracklet #Delta#eta', 'colz', plotpath+'RecoTracklet_dPhi_dEta_altrange_LinZ_layer12')
    Draw_2Dhist(hM_NTracklet_NClusLayer1_proto[0], False, False, 0.17, 'Number of proto-tracklets', 'Number of clusters (MVTX 1st layer)', 'colz', plotpath+'ProtoTracklet_NTracklet_NClusLayer1')
    Draw_2Dhist(hM_NTracklet_NClusLayer1_reco[0], False, False, 0.17, 'Number of reco-tracklets', 'Number of clusters (MVTX 1st layer)', 'colz', plotpath+'RecoTracklet_NTracklet_NClusLayer1')

    # Draw_2Dhist(hist, logz, norm1, rmargin, XaxisName, YaxisName, drawopt, outname):
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt[0], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_layer12')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt[1], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_layer23')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt[2], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_layer13')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt_altrange2[0], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_altrange2_layer12')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt_altrange2[1], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_altrange2_layer23')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt_altrange2[2], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_altrange2_layer13')
    Draw_2Dhist(hM_NTracklet_NClusLayer1_recogmatch[0], True, False, 0.17, 'Number of gen-matched reco-tracklets', 'Number of clusters (1st layer)', 'colz', plotpath+'RecoTracklet_genmatched_NTracklet_NClusLayer1_layer12')
    # Draw_2Dhist(hM_NTracklet_NClusLayer1_recogmatch[1], True, False, 0.17, 'Number of gen-matched reco-tracklets', 'Number of clusters (1st layer)', 'colz', plotpath+'RecoTracklet_genmatched_NTracklet_NClusLayer1_layer23')
    # Draw_2Dhist(hM_NTracklet_NClusLayer1_recogmatch[2], True, False, 0.17, 'Number of gen-matched reco-tracklets', 'Number of clusters (1st layer)', 'colz', plotpath+'RecoTracklet_genmatched_NTracklet_NClusLayer1_layer13')

    Draw_2Dhist(hM_Eta_vtxZ_proto_incl[0], False, False, 0.16, 'Proto-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'ProtoTracklet_Eta_vtxZ_layer12')
    Draw_2Dhist(hM_Eta_vtxZ_proto_incl[1], False, False, 0.16, 'Proto-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'ProtoTracklet_Eta_vtxZ_layer23')
    Draw_2Dhist(hM_Eta_vtxZ_proto_incl[2], False, False, 0.16, 'Proto-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'ProtoTracklet_Eta_vtxZ_layer13')
    Draw_2Dhist(hM_Eta_vtxZ_reco_incl[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_layer12')
    Draw_2Dhist(hM_Eta_vtxZ_reco_incl[1], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_layer23')
    Draw_2Dhist(hM_Eta_vtxZ_reco_incl[2], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_layer13')
    Draw_2Dhist(hM_Eta_vtxZ_recogmatch_incl[0], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_genmatched_layer12')
    Draw_2Dhist(hM_Eta_vtxZ_recogmatch_incl[1], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_genmatched_layer23')
    Draw_2Dhist(hM_Eta_vtxZ_recogmatch_incl[2], False, False, 0.14, 'Reco-tracklet #eta', 'Primary vertex V_{z} (cm)', 'colz', plotpath+'RecoTracklet_Eta_vtxZ_genmatched_layer13')

    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt[0], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_layer12')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt[1], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_layer23')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt[2], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_layer13')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt_altrange2[0], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_altrange2_layer12')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt_altrange2[1], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_altrange2_layer23')
    Draw_2Dhist(hM_GenmatchedRecotkldR_GenhadPt_altrange2[2], True, False, 0.17, 'Gen-matched reco-tracklet #Delta R', 'Gen-hadron p_{T}', 'colz', plotpath+'RecoTracklet_genmatched_dR_GenhadronPt_altrange2_layer13')
