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
gStyle.SetPalette(kLightTemperature)

if __name__ == '__main__':
    plotpath = './MVTXTrkrHits/MVTXHits_ana325private_singleEvtDst_Nevt2000_sepCondor_20221021/'
    os.makedirs(plotpath, exist_ok=True)

    df = ROOT.RDataFrame('EventTree', '/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters_20221021_Nevt2000.root')
    hM_G4HitT = df.Define("G4HitT_new", "G4HitT[G4Hit_IsSameClus==0]").Histo1D(("hM_G4HitT", " ", 100, 0, 100), "G4HitT_new")

    hM_G4HitXY = df.Define("G4HitX_new", "G4HitX[G4Hit_IsSameClus==0]").Define("G4HitY_new", "G4HitY[G4Hit_IsSameClus==0]").Histo2D(("hM_G4HitXY", " ", 10000, -100, 100, 10000, -100, 100), "G4HitX_new", "G4HitY_new")
    hM_G4HitXY_layerNot99 = df.Define("G4HitX_new2", "G4HitX[G4Hit_IsSameClus==0&&G4HitLayer!=99]").Define("G4HitY_new2", "G4HitY[G4Hit_IsSameClus==0&&G4HitLayer!=99]").Histo2D(("hM_G4HitXY_layerNot99", " ", 10000, -100, 100, 10000, -100, 100), "G4HitX_new2", "G4HitY_new2")
    hM_G4HitXY_zoomin = df.Define("G4HitX_new", "G4HitX[G4Hit_IsSameClus==0]").Define("G4HitY_new", "G4HitY[G4Hit_IsSameClus==0]").Histo2D(("hM_G4HitXY_zoomin", " ", 200, -5, 5, 200, -5, 5), "G4HitX_new", "G4HitY_new")
    hM_G4HitR_G4HitT = df.Define("G4HitX_new", "G4HitX[G4Hit_IsSameClus==0]").Define("G4HitY_new", "G4HitY[G4Hit_IsSameClus==0]").Define("G4HitR", "sqrt(G4HitX_new*G4HitX_new+G4HitY_new*G4HitY_new)").Define("G4HitT_new", "G4HitT[G4Hit_IsSameClus==0]").Histo2D(("hM_G4HitR_G4HitT", " ", 100, 0, 100, 100, 0, 20000), "G4HitR", "G4HitT_new")
    hM_G4HitR_G4HitT_range2 = df.Define("G4HitX_new", "G4HitX[G4Hit_IsSameClus==0]").Define("G4HitY_new", "G4HitY[G4Hit_IsSameClus==0]").Define("G4HitR", "sqrt(G4HitX_new*G4HitX_new+G4HitY_new*G4HitY_new)").Define("G4HitT_new", "G4HitT[G4Hit_IsSameClus==0]").Histo2D(("hM_G4HitR_G4HitT_range2", " ", 100, 2, 5, 100, 0, 50), "G4HitR", "G4HitT_new")

    hM_TrkrHitRow_TrkrHitColumn = df.Histo2D(("hM_TrkrHitRow_TrkrHitColumn", " ", 520, 0, 520, 1050, 0, 1050), "TrkrHitRow", "TrkrHitColumn")
    hM_TrkrHitStaveID_TrkrHitLayerID = df.Histo2D(("hM_TrkrHitStaveID_TrkrHitLayerID", " ", 20, 0, 20, 3, 0, 3), "TrkrHitStaveID", "TrkrHitLayerID")
    hM_TrkrHitChipID_TrkrHitLayerID = df.Histo2D(("hM_TrkrHitChipID_TrkrHitLayerID", " ", 9, 0, 9, 3, 0, 3), "TrkrHitChipID", "TrkrHitLayerID")
    hM_TrkrHitChipID_TrkrHitStaveID = df.Histo2D(("hM_TrkrHitChipID_TrkrHitStaveID", " ", 9, 0, 9, 20, 0, 20), "TrkrHitChipID", "TrkrHitStaveID")
    hM_TrkrHitChipID_G4HitZ = df.Histo2D(("hM_TrkrHitChipID_G4HitZ", " ", 9, 0, 9, 60, -15, 15), "TrkrHitChipID", "G4HitZ")


    # Draw_1Dhist(hist, norm1, logy, ymaxscale, XaxisName, Ytitle_unit, outname)
    Draw_1Dhist(hM_G4HitT, False, True, 50, 'G4Hit time (ns)', 'ns', plotpath+"G4HitT")
    # Draw_2Dhist(hist, logz, norm1, rmargin, XaxisName, YaxisName, outname)
    Draw_2Dhist(hM_G4HitXY, False, False, 0.15, 'G4Hit x pos (cm)', 'G4Hit y pos (cm)', plotpath+"G4HitXY")
    Draw_2Dhist(hM_G4HitXY_layerNot99, False, False, 0.15, 'G4Hit x pos (cm)', 'G4Hit y pos (cm)', plotpath+"G4HitXY_layerNot99")
    Draw_2Dhist(hM_G4HitXY_zoomin, False, False, 0.15, 'G4Hit x pos (cm)', 'G4Hit y pos (cm)', plotpath+"G4HitXY_zoomin")
    Draw_2Dhist(hM_G4HitR_G4HitT, True, False, 0.15, 'G4Hit radius (cm)', 'G4Hit time (ns)', plotpath+"G4HitR_G4HitT")
    Draw_2Dhist(hM_G4HitR_G4HitT_range2, True, False, 0.15, 'G4Hit radius (cm)', 'G4Hit time (ns)', plotpath+"G4HitR_G4HitT_range2")
    Draw_2Dhist(hM_TrkrHitRow_TrkrHitColumn, False, False, 0.13, 'TrkrHit row', 'TrkrHit column', plotpath+"TrkrHitRow_TrkrHitColumn")
    Draw_2Dhist(hM_TrkrHitStaveID_TrkrHitLayerID, False, False, 0.15, 'TrkrHit Stave ID', 'TrkrHit Layer ID', plotpath+"TrkrHitStaveID_TrkrHitLayerID")
    Draw_2Dhist(hM_TrkrHitChipID_TrkrHitLayerID, False, False, 0.15, 'TrkrHit Chip ID', 'TrkrHit Layer ID', plotpath+"TrkrHitChipID_TrkrHitLayerID")
    Draw_2Dhist(hM_TrkrHitChipID_G4HitZ, False, False, 0.15, 'TrkrHit Chip ID', 'G4Hit z pos (cm)', plotpath+"TrkrHitChipID_G4HitZ")