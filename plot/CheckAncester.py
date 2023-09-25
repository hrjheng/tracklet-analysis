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
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)

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

if __name__ == '__main__':
    plotpath = './GenMatching_checkAncestor/'
    os.makedirs(plotpath, exist_ok=True)
    fname = '/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters.root'

    G4PAncPID_all = np.array([])

    hM_G4PAnc_Pt = TH1F('hM_G4PAnc_Pt', 'hM_G4PAnc_Pt', 100, 0, 10)
    hM_G4PAnc_Eta = TH1F('hM_G4PAnc_Eta', 'hM_G4PAnc_Eta', 80, -8, 8)
    hM_G4PAnc_Phi = TH1F('hM_G4PAnc_Phi', 'hM_G4PAnc_Phi', 35, -3.5, 3.5)
    hM_G4PAnc_E = TH1F('hM_G4PAnc_E', 'hM_G4PAnc_E', 100, 0, 50)


    f = TFile(fname, 'r')
    tree = f.Get('EventTree')
    for idx in range(tree.GetEntries()):
        tree.GetEntry(idx)
        
        G4PAncPID_all = np.concatenate((G4PAncPID_all,np.array(tree.G4PfromClus_AncPID)))

        for i in range(len(tree.UniqueAncG4P_Eta)):
            hM_G4PAnc_Pt.Fill(tree.UniqueAncG4P_Pt[i])
            hM_G4PAnc_Eta.Fill(tree.UniqueAncG4P_Eta[i])
            hM_G4PAnc_Phi.Fill(tree.UniqueAncG4P_Phi[i])
            hM_G4PAnc_E.Fill(tree.UniqueAncG4P_E[i])

    unique, counts = np.unique(G4PAncPID_all, return_counts=True)
    print(np.asarray((unique, counts)).T)

    Draw_1Dhist(hM_G4PAnc_Pt, False, True, 150, 'G4Particle p_{T} (GeV)', 'GeV', plotpath+'G4PAnc_Pt')
    Draw_1Dhist(hM_G4PAnc_Eta, False, True, 150, 'G4Particle #eta', '', plotpath+'G4PAnc_Eta')
    Draw_1Dhist(hM_G4PAnc_Phi, False, True, 150, 'G4Particle #phi', '', plotpath+'G4PAnc_Phi')
    Draw_1Dhist(hM_G4PAnc_E, False, True, 150, 'G4Particle E (GeV)', 'GeV', plotpath+'G4PAnc_E')

    fig = plt.figure(figsize = (10, 5))
    plt.bar(unique, counts, color='black', width = 1.2)
    plt.xlabel("Unique PID")
    plt.ylabel("Counts")
    plt.yscale("log")
    plt.savefig(plotpath+"G4PAncPID.pdf")