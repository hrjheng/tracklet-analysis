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

if __name__ == '__main__':
    plotpath = './ClusRecoIneff/'
    os.makedirs(plotpath, exist_ok=True)

    histfilepath = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/'
    filelist = [histfilepath + 'Hists_ClusterRecoIneff_SimplePions.root']
    # histfilepath = '/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/'
    # filelist = [histfilepath + 'Hists_ClusterRecoIneff_HijingAuAuNoPileup_ana325private.root']
    filetype = ['sim']

    for i, fname in enumerate(filelist):
        hM_clusterl2_eta = GetHistogram(fname, 'hM_recotkll2_eta_{}'.format(filetype[i]))
        hM_clusterl2_eta_nonmatched = GetHistogram(fname, 'hM_recotkll2_eta_nonmatched_{}'.format(filetype[i]))

        err_ClusRecoIneff = TGraphAsymmErrors()
        err_ClusRecoIneff.BayesDivide(hM_clusterl2_eta_nonmatched, hM_clusterl2_eta)

        Draw_simpleEff(err_ClusRecoIneff, [-3,3.], [0,1.05], '#eta', 'Inefficiency', plotpath + 'ClusRecoIneff_{}'.format(filetype[i])+'_eta')