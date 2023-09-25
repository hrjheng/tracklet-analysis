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
from plot.plotUtil import *

gROOT.LoadMacro('./plot/sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)

if __name__ == '__main__':
    plotpath = './plot/systematics/'
    os.makedirs(plotpath, exist_ok=True)

    for cb in range(20):
        hM_dNdEtafinal_layer12 = GetHistogram('./plot/corrections/layer12_CentBin{}to{}/correction_hists.root'.format(cb,cb), 'h1WEfinal')
        hM_dNdEtafinal_layer23 = GetHistogram('./plot/corrections/layer23_CentBin{}to{}/correction_hists.root'.format(cb,cb), 'h1WEfinal')
        hM_dNdEtafinal_layer13 = GetHistogram('./plot/corrections/layer13_CentBin{}to{}/correction_hists.root'.format(cb,cb), 'h1WEfinal')

        nbins = hM_dNdEtafinal_layer12.GetNbinsX()
        dNdEta_eta0_l12 = hM_dNdEtafinal_layer12.GetBinContent(int((nbins + 1) / 2))
        dNdEta_eta0_l23 = hM_dNdEtafinal_layer23.GetBinContent(int((nbins + 1) / 2))
        dNdEta_eta0_l13 = hM_dNdEtafinal_layer13.GetBinContent(int((nbins + 1) / 2))
        print(dNdEta_eta0_l12, dNdEta_eta0_l23, dNdEta_eta0_l13)
        # lst_hM_dNdEtafinal_layer = [hM_dNdEtafinal_layer12, hM_dNdEtafinal_layer23, hM_dNdEtafinal_layer13]
        Draw_1DhistsComp([hM_dNdEtafinal_layer12, hM_dNdEtafinal_layer23, hM_dNdEtafinal_layer13], False, False, False, 1.3, '#eta', '', plotpath+'dNdEta_layervariation_CentBin{}to{}'.format(cb,cb))
