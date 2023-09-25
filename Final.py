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
from array import array
import json
from plot.plotUtil import GetHistogram
from measurements.measurements import *

gROOT.LoadMacro('./plot/sPHENIXStyle/sPhenixStyle.C')
gROOT.ProcessLine('SetsPhenixStyle()')
gROOT.SetBatch(True)

TickSize = 0.03
AxisTitleSize = 0.05
AxisLabelSize = 0.04
LeftMargin = 0.15
RightMargin = 0.08
TopMargin = 0.08
BottomMargin = 0.13
FillAlpha = 0.5

def msrmnt_sphenix_auau_200gev():
    dNdEta_eta0, Centrality, dNdEta_eta0_err, Centrality_err = [], [], [], []

    # ToDo: modify when the systematic uncertainties are available
    for cb in range(20):
        hM_dNdEtafinal_layer12 = GetHistogram('./plot/corrections/layer12_CentBin{}to{}/correction_hists.root'.format(cb,cb), 'h1WEfinal')
        hM_dNdEtafinal_layer23 = GetHistogram('./plot/corrections/layer23_CentBin{}to{}/correction_hists.root'.format(cb,cb), 'h1WEfinal')
        hM_dNdEtafinal_layer13 = GetHistogram('./plot/corrections/layer13_CentBin{}to{}/correction_hists.root'.format(cb,cb), 'h1WEfinal')

        nbins = hM_dNdEtafinal_layer12.GetNbinsX()
        dNdEta_eta0_l12 = hM_dNdEtafinal_layer12.GetBinContent(int((nbins + 1) / 2))
        dNdEta_eta0_l23 = hM_dNdEtafinal_layer23.GetBinContent(int((nbins + 1) / 2))
        dNdEta_eta0_l13 = hM_dNdEtafinal_layer13.GetBinContent(int((nbins + 1) / 2))
        
        dNdEta_eta0.append(dNdEta_eta0_l12)
        Centrality.append((2 * cb + 1) / 2 * 5)
        dNdEta_eta0_err.append(0)
        Centrality_err.append(0)

    dNdEta_eta0 = array('f', dNdEta_eta0)
    Centrality = array('f', Centrality)
    dNdEta_eta0_err = array('f', dNdEta_eta0_err)
    Centrality_err = array('f', Centrality_err)
    gsphenix_auau_200gev = TGraphErrors(len(Centrality), Centrality, dNdEta_eta0, Centrality_err, dNdEta_eta0_err)

    return gsphenix_auau_200gev


def msrmnt_dict():
    msrmntdict = {}

    msrmntdict['cms_pbpb_2p76'] = cms_pbpb_2p76()
    msrmntdict['cms_xexe_5p44'] = cms_xexe_5p44()
    msrmntdict['atlas_pbpb_2p76'] = atlas_pbpb_2p76()
    msrmntdict['alice_pbpb_2p76'] = alice_pbpb_2p76()
    msrmntdict['alice_pbpb_5p02'] = alice_pbpb_5p02()
    msrmntdict['alice_xexe_5p44'] = alice_xexe_5p44()
    msrmntdict['phobos_auau_0p2'] = phobos_auau_0p2()
    msrmntdict['phobos_auau_0p13'] = phobos_auau_0p13()
    msrmntdict['phobos_auau_0p0624'] = phobos_auau_0p0624()
    msrmntdict['phobos_auau_0p0196'] = phobos_auau_0p0196()
    msrmntdict['phenix_auau_0p2'] = phenix_auau_0p2()
    msrmntdict['phenix_auau_0p13'] = phenix_auau_0p13()
    msrmntdict['phenix_auau_0p0624'] = phenix_auau_0p0624()
    msrmntdict['phenix_auau_0p039'] = phenix_auau_0p039()
    msrmntdict['phenix_auau_0p027'] = phenix_auau_0p027()
    msrmntdict['phenix_auau_0p0196'] = phenix_auau_0p0196()
    msrmntdict['phenix_auau_0p0145'] = phenix_auau_0p0145()
    msrmntdict['phenix_auau_0p0077'] = phenix_auau_0p0077()
    msrmntdict['brahms_auau_0p2'] = brahms_auau_0p2()
    msrmntdict['sphenix_auau_200gev'] = msrmnt_sphenix_auau_200gev()

    return msrmntdict


def gstyle(g, mstyle, msize, color, lwidth, alpha):
    g.SetMarkerStyle(mstyle)
    g.SetMarkerSize(msize)
    g.SetMarkerColor(color)
    g.SetLineColor(color)
    g.SetLineWidth(lwidth)
    g.SetFillColorAlpha(color, alpha)


def dNdeta_eta0_Summary(jsonpath, canvname, axisrange = [-1, 101, 1, 3000], canvsize = [1000, 500], grlegpos = [0.68, 0.18, 0.92, 0.45], prelimtext='Work-in-progress'):
    xmin = axisrange[0]
    xmax = axisrange[1]
    ymin = axisrange[2]
    ymax = axisrange[3]
    # Create a canvas
    c = TCanvas('c_'+canvname, 'c_'+canvname, canvsize[0], canvsize[1])
    c.cd()
    # Set left margin and logarithmic y-axis
    gPad.SetLeftMargin(LeftMargin)
    gPad.SetLogy()
    # Create a frame for the graph
    gframe = TH1F('gframe_' + canvname, 'gframe_' + canvname, 1, xmin, xmax)
    gframe.SetLabelOffset(999, 'X')
    gframe.SetTickLength(0, 'X')
    gframe.GetXaxis().SetTitle('Centrality (%)')
    gframe.GetYaxis().SetTitle('#frac{dN_{ch}}{d#eta}#lbar_{#eta=0}')
    gframe.GetYaxis().SetMoreLogLabels()
    gframe.GetYaxis().SetTitleOffset(1.35)
    gframe.SetAxisRange(ymin, ymax, 'Y')
    gframe.Draw()
    # Create and draw the x-axis on top
    axis = TGaxis(xmax, ymin, xmin, ymin, xmin, xmax, 511, '-')
    axis.SetLabelOffset(-0.032)
    axis.SetLabelFont(43)
    axis.SetLabelSize(22)
    axis.Draw()
    # Create and draw the y-axis on top
    axisup = TGaxis(xmax, ymax, xmin, ymax, xmin, xmax, 511, '+')
    axisup.SetLabelOffset(1)
    axisup.Draw()
    # Create the legend for the TGraphs
    gr_leg = TLegend(grlegpos[0], grlegpos[1], grlegpos[2], grlegpos[3])
    gr_leg.SetTextSize(0.035)
    gr_leg.SetFillStyle(0)
    # Load measurement dictionaries
    dict_msrmnt= msrmnt_dict()
    # Load the json file
    with open(jsonpath, 'r') as fp:
        grtags = json.load(fp)
        for tagName, par in grtags.items():
            gr = dict_msrmnt[tagName]
            gstyle(gr, par['markerstyle'], par['markersize'], TColor.GetColor(par['markercolor']), par['linewidth'], par['alpha'])
            gr.Draw(par['DrawOption'][0])
            if par['DrawOption'][1] != '':
                gr.Draw(par['DrawOption'][1])
            gr_leg.AddEntry(gr, par['Legendtext'], par['Legendstyle'])

    gr_leg.Draw()
    sphnxleg = TLegend((1 - RightMargin) - 0.83, (1 - TopMargin) - 0.08, (1 - RightMargin) - 0.4, (1 - TopMargin) - 0.03)
    sphnxleg.SetTextSize(0.045)
    sphnxleg.SetFillStyle(0)
    sphnxleg.AddEntry('', '#it{#bf{sPHENIX}} ' + prelimtext, '')
    # sphnxleg/AddEntry('', 'Au+Au #sqrt{s_{NN}}=200 GeV', '')
    sphnxleg.Draw()
    c.SaveAs(plotpath + 'dNdEta_eta0_{}.png'.format(canvname));
    c.SaveAs(plotpath + 'dNdEta_eta0_{}.pdf'.format(canvname));


if __name__ == '__main__':
    plotpath = './plot/dNdEtaFinal/'
    os.makedirs(plotpath, exist_ok=True)

    dNdeta_eta0_Summary('./measurements/RHIC_plotparam.json', 'RHIC')
    dNdeta_eta0_Summary('./measurements/LHC_plotparam.json', 'LHC')
    