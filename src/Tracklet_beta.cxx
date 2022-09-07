#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include <TObjString.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeIndex.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include "GenParticleData.h"
#include "Utilities.h"

int main(int argc, char *argv[])
{
    int NevtToRun_ = TString(argv[1]).Atoi();
    int skip = TString(argv[2]).Atoi();
    int layer_ = TString(argv[3]).Atoi();
    int l1 = ((int)layer_ / 10) - 1;
    int l2 = (layer_ % 10) - 1;
    int iniEvt = skip;



    TString EvtGenHadron_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_RecoClusters_GenHadron_Nevt1000.root";
    TString infilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal/TrackletAna_minitree_layer%d_Evt%dto%d_RandhitCase0_MisAlignNum0.root", layer_, iniEvt, iniEvt + NevtToRun_);
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal/TrackletAna_minitree_layer%d_Evt%dto%d_RandhitCase%d_MisAlignNum%d.root", layer_, iniEvt, iniEvt + NevtToRun_, randhit_case, misalignment_num);
}