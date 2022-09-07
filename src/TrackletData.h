#ifndef TRACKLETDATA_H
#define TRACKLETDATA_H

#include <TTree.h>

#include "Utilities.h"

class TrackletData
{
public:
    int event, NhitsLayer1, NRecotkl_Total, NRecotkl_GenMatched;
    // int update_l1, update_l2;
    float PV_z, TruthPV_trig_z, TruthPV_mostNpart_z;
    // vector<vector<Hit *>> MVTXlayers;
    // vector<Tracklet *> ProtoTkls, RecoTkls;
    // vector<GenHadron *> GenHadrons;
    vector<float> recotkl_eta1, recotkl_phi1, recotkl_eta2, recotkl_phi2, recotkl_deta, recotkl_dphi, recotkl_dR;
    vector<float> recotkl_genmatched_eta1, recotkl_genmatched_phi1, recotkl_genmatched_eta2, recotkl_genmatched_phi2, recotkl_genmatched_deta, recotkl_genmatched_dphi, recotkl_genmatched_dR;
};

void SetMinitree(TTree *outTree, TrackletData &tkldata)
{
    outTree->Branch("event", &tkldata.event);
    outTree->Branch("NhitsLayer1", &tkldata.NhitsLayer1);
    outTree->Branch("PV_z", &tkldata.PV_z);
    outTree->Branch("NRecotkl_Total", &tkldata.NRecotkl_Total);
    outTree->Branch("recotkl_eta1", &tkldata.recotkl_eta1);
    outTree->Branch("recotkl_phi1", &tkldata.recotkl_phi1);
    outTree->Branch("recotkl_eta2", &tkldata.recotkl_eta2);
    outTree->Branch("recotkl_phi2", &tkldata.recotkl_phi2);
    outTree->Branch("recotkl_deta", &tkldata.recotkl_deta);
    outTree->Branch("recotkl_dphi", &tkldata.recotkl_dphi);
    outTree->Branch("recotkl_dR", &tkldata.recotkl_dR);
    outTree->Branch("NRecotkl_GenMatched", &tkldata.NRecotkl_GenMatched);
    outTree->Branch("recotkl_genmatched_eta1", &tkldata.recotkl_genmatched_eta1);
    outTree->Branch("recotkl_genmatched_phi1", &tkldata.recotkl_genmatched_phi1);
    outTree->Branch("recotkl_genmatched_eta2", &tkldata.recotkl_genmatched_eta2);
    outTree->Branch("recotkl_genmatched_phi2", &tkldata.recotkl_genmatched_phi2);
    outTree->Branch("recotkl_genmatched_deta", &tkldata.recotkl_genmatched_deta);
    outTree->Branch("recotkl_genmatched_dphi", &tkldata.recotkl_genmatched_dphi);
    outTree->Branch("recotkl_genmatched_dR", &tkldata.recotkl_genmatched_dR);
}

#endif