#ifndef VTXDATA_H
#define VTXDATA_H

#include <vector>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <stdlib.h>
#include <map>

#include <TFile.h>
#include <TTree.h>

using namespace std;

struct VtxData
{
    int event, NhitsLayer1, NTruthVtx;
    float PV_z, TruthPV_trig_z, TruthPV_mostNpart_z;
};

void SetVtxMinitree(TTree *outTree, VtxData &vtxdata)
{
    outTree->Branch("event", &vtxdata.event);
    outTree->Branch("NhitsLayer1", &vtxdata.NhitsLayer1);
    outTree->Branch("NTruthVtx", &vtxdata.NTruthVtx);
    outTree->Branch("PV_z", &vtxdata.PV_z);
    outTree->Branch("TruthPV_trig_z", &vtxdata.TruthPV_trig_z);
    outTree->Branch("TruthPV_mostNpart_z", &vtxdata.TruthPV_mostNpart_z);
}

std::map<int, float> EvtVtx_map_tklcluster(const char *vtxfname)
{
    std::map<int, float> EvtVtx_map;

    TFile *f = new TFile(vtxfname, "READ");
    TTree *t = (TTree *)f->Get("minitree");
    int event;
    float PV_z;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("PV_z", &PV_z);
    for (int ev = 0; ev < t->GetEntriesFast(); ev++)
    {
        t->GetEntry(ev);
        EvtVtx_map[ev] = PV_z;
    }

    return EvtVtx_map;
}

#endif