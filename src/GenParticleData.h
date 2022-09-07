#ifndef GENPARTICLEDATA_H
#define GENPARTICLEDATA_H

#include <vector>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <stdlib.h>
#include <map>

#include <TFile.h>
#include <TTree.h>

#include "Utilities.h"

using namespace std;

struct GenParticleData
{
    int event, NGenHadron;
    vector<int> GenHadron_PID, GenHadron_PrimaryID, GenHadron_ParentID;
    vector<float> GenHadron_pT, GenHadron_Eta, GenHadron_Phi, GenHadron_E;
};

void SetGenParticleMinitree(TTree *outTree, GenParticleData &genparticledata)
{
    outTree->Branch("event", &genparticledata.event);
    outTree->Branch("NGenHadron", &genparticledata.NGenHadron);
    outTree->Branch("GenHadron_PID", &genparticledata.GenHadron_PID);
    outTree->Branch("GenHadron_PrimaryID", &genparticledata.GenHadron_PrimaryID);
    outTree->Branch("GenHadron_ParentID", &genparticledata.GenHadron_ParentID);
    outTree->Branch("GenHadron_pT", &genparticledata.GenHadron_pT);
    outTree->Branch("GenHadron_Eta", &genparticledata.GenHadron_Eta);
    outTree->Branch("GenHadron_Phi", &genparticledata.GenHadron_Phi);
    outTree->Branch("GenHadron_E", &genparticledata.GenHadron_E);
}

std::map<int, vector<float>> Map_EvtGenHadron(const char *fname)
{
    std::map<int, vector<float>> EvtGenHadron_map;

    vector<float> vEtaPhi;

    TFile *f = new TFile(fname, "READ");
    TTree *t = (TTree *)f->Get("minitree");
    int event;
    vector<float> *GenHadron_Eta = 0, *GenHadron_Phi = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("GenHadron_Eta", &GenHadron_Eta);
    t->SetBranchAddress("GenHadron_Phi", &GenHadron_Phi);
    for (int ev = 0; ev < t->GetEntriesFast(); ev++)
    {
        t->GetEntry(ev);

        for (size_t ipart = 0; ipart < GenHadron_Eta->size(); ipart++)
        {
            vEtaPhi.push_back(GenHadron_Eta->at(ipart));
            vEtaPhi.push_back(GenHadron_Phi->at(ipart));
        }

        EvtGenHadron_map[ev] = vEtaPhi;
        CleanVec(vEtaPhi);
    }

    return EvtGenHadron_map;
}

#endif