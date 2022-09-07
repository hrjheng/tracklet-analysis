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
    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/MVTXRecoClusters_Nevt1000.root";
    TString outfilename = "/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_RecoClusters_GenHadron_Nevt1000.root";

    GenParticleData genparticledata = {};

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("EventTree");
    t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    int event;
    vector<float> *G4PartPt = 0, *G4PartEta = 0, *G4PartPhi = 0, *G4PartE = 0;
    vector<int> *G4PartPID = 0, *G4PartPrimaryID = 0, *G4PartParentID = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("G4PartPID", &G4PartPID);
    t->SetBranchAddress("G4PartPrimaryID", &G4PartPrimaryID);
    t->SetBranchAddress("G4PartParentID", &G4PartParentID);
    t->SetBranchAddress("G4PartPt", &G4PartPt);
    t->SetBranchAddress("G4PartEta", &G4PartEta);
    t->SetBranchAddress("G4PartPhi", &G4PartPhi);
    t->SetBranchAddress("G4PartE", &G4PartE);

    TFile *outfile = new TFile(outfilename, "RECREATE");
    TTree *minitree = new TTree("minitree", "Minitree of generated charged hadrons");
    SetGenParticleMinitree(minitree, genparticledata);
    for (int i = 0; i < index->GetN(); i++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);

        printf("event=%d has a total of %d generated particles\n", event, G4PartPID->size());
        for (size_t ipart = 0; ipart < G4PartPID->size(); ipart++)
        {
            if (abs(G4PartPID->at(ipart)) < 9 ||
                (abs(G4PartPID->at(ipart)) > 10 && abs(G4PartPID->at(ipart)) < 19) ||
                (abs(G4PartPID->at(ipart)) > 20 && abs(G4PartPID->at(ipart)) < 26) ||
                (abs(G4PartPID->at(ipart)) > 31 && abs(G4PartPID->at(ipart)) < 38))
                continue;

            if (G4PartPt->at(ipart) < 0.03) continue;
            if (fabs(G4PartEta->at(ipart)) > 2.5) continue;
                
            genparticledata.GenHadron_PID.push_back(G4PartPID->at(ipart));
            genparticledata.GenHadron_PrimaryID.push_back(G4PartPrimaryID->at(ipart));
            genparticledata.GenHadron_ParentID.push_back(G4PartParentID->at(ipart));
            genparticledata.GenHadron_pT.push_back(G4PartPt->at(ipart));
            genparticledata.GenHadron_Eta.push_back(G4PartEta->at(ipart));
            genparticledata.GenHadron_Phi.push_back(G4PartPhi->at(ipart));
            genparticledata.GenHadron_E.push_back(G4PartE->at(ipart));

            // printf("ipart=%d, PID=%d, Primary ID=%d, Parent ID=%d, (pT,Eta,Phi,E)=(%f,%f,%f,%f)\n", ipart, G4PartPID->at(ipart), G4PartPrimaryID->at(ipart), G4PartParentID->at(ipart), G4PartPt->at(ipart), G4PartEta->at(ipart), G4PartPhi->at(ipart), G4PartE->at(ipart));
        }

        genparticledata.NGenHadron = genparticledata.GenHadron_PID.size();
        minitree->Fill();

        CleanVec(genparticledata.GenHadron_PID);
        CleanVec(genparticledata.GenHadron_PrimaryID);
        CleanVec(genparticledata.GenHadron_ParentID);
        CleanVec(genparticledata.GenHadron_pT);
        CleanVec(genparticledata.GenHadron_Eta);
        CleanVec(genparticledata.GenHadron_Phi);
        CleanVec(genparticledata.GenHadron_E);
    }

    return 0;
}