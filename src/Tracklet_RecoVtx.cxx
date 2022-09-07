#include "Tracklet.h"

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

#include "VtxData.h"

int main(int argc, char *argv[])
{
    int dPhiCutbin = TString(argv[1]).Atoi();
    int dZCutbin = TString(argv[2]).Atoi();
    // TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/MVTXRecoClusters_Ntuple.root";
    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/MVTXRecoClusters_NoPileup_Nevt2000.root";
    TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin%d_dZCutbin%d_NoPileup_Nevt2000.root", dPhiCutbin, dZCutbin);
    // int vtxzmean = -3;
    // int vtxzwidth = 8;
    // TString infilename = Form("/sphenix/user/hjheng/TrackletAna/data/myMVTXHits_SIMPLE1000Pions_VtxZmean%dcm_VtxZwidth%dcm.root", vtxzmean, vtxzwidth);
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_Simple1000Pions_VtxZmean%dcm_VtxZwidth%dcm_RecoVtx_TklCluster.root", vtxzmean, vtxzwidth);

    vector<Hit *> MVTXlayer1;
    vector<Hit *> MVTXlayer2;

    VtxData vtxdata = {};

    TH1F *hM_NTruthVtx = new TH1F("hM_NTruthVtx", "hM_NTruthVtx", 11, 0, 11);
    TH1F *hM_DiffVtxZ_trig = new TH1F("hM_DiffVtxZ_trig", "hM_DiffVtxZ_trig", 50, -5, 5);
    TH1F *hM_DiffVtxZ_mostNpart = new TH1F("hM_DiffVtxZ_mostNpart", "hM_DiffVtxZ_mostNpart", 50, -5, 5);
    TH1F *hM_DiffVtxZ_Nvtx1 = new TH1F("hM_DiffVtxZ_Nvtx1", "hM_DiffVtxZ_Nvtx1", 50, -5, 5);
    TH2F *hM_TruthPVzTrig_RecoPVz = new TH2F("hM_TruthPVzTrig_RecoPVz", "hM_TruthPVzTrig_RecoPVz", 100, -25, 25, 100, -25, 25);
    TH2F *hM_TruthPVzMostNpart_RecoPVz = new TH2F("hM_TruthPVzMostNpart_RecoPVz", "hM_TruthPVzMostNpart_RecoPVz", 100, -25, 25, 100, -25, 25);
    TH2F *hM_TruthPVz_RecoPVz_Nvtx1 = new TH2F("hM_TruthPVz_RecoPVz_Nvtx1", "hM_TruthPVz_RecoPVz_Nvtx1", 100, -25, 25, 100, -25, 25);
    TH2F *hM_DiffVtxZTrig_NhitsLayer1 = new TH2F("hM_DiffVtxZTrig_NhitsLayer1", "hM_DiffVtxZTrig_NhitsLayer1", 100, -5, 5, 200, 0, 20000);
    TH2F *hM_DiffVtxZMostNpart_NhitsLayer1 = new TH2F("hM_DiffVtxZMostNpart_NhitsLayer1", "hM_DiffVtxZMostNpart_NhitsLayer1", 100, -5, 5, 200, 0, 20000);
    TH2F *hM_DiffVtxZ_NhitsLayer1_Nvtx1 = new TH2F("hM_DiffVtxZ_NhitsLayer1_Nvtx1", "hM_DiffVtxZ_NhitsLayer1_Nvtx1", 100, -5, 5, 200, 0, 20000);

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("EventTree");
    t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    int event, NTruthVtx;
    float TruthPV_trig_z, TruthPV_mostNpart_z;
    vector<int> *ClusLayer = 0;
    vector<float> *ClusX = 0, *ClusY = 0, *ClusZ = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("NTruthVtx", &NTruthVtx);
    t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
    t->SetBranchAddress("TruthPV_mostNpart_z", &TruthPV_mostNpart_z);
    t->SetBranchAddress("ClusLayer", &ClusLayer);
    t->SetBranchAddress("ClusX", &ClusX);
    t->SetBranchAddress("ClusY", &ClusY);
    t->SetBranchAddress("ClusZ", &ClusZ);

    TFile *outfile = new TFile(outfilename, "RECREATE");
    TTree *minitree = new TTree("minitree", "Minitree of reconstructed vertices");
    SetVtxMinitree(minitree, vtxdata);

    for (int i = 0; i < index->GetN(); i++)
    // for (Long64_t ev = 0; ev < t->GetEntriesFast(); ev++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);
        // cout << "event = " << event << " local = " << local << endl;
        cout << "event=" << event << " has a total of " << ClusLayer->size() << " clusters" << endl;

        CleanVec(MVTXlayer1);
        CleanVec(MVTXlayer2);

        for (size_t ihit = 0; ihit < ClusLayer->size(); ihit++)
        {
            if (ClusLayer->at(ihit) == 0)
            {
                Hit *hit = new Hit(ClusX->at(ihit), ClusY->at(ihit), ClusZ->at(ihit), 0., 0., 0., 0);
                MVTXlayer1.push_back(hit);
            }
            else if (ClusLayer->at(ihit) == 1)
            {
                Hit *hit = new Hit(ClusX->at(ihit), ClusY->at(ihit), ClusZ->at(ihit), 0., 0., 0., 1);
                MVTXlayer2.push_back(hit);
            }
            else
                continue;
        }

        cout << "# of clusters in Layer 1 = " << MVTXlayer1.size() << ", Layer 2 = " << MVTXlayer2.size() << endl;

        // Values for CMS Pb-Pb analysis: float dPhi_cut = 0.08; float dZ_cut = 0.14;
        // float dPhi_cut = dPhiCutbin * 0.01;
        // float dZ_cut = dZCutbin * 0.01;
        float PV_z = TrackletPV_cluster(event, MVTXlayer1, MVTXlayer2, dPhiCutbin, dZCutbin, false);

        cout << "NTruthVtx = " << NTruthVtx << endl
             << "Truth primary (most N particles) Z position = " << TruthPV_mostNpart_z << endl
             << "Recontructed primary vertex Z position = " << PV_z << endl;

        hM_NTruthVtx->Fill(NTruthVtx);
        hM_DiffVtxZ_trig->Fill((PV_z - TruthPV_trig_z));
        hM_DiffVtxZ_mostNpart->Fill((PV_z - TruthPV_mostNpart_z));
        hM_TruthPVzTrig_RecoPVz->Fill(TruthPV_trig_z, PV_z);
        hM_TruthPVzMostNpart_RecoPVz->Fill(TruthPV_mostNpart_z, PV_z);
        hM_DiffVtxZTrig_NhitsLayer1->Fill((PV_z - TruthPV_trig_z), (int)MVTXlayer1.size());
        hM_DiffVtxZMostNpart_NhitsLayer1->Fill((PV_z - TruthPV_mostNpart_z), (int)MVTXlayer1.size());
        if (NTruthVtx == 1)
        {
            if (TruthPV_trig_z != TruthPV_mostNpart_z)
                cout << "For events with only 1 truth vertex, (TruthPV_trig_z,TruthPV_mostNpart_z) = (" << TruthPV_trig_z << "," << TruthPV_mostNpart_z << "). Check!" << endl;

            hM_DiffVtxZ_Nvtx1->Fill((PV_z - TruthPV_trig_z));
            hM_TruthPVz_RecoPVz_Nvtx1->Fill(TruthPV_trig_z, PV_z);
            hM_DiffVtxZ_NhitsLayer1_Nvtx1->Fill((PV_z - TruthPV_mostNpart_z), (int)MVTXlayer1.size());
        }

        vtxdata.event = event;
        vtxdata.NTruthVtx = NTruthVtx;
        vtxdata.NhitsLayer1 = (int)MVTXlayer1.size();
        vtxdata.PV_z = PV_z;
        vtxdata.TruthPV_trig_z = TruthPV_trig_z;
        vtxdata.TruthPV_mostNpart_z = TruthPV_mostNpart_z;
        minitree->Fill();
    }

    outfile->cd();
    minitree->Write();
    hM_NTruthVtx->Write();
    hM_DiffVtxZ_trig->Write();
    hM_DiffVtxZ_mostNpart->Write();
    hM_DiffVtxZ_Nvtx1->Write();
    hM_TruthPVzTrig_RecoPVz->Write();
    hM_TruthPVzMostNpart_RecoPVz->Write();
    hM_TruthPVz_RecoPVz_Nvtx1->Write();
    hM_DiffVtxZTrig_NhitsLayer1->Write();
    hM_DiffVtxZMostNpart_NhitsLayer1->Write();
    hM_DiffVtxZ_NhitsLayer1_Nvtx1->Write();
    outfile->Close();

    f->Close();

    return 0;
}
