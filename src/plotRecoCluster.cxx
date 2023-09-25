#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TTreeIndex.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "GenHadron.h"
#include "Tracklet.h"
#include "Vertex.h"
#include "misalignment.h"

void Fillhist_CluspairMinDr(vector<Hit *> hits, TH1F *hM_cluspairDr)
{
    vector<bool> paired;
    paired.resize(hits.size());
    std::fill(paired.begin(), paired.end(), false);

    for (int i = 0; i < hits.size(); i++)
    {
        if (paired[i])
            continue;
        float minDr = 999;
        int minDr_j = -1;
        for (int j = i + 1; j < hits.size(); j++)
        {
            if (paired[j])
                continue;
            float dr = deltaR(hits[i]->Eta(), hits[i]->Phi(), hits[j]->Eta(), hits[j]->Phi());
            if (dr < minDr)
            {
                minDr = dr;
                minDr_j = j;
            }
        }
        if (minDr_j == -1)
            continue;
        else
        {
            if (minDr < 0.0001)
                cout << "i = " << i << " (eta,phi)=(" << hits[i]->Eta() << "," << hits[i]->Phi() << "), minDr_j = " << minDr_j << " (eta,phi)=(" << hits[minDr_j]->Eta() << "," << hits[minDr_j]->Phi()
                     << "), minDr = " << minDr << endl;
            paired[i] = true;
            paired[minDr_j] = true;
            hM_cluspairDr->Fill(minDr);
        }
    }
}

void Fillhist_CluspairDr(vector<Hit *> hits, TH1F *hM_cluspairDr)
{
    for (int i = 0; i < hits.size(); i++)
    {
        for (int j = i + 1; j < hits.size(); j++)
        {
            float dr = deltaR(hits[i]->Eta(), hits[i]->Phi(), hits[j]->Eta(), hits[j]->Phi());
            hM_cluspairDr->Fill(dr);
        }
    }
}

int main(int argc, char *argv[])
{
    float gap_north = TString(argv[1]).Atof(); // Nominal: 2.9mm
    float gap_upper = TString(argv[2]).Atof(); // Nominal: (gap_north / 2.) mm
    float cent_shift = TString(argv[3]).Atof(); // Nominal: 0.0mm
    TString gap_north_str = TString(argv[1]).Replace(TString(argv[1]).Index("."), 1, "p");
    TString gap_upper_str = TString(argv[2]).Replace(TString(argv[2]).Index("."), 1, "p");
    TString cent_shift_str = TString(argv[3]).Replace(TString(argv[3]).Index("."), 1, "p");

    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters.root";
    // TString outfilename = "/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_RecoClusters_twohalfangle_extreme.root";
    TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_RecoClusters_twohalfangle_GapNorth%s_GapUpper%s_CentShift%s_updated.root", gap_north_str.Data(), gap_upper_str.Data(), cent_shift_str.Data());
    // TString EvtVtx_map_filename =
    // "/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_twohalfangle_extreme.root";
    TString EvtVtx_map_filename = Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth%s_GapUpper%s_CentShift%s_3DVertex_twohalves.root", gap_north_str.Data(), gap_upper_str.Data(), cent_shift_str.Data());

    // TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters.root";
    // TString outfilename = "/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_RecoClusters.root";
    // TString EvtVtx_map_filename =
    // "/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root";

    // TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/Pion/Ntuple_Pions_Npart500_VtxZmean0cm_VtxZwidth0cm_pT0p001to5GeV.root";
    // TString outfilename = "/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/Hists_Pions_Npart500_VtxZmean0cm_VtxZwidth0cm_pT0p001to5GeV.root";
    // TString EvtVtx_map_filename =
    // "/sphenix/user/hjheng/TrackletAna/minitree/Pion_RecoVtx_Optimization/SimplePion_Npart500_VtxZmean0cm_VtxZwidth0cm_pT0p001to5GeV/TrackletAna_RecoClusters_SimplePion_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root";
    std::map<int, float> EvtVtx_map = EvtVtx_map_tklcluster(EvtVtx_map_filename.Data());

    TH1F *hM_ClusX_all = new TH1F("hM_ClusX_all", "hM_ClusX_all", 120, -6, 6);
    TH1F *hM_ClusX_layer1 = new TH1F("hM_ClusX_layer1", "hM_ClusX_layer1", 120, -6, 6);
    TH1F *hM_ClusX_layer2 = new TH1F("hM_ClusX_layer2", "hM_ClusX_layer2", 120, -6, 6);
    TH1F *hM_ClusX_layer3 = new TH1F("hM_ClusX_layer3", "hM_ClusX_layer3", 120, -6, 6);
    TH1F *hM_ClusY_all = new TH1F("hM_ClusY_all", "hM_ClusY_all", 120, -6, 6);
    TH1F *hM_ClusY_layer1 = new TH1F("hM_ClusY_layer1", "hM_ClusY_layer1", 120, -6, 6);
    TH1F *hM_ClusY_layer2 = new TH1F("hM_ClusY_layer2", "hM_ClusY_layer2", 120, -6, 6);
    TH1F *hM_ClusY_layer3 = new TH1F("hM_ClusY_layer3", "hM_ClusY_layer3", 120, -6, 6);
    TH1F *hM_ClusZ_all = new TH1F("hM_ClusZ_all", "hM_ClusZ_all", 150, -15, 15);
    TH1F *hM_ClusZ_layer1 = new TH1F("hM_ClusZ_layer1", "hM_ClusZ_layer1", 150, -15, 15);
    TH1F *hM_ClusZ_layer2 = new TH1F("hM_ClusZ_layer2", "hM_ClusZ_layer2", 150, -15, 15);
    TH1F *hM_ClusZ_layer3 = new TH1F("hM_ClusZ_layer3", "hM_ClusZ_layer3", 150, -15, 15);
    TH1F *hM_ClusR_all = new TH1F("hM_ClusR_all", "hM_ClusR_all", 120, 2, 5);
    TH1F *hM_ClusR_layer1 = new TH1F("hM_ClusR_layer1", "hM_ClusR_layer1", 120, 2, 5);
    TH1F *hM_ClusR_layer2 = new TH1F("hM_ClusR_layer2", "hM_ClusR_layer2", 120, 2, 5);
    TH1F *hM_ClusR_layer3 = new TH1F("hM_ClusR_layer3", "hM_ClusR_layer3", 120, 2, 5);
    TH1F *hM_ClusEtaPV_all = new TH1F("hM_ClusEtaPV_all", "hM_ClusEtaPV_all", 160, -4, 4);
    TH1F *hM_ClusEtaPV_layer1 = new TH1F("hM_ClusEtaPV_layer1", "hM_ClusEtaPV_layer1", 160, -4, 4);
    TH1F *hM_ClusEtaPV_layer2 = new TH1F("hM_ClusEtaPV_layer2", "hM_ClusEtaPV_layer2", 160, -4, 4);
    TH1F *hM_ClusEtaPV_layer3 = new TH1F("hM_ClusEtaPV_layer3", "hM_ClusEtaPV_layer3", 160, -4, 4);
    TH1F *hM_ClusEtaOri_all = new TH1F("hM_ClusEtaOri_all", "hM_ClusEtaOri_all", 120, -3, 3);
    TH1F *hM_ClusEtaOri_layer1 = new TH1F("hM_ClusEtaOri_layer1", "hM_ClusEtaOri_layer1", 120, -3, 3);
    TH1F *hM_ClusEtaOri_layer2 = new TH1F("hM_ClusEtaOri_layer2", "hM_ClusEtaOri_layer2", 120, -3, 3);
    TH1F *hM_ClusEtaOri_layer3 = new TH1F("hM_ClusEtaOri_layer3", "hM_ClusEtaOri_layer3", 120, -3, 3);
    TH1F *hM_ClusPhi_all = new TH1F("hM_ClusPhi_all", "hM_ClusPhi_all", 140, -3.5, 3.5);
    TH1F *hM_ClusPhi_layer1 = new TH1F("hM_ClusPhi_layer1", "hM_ClusPhi_layer1", 140, -3.5, 3.5);
    TH1F *hM_ClusPhi_layer2 = new TH1F("hM_ClusPhi_layer2", "hM_ClusPhi_layer2", 140, -3.5, 3.5);
    TH1F *hM_ClusPhi_layer3 = new TH1F("hM_ClusPhi_layer3", "hM_ClusPhi_layer3", 140, -3.5, 3.5);
    TH1F *hM_ClusADC_all = new TH1F("hM_ClusADC_all", "hM_ClusADC_all", 50, 0, 50);
    TH1F *hM_ClusADC_layer1 = new TH1F("hM_ClusADC_layer1", "hM_ClusADC_layer1", 50, 0, 50);
    TH1F *hM_ClusADC_layer2 = new TH1F("hM_ClusADC_layer2", "hM_ClusADC_layer2", 50, 0, 50);
    TH1F *hM_ClusADC_layer3 = new TH1F("hM_ClusADC_layer3", "hM_ClusADC_layer3", 50, 0, 50);
    TH1F *hM_ClusZSize_all = new TH1F("hM_ClusZSize_all", "hM_ClusZSize_all", 20, 0, 20);
    TH1F *hM_ClusZSize_layer1 = new TH1F("hM_ClusZSize_layer1", "hM_ClusZSize_layer1", 20, 0, 20);
    TH1F *hM_ClusZSize_layer2 = new TH1F("hM_ClusZSize_layer2", "hM_ClusZSize_layer2", 20, 0, 20);
    TH1F *hM_ClusZSize_layer3 = new TH1F("hM_ClusZSize_layer3", "hM_ClusZSize_layer3", 20, 0, 20);
    TH1F *hM_ClusPhiSize_all = new TH1F("hM_ClusPhiSize_all", "hM_ClusPhiSize_all", 20, 0, 20);
    TH1F *hM_ClusPhiSize_layer1 = new TH1F("hM_ClusPhiSize_layer1", "hM_ClusPhiSize_layer1", 20, 0, 20);
    TH1F *hM_ClusPhiSize_layer2 = new TH1F("hM_ClusPhiSize_layer2", "hM_ClusPhiSize_layer2", 20, 0, 20);
    TH1F *hM_ClusPhiSize_layer3 = new TH1F("hM_ClusPhiSize_layer3", "hM_ClusPhiSize_layer3", 20, 0, 20);
    TH1F *hM_ClusPairDr_layer1 = new TH1F("hM_ClusPairDr_layer1", "hM_ClusPairDr_layer1", 50, 0, 0.1);
    TH1F *hM_ClusPairDr_layer2 = new TH1F("hM_ClusPairDr_layer2", "hM_ClusPairDr_layer2", 50, 0, 0.1);
    TH1F *hM_ClusPairDr_layer3 = new TH1F("hM_ClusPairDr_layer3", "hM_ClusPairDr_layer3", 50, 0, 0.1);
    TH2F *hM_ClusX_ClusY_all = new TH2F("hM_ClusX_ClusY_all", "hM_ClusX_ClusY_all", 240, -6, 6, 240, -6, 6);
    TH2F *hM_ClusZ_ClusPhi_all = new TH2F("hM_ClusZ_ClusPhi_all", "hM_ClusZ_ClusPhi_all", 900, -15, 15, 350, -3.5, 3.5);
    TH2F *hM_ClusZ_ClusPhi_layer1 = new TH2F("hM_ClusZ_ClusPhi_layer1", "hM_ClusZ_ClusPhi_layer1", 900, -15, 15, 350, -3.5, 3.5);
    TH2F *hM_ClusZ_ClusPhi_layer2 = new TH2F("hM_ClusZ_ClusPhi_layer2", "hM_ClusZ_ClusPhi_layer2", 900, -15, 15, 350, -3.5, 3.5);
    TH2F *hM_ClusZ_ClusPhi_layer3 = new TH2F("hM_ClusZ_ClusPhi_layer3", "hM_ClusZ_ClusPhi_layer3", 900, -15, 15, 350, -3.5, 3.5);
    TH2F *hM_ClusEta_ClusZSize_all = new TH2F("hM_ClusEta_ClusZSize_all", "hM_ClusEta_ClusZSize_all", 160, -4, 4, 20, 0, 20);
    TH2F *hM_ClusEta_ClusZSize_layer1 = new TH2F("hM_ClusEta_ClusZSize_layer1", "hM_ClusEta_ClusZSize_layer1", 160, -4, 4, 20, 0, 20);
    TH2F *hM_ClusEta_ClusZSize_layer2 = new TH2F("hM_ClusEta_ClusZSize_layer2", "hM_ClusEta_ClusZSize_layer2", 160, -4, 4, 20, 0, 20);
    TH2F *hM_ClusEta_ClusZSize_layer3 = new TH2F("hM_ClusEta_ClusZSize_layer3", "hM_ClusEta_ClusZSize_layer3", 160, -4, 4, 20, 0, 20);
    TH2F *hM_ClusPhi_ClusPhiSize_all = new TH2F("hM_ClusPhi_ClusPhiSize_all", "hM_ClusPhi_ClusPhiSize_all", 140, -3.5, 3.5, 20, 0, 20);
    TH2F *hM_ClusPhi_ClusPhiSize_layer1 = new TH2F("hM_ClusPhi_ClusPhiSize_layer1", "hM_ClusPhi_ClusPhiSize_layer1", 140, -3.5, 3.5, 20, 0, 20);
    TH2F *hM_ClusPhi_ClusPhiSize_layer2 = new TH2F("hM_ClusPhi_ClusPhiSize_layer2", "hM_ClusPhi_ClusPhiSize_layer2", 140, -3.5, 3.5, 20, 0, 20);
    TH2F *hM_ClusPhi_ClusPhiSize_layer3 = new TH2F("hM_ClusPhi_ClusPhiSize_layer3", "hM_ClusPhi_ClusPhiSize_layer3", 140, -3.5, 3.5, 20, 0, 20);
    TH2F *hM_ClusZSize_ClusPhiSize_all = new TH2F("hM_ClusZSize_ClusPhiSize_all", "hM_ClusZSize_ClusPhiSize_all", 20, 0, 20, 20, 0, 20);
    TH2F *hM_ClusZSize_ClusPhiSize_layer1 = new TH2F("hM_ClusZSize_ClusPhiSize_layer1", "hM_ClusZSize_ClusPhiSize_layer1", 20, 0, 20, 20, 0, 20);
    TH2F *hM_ClusZSize_ClusPhiSize_layer2 = new TH2F("hM_ClusZSize_ClusPhiSize_layer2", "hM_ClusZSize_ClusPhiSize_layer2", 20, 0, 20, 20, 0, 20);
    TH2F *hM_ClusZSize_ClusPhiSize_layer3 = new TH2F("hM_ClusZSize_ClusPhiSize_layer3", "hM_ClusZSize_ClusPhiSize_layer3", 20, 0, 20, 20, 0, 20);
    TH2F *hM_ClusEta_ClusADC_all = new TH2F("hM_ClusEta_ClusADC_all", "hM_ClusEta_ClusADC_all", 160, -4, 4, 50, 0, 50);
    TH2F *hM_ClusEta_ClusADC_layer1 = new TH2F("hM_ClusEta_ClusADC_layer1", "hM_ClusEta_ClusADC_layer1", 160, -4, 4, 50, 0, 50);
    TH2F *hM_ClusEta_ClusADC_layer2 = new TH2F("hM_ClusEta_ClusADC_layer2", "hM_ClusEta_ClusADC_layer2", 160, -4, 4, 50, 0, 50);
    TH2F *hM_ClusEta_ClusADC_layer3 = new TH2F("hM_ClusEta_ClusADC_layer3", "hM_ClusEta_ClusADC_layer3", 160, -4, 4, 50, 0, 50);
    TH2F *hM_ClusZSize_ClusADC_all = new TH2F("hM_ClusZSize_ClusADC_all", "hM_ClusZSize_ClusADC_all", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusZSize_ClusADC_layer1 = new TH2F("hM_ClusZSize_ClusADC_layer1", "hM_ClusZSize_ClusADC_layer1", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusZSize_ClusADC_layer2 = new TH2F("hM_ClusZSize_ClusADC_layer2", "hM_ClusZSize_ClusADC_layer2", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusZSize_ClusADC_layer3 = new TH2F("hM_ClusZSize_ClusADC_layer3", "hM_ClusZSize_ClusADC_layer3", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusPhiSize_ClusADC_all = new TH2F("hM_ClusPhiSize_ClusADC_all", "hM_ClusPhiSize_ClusADC_all", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusPhiSize_ClusADC_layer1 = new TH2F("hM_ClusPhiSize_ClusADC_layer1", "hM_ClusPhiSize_ClusADC_layer1", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusPhiSize_ClusADC_layer2 = new TH2F("hM_ClusPhiSize_ClusADC_layer2", "hM_ClusPhiSize_ClusADC_layer2", 20, 0, 20, 50, 0, 50);
    TH2F *hM_ClusPhiSize_ClusADC_layer3 = new TH2F("hM_ClusPhiSize_ClusADC_layer3", "hM_ClusPhiSize_ClusADC_layer3", 20, 0, 20, 50, 0, 50);

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("EventTree");
    t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    int event, NTruthVtx;
    float TruthPV_trig_z, TruthPV_mostNpart_z;
    vector<int> *ClusLayer = 0;
    vector<float> *ClusX = 0, *ClusY = 0, *ClusZ = 0, *ClusR = 0, *ClusPhi = 0, *ClusEta = 0, *ClusPhiSize = 0, *ClusZSize = 0;
    vector<unsigned int> *ClusAdc = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("NTruthVtx", &NTruthVtx);
    t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
    t->SetBranchAddress("TruthPV_mostNpart_z", &TruthPV_mostNpart_z);
    t->SetBranchAddress("ClusLayer", &ClusLayer);
    t->SetBranchAddress("ClusX", &ClusX);
    t->SetBranchAddress("ClusY", &ClusY);
    t->SetBranchAddress("ClusZ", &ClusZ);
    t->SetBranchAddress("ClusR", &ClusR);
    t->SetBranchAddress("ClusPhi", &ClusPhi);
    t->SetBranchAddress("ClusEta", &ClusEta);
    t->SetBranchAddress("ClusAdc", &ClusAdc);
    t->SetBranchAddress("ClusPhiSize", &ClusPhiSize);
    t->SetBranchAddress("ClusZSize", &ClusZSize);
    for (int i = 0; i < t->GetEntriesFast(); i++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);

        float PV_z = EvtVtx_map[event];

        vector<Hit *> l1hits, l2hits, l3hits;
        l1hits.clear();
        l2hits.clear();
        l3hits.clear();

        for (size_t i = 0; i < ClusLayer->size(); i++)
        {   
            float x0 = ClusX->at(i);
            float y0 = ClusY->at(i);
            float z0 = ClusZ->at(i);
            vector<float> cpos = {x0, y0, z0};
            UpdatePos_GapTwoHalves(cpos, gap_north, cent_shift, gap_upper);
            Hit *hit = new Hit(cpos[0], cpos[1], cpos[2], 0., 0., PV_z, ClusLayer->at(i));

            hM_ClusX_all->Fill(hit->posX());
            hM_ClusY_all->Fill(hit->posY());
            hM_ClusZ_all->Fill(hit->posZ());
            hM_ClusR_all->Fill(hit->rho());
            hM_ClusEtaPV_all->Fill(hit->Eta());
            hM_ClusEtaOri_all->Fill(ClusEta->at(i));
            hM_ClusPhi_all->Fill(hit->Phi());
            hM_ClusADC_all->Fill(ClusAdc->at(i));
            hM_ClusZSize_all->Fill(ClusZSize->at(i));
            hM_ClusPhiSize_all->Fill(ClusPhiSize->at(i));
            hM_ClusX_ClusY_all->Fill(hit->posX(), hit->posY());
            hM_ClusZ_ClusPhi_all->Fill(hit->posZ(), hit->Phi());
            hM_ClusEta_ClusZSize_all->Fill(hit->Eta(), ClusZSize->at(i));
            hM_ClusPhi_ClusPhiSize_all->Fill(hit->Phi(), ClusPhiSize->at(i));
            hM_ClusZSize_ClusPhiSize_all->Fill(ClusZSize->at(i), ClusPhiSize->at(i));
            hM_ClusEta_ClusADC_all->Fill(hit->Eta(), ClusAdc->at(i));
            hM_ClusZSize_ClusADC_all->Fill(ClusZSize->at(i), ClusAdc->at(i));
            hM_ClusPhiSize_ClusADC_all->Fill(ClusPhiSize->at(i), ClusAdc->at(i));
            if (ClusLayer->at(i) == 0)
            {
                l1hits.push_back(hit);
                hM_ClusX_layer1->Fill(hit->posX());
                hM_ClusY_layer1->Fill(hit->posY());
                hM_ClusZ_layer1->Fill(hit->posZ());
                hM_ClusR_layer1->Fill(hit->rho());
                hM_ClusEtaPV_layer1->Fill(hit->Eta());
                hM_ClusEtaOri_layer1->Fill(ClusEta->at(i));
                hM_ClusPhi_layer1->Fill(hit->Phi());
                hM_ClusADC_layer1->Fill(ClusAdc->at(i));
                hM_ClusZSize_layer1->Fill(ClusZSize->at(i));
                hM_ClusPhiSize_layer1->Fill(ClusPhiSize->at(i));
                hM_ClusZ_ClusPhi_layer1->Fill(hit->posZ(), hit->Phi());
                hM_ClusEta_ClusZSize_layer1->Fill(hit->Eta(), ClusZSize->at(i));
                hM_ClusPhi_ClusPhiSize_layer1->Fill(hit->Phi(), ClusPhiSize->at(i));
                hM_ClusZSize_ClusPhiSize_layer1->Fill(ClusZSize->at(i), ClusPhiSize->at(i));
                hM_ClusEta_ClusADC_layer1->Fill(hit->Eta(), ClusAdc->at(i));
                hM_ClusZSize_ClusADC_layer1->Fill(ClusZSize->at(i), ClusAdc->at(i));
                hM_ClusPhiSize_ClusADC_layer1->Fill(ClusPhiSize->at(i), ClusAdc->at(i));
            }
            else if (ClusLayer->at(i) == 1)
            {
                l2hits.push_back(hit);
                hM_ClusX_layer2->Fill(hit->posX());
                hM_ClusY_layer2->Fill(hit->posY());
                hM_ClusZ_layer2->Fill(hit->posZ());
                hM_ClusR_layer2->Fill(hit->rho());
                hM_ClusEtaPV_layer2->Fill(hit->Eta());
                hM_ClusEtaOri_layer2->Fill(ClusEta->at(i));
                hM_ClusPhi_layer2->Fill(hit->Phi());
                hM_ClusADC_layer2->Fill(ClusAdc->at(i));
                hM_ClusZSize_layer2->Fill(ClusZSize->at(i));
                hM_ClusPhiSize_layer2->Fill(ClusPhiSize->at(i));
                hM_ClusZ_ClusPhi_layer2->Fill(hit->posZ(), hit->Phi());
                hM_ClusEta_ClusZSize_layer2->Fill(hit->Eta(), ClusZSize->at(i));
                hM_ClusPhi_ClusPhiSize_layer2->Fill(hit->Phi(), ClusPhiSize->at(i));
                hM_ClusZSize_ClusPhiSize_layer2->Fill(ClusZSize->at(i), ClusPhiSize->at(i));
                hM_ClusEta_ClusADC_layer2->Fill(hit->Eta(), ClusAdc->at(i));
                hM_ClusZSize_ClusADC_layer2->Fill(ClusZSize->at(i), ClusAdc->at(i));
                hM_ClusPhiSize_ClusADC_layer2->Fill(ClusPhiSize->at(i), ClusAdc->at(i));
            }
            else if (ClusLayer->at(i) == 2)
            {
                l3hits.push_back(hit);
                hM_ClusX_layer3->Fill(hit->posX());
                hM_ClusY_layer3->Fill(hit->posY());
                hM_ClusZ_layer3->Fill(hit->posZ());
                hM_ClusR_layer3->Fill(hit->rho());
                hM_ClusEtaPV_layer3->Fill(hit->Eta());
                hM_ClusEtaOri_layer3->Fill(ClusEta->at(i));
                hM_ClusPhi_layer3->Fill(hit->Phi());
                hM_ClusADC_layer3->Fill(ClusAdc->at(i));
                hM_ClusZSize_layer3->Fill(ClusZSize->at(i));
                hM_ClusPhiSize_layer3->Fill(ClusPhiSize->at(i));
                hM_ClusZ_ClusPhi_layer3->Fill(hit->posZ(), hit->Phi());
                hM_ClusEta_ClusZSize_layer3->Fill(hit->Eta(), ClusZSize->at(i));
                hM_ClusPhi_ClusPhiSize_layer3->Fill(hit->Phi(), ClusPhiSize->at(i));
                hM_ClusZSize_ClusPhiSize_layer3->Fill(ClusZSize->at(i), ClusPhiSize->at(i));
                hM_ClusEta_ClusADC_layer3->Fill(hit->Eta(), ClusAdc->at(i));
                hM_ClusZSize_ClusADC_layer3->Fill(ClusZSize->at(i), ClusAdc->at(i));
                hM_ClusPhiSize_ClusADC_layer3->Fill(ClusPhiSize->at(i), ClusAdc->at(i));
            }
            else
            {
                cout << "[WARNING] ClusLayer = " << ClusLayer->at(i) << ", which is not MVTX. Check!" << endl;
                continue;
            }
        }

        // Fillhist_CluspairMinDr(l1hits, hM_ClusPairDr_layer1);
        // Fillhist_CluspairMinDr(l2hits, hM_ClusPairDr_layer2);
        // Fillhist_CluspairMinDr(l3hits, hM_ClusPairDr_layer3);
        Fillhist_CluspairDr(l1hits, hM_ClusPairDr_layer1);
        Fillhist_CluspairDr(l2hits, hM_ClusPairDr_layer2);
        Fillhist_CluspairDr(l3hits, hM_ClusPairDr_layer3);
    }

    TFile *fout = new TFile(outfilename, "RECREATE");
    fout->cd();
    hM_ClusX_all->Write();
    hM_ClusX_layer1->Write();
    hM_ClusX_layer2->Write();
    hM_ClusX_layer3->Write();
    hM_ClusY_all->Write();
    hM_ClusY_layer1->Write();
    hM_ClusY_layer2->Write();
    hM_ClusY_layer3->Write();
    hM_ClusZ_all->Write();
    hM_ClusZ_layer1->Write();
    hM_ClusZ_layer2->Write();
    hM_ClusZ_layer3->Write();
    hM_ClusR_all->Write();
    hM_ClusR_layer1->Write();
    hM_ClusR_layer2->Write();
    hM_ClusR_layer3->Write();
    hM_ClusEtaPV_all->Write();
    hM_ClusEtaPV_layer1->Write();
    hM_ClusEtaPV_layer2->Write();
    hM_ClusEtaPV_layer3->Write();
    hM_ClusEtaOri_all->Write();
    hM_ClusEtaOri_layer1->Write();
    hM_ClusEtaOri_layer2->Write();
    hM_ClusEtaOri_layer3->Write();
    hM_ClusPhi_all->Write();
    hM_ClusPhi_layer1->Write();
    hM_ClusPhi_layer2->Write();
    hM_ClusPhi_layer3->Write();
    hM_ClusADC_all->Write();
    hM_ClusADC_layer1->Write();
    hM_ClusADC_layer2->Write();
    hM_ClusADC_layer3->Write();
    hM_ClusZSize_all->Write();
    hM_ClusZSize_layer1->Write();
    hM_ClusZSize_layer2->Write();
    hM_ClusZSize_layer3->Write();
    hM_ClusPhiSize_all->Write();
    hM_ClusPhiSize_layer1->Write();
    hM_ClusPhiSize_layer2->Write();
    hM_ClusPhiSize_layer3->Write();
    hM_ClusPairDr_layer1->Write();
    hM_ClusPairDr_layer2->Write();
    hM_ClusPairDr_layer3->Write();
    hM_ClusX_ClusY_all->Write();
    hM_ClusZ_ClusPhi_all->Write();
    hM_ClusZ_ClusPhi_layer1->Write();
    hM_ClusZ_ClusPhi_layer2->Write();
    hM_ClusZ_ClusPhi_layer3->Write();
    hM_ClusEta_ClusZSize_all->Write();
    hM_ClusEta_ClusZSize_layer1->Write();
    hM_ClusEta_ClusZSize_layer2->Write();
    hM_ClusEta_ClusZSize_layer3->Write();
    hM_ClusPhi_ClusPhiSize_all->Write();
    hM_ClusPhi_ClusPhiSize_layer1->Write();
    hM_ClusPhi_ClusPhiSize_layer2->Write();
    hM_ClusPhi_ClusPhiSize_layer3->Write();
    hM_ClusEta_ClusADC_all->Write();
    hM_ClusEta_ClusADC_layer1->Write();
    hM_ClusEta_ClusADC_layer2->Write();
    hM_ClusEta_ClusADC_layer3->Write();
    hM_ClusZSize_ClusPhiSize_all->Write();
    hM_ClusZSize_ClusPhiSize_layer1->Write();
    hM_ClusZSize_ClusPhiSize_layer2->Write();
    hM_ClusZSize_ClusPhiSize_layer3->Write();
    hM_ClusZSize_ClusADC_all->Write();
    hM_ClusZSize_ClusADC_layer1->Write();
    hM_ClusZSize_ClusADC_layer2->Write();
    hM_ClusZSize_ClusADC_layer3->Write();
    hM_ClusPhiSize_ClusADC_all->Write();
    hM_ClusPhiSize_ClusADC_layer1->Write();
    hM_ClusPhiSize_ClusADC_layer2->Write();
    hM_ClusPhiSize_ClusADC_layer3->Write();
    fout->Close();

    system(Form("cd ./plot/; python plotRecoCluster.py -f %s", outfilename.Data()));
}