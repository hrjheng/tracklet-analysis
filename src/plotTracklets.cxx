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

const int NBins = 50;
float edges[NBins + 1] = {1.00000000e-04, 1.18571250e-04, 1.40591414e-04, 1.66700998e-04, 1.97659458e-04, 2.34367291e-04, 2.77892228e-04, 3.29500290e-04, 3.90692614e-04, 4.63249118e-04, 5.49280272e-04, 6.51288487e-04, 7.72240903e-04, 9.15655696e-04, 1.08570441e-03, 1.28733329e-03, 1.52640718e-03,
                          1.80988009e-03, 2.14599745e-03, 2.54453601e-03, 3.01708817e-03, 3.57739917e-03, 4.24176693e-03, 5.02951609e-03, 5.96356012e-03, 7.07106781e-03, 8.38425353e-03, 9.94131425e-03, 1.17875406e-02, 1.39766343e-02, 1.65722701e-02, 1.96499479e-02, 2.32991889e-02, 2.76261397e-02,
                          3.27566592e-02, 3.88399805e-02, 4.60530506e-02, 5.46056779e-02, 6.47466352e-02, 7.67708949e-02, 9.10282102e-02, 1.07933287e-01, 1.27977848e-01, 1.51744935e-01, 1.79925867e-01, 2.13340350e-01, 2.52960321e-01, 2.99938216e-01, 3.55640493e-01, 4.21687380e-01, 5.00000000e-01};
vector<int> layer = {12, 23, 13};

void hist_CluspairDr(TString tkl12treename, TString tkl23treename, TFile *histfile)
{
    TH1F *hM_ClusPairDr_layer1 = new TH1F("hM_ClusPairDr_layer1", "hM_ClusPairDr_layer1", 100, 0, 0.1);
    TH1F *hM_ClusPairDr_layer2 = new TH1F("hM_ClusPairDr_layer2", "hM_ClusPairDr_layer2", 100, 0, 0.1);
    TH1F *hM_ClusPairDr_layer3 = new TH1F("hM_ClusPairDr_layer3", "hM_ClusPairDr_layer3", 100, 0, 0.1);

    TFile *f_tkl12minitree = new TFile(tkl12treename, "READ");
    TTree *t12 = (TTree *)f_tkl12minitree->Get("minitree_12");
    t12->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index12 = (TTreeIndex *)t12->GetTreeIndex();
    vector<float> *recoclus_eta_l1 = 0, *recoclus_phi_l1 = 0, *recoclus_eta_l2 = 0, *recoclus_phi_l2 = 0;
    t12->SetBranchAddress("recoclus_eta_l1", &recoclus_eta_l1);
    t12->SetBranchAddress("recoclus_phi_l1", &recoclus_phi_l1);
    t12->SetBranchAddress("recoclus_eta_l2", &recoclus_eta_l2);
    t12->SetBranchAddress("recoclus_phi_l2", &recoclus_phi_l2);
    for (int ev = 0; ev < t12->GetEntriesFast(); ev++)
    {
        Long64_t local = t12->LoadTree(index12->GetIndex()[ev]);
        t12->GetEntry(local);
        for (size_t i = 0; i < recoclus_eta_l1->size(); i++)
        {
            for (int j = i + 1; j < recoclus_eta_l1->size(); j++)
            {
                float dr = deltaR(recoclus_eta_l1->at(i), recoclus_phi_l1->at(i), recoclus_eta_l1->at(j), recoclus_phi_l1->at(j));
                hM_ClusPairDr_layer1->Fill(dr);
            }
        }
        for (size_t i = 0; i < recoclus_eta_l2->size(); i++)
        {
            for (int j = i + 1; j < recoclus_eta_l2->size(); j++)
            {
                float dr = deltaR(recoclus_eta_l2->at(i), recoclus_phi_l2->at(i), recoclus_eta_l2->at(j), recoclus_phi_l2->at(j));
                hM_ClusPairDr_layer2->Fill(dr);
            }
        }
    }
    f_tkl12minitree->Close();

    TFile *f_tkl23minitree = new TFile(tkl23treename, "READ");
    TTree *t23 = (TTree *)f_tkl23minitree->Get("minitree_23");
    t23->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index23 = (TTreeIndex *)t23->GetTreeIndex();
    vector<float> *recoclus_eta_l3 = 0, *recoclus_phi_l3 = 0;
    t23->SetBranchAddress("recoclus_eta_l3", &recoclus_eta_l3);
    t23->SetBranchAddress("recoclus_phi_l3", &recoclus_phi_l3);
    for (int ev = 0; ev < t23->GetEntriesFast(); ev++)
    {
        Long64_t local = t23->LoadTree(index23->GetIndex()[ev]);
        t23->GetEntry(local);
        for (size_t i = 0; i < recoclus_eta_l3->size(); i++)
        {
            for (int j = i + 1; j < recoclus_eta_l3->size(); j++)
            {
                float dr = deltaR(recoclus_eta_l3->at(i), recoclus_phi_l3->at(i), recoclus_eta_l3->at(j), recoclus_phi_l3->at(j));
                hM_ClusPairDr_layer3->Fill(dr);
            }
        }
    }
    f_tkl23minitree->Close();

    histfile->cd();
    hM_ClusPairDr_layer1->Write();
    hM_ClusPairDr_layer2->Write();
    hM_ClusPairDr_layer3->Write();
}

void makehist(int randhitcase, int clussplitcase, int misalignnum, TString drcut)
{
    // vector<TString> infilename = {Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer12_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
    //                               Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer23_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
    //                               Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer13_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data())};
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data());
    vector<TString> infilename = {Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer12_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth2p9.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
                                  Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer23_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth2p9.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
                                  Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer13_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth2p9.root", randhitcase, clussplitcase, misalignnum, drcut.Data())};
    TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth2p9.root", randhitcase, clussplitcase, misalignnum, drcut.Data());
    // vector<TString> infilename = {Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer12_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth3p5.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
    //                               Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer23_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth3p5.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
    //                               Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer13_Evt0to2000_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth3p5.root", randhitcase, clussplitcase, misalignnum, drcut.Data())};
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_Tracklets_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth3p5.root", randhitcase, clussplitcase, misalignnum, drcut.Data());
    // vector<TString> infilename = {Form("/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer12_Evt0to500_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
    //                               Form("/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer23_Evt0to500_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data()),
    //                               Form("/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer13_Evt0to500_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data())};
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/Hists_Tracklets_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", randhitcase, clussplitcase, misalignnum, drcut.Data());

    TFile *fout = new TFile(outfilename, "RECREATE");
    hist_CluspairDr(infilename[0], infilename[1], fout);
    for (size_t i = 0; i < layer.size(); i++)
    {
        TH1F *hM_dEta_proto = new TH1F(Form("hM_dEta_proto_layer%d", layer[i]), Form("hM_dEta_proto_layer%d", layer[i]), 100, -0.5, 0.5);
        TH1F *hM_dEta_proto_altrange = new TH1F(Form("hM_dEta_proto_altrange_layer%d", layer[i]), Form("hM_dEta_proto_altrange_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_proto_altrange2 = new TH1F(Form("hM_dEta_proto_altrange2_layer%d", layer[i]), Form("hM_dEta_proto_altrange2_layer%d", layer[i]), 100, -0.01, 0.01);
        TH1F *hM_dPhi_proto = new TH1F(Form("hM_dPhi_proto_layer%d", layer[i]), Form("hM_dPhi_proto_layer%d", layer[i]), 100, -0.5, 0.5);
        TH1F *hM_dPhi_proto_altrange = new TH1F(Form("hM_dPhi_proto_altrange_layer%d", layer[i]), Form("hM_dPhi_proto_altrange_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dR_proto = new TH1F(Form("hM_dR_proto_layer%d", layer[i]), Form("hM_dR_proto_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_proto_altrange = new TH1F(Form("hM_dR_proto_altrange_layer%d", layer[i]), Form("hM_dR_proto_altrange_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_proto_LogX = new TH1F(Form("hM_dR_proto_LogX_layer%d", layer[i]), Form("hM_dR_proto_LogX_layer%d", layer[i]), NBins, edges);
        TH1F *hM_Eta_proto = new TH1F(Form("hM_Eta_proto_layer%d", layer[i]), Form("hM_Eta_proto_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Phi_proto = new TH1F(Form("hM_Phi_proto_layer%d", layer[i]), Form("hM_Phi_proto_layer%d", layer[i]), 140, -3.5, 3.5);
        TH1F *hM_NTracklet_proto = new TH1F(Form("hM_NTracklet_proto_layer%d", layer[i]), Form("hM_NTracklet_proto_layer%d", layer[i]), 100, 0, 5000);
        TH1F *hM_Eta_proto_PVzRange1 = new TH1F(Form("hM_Eta_proto_PVzRange1_layer%d", layer[i]), Form("hM_Eta_proto_PVzRange1_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_proto_PVzRange2 = new TH1F(Form("hM_Eta_proto_PVzRange2_layer%d", layer[i]), Form("hM_Eta_proto_PVzRange2_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_proto_PVzRange3 = new TH1F(Form("hM_Eta_proto_PVzRange3_layer%d", layer[i]), Form("hM_Eta_proto_PVzRange3_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_proto_PVzRange4 = new TH1F(Form("hM_Eta_proto_PVzRange4_layer%d", layer[i]), Form("hM_Eta_proto_PVzRange4_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_proto_PVzRange5 = new TH1F(Form("hM_Eta_proto_PVzRange5_layer%d", layer[i]), Form("hM_Eta_proto_PVzRange5_layer%d", layer[i]), 160, -4, 4);

        TH1F *hM_dEta_reco = new TH1F(Form("hM_dEta_reco_layer%d", layer[i]), Form("hM_dEta_reco_layer%d", layer[i]), 100, -0.5, 0.5);
        TH1F *hM_dEta_reco_altrange = new TH1F(Form("hM_dEta_reco_altrange_layer%d", layer[i]), Form("hM_dEta_reco_altrange_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_reco_altrange2 = new TH1F(Form("hM_dEta_reco_altrange2_layer%d", layer[i]), Form("hM_dEta_reco_altrange2_layer%d", layer[i]), 100, -0.01, 0.01);
        TH1F *hM_dPhi_reco = new TH1F(Form("hM_dPhi_reco_layer%d", layer[i]), Form("hM_dPhi_reco_layer%d", layer[i]), 100, -0.5, 0.5);
        TH1F *hM_dPhi_reco_altrange = new TH1F(Form("hM_dPhi_reco_altrange_layer%d", layer[i]), Form("hM_dPhi_reco_altrange_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dR_reco = new TH1F(Form("hM_dR_reco_layer%d", layer[i]), Form("hM_dR_reco_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_reco_altrange = new TH1F(Form("hM_dR_reco_altrange_layer%d", layer[i]), Form("hM_dR_reco_altrange_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_reco_LogX = new TH1F(Form("hM_dR_reco_LogX_layer%d", layer[i]), Form("hM_dR_reco_LogX_layer%d", layer[i]), NBins, edges);
        TH1F *hM_Eta_reco = new TH1F(Form("hM_Eta_reco_layer%d", layer[i]), Form("hM_Eta_reco_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Phi_reco = new TH1F(Form("hM_Phi_reco_layer%d", layer[i]), Form("hM_Phi_reco_layer%d", layer[i]), 140, -3.5, 3.5);
        TH1F *hM_NTracklet_reco = new TH1F(Form("hM_NTracklet_reco_layer%d", layer[i]), Form("hM_NTracklet_reco_layer%d", layer[i]), 100, 0, 5000);
        TH1F *hM_Eta_reco_PVzRange1 = new TH1F(Form("hM_Eta_reco_PVzRange1_layer%d", layer[i]), Form("hM_Eta_reco_PVzRange1_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_reco_PVzRange2 = new TH1F(Form("hM_Eta_reco_PVzRange2_layer%d", layer[i]), Form("hM_Eta_reco_PVzRange2_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_reco_PVzRange3 = new TH1F(Form("hM_Eta_reco_PVzRange3_layer%d", layer[i]), Form("hM_Eta_reco_PVzRange3_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_reco_PVzRange4 = new TH1F(Form("hM_Eta_reco_PVzRange4_layer%d", layer[i]), Form("hM_Eta_reco_PVzRange4_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Eta_reco_PVzRange5 = new TH1F(Form("hM_Eta_reco_PVzRange5_layer%d", layer[i]), Form("hM_Eta_reco_PVzRange5_layer%d", layer[i]), 160, -4, 4);

        TH1F *hM_dEta_recogmatch_ghadPtIncl = new TH1F(Form("hM_dEta_recogmatch_ghadPtIncl_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPtIncl_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPt0to0p5 = new TH1F(Form("hM_dEta_recogmatch_ghadPt0to0p5_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPt0to0p5_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPt0p5to1 = new TH1F(Form("hM_dEta_recogmatch_ghadPt0p5to1_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPt0p5to1_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPt1to2 = new TH1F(Form("hM_dEta_recogmatch_ghadPt1to2_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPt1to2_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPt2to3 = new TH1F(Form("hM_dEta_recogmatch_ghadPt2to3_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPt2to3_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPtg3 = new TH1F(Form("hM_dEta_recogmatch_ghadPtg3_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPtg3_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPtIncl_altrange = new TH1F(Form("hM_dEta_recogmatch_ghadPtIncl_altrange_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPtIncl_altrange_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dEta_recogmatch_ghadPtIncl_altrange2 = new TH1F(Form("hM_dEta_recogmatch_ghadPtIncl_altrange2_layer%d", layer[i]), Form("hM_dEta_recogmatch_ghadPtIncl_altrange2_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPtIncl = new TH1F(Form("hM_dPhi_recogmatch_ghadPtIncl_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPtIncl_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPt0to0p5 = new TH1F(Form("hM_dPhi_recogmatch_ghadPt0to0p5_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPt0to0p5_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPt0p5to1 = new TH1F(Form("hM_dPhi_recogmatch_ghadPt0p5to1_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPt0p5to1_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPt1to2 = new TH1F(Form("hM_dPhi_recogmatch_ghadPt1to2_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPt1to2_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPt2to3 = new TH1F(Form("hM_dPhi_recogmatch_ghadPt2to3_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPt2to3_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPtg3 = new TH1F(Form("hM_dPhi_recogmatch_ghadPtg3_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPtg3_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dPhi_recogmatch_ghadPtIncl_altrange = new TH1F(Form("hM_dPhi_recogmatch_ghadPtIncl_altrange_layer%d", layer[i]), Form("hM_dPhi_recogmatch_ghadPtIncl_altrange_layer%d", layer[i]), 100, -0.05, 0.05);
        TH1F *hM_dR_recogmatch_ghadPtIncl = new TH1F(Form("hM_dR_recogmatch_ghadPtIncl_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPtIncl_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_recogmatch_ghadPt0to0p5 = new TH1F(Form("hM_dR_recogmatch_ghadPt0to0p5_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt0to0p5_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_recogmatch_ghadPt0p5to1 = new TH1F(Form("hM_dR_recogmatch_ghadPt0p5to1_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt0p5to1_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_recogmatch_ghadPt1to2 = new TH1F(Form("hM_dR_recogmatch_ghadPt1to2_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt1to2_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_recogmatch_ghadPt2to3 = new TH1F(Form("hM_dR_recogmatch_ghadPt2to3_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt2to3_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_recogmatch_ghadPtg3 = new TH1F(Form("hM_dR_recogmatch_ghadPtg3_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPtg3_layer%d", layer[i]), 50, 0, 0.05);
        TH1F *hM_dR_recogmatch_ghadPtIncl_LogX = new TH1F(Form("hM_dR_recogmatch_ghadPtIncl_LogX_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPtIncl_LogX_layer%d", layer[i]), NBins, edges);
        TH1F *hM_dR_recogmatch_ghadPtIncl_fr = new TH1F(Form("hM_dR_recogmatch_ghadPtIncl_fr_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPtIncl_fr_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_recogmatch_ghadPt0to0p5_fr = new TH1F(Form("hM_dR_recogmatch_ghadPt0to0p5_fr_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt0to0p5_fr_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_recogmatch_ghadPt0p5to1_fr = new TH1F(Form("hM_dR_recogmatch_ghadPt0p5to1_fr_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt0p5to1_fr_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_recogmatch_ghadPt1to2_fr = new TH1F(Form("hM_dR_recogmatch_ghadPt1to2_fr_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt1to2_fr_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_recogmatch_ghadPt2to3_fr = new TH1F(Form("hM_dR_recogmatch_ghadPt2to3_fr_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPt2to3_fr_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_dR_recogmatch_ghadPtg3_fr = new TH1F(Form("hM_dR_recogmatch_ghadPtg3_fr_layer%d", layer[i]), Form("hM_dR_recogmatch_ghadPtg3_fr_layer%d", layer[i]), 100, 0, 0.5);
        TH1F *hM_Eta_recogmatch = new TH1F(Form("hM_Eta_recogmatch_layer%d", layer[i]), Form("hM_Eta_recogmatch_layer%d", layer[i]), 160, -4, 4);
        TH1F *hM_Phi_recogmatch = new TH1F(Form("hM_Phi_recogmatch_layer%d", layer[i]), Form("hM_Phi_recogmatch_layer%d", layer[i]), 140, -3.5, 3.5);
        TH1F *hM_NTracklet_recogmatch = new TH1F(Form("hM_NTracklet_recogmatch_layer%d", layer[i]), Form("hM_NTracklet_recogmatch_layer%d", layer[i]), 100, 0, 5000);

        TH2F *hM_hit1Eta_hit2Eta_proto = new TH2F(Form("hM_hit1Eta_hit2Eta_proto_layer%d", layer[i]), Form("hM_hit1Eta_hit2Eta_proto_layer%d", layer[i]), 160, -4, 4, 160, -4, 4);
        TH2F *hM_hit1Phi_hit2Phi_proto = new TH2F(Form("hM_hit1Phi_hit2Phi_proto_layer%d", layer[i]), Form("hM_hit1Phi_hit2Phi_proto_layer%d", layer[i]), 140, -3.5, 3.5, 140, -3.5, 3.5);
        TH2F *hM_dPhi_dEta_proto = new TH2F(Form("hM_dPhi_dEta_proto_layer%d", layer[i]), Form("hM_dPhi_dEta_proto_layer%d", layer[i]), 120, -0.6, 0.6, 120, -0.6, 0.6);
        TH2F *hM_dPhi_dEta_proto_altrange = new TH2F(Form("hM_dPhi_dEta_proto_altrange_layer%d", layer[i]), Form("hM_dPhi_dEta_proto_altrange_layer%d", layer[i]), 100, -0.05, 0.05, 100, -0.05, 0.05);
        TH2F *hM_NTracklet_NClusLayer1_proto = new TH2F(Form("hM_NTracklet_NClusLayer1_proto_layer%d", layer[i]), Form("hM_NTracklet_NClusLayer1_proto_layer%d", layer[i]), 250, 0, 5000, 250, 0, 5000);

        TH2F *hM_hit1Eta_hit2Eta_reco = new TH2F(Form("hM_hit1Eta_hit2Eta_reco_layer%d", layer[i]), Form("hM_hit1Eta_hit2Eta_reco_layer%d", layer[i]), 160, -4, 4, 160, -4, 4);
        TH2F *hM_hit1Phi_hit2Phi_reco = new TH2F(Form("hM_hit1Phi_hit2Phi_reco_layer%d", layer[i]), Form("hM_hit1Phi_hit2Phi_reco_layer%d", layer[i]), 140, -3.5, 3.5, 140, -3.5, 3.5);
        TH2F *hM_dPhi_dEta_reco = new TH2F(Form("hM_dPhi_dEta_reco_layer%d", layer[i]), Form("hM_dPhi_dEta_reco_layer%d", layer[i]), 120, -0.6, 0.6, 120, -0.6, 0.6);
        TH2F *hM_dPhi_dEta_reco_altrange = new TH2F(Form("hM_dPhi_dEta_reco_altrange_layer%d", layer[i]), Form("hM_dPhi_dEta_reco_altrange_layer%d", layer[i]), 100, -0.05, 0.05, 100, -0.05, 0.05);
        TH2F *hM_NTracklet_NClusLayer1_reco = new TH2F(Form("hM_NTracklet_NClusLayer1_reco_layer%d", layer[i]), Form("hM_NTracklet_NClusLayer1_reco_layer%d", layer[i]), 250, 0, 5000, 250, 0, 5000);

        TH2F *hM_NTracklet_NClusLayer1_recogmatch = new TH2F(Form("hM_NTracklet_NClusLayer1_recogmatch_layer%d", layer[i]), Form("hM_NTracklet_NClusLayer1_recogmatch_layer%d", layer[i]), 250, 0, 5000, 250, 0, 5000);

        TH2F *hM_Eta_vtxZ_proto_incl = new TH2F(Form("hM_Eta_vtxZ_proto_incl_layer%d", layer[i]), Form("hM_Eta_vtxZ_proto_incl_layer%d", layer[i]), 280, -3.5, 3.5, 240, -12, 12);
        TH2F *hM_Eta_vtxZ_reco_incl = new TH2F(Form("hM_Eta_vtxZ_reco_incl_layer%d", layer[i]), Form("hM_Eta_vtxZ_reco_incl_layer%d", layer[i]), 280, -3.5, 3.5, 240, -12, 12);
        TH2F *hM_Eta_vtxZ_recogmatch_incl = new TH2F(Form("hM_Eta_vtxZ_recogmatch_incl_layer%d", layer[i]), Form("hM_Eta_vtxZ_recogmatch_incl_layer%d", layer[i]), 280, -3.5, 3.5, 240, -12, 12);

        TH2F *hM_GenmatchedRecotkldR_GenhadPt = new TH2F(Form("hM_GenmatchedRecotkldR_GenhadPt_layer%d", layer[i]), Form("hM_GenmatchedRecotkldR_GenhadPt_layer%d", layer[i]), 100, 0, 0.5, 100, 0, 5);
        TH2F *hM_GenmatchedRecotkldR_GenhadPt_altrange2 = new TH2F(Form("hM_GenmatchedRecotkldR_GenhadPt_altrange2_layer%d", layer[i]), Form("hM_GenmatchedRecotkldR_GenhadPt_altrange2_layer%d", layer[i]), 100, 0, 0.05, 100, 0, 5);

        TFile *f = new TFile(infilename[i], "READ");
        TTree *t = (TTree *)f->Get(Form("minitree_%d", layer[i]));
        t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
        TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
        int event, NhitsLayer1, NPrototkl, NRecotkl_Raw, NRecotkl_GenMatched, NGenHadron;
        float PV_z, TruthPV_trig_z;
        vector<float> *prototkl_eta = 0, *prototkl_phi = 0, *prototkl_clus2eta = 0, *prototkl_clus2phi = 0, *prototkl_deta = 0, *prototkl_dphi = 0, *prototkl_dR = 0;
        vector<float> *recotklraw_eta = 0, *recotklraw_phi = 0, *recotklraw_clus2eta = 0, *recotklraw_clus2phi = 0, *recotklraw_deta = 0, *recotklraw_dphi = 0, *recotklraw_dR = 0;
        vector<float> *recotklgm_eta = 0, *recotklgm_phi = 0, *recotklgm_clus2eta = 0, *recotklgm_clus2phi = 0, *recotklgm_deta = 0, *recotklgm_dphi = 0, *recotklgm_dR = 0;
        vector<float> *GenHadron_matched_Pt = 0, *GenHadron_matched_eta = 0, *GenHadron_matched_phi = 0, *GenHadron_matched_E = 0;
        t->SetBranchAddress("event", &event);
        t->SetBranchAddress("NhitsLayer1", &NhitsLayer1);
        t->SetBranchAddress("NPrototkl", &NPrototkl);
        t->SetBranchAddress("NRecotkl_Raw", &NRecotkl_Raw);
        t->SetBranchAddress("NRecotkl_GenMatched", &NRecotkl_GenMatched);
        t->SetBranchAddress("NGenHadron", &NGenHadron);
        t->SetBranchAddress("PV_z", &PV_z);
        t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
        t->SetBranchAddress("prototkl_eta", &prototkl_eta);
        t->SetBranchAddress("prototkl_phi", &prototkl_phi);
        t->SetBranchAddress("prototkl_clus2eta", &prototkl_clus2eta);
        t->SetBranchAddress("prototkl_clus2phi", &prototkl_clus2phi);
        t->SetBranchAddress("prototkl_deta", &prototkl_deta);
        t->SetBranchAddress("prototkl_dphi", &prototkl_dphi);
        t->SetBranchAddress("prototkl_dR", &prototkl_dR);
        t->SetBranchAddress("recotklraw_eta", &recotklraw_eta);
        t->SetBranchAddress("recotklraw_phi", &recotklraw_phi);
        t->SetBranchAddress("recotklraw_clus2eta", &recotklraw_clus2eta);
        t->SetBranchAddress("recotklraw_clus2phi", &recotklraw_clus2phi);
        t->SetBranchAddress("recotklraw_deta", &recotklraw_deta);
        t->SetBranchAddress("recotklraw_dphi", &recotklraw_dphi);
        t->SetBranchAddress("recotklraw_dR", &recotklraw_dR);
        t->SetBranchAddress("recotklgm_eta", &recotklgm_eta);
        t->SetBranchAddress("recotklgm_phi", &recotklgm_phi);
        t->SetBranchAddress("recotklgm_clus2eta", &recotklgm_clus2eta);
        t->SetBranchAddress("recotklgm_clus2phi", &recotklgm_clus2phi);
        t->SetBranchAddress("recotklgm_deta", &recotklgm_deta);
        t->SetBranchAddress("recotklgm_dphi", &recotklgm_dphi);
        t->SetBranchAddress("recotklgm_dR", &recotklgm_dR);
        t->SetBranchAddress("GenHadron_matched_Pt", &GenHadron_matched_Pt);
        t->SetBranchAddress("GenHadron_matched_eta", &GenHadron_matched_eta);
        t->SetBranchAddress("GenHadron_matched_phi", &GenHadron_matched_phi);
        t->SetBranchAddress("GenHadron_matched_E", &GenHadron_matched_E);
        for (int ev = 0; ev < t->GetEntriesFast(); ev++)
        {
            Long64_t local = t->LoadTree(index->GetIndex()[ev]);
            t->GetEntry(local);

            if (PV_z < -99.)
                continue;

            cout << "layer=" << layer[i] << "; event=" << event << "; NhitsLayer1=" << NhitsLayer1 << "; NPrototkl=" << NPrototkl << "; NRecotkl_Raw=" << NRecotkl_Raw << "; NRecotkl_GenMatched=" << NRecotkl_GenMatched << "; NGenHadron=" << NGenHadron << endl;

            hM_NTracklet_proto->Fill(NPrototkl);
            hM_NTracklet_NClusLayer1_proto->Fill(NPrototkl, NhitsLayer1);
            for (size_t j = 0; j < prototkl_eta->size(); j++)
            {
                hM_dEta_proto->Fill(prototkl_deta->at(j));
                hM_dEta_proto_altrange->Fill(prototkl_deta->at(j));
                hM_dEta_proto_altrange2->Fill(prototkl_deta->at(j));
                hM_dPhi_proto->Fill(prototkl_dphi->at(j));
                hM_dPhi_proto_altrange->Fill(prototkl_dphi->at(j));
                hM_dR_proto->Fill(prototkl_dR->at(j));
                hM_dR_proto_altrange->Fill(prototkl_dR->at(j));
                hM_dR_proto_LogX->Fill(prototkl_dR->at(j));
                hM_Eta_proto->Fill(prototkl_eta->at(j));
                hM_Phi_proto->Fill(prototkl_phi->at(j));

                hM_hit1Eta_hit2Eta_proto->Fill(prototkl_eta->at(j), prototkl_clus2eta->at(j));
                hM_hit1Phi_hit2Phi_proto->Fill(prototkl_phi->at(j), prototkl_clus2phi->at(j));
                hM_dPhi_dEta_proto->Fill(prototkl_dphi->at(j), prototkl_deta->at(j));
                hM_dPhi_dEta_proto_altrange->Fill(prototkl_dphi->at(j), prototkl_deta->at(j));
                hM_Eta_vtxZ_proto_incl->Fill(prototkl_eta->at(j), PV_z);

                if (PV_z < -7)
                    hM_Eta_proto_PVzRange1->Fill(prototkl_eta->at(j));
                else if (PV_z >= -7 && PV_z < -2.5)
                    hM_Eta_proto_PVzRange2->Fill(prototkl_eta->at(j));
                else if (PV_z >= -2.5 && PV_z < 2.5)
                    hM_Eta_proto_PVzRange3->Fill(prototkl_eta->at(j));
                else if (PV_z >= 2.5 && PV_z < 7)
                    hM_Eta_proto_PVzRange4->Fill(prototkl_eta->at(j));
                else if (PV_z >= 7)
                    hM_Eta_proto_PVzRange5->Fill(prototkl_eta->at(j));
                else
                    cout << "Weird! PVz = " << PV_z << endl;
            }

            hM_NTracklet_reco->Fill(NRecotkl_Raw);
            hM_NTracklet_NClusLayer1_reco->Fill(NRecotkl_Raw, NhitsLayer1);
            for (size_t j = 0; j < recotklraw_eta->size(); j++)
            {
                hM_dEta_reco->Fill(recotklraw_deta->at(j));
                hM_dEta_reco_altrange->Fill(recotklraw_deta->at(j));
                hM_dEta_reco_altrange2->Fill(recotklraw_deta->at(j));
                hM_dPhi_reco->Fill(recotklraw_dphi->at(j));
                hM_dPhi_reco_altrange->Fill(recotklraw_dphi->at(j));
                hM_dR_reco->Fill(recotklraw_dR->at(j));
                hM_dR_reco_altrange->Fill(recotklraw_dR->at(j));
                hM_dR_reco_LogX->Fill(recotklraw_dR->at(j));
                hM_Eta_reco->Fill(recotklraw_eta->at(j));
                hM_Phi_reco->Fill(recotklraw_phi->at(j));
                hM_hit1Eta_hit2Eta_reco->Fill(recotklraw_eta->at(j), recotklraw_clus2eta->at(j));
                hM_hit1Phi_hit2Phi_reco->Fill(recotklraw_phi->at(j), recotklraw_clus2phi->at(j));
                hM_dPhi_dEta_reco->Fill(recotklraw_dphi->at(j), recotklraw_deta->at(j));
                hM_dPhi_dEta_reco_altrange->Fill(recotklraw_dphi->at(j), recotklraw_deta->at(j));
                hM_Eta_vtxZ_reco_incl->Fill(recotklraw_eta->at(j), PV_z);

                if (PV_z < -7)
                    hM_Eta_reco_PVzRange1->Fill(recotklraw_eta->at(j));
                else if (PV_z >= -7 && PV_z < -2.5)
                    hM_Eta_reco_PVzRange2->Fill(recotklraw_eta->at(j));
                else if (PV_z >= -2.5 && PV_z < 2.5)
                    hM_Eta_reco_PVzRange3->Fill(recotklraw_eta->at(j));
                else if (PV_z >= 2.5 && PV_z < 7)
                    hM_Eta_reco_PVzRange4->Fill(recotklraw_eta->at(j));
                else if (PV_z >= 7)
                    hM_Eta_reco_PVzRange5->Fill(recotklraw_eta->at(j));
                else
                    cout << "Weird! PVz = " << PV_z << endl;
            }

            if (recotklgm_eta->size() != GenHadron_matched_Pt->size())
                cout << "Size of gmatched recotkl vector != size of matched genhadron. Check!" << endl;

            hM_NTracklet_recogmatch->Fill(NRecotkl_GenMatched);
            hM_NTracklet_NClusLayer1_recogmatch->Fill(NRecotkl_GenMatched, NhitsLayer1);
            for (size_t j = 0; j < recotklgm_eta->size(); j++)
            {
                hM_dEta_recogmatch_ghadPtIncl->Fill(recotklgm_deta->at(j));
                hM_dEta_recogmatch_ghadPtIncl_altrange->Fill(recotklgm_deta->at(j));
                hM_dEta_recogmatch_ghadPtIncl_altrange2->Fill(recotklgm_deta->at(j));
                hM_dPhi_recogmatch_ghadPtIncl->Fill(recotklgm_dphi->at(j));
                hM_dPhi_recogmatch_ghadPtIncl_altrange->Fill(recotklgm_dphi->at(j));
                hM_dR_recogmatch_ghadPtIncl->Fill(recotklgm_dR->at(j));
                hM_dR_recogmatch_ghadPtIncl_LogX->Fill(recotklgm_dR->at(j));
                hM_dR_recogmatch_ghadPtIncl_fr->Fill(recotklgm_dR->at(j));
                hM_Eta_recogmatch->Fill(recotklgm_eta->at(j));
                hM_Phi_recogmatch->Fill(recotklgm_phi->at(j));
                hM_Eta_vtxZ_recogmatch_incl->Fill(recotklgm_eta->at(j), PV_z);
                hM_GenmatchedRecotkldR_GenhadPt->Fill(recotklgm_dR->at(j), GenHadron_matched_Pt->at(j));
                hM_GenmatchedRecotkldR_GenhadPt_altrange2->Fill(recotklgm_dR->at(j), GenHadron_matched_Pt->at(j));

                if (GenHadron_matched_Pt->at(j) > 0 && GenHadron_matched_Pt->at(j) <= 0.5)
                {
                    hM_dR_recogmatch_ghadPt0to0p5->Fill(recotklgm_dR->at(j));
                    hM_dR_recogmatch_ghadPt0to0p5_fr->Fill(recotklgm_dR->at(j));
                    hM_dEta_recogmatch_ghadPt0to0p5->Fill(recotklgm_deta->at(j));
                    hM_dPhi_recogmatch_ghadPt0to0p5->Fill(recotklgm_dphi->at(j));
                }
                else if (GenHadron_matched_Pt->at(j) > 0.5 && GenHadron_matched_Pt->at(j) <= 1)
                {
                    hM_dR_recogmatch_ghadPt0p5to1->Fill(recotklgm_dR->at(j));
                    hM_dR_recogmatch_ghadPt0p5to1_fr->Fill(recotklgm_dR->at(j));
                    hM_dEta_recogmatch_ghadPt0p5to1->Fill(recotklgm_deta->at(j));
                    hM_dPhi_recogmatch_ghadPt0p5to1->Fill(recotklgm_dphi->at(j));
                }
                else if (GenHadron_matched_Pt->at(j) > 1 && GenHadron_matched_Pt->at(j) <= 2)
                {
                    hM_dR_recogmatch_ghadPt1to2->Fill(recotklgm_dR->at(j));
                    hM_dR_recogmatch_ghadPt1to2_fr->Fill(recotklgm_dR->at(j));
                    hM_dEta_recogmatch_ghadPt1to2->Fill(recotklgm_deta->at(j));
                    hM_dPhi_recogmatch_ghadPt1to2->Fill(recotklgm_dphi->at(j));
                }
                else if (GenHadron_matched_Pt->at(j) > 2 && GenHadron_matched_Pt->at(j) <= 3)
                {
                    hM_dR_recogmatch_ghadPt2to3->Fill(recotklgm_dR->at(j));
                    hM_dR_recogmatch_ghadPt2to3_fr->Fill(recotklgm_dR->at(j));
                    hM_dEta_recogmatch_ghadPt2to3->Fill(recotklgm_deta->at(j));
                    hM_dPhi_recogmatch_ghadPt2to3->Fill(recotklgm_dphi->at(j));
                }
                else if (GenHadron_matched_Pt->at(j) >= 3)
                {
                    hM_dR_recogmatch_ghadPtg3->Fill(recotklgm_dR->at(j));
                    hM_dR_recogmatch_ghadPtg3_fr->Fill(recotklgm_dR->at(j));
                    hM_dEta_recogmatch_ghadPtg3->Fill(recotklgm_deta->at(j));
                    hM_dPhi_recogmatch_ghadPtg3->Fill(recotklgm_dphi->at(j));
                }
                else
                {
                    cout << "GenHadron Pt = " << GenHadron_matched_Pt->at(j) << ", weird! Check!" << endl;
                }
            }
        }

        f->Close();

        fout->cd();
        hM_dEta_proto->Write();
        hM_dEta_reco->Write();
        hM_dEta_recogmatch_ghadPtIncl->Write();
        hM_dEta_proto->Write();
        hM_dEta_proto_altrange->Write();
        hM_dEta_proto_altrange2->Write();
        hM_dPhi_proto->Write();
        hM_dPhi_proto_altrange->Write();
        hM_dR_proto->Write();
        hM_dR_proto_altrange->Write();
        hM_dR_proto_LogX->Write();
        hM_Eta_proto->Write();
        hM_Phi_proto->Write();
        hM_NTracklet_proto->Write();
        hM_Eta_proto_PVzRange1->Write();
        hM_Eta_proto_PVzRange2->Write();
        hM_Eta_proto_PVzRange3->Write();
        hM_Eta_proto_PVzRange4->Write();
        hM_Eta_proto_PVzRange5->Write();
        hM_dEta_reco->Write();
        hM_dEta_reco_altrange->Write();
        hM_dEta_reco_altrange2->Write();
        hM_dPhi_reco->Write();
        hM_dPhi_reco_altrange->Write();
        hM_dR_reco->Write();
        hM_dR_reco_altrange->Write();
        hM_dR_reco_LogX->Write();
        hM_Eta_reco->Write();
        hM_Phi_reco->Write();
        hM_NTracklet_reco->Write();
        hM_Eta_reco_PVzRange1->Write();
        hM_Eta_reco_PVzRange2->Write();
        hM_Eta_reco_PVzRange3->Write();
        hM_Eta_reco_PVzRange4->Write();
        hM_Eta_reco_PVzRange5->Write();
        hM_dEta_recogmatch_ghadPtIncl->Write();
        hM_dEta_recogmatch_ghadPt0to0p5->Write();
        hM_dEta_recogmatch_ghadPt0p5to1->Write();
        hM_dEta_recogmatch_ghadPt1to2->Write();
        hM_dEta_recogmatch_ghadPt2to3->Write();
        hM_dEta_recogmatch_ghadPtg3->Write();
        hM_dEta_recogmatch_ghadPtIncl_altrange->Write();
        hM_dEta_recogmatch_ghadPtIncl_altrange2->Write();
        hM_dPhi_recogmatch_ghadPtIncl->Write();
        hM_dPhi_recogmatch_ghadPt0to0p5->Write();
        hM_dPhi_recogmatch_ghadPt0p5to1->Write();
        hM_dPhi_recogmatch_ghadPt1to2->Write();
        hM_dPhi_recogmatch_ghadPt2to3->Write();
        hM_dPhi_recogmatch_ghadPtg3->Write();
        hM_dPhi_recogmatch_ghadPtIncl_altrange->Write();
        hM_dR_recogmatch_ghadPtIncl->Write();
        hM_dR_recogmatch_ghadPt0to0p5->Write();
        hM_dR_recogmatch_ghadPt0p5to1->Write();
        hM_dR_recogmatch_ghadPt1to2->Write();
        hM_dR_recogmatch_ghadPt2to3->Write();
        hM_dR_recogmatch_ghadPtg3->Write();
        hM_dR_recogmatch_ghadPtIncl_LogX->Write();
        hM_dR_recogmatch_ghadPtIncl_fr->Write();
        hM_dR_recogmatch_ghadPt0to0p5_fr->Write();
        hM_dR_recogmatch_ghadPt0p5to1_fr->Write();
        hM_dR_recogmatch_ghadPt1to2_fr->Write();
        hM_dR_recogmatch_ghadPt2to3_fr->Write();
        hM_dR_recogmatch_ghadPtg3_fr->Write();
        hM_Eta_recogmatch->Write();
        hM_Phi_recogmatch->Write();
        hM_NTracklet_recogmatch->Write();
        hM_hit1Eta_hit2Eta_proto->Write();
        hM_hit1Phi_hit2Phi_proto->Write();
        hM_dPhi_dEta_proto->Write();
        hM_dPhi_dEta_proto_altrange->Write();
        hM_NTracklet_NClusLayer1_proto->Write();
        hM_hit1Eta_hit2Eta_reco->Write();
        hM_hit1Phi_hit2Phi_reco->Write();
        hM_dPhi_dEta_reco->Write();
        hM_dPhi_dEta_reco_altrange->Write();
        hM_NTracklet_NClusLayer1_reco->Write();
        hM_NTracklet_NClusLayer1_recogmatch->Write();
        hM_Eta_vtxZ_proto_incl->Write();
        hM_Eta_vtxZ_reco_incl->Write();
        hM_Eta_vtxZ_recogmatch_incl->Write();
        hM_GenmatchedRecotkldR_GenhadPt->Write();
        hM_GenmatchedRecotkldR_GenhadPt_altrange2->Write();
    }

    fout->Close();

    // int randhitcase, int clussplitcase, int misalignnum, TString drcut
    system(Form("cd ./plot/; python plotTracklet.py -f %s -r %d -s %d -m %d -d %s", outfilename.Data(), randhitcase, clussplitcase, misalignnum, drcut.Data()));
}

int main(int argc, char *argv[])
{
    // makehist(int randhitcase, int clussplitcase, int misalignnum, TString drcut)
    // Nominal
    // makehist(0, 0, 0, "0p5");
    // randhitcase
    // for (int i = 1; i < 4; i++)
    //     makehist(i, 0, 0, "0p5");
    // clussplitcase
    // for (int i = 1; i < 4; i++)
    //     makehist(0, i, 0, "0p5");
    // misalignnum
    // for (int i = 1; i < 40; i++)
        // makehist(0, 0, i, "0p5");
    // makehist(0, 0, 100, "0p5");
    makehist(0, 0, 100, "0p5");
    // drcut
    // vector<TString> drcutlist = {"0p4", "0p6", "999"};
    // for (size_t i = 0; i < drcutlist.size(); i++)
    //     makehist(0, 0, 0, drcutlist[i]);

    return 0;
}