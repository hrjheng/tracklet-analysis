#ifndef HISTS_H
#define HISTS_H

#include <TString.h>
#include <TFile.h>

#include "Tracklet.h"
#include "TrackletData.h"

struct Hists
{
    TH1F *hM_PVz;
    TH1F *hM_NhitsLayer1;

    // Clusters
    TH1F *hM_ClusX_all;
    TH1F *hM_ClusX_layer1;
    TH1F *hM_ClusX_layer2;
    TH1F *hM_ClusX_layer3;
    TH1F *hM_ClusY_all;
    TH1F *hM_ClusY_layer1;
    TH1F *hM_ClusY_layer2;
    TH1F *hM_ClusY_layer3;
    TH1F *hM_ClusZ_all;
    TH1F *hM_ClusZ_layer1;
    TH1F *hM_ClusZ_layer2;
    TH1F *hM_ClusZ_layer3;
    TH1F *hM_ClusEta_all;
    TH1F *hM_ClusEta_layer1;
    TH1F *hM_ClusEta_layer2;
    TH1F *hM_ClusEta_layer3;
    TH1F *hM_ClusPhi_all;
    TH1F *hM_ClusPhi_layer1;
    TH1F *hM_ClusPhi_layer2;
    TH1F *hM_ClusPhi_layer3;
    TH2F *hM_ClusX_ClusY_all;
    TH2F *hM_ClusZ_ClusPhi_all;
    TH2F *hM_ClusZ_ClusPhi_layer1;
    TH2F *hM_ClusZ_ClusPhi_layer2;
    TH2F *hM_ClusZ_ClusPhi_layer3;
    TH2F *hM_ClusZ_TruthPVz_all;
    TH2F *hM_ClusZ_TruthPVz_layer1;
    TH2F *hM_ClusZ_TruthPVz_layer2;
    TH2F *hM_ClusZ_TruthPVz_layer3;
    TH2F *hM_ClusPhi_TruthPVz_all;
    TH2F *hM_ClusPhi_TruthPVz_layer1;
    TH2F *hM_ClusPhi_TruthPVz_layer2;
    TH2F *hM_ClusPhi_TruthPVz_layer3;

    // Tracklets
    TH2F *hM_dPhi_dEta_woCuts;
    TH1F *hM_Eta_hit1_proto;
    TH1F *hM_Phi_hit1_proto;
    TH1F *hM_Eta_hit2_proto;
    TH1F *hM_Phi_hit2_proto;
    TH1F *hM_dEta_proto;
    TH1F *hM_dEta_proto_altrange;
    TH1F *hM_dEta_proto_altrange2;
    TH1F *hM_dPhi_proto;
    TH1F *hM_dPhi_proto_altrange;
    TH1F *hM_dR_proto;
    TH1F *hM_dR_proto_altrange;
    TH1F *hM_dR_proto_LogX;
    TH1F *hM_Eta_proto;
    TH1F *hM_Eta_proto_PVzRange1; // PVz < -7 cm
    TH1F *hM_Eta_proto_PVzRange2; // -7 < PVz < -2.5 cm
    TH1F *hM_Eta_proto_PVzRange3; // -2.5 < PVz < 2.5 cm
    TH1F *hM_Eta_proto_PVzRange4; // 2.5 < PVz < 7 cm
    TH1F *hM_Eta_proto_PVzRange5; // PVz > 7 cm
    TH1F *hM_Phi_proto;
    TH2F *hM_hit1Eta_hit2Eta_proto;
    TH2F *hM_hit1Phi_hit2Phi_proto;
    TH2F *hM_Eta_vtxZ_proto;
    TH2F *hM_Eta_vtxZ_anabin_proto;

    TH2F *hM_dPhi_dEta_proto;
    TH2F *hM_dPhi_dEta_proto_altrange;
    TH1F *hM_Eta_hit1_reco;
    TH1F *hM_Phi_hit1_reco;
    TH1F *hM_Eta_hit2_reco;
    TH1F *hM_Phi_hit2_reco;
    TH1F *hM_dEta_reco;
    TH1F *hM_dEta_reco_altrange;
    TH1F *hM_dEta_reco_altrange2;
    TH1F *hM_dPhi_reco;
    TH1F *hM_dPhi_reco_altrange;
    TH1F *hM_dR_reco;
    TH1F *hM_dR_reco_altrange;
    TH1F *hM_dR_reco_LogX;
    TH1F *hM_Eta_reco;
    TH1F *hM_Phi_reco;
    TH2F *hM_hit1Eta_hit2Eta_reco;
    TH2F *hM_hit1Phi_hit2Phi_reco;
    TH2F *hM_Eta_vtxZ_reco;
    TH2F *hM_Eta_vtxZ_anabin_reco;
    TH2F *hM_dPhi_dEta_reco;
    TH2F *hM_dPhi_dEta_reco_altrange;

    // After gen hadron matching
    TH2F *hM_Eta_vtxZ_reco_genmatched;
    TH2F *hM_Eta_vtxZ_anabin_reco_genmatched;
};

void SetupHists(Hists &hists)
{
    const int NBins = 50;
    float edges[NBins + 1] = {
        1.00000000e-04, 1.18571250e-04, 1.40591414e-04, 1.66700998e-04,
        1.97659458e-04, 2.34367291e-04, 2.77892228e-04, 3.29500290e-04,
        3.90692614e-04, 4.63249118e-04, 5.49280272e-04, 6.51288487e-04,
        7.72240903e-04, 9.15655696e-04, 1.08570441e-03, 1.28733329e-03,
        1.52640718e-03, 1.80988009e-03, 2.14599745e-03, 2.54453601e-03,
        3.01708817e-03, 3.57739917e-03, 4.24176693e-03, 5.02951609e-03,
        5.96356012e-03, 7.07106781e-03, 8.38425353e-03, 9.94131425e-03,
        1.17875406e-02, 1.39766343e-02, 1.65722701e-02, 1.96499479e-02,
        2.32991889e-02, 2.76261397e-02, 3.27566592e-02, 3.88399805e-02,
        4.60530506e-02, 5.46056779e-02, 6.47466352e-02, 7.67708949e-02,
        9.10282102e-02, 1.07933287e-01, 1.27977848e-01, 1.51744935e-01,
        1.79925867e-01, 2.13340350e-01, 2.52960321e-01, 2.99938216e-01,
        3.55640493e-01, 4.21687380e-01, 5.00000000e-01};

    hists.hM_PVz = new TH1F("hM_PVz", "hM_PVz", 60, -12, 12);
    hists.hM_NhitsLayer1 = new TH1F("hM_NhitsLayer1", "hM_NhitsLayer1", 200, 0, 20000);

    hists.hM_ClusX_all = new TH1F("hM_ClusX_all", "hM_ClusX_all", 120, -6, 6);
    hists.hM_ClusX_layer1 = new TH1F("hM_ClusX_layer1", "hM_ClusX_layer1", 120, -6, 6);
    hists.hM_ClusX_layer2 = new TH1F("hM_ClusX_layer2", "hM_ClusX_layer2", 120, -6, 6);
    hists.hM_ClusX_layer3 = new TH1F("hM_ClusX_layer3", "hM_ClusX_layer3", 120, -6, 6);
    hists.hM_ClusY_all = new TH1F("hM_ClusY_all", "hM_ClusY_all", 120, -6, 6);
    hists.hM_ClusY_layer1 = new TH1F("hM_ClusY_layer1", "hM_ClusY_layer1", 120, -6, 6);
    hists.hM_ClusY_layer2 = new TH1F("hM_ClusY_layer2", "hM_ClusY_layer2", 120, -6, 6);
    hists.hM_ClusY_layer3 = new TH1F("hM_ClusY_layer3", "hM_ClusY_layer3", 120, -6, 6);
    hists.hM_ClusZ_all = new TH1F("hM_ClusZ_all", "hM_ClusZ_all", 150, -15, 15);
    hists.hM_ClusZ_layer1 = new TH1F("hM_ClusZ_layer1", "hM_ClusZ_layer1", 150, -15, 15);
    hists.hM_ClusZ_layer2 = new TH1F("hM_ClusZ_layer2", "hM_ClusZ_layer2", 150, -15, 15);
    hists.hM_ClusZ_layer3 = new TH1F("hM_ClusZ_layer3", "hM_ClusZ_layer3", 150, -15, 15);
    hists.hM_ClusEta_all = new TH1F("hM_ClusEta_all", "hM_ClusEta_all", 120, -3, 3);
    hists.hM_ClusEta_layer1 = new TH1F("hM_ClusEta_layer1", "hM_ClusEta_layer1", 120, -3, 3);
    hists.hM_ClusEta_layer2 = new TH1F("hM_ClusEta_layer2", "hM_ClusEta_layer2", 120, -3, 3);
    hists.hM_ClusEta_layer3 = new TH1F("hM_ClusEta_layer3", "hM_ClusEta_layer3", 120, -3, 3);
    hists.hM_ClusPhi_all = new TH1F("hM_ClusPhi_all", "hM_ClusPhi_all", 70, -3.5, 3.5);
    hists.hM_ClusPhi_layer1 = new TH1F("hM_ClusPhi_layer1", "hM_ClusPhi_layer1", 70, -3.5, 3.5);
    hists.hM_ClusPhi_layer2 = new TH1F("hM_ClusPhi_layer2", "hM_ClusPhi_layer2", 70, -3.5, 3.5);
    hists.hM_ClusPhi_layer3 = new TH1F("hM_ClusPhi_layer3", "hM_ClusPhi_layer3", 70, -3.5, 3.5);
    hists.hM_ClusX_ClusY_all = new TH2F("hM_ClusX_ClusY_all", "hM_ClusX_ClusY_all", 240, -6, 6, 240, -6, 6);
    hists.hM_ClusZ_ClusPhi_all = new TH2F("hM_ClusZ_ClusPhi_all", "hM_ClusZ_ClusPhi_all", 150, -15, 15, 70, -3.5, 3.5);
    hists.hM_ClusZ_ClusPhi_layer1 = new TH2F("hM_ClusZ_ClusPhi_layer1", "hM_ClusZ_ClusPhi_layer1", 150, -15, 15, 70, -3.5, 3.5);
    hists.hM_ClusZ_ClusPhi_layer2 = new TH2F("hM_ClusZ_ClusPhi_layer2", "hM_ClusZ_ClusPhi_layer2", 150, -15, 15, 70, -3.5, 3.5);
    hists.hM_ClusZ_ClusPhi_layer3 = new TH2F("hM_ClusZ_ClusPhi_layer3", "hM_ClusZ_ClusPhi_layer3", 150, -15, 15, 70, -3.5, 3.5);
    hists.hM_ClusZ_TruthPVz_all = new TH2F("hM_ClusZ_TruthPVz_all", "hM_ClusZ_TruthPVz_all", 150, -15, 15, 120, -12, 12);
    hists.hM_ClusZ_TruthPVz_layer1 = new TH2F("hM_ClusZ_TruthPVz_layer1", "hM_ClusZ_TruthPVz_layer1", 150, -15, 15, 120, -12, 12);
    hists.hM_ClusZ_TruthPVz_layer2 = new TH2F("hM_ClusZ_TruthPVz_layer2", "hM_ClusZ_TruthPVz_layer2", 150, -15, 15, 120, -12, 12);
    hists.hM_ClusZ_TruthPVz_layer3 = new TH2F("hM_ClusZ_TruthPVz_layer3", "hM_ClusZ_TruthPVz_layer3", 150, -15, 15, 120, -12, 12);
    hists.hM_ClusPhi_TruthPVz_all = new TH2F("hM_ClusPhi_TruthPVz_all", "hM_ClusPhi_TruthPVz_all", 70, -3.5, 3.5, 120, -12, 12);
    hists.hM_ClusPhi_TruthPVz_layer1 = new TH2F("hM_ClusPhi_TruthPVz_layer1", "hM_ClusPhi_TruthPVz_layer1", 70, -3.5, 3.5, 120, -12, 12);
    hists.hM_ClusPhi_TruthPVz_layer2 = new TH2F("hM_ClusPhi_TruthPVz_layer2", "hM_ClusPhi_TruthPVz_layer2", 70, -3.5, 3.5, 120, -12, 12);
    hists.hM_ClusPhi_TruthPVz_layer3 = new TH2F("hM_ClusPhi_TruthPVz_layer3", "hM_ClusPhi_TruthPVz_layer3", 70, -3.5, 3.5, 120, -12, 12);

    hists.hM_dPhi_dEta_woCuts = new TH2F("hM_dPhi_dEta_woCuts", "hM_dPhi_dEta_woCuts", 140, -3.5, 3.5, 200, -5, 5);
    hists.hM_Eta_hit1_proto = new TH1F("hM_Eta_hit1_proto", "hM_Eta_hit1_proto", 80, -4, 4);
    hists.hM_Phi_hit1_proto = new TH1F("hM_Phi_hit1_proto", "hM_Phi_hit1_proto", 70, -3.5, 3.5);
    hists.hM_Eta_hit2_proto = new TH1F("hM_Eta_hit2_proto", "hM_Eta_hit2_proto", 80, -4, 4);
    hists.hM_Phi_hit2_proto = new TH1F("hM_Phi_hit2_proto", "hM_Phi_hit2_proto", 70, -3.5, 3.5);
    hists.hM_dEta_proto = new TH1F("hM_dEta_proto", "hM_dEta_proto", 100, -0.5, 0.5);
    hists.hM_dEta_proto_altrange = new TH1F("hM_dEta_proto_altrange", "hM_dEta_proto_altrange", 100, -0.05, 0.05);
    hists.hM_dEta_proto_altrange2 = new TH1F("hM_dEta_proto_altrange2", "hM_dEta_proto_altrange2", 100, -0.01, 0.01);
    hists.hM_dPhi_proto = new TH1F("hM_dPhi_proto", "hM_dPhi_proto", 100, -0.5, 0.5);
    hists.hM_dPhi_proto_altrange = new TH1F("hM_dPhi_proto_altrange", "hM_dPhi_proto_altrange", 100, -0.05, 0.05);
    hists.hM_dR_proto = new TH1F("hM_dR_proto", "hM_dR_proto", 100, 0, 0.5);
    hists.hM_dR_proto_altrange = new TH1F("hM_dR_proto_altrange", "hM_dR_proto_altrange", 50, 0, 0.05);
    hists.hM_dR_proto_LogX = new TH1F("hM_dR_proto_LogX", "hM_dR_proto_LogX", NBins, edges);
    hists.hM_Eta_proto = new TH1F("hM_Eta_proto", "hM_Eta_proto", 80, -4, 4);
    hists.hM_Eta_proto_PVzRange1 = new TH1F("hM_Eta_proto_PVzRange1", "hM_Eta_proto_PVzRange1", 80, -4, 4);
    hists.hM_Eta_proto_PVzRange2 = new TH1F("hM_Eta_proto_PVzRange2", "hM_Eta_proto_PVzRange2", 80, -4, 4);
    hists.hM_Eta_proto_PVzRange3 = new TH1F("hM_Eta_proto_PVzRange3", "hM_Eta_proto_PVzRange3", 80, -4, 4);
    hists.hM_Eta_proto_PVzRange4 = new TH1F("hM_Eta_proto_PVzRange4", "hM_Eta_proto_PVzRange4", 80, -4, 4);
    hists.hM_Eta_proto_PVzRange5 = new TH1F("hM_Eta_proto_PVzRange5", "hM_Eta_proto_PVzRange5", 80, -4, 4);
    hists.hM_Phi_proto = new TH1F("hM_Phi_proto", "hM_Phi_proto", 70, -3.5, 3.5);
    hists.hM_hit1Eta_hit2Eta_proto = new TH2F("hM_hit1Eta_hit2Eta_proto", "hM_hit1Eta_hit2Eta_proto", 80, -4, 4, 80, -4, 4);
    hists.hM_hit1Phi_hit2Phi_proto = new TH2F("hM_hit1Phi_hit2Phi_proto", "hM_hit1Phi_hit2Phi_proto", 70, -3.5, 3.5, 70, -3.5, 3.5);
    hists.hM_dPhi_dEta_proto = new TH2F("hM_dPhi_dEta_proto", "hM_dPhi_dEta_proto", 120, -0.6, 0.6, 120, -0.6, 0.6);
    hists.hM_dPhi_dEta_proto_altrange = new TH2F("hM_dPhi_dEta_proto_altrange", "hM_dPhi_dEta_proto_altrange", 100, -0.05, 0.05, 100, -0.05, 0.05);
    hists.hM_Eta_vtxZ_proto = new TH2F("hM_Eta_vtxZ_proto", "hM_Eta_vtxZ_proto", 140, -3.5, 3.5, 120, -12, 12);
    hists.hM_Eta_vtxZ_anabin_proto = new TH2F("hM_Eta_vtxZ_anabin_proto", "hM_Eta_vtxZ_anabin_proto", 12, -3.0, 3.0, 4, -10, 10);

    hists.hM_Eta_hit1_reco = new TH1F("hM_Eta_hit1_reco", "hM_Eta_hit1_reco", 80, -4, 4);
    hists.hM_Phi_hit1_reco = new TH1F("hM_Phi_hit1_reco", "hM_Phi_hit1_reco", 100, -0.5, 0.5);
    hists.hM_Eta_hit2_reco = new TH1F("hM_Eta_hit2_reco", "hM_Eta_hit2_reco", 80, -4, 4);
    hists.hM_Phi_hit2_reco = new TH1F("hM_Phi_hit2_reco", "hM_Phi_hit2_reco", 100, -0.5, 0.5);
    hists.hM_dEta_reco = new TH1F("hM_dEta_reco", "hM_dEta_reco", 100, -0.5, 0.5);
    hists.hM_dEta_reco_altrange = new TH1F("hM_dEta_reco_altrange", "hM_dEta_reco_altrange", 100, -0.05, 0.05);
    hists.hM_dEta_reco_altrange2 = new TH1F("hM_dEta_reco_altrange2", "hM_dEta_reco_altrange2", 100, -0.01, 0.01);
    hists.hM_dPhi_reco = new TH1F("hM_dPhi_reco", "hM_dPhi_reco", 100, -0.5, 0.5);
    hists.hM_dPhi_reco_altrange = new TH1F("hM_dPhi_reco_altrange", "hM_dPhi_reco_altrange", 100, -0.05, 0.05);
    hists.hM_dR_reco = new TH1F("hM_dR_reco", "hM_dR_reco", 100, 0, 0.5);
    hists.hM_dR_reco_altrange = new TH1F("hM_dR_reco_altrange", "hM_dR_reco_altrange", 50, 0, 0.05);
    hists.hM_dR_reco_LogX = new TH1F("hM_dR_reco_LogX", "hM_dR_reco_LogX", NBins, edges);
    hists.hM_Eta_reco = new TH1F("hM_Eta_reco", "hM_Eta_reco", 80, -4, 4);
    hists.hM_Phi_reco = new TH1F("hM_Phi_reco", "hM_Phi_reco", 70, -3.5, 3.5);
    hists.hM_hit1Eta_hit2Eta_reco = new TH2F("hM_hit1Eta_hit2Eta_reco", "hM_hit1Eta_hit2Eta_reco", 80, -4, 4, 80, -4, 4);
    hists.hM_hit1Phi_hit2Phi_reco = new TH2F("hM_hit1Phi_hit2Phi_reco", "hM_hit1Phi_hit2Phi_reco", 70, -3.5, 3.5, 70, -3.5, 3.5);
    hists.hM_dPhi_dEta_reco = new TH2F("hM_dPhi_dEta_reco", "hM_dPhi_dEta_reco", 120, -0.6, 0.6, 120, -0.6, 0.6);
    hists.hM_dPhi_dEta_reco_altrange = new TH2F("hM_dPhi_dEta_reco_altrange", "hM_dPhi_dEta_reco_altrange", 100, -0.05, 0.05, 100, -0.05, 0.05);
    hists.hM_Eta_vtxZ_reco = new TH2F("hM_Eta_vtxZ_reco", "hM_Eta_vtxZ_reco", 140, -3.5, 3.5, 120, -12, 12);
    hists.hM_Eta_vtxZ_anabin_reco = new TH2F("hM_Eta_vtxZ_anabin_reco", "hM_Eta_vtxZ_anabin_reco", 12, -3.0, 3.0, 4, -10, 10);

    hists.hM_Eta_vtxZ_reco_genmatched = new TH2F("hM_Eta_vtxZ_reco_genmatched", "hM_Eta_vtxZ_reco_genmatched", 140, -3.5, 3.5, 120, -12, 12);
    hists.hM_Eta_vtxZ_anabin_reco_genmatched = new TH2F("hM_Eta_vtxZ_anabin_reco_genmatched", "hM_Eta_vtxZ_anabin_reco_genmatched", 140, -3.5, 3.5, 120, -12, 12);
}

void FillHist_Cluster(vector<vector<Hit *>> &MVTXlayers, float truthPVz, Hists &hists)
{
    for (size_t i = 0; i < MVTXlayers.size(); i++)
    {
        for (auto h : MVTXlayers[i])
        {
            hists.hM_ClusX_all->Fill(h->posX());
            hists.hM_ClusY_all->Fill(h->posY());
            hists.hM_ClusZ_all->Fill(h->posZ());
            hists.hM_ClusEta_all->Fill(h->Eta());
            hists.hM_ClusPhi_all->Fill(h->Phi());
            hists.hM_ClusX_ClusY_all->Fill(h->posX(), h->posY());
            hists.hM_ClusZ_ClusPhi_all->Fill(h->posZ(), h->Phi());
            hists.hM_ClusZ_TruthPVz_all->Fill(h->posZ(), truthPVz);
            hists.hM_ClusPhi_TruthPVz_all->Fill(h->Phi(), truthPVz);
            if (h->Layer() == 0)
            {
                hists.hM_ClusX_layer1->Fill(h->posX());
                hists.hM_ClusY_layer1->Fill(h->posY());
                hists.hM_ClusZ_layer1->Fill(h->posZ());
                hists.hM_ClusEta_layer1->Fill(h->Eta());
                hists.hM_ClusPhi_layer1->Fill(h->Phi());
                hists.hM_ClusZ_ClusPhi_layer1->Fill(h->posZ(), h->Phi());
                hists.hM_ClusZ_TruthPVz_layer1->Fill(h->posZ(), truthPVz);
                hists.hM_ClusPhi_TruthPVz_layer1->Fill(h->Phi(), truthPVz);
            }
            else if (h->Layer() == 1)
            {
                hists.hM_ClusX_layer2->Fill(h->posX());
                hists.hM_ClusY_layer2->Fill(h->posY());
                hists.hM_ClusZ_layer2->Fill(h->posZ());
                hists.hM_ClusEta_layer2->Fill(h->Eta());
                hists.hM_ClusPhi_layer2->Fill(h->Phi());
                hists.hM_ClusZ_ClusPhi_layer2->Fill(h->posZ(), h->Phi());
                hists.hM_ClusZ_TruthPVz_layer2->Fill(h->posZ(), truthPVz);
                hists.hM_ClusPhi_TruthPVz_layer2->Fill(h->Phi(), truthPVz);
            }
            else if (h->Layer() == 2)
            {
                hists.hM_ClusX_layer3->Fill(h->posX());
                hists.hM_ClusY_layer3->Fill(h->posY());
                hists.hM_ClusZ_layer3->Fill(h->posZ());
                hists.hM_ClusEta_layer3->Fill(h->Eta());
                hists.hM_ClusPhi_layer3->Fill(h->Phi());
                hists.hM_ClusZ_ClusPhi_layer3->Fill(h->posZ(), h->Phi());
                hists.hM_ClusZ_TruthPVz_layer3->Fill(h->posZ(), truthPVz);
                hists.hM_ClusPhi_TruthPVz_layer3->Fill(h->Phi(), truthPVz);
            }
            else
            {
                printf("[WARNING] Hit layer = %d, which is not correct. Should check! \n", h->Layer());
            }
        }
    }
}

void SaveToFile(Hists &hists, TString outname)
{
    TFile *fout = new TFile(outname, "RECREATE");
    fout->cd();
    hists.hM_PVz->Write();
    hists.hM_ClusX_all->Write();
    hists.hM_ClusY_all->Write();
    hists.hM_ClusZ_all->Write();
    hists.hM_ClusEta_all->Write();
    hists.hM_ClusPhi_all->Write();
    hists.hM_ClusX_ClusY_all->Write();
    hists.hM_ClusZ_ClusPhi_all->Write();
    hists.hM_ClusZ_TruthPVz_all->Write();
    hists.hM_ClusPhi_TruthPVz_all->Write();
    hists.hM_ClusX_layer1->Write();
    hists.hM_ClusY_layer1->Write();
    hists.hM_ClusZ_layer1->Write();
    hists.hM_ClusEta_layer1->Write();
    hists.hM_ClusPhi_layer1->Write();
    hists.hM_ClusZ_ClusPhi_layer1->Write();
    hists.hM_ClusZ_TruthPVz_layer1->Write();
    hists.hM_ClusPhi_TruthPVz_layer1->Write();
    hists.hM_ClusX_layer2->Write();
    hists.hM_ClusY_layer2->Write();
    hists.hM_ClusZ_layer2->Write();
    hists.hM_ClusEta_layer2->Write();
    hists.hM_ClusPhi_layer2->Write();
    hists.hM_ClusZ_ClusPhi_layer2->Write();
    hists.hM_ClusZ_TruthPVz_layer2->Write();
    hists.hM_ClusPhi_TruthPVz_layer2->Write();
    hists.hM_ClusX_layer3->Write();
    hists.hM_ClusY_layer3->Write();
    hists.hM_ClusZ_layer3->Write();
    hists.hM_ClusEta_layer3->Write();
    hists.hM_ClusPhi_layer3->Write();
    hists.hM_ClusZ_ClusPhi_layer3->Write();
    hists.hM_ClusZ_TruthPVz_layer3->Write();
    hists.hM_ClusPhi_TruthPVz_layer3->Write();

    hists.hM_dPhi_dEta_woCuts->Write();
    hists.hM_Eta_hit1_proto->Write();
    hists.hM_Phi_hit1_proto->Write();
    hists.hM_Eta_hit2_proto->Write();
    hists.hM_Phi_hit2_proto->Write();
    hists.hM_dEta_proto->Write();
    hists.hM_dEta_proto_altrange->Write();
    hists.hM_dEta_proto_altrange2->Write();
    hists.hM_dPhi_proto->Write();
    hists.hM_dPhi_proto_altrange->Write();
    hists.hM_Eta_proto->Write();
    hists.hM_Eta_proto_PVzRange1->Write();
    hists.hM_Eta_proto_PVzRange2->Write();
    hists.hM_Eta_proto_PVzRange3->Write();
    hists.hM_Eta_proto_PVzRange4->Write();
    hists.hM_Eta_proto_PVzRange5->Write();
    hists.hM_Phi_proto->Write();
    hists.hM_hit1Eta_hit2Eta_proto->Write();
    hists.hM_hit1Phi_hit2Phi_proto->Write();
    hists.hM_dR_proto->Write();
    hists.hM_dR_proto_altrange->Write();
    hists.hM_dR_proto_LogX->Write();
    hists.hM_Eta_vtxZ_proto->Write();
    hists.hM_Eta_vtxZ_anabin_proto->Write();
    hists.hM_dPhi_dEta_proto->Write();
    hists.hM_dPhi_dEta_proto_altrange->Write();

    hists.hM_Eta_hit1_reco->Write();
    hists.hM_Phi_hit1_reco->Write();
    hists.hM_Eta_hit2_reco->Write();
    hists.hM_Phi_hit2_reco->Write();
    hists.hM_dEta_reco->Write();
    hists.hM_dEta_reco_altrange->Write();
    hists.hM_dEta_reco_altrange2->Write();
    hists.hM_dPhi_reco->Write();
    hists.hM_dPhi_reco_altrange->Write();
    hists.hM_Eta_reco->Write();
    hists.hM_Phi_reco->Write();
    hists.hM_hit1Eta_hit2Eta_reco->Write();
    hists.hM_hit1Phi_hit2Phi_reco->Write();
    hists.hM_dR_reco->Write();
    hists.hM_dR_reco_altrange->Write();
    hists.hM_dR_reco_LogX->Write();
    hists.hM_Eta_vtxZ_reco->Write();
    hists.hM_Eta_vtxZ_anabin_reco->Write();
    hists.hM_dPhi_dEta_reco->Write();
    hists.hM_dPhi_dEta_reco_altrange->Write();

    hists.hM_Eta_vtxZ_reco_genmatched->Write();
    hists.hM_Eta_vtxZ_anabin_reco_genmatched->Write();

    fout->Close();
}

#endif