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
#include <TMath.h>

#include "Tracklet.h"
#include "GenHadron.h"
#include "TrackletData.h"
#include "VtxData.h"
#include "misalignment.h"

int main(int argc, char *argv[])
{
    int NevtToRun_ = TString(argv[1]).Atoi();
    int skip = TString(argv[2]).Atoi();
    int layer_ = TString(argv[3]).Atoi();
    int l1 = ((int)layer_ / 10) - 1;
    int l2 = (layer_ % 10) - 1;
    int randhit_case = TString(argv[4]).Atoi();     // 0: nominal, 1: 1%, 2: 5%, 3: 10%
    int misalignment_num = TString(argv[5]).Atoi(); // 0: nominal (no mis-alignment)
    vector<float> randhit_prec = {0., 1., 5., 10.};
    // fprintf(stderr, "Random hit case %d: + %.2f %% random hits \n", randhit_case, randhit_prec[randhit_case]);
    int iniEvt = skip;
    // Mis-alignment study
    int misalign_case = ((int)misalignment_num / 10);
    int misalign_layer = (misalignment_num % 10) - 1;
    vector<float> misalign_param = misalignment(misalign_case);
    float shift_z = misalign_param[0];
    float shift_r = misalign_param[1];
    float shift_angle = misalign_param[2];

    TString EvtVtx_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_RecoClusters_RecoVtx_TklCluster_Nevt4000.root";
    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/MVTXRecoClusters_Nevt4000.root";
    TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal/TrackletAna_minitree_layer%d_Evt%dto%d_RandhitCase%d_MisAlignNum%d.root", layer_, iniEvt, iniEvt + NevtToRun_, randhit_case, misalignment_num);
    TString histfilename = Form("./plot/hists/AuAu_Nominal/TrackletAna_histograms_layer%d_Evt%dto%d_RandhitCase%d_MisAlignNum%d.root", layer_, iniEvt, iniEvt + NevtToRun_, randhit_case, misalignment_num);

    // For simple pions event
    // int vtxzmean = -3;
    // int vtxzwidth = 8;
    // TString EvtVtx_map_filename = Form("/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_Simple1000Pions_VtxZmean%dcm_VtxZwidth%dcm_RecoVtx_TklCluster.root", vtxzmean, vtxzwidth);
    // TString infilename = Form("/sphenix/user/hjheng/TrackletAna/data/myMVTXHits_SIMPLE1000Pions_VtxZmean%dcm_VtxZwidth%dcm.root", vtxzmean, vtxzwidth);
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/TrackletAna_Simple1000Pions_VtxZmean%dcm_VtxZwidth%dcm_minitree_layer%d_Evt%dto%d.root", vtxzmean, vtxzwidth, layer_, iniEvt, iniEvt + NevtToRun_);
    // TString histfilename = Form("./plot/hists/TrackletAna_histograms_Simple1000Pions_VtxZmean%dcm_VtxZwidth%dcm_layer%d_Evt%dto%d.root", vtxzmean, vtxzwidth, layer_, iniEvt, iniEvt + NevtToRun_);

    vector<vector<Hit *>> MVTXlayers = {{}, {}, {}};
    vector<Tracklet *> ProtoTkls, RecoTkls, RecoTkls_GenMatched;
    vector<GenHadron *> GenHadrons;
    Hists hists = {};
    SetupHists(hists);
    TrackletData tkldata = {};

    std::map<int, float> EvtVtx_map = EvtVtx_map_tklcluster(EvtVtx_map_filename.Data());

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("EventTree");
    t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    int event, NTruthVtx;
    float TruthPV_trig_z, TruthPV_mostNpart_z;
    vector<int> *ClusLayer = 0;
    vector<float> *ClusX = 0, *ClusY = 0, *ClusZ = 0;
    vector<int> *G4PartPID = 0, *G4PartPrimaryID = 0, *G4PartParentID = 0;
    vector<float> *G4PartPt = 0, *G4PartEta = 0, *G4PartPhi = 0, *G4PartE = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("NTruthVtx", &NTruthVtx);
    t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
    t->SetBranchAddress("TruthPV_mostNpart_z", &TruthPV_mostNpart_z);
    t->SetBranchAddress("ClusLayer", &ClusLayer);
    t->SetBranchAddress("ClusX", &ClusX);
    t->SetBranchAddress("ClusY", &ClusY);
    t->SetBranchAddress("ClusZ", &ClusZ);
    t->SetBranchAddress("G4PartPID", &G4PartPID);
    t->SetBranchAddress("G4PartPrimaryID", &G4PartPrimaryID);
    t->SetBranchAddress("G4PartParentID", &G4PartParentID);
    t->SetBranchAddress("G4PartPt", &G4PartPt);
    t->SetBranchAddress("G4PartEta", &G4PartEta);
    t->SetBranchAddress("G4PartPhi", &G4PartPhi);
    t->SetBranchAddress("G4PartE", &G4PartE);

    TFile *outfile = new TFile(outfilename, "RECREATE");
    TTree *minitree = new TTree(Form("minitree_%d", layer_), "Minitree of Reconstructed Tracklets");
    SetMinitree(minitree, tkldata);

    for (int i = skip; i < skip + NevtToRun_; i++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);

        // cout << "event = " << event << " local = " << local << endl;
        float PV_z = EvtVtx_map[event];
        vector<float> PV = {0, 0, PV_z};
        cout << "event = " << event << " NTruthVtx = " << NTruthVtx << " PV_z = " << PV_z << endl;

        if (NTruthVtx != 1)
            continue;

        if (fabs(PV_z) > 10.)
            continue;

        for (size_t ihit = 0; ihit < ClusLayer->size(); ihit++)
        {
            // Mis-alignment study
            float clusposx = (ClusLayer->at(ihit) == misalign_layer) ? ClusX->at(ihit) + shift_r * TMath::Cos(shift_angle) : ClusX->at(ihit);
            float clusposy = (ClusLayer->at(ihit) == misalign_layer) ? ClusY->at(ihit) + shift_r * TMath::Sin(shift_angle) : ClusY->at(ihit);
            float clusposz = (ClusLayer->at(ihit) == misalign_layer) ? ClusZ->at(ihit) + shift_z : ClusZ->at(ihit);
            Hit *hit = new Hit(clusposx, clusposy, clusposz, 0., 0., 0., ClusLayer->at(ihit));
            MVTXlayers[ClusLayer->at(ihit)].push_back(hit);
        }

        tkldata.NhitsLayer1 = MVTXlayers[0].size();
        FillHist_Cluster(MVTXlayers, TruthPV_trig_z, hists);

        // fprintf(stderr, "event=%d has a total of %d clusters.\n", event, ClusLayer->size());
        // fprintf(stderr, "# of clusters in Layer %d = %d, Layer %d = %d \n", l1, MVTXlayers[l1].size(), l2, MVTXlayers[l2].size());

        // for (size_t i = 0; i < 5; i++)
        // {
        //     cout << "Before recalculating: hit" << i << " " << MVTXlayers[l1][i]->Eta() << " " << MVTXlayers[l1][i]->Phi() << " " << endl;
        // }
        UpdateHits(MVTXlayers[l1], PV);
        UpdateHits(MVTXlayers[l2], PV);
        // for (size_t i = 0; i < 5; i++)
        // {
        //     cout << "After recalculating: hit" << i << " " << MVTXlayers[l1][i]->Eta() << " " << MVTXlayers[l1][i]->Phi() << " " << endl;
        // }

        // Random hits
        if (randhit_case != 0)
        {
            int Nrandhit_l1 = (int)randhit_prec[randhit_case] * 0.01 * MVTXlayers[l1].size();
            int Nrandhit_l2 = (int)randhit_prec[randhit_case] * 0.01 * MVTXlayers[l2].size();

            for (int i = 0; i < Nrandhit_l1; i++)
            {
                MVTXlayers[l1].push_back(RandomHit(-3.0, 3.0, -TMath::Pi(), TMath::Pi()));
            }

            for (int i = 0; i < Nrandhit_l2; i++)
            {
                MVTXlayers[l2].push_back(RandomHit(-3.0, 3.0, -TMath::Pi(), TMath::Pi()));
            }
        }

        ProtoTkls = ProtoTracklets(MVTXlayers[l1], MVTXlayers[l2], hists);
        RecoTkls = RecoTracklets(ProtoTkls, hists);

        // Generated hadrons
        for (size_t ihad = 0; ihad < G4PartPID->size(); ihad++)
        {
            if (abs(G4PartPID->at(ihad)) < 9 ||
                (abs(G4PartPID->at(ihad)) > 10 && abs(G4PartPID->at(ihad)) < 19) ||
                (abs(G4PartPID->at(ihad)) > 20 && abs(G4PartPID->at(ihad)) < 26) ||
                (abs(G4PartPID->at(ihad)) > 31 && abs(G4PartPID->at(ihad)) < 38))
                continue;

            if (G4PartPt->at(ihad) < 0.03)
                continue;
            if (fabs(G4PartEta->at(ihad)) > 2.5)
                continue;

            GenHadron *genhadron = new GenHadron(G4PartEta->at(ihad), G4PartPhi->at(ihad));
            GenHadrons.push_back(genhadron);
        }

        RecoTkls_GenMatched = GenMatch_Recotkl(RecoTkls, GenHadrons, hists);

        hists.hM_PVz->Fill(PV_z);
        hists.hM_NhitsLayer1->Fill(MVTXlayers[0].size());

        tkldata.event = event;
        tkldata.PV_z = PV_z;
        tkldata.TruthPV_trig_z = TruthPV_trig_z;
        tkldata.TruthPV_mostNpart_z = TruthPV_mostNpart_z;
        tkldata.NRecotkl_Total = RecoTkls.size();
        for (auto &tkl : RecoTkls)
        {
            tkldata.recotkl_eta1.push_back(tkl->Hit1()->Eta());
            tkldata.recotkl_phi1.push_back(tkl->Hit1()->Phi());
            tkldata.recotkl_eta2.push_back(tkl->Hit2()->Eta());
            tkldata.recotkl_phi2.push_back(tkl->Hit2()->Phi());
            tkldata.recotkl_deta.push_back(tkl->dEta());
            tkldata.recotkl_dphi.push_back(tkl->dPhi());
            tkldata.recotkl_dR.push_back(tkl->dR());
        }

        tkldata.NRecotkl_GenMatched = RecoTkls_GenMatched.size();
        for (auto &tkl : RecoTkls_GenMatched)
        {
            tkldata.recotkl_genmatched_eta1.push_back(tkl->Hit1()->Eta());
            tkldata.recotkl_genmatched_phi1.push_back(tkl->Hit1()->Phi());
            tkldata.recotkl_genmatched_eta2.push_back(tkl->Hit2()->Eta());
            tkldata.recotkl_genmatched_phi2.push_back(tkl->Hit2()->Phi());
            tkldata.recotkl_genmatched_deta.push_back(tkl->dEta());
            tkldata.recotkl_genmatched_dphi.push_back(tkl->dPhi());
            tkldata.recotkl_genmatched_dR.push_back(tkl->dR());
        }

        minitree->Fill();

        /*
        int irecotkl = 0;
        for (auto &recotkl : v_recotkl_L12)
        {
            printf("ith recotkl = %d,\t dR = %f,\t hit1 eta = %f,\t hit1 phi = %f,\t IsMatchedTkl = %s,\t hit2 eta = %f,\t hit2 phi = %f,\t IsMatchedTkl = %s,\t # of matched hits @ L1 = %d\n", irecotkl, recotkl->dR(), recotkl->Hit1()->Eta(), recotkl->Hit1()->Phi(), recotkl->Hit1()->IsMatchedTkl() ? "true" : "false", recotkl->Hit2()->Eta(), recotkl->Hit2()->Phi(), recotkl->Hit2()->IsMatchedTkl() ? "true" : "false", recotkl->Hit2()->MatchedL1Hits().size());
            irecotkl++;
        }
        */

        for (size_t i = 0; i < MVTXlayers.size(); i++)
        {
            CleanVec(MVTXlayers[i]);
        }
        CleanVec(ProtoTkls);
        CleanVec(RecoTkls);
        CleanVec(RecoTkls_GenMatched);
        CleanVec(GenHadrons);
        CleanVec(tkldata.recotkl_eta1);
        CleanVec(tkldata.recotkl_phi1);
        CleanVec(tkldata.recotkl_eta2);
        CleanVec(tkldata.recotkl_phi2);
        CleanVec(tkldata.recotkl_deta);
        CleanVec(tkldata.recotkl_dphi);
        CleanVec(tkldata.recotkl_dR);
        CleanVec(tkldata.recotkl_genmatched_eta1);
        CleanVec(tkldata.recotkl_genmatched_phi1);
        CleanVec(tkldata.recotkl_genmatched_eta2);
        CleanVec(tkldata.recotkl_genmatched_phi2);
        CleanVec(tkldata.recotkl_genmatched_deta);
        CleanVec(tkldata.recotkl_genmatched_dphi);
        CleanVec(tkldata.recotkl_genmatched_dR);
    }

    SaveToFile(hists, histfilename);

    outfile->cd();
    minitree->Write();
    outfile->Close();

    f->Close();

    return 0;
}
