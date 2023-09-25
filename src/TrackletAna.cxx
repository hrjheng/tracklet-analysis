#include <TFile.h>
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
#include "pdgidfunc.h"

int main(int argc, char *argv[])
{
    // Usage: ./TrackletAna [NevtToRun] [skip] [layer] [randhit_case] [clusplit_case] [misalignment_num] [dRCut]
    // Example: ./TrackletAna 2000 0 12 0 0 0 0.5
    int NevtToRun_ = TString(argv[1]).Atoi();
    int skip = TString(argv[2]).Atoi();
    int layer_ = TString(argv[3]).Atoi();
    int randhit_case = TString(argv[4]).Atoi();     // 0: nominal, 1: 1%, 2: 5%, 3: 10%
    int clusplit_case = TString(argv[5]).Atoi();    // 0: nominal, 1: 1%, 2: 2%, 3: 5%
    int misalignment_num = TString(argv[6]).Atoi(); // 0: nominal (no mis-alignment)
    float dRCut = TString(argv[7]).Atof();          // Nominal: 0.5, variation \pm 0.1

    int l1 = ((int)layer_ / 10) - 1;
    int l2 = (layer_ % 10) - 1;
    vector<float> randhit_prec = {0., 1., 5., 10.};
    vector<float> clusplit_prec = {0., 1., 2., 5.};
    int iniEvt = skip;
    // mis-alignment study
    int misalign_case, misalign_layer;
    float shift_z, shift_r, shift_angle, twohalf_angle;
    vector<float> misalign_param;
    misalign_param.clear();
    if (misalignment_num < 40)
    {
        misalign_case = ((int)misalignment_num / 10);
        misalign_layer = (misalignment_num % 10) - 1;
        misalign_param = misalignment(misalign_case);
        shift_z = misalign_param[0];
        shift_r = misalign_param[1];
        shift_angle = misalign_param[2];
        twohalf_angle = 0;
    }
    else
    {
        misalign_case = 0;
        misalign_layer = 0;
        misalign_param = misalignment(misalignment_num);
        shift_z = 0;
        shift_r = 0;
        shift_angle = 0;
        twohalf_angle = misalign_param[0];
    }

    cout << "[Run Info] NevtToRun: " << NevtToRun_ << ", skip: " << skip << ", layer: " << layer_ << ", randhit_case: " << randhit_case << ", clusplit_case: " << clusplit_case << ", misalignment_num: " << misalignment_num << ", dRCut: " << dRCut << endl;

    // Optimized cut values for vertex finding algorithm: dPhi = 0.03 & |dZ| = 0.10
    // Feild on
    // TString EvtVtx_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth2p9_3DVertex_twohalves.root";
    // TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters.root";
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer%d_Evt%dto%d_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_twohalfangle_GapNorth2p9.root", layer_, iniEvt, iniEvt + NevtToRun_, randhit_case, clusplit_case, misalignment_num, number_to_string(dRCut).c_str());
    // Feild off
    TString EvtVtx_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/HijingAuAuMB_NoPileup_0T_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth2p9_GapUpper1p45_CentShift0p0_3DVertex_twohalves.root";
    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/HijingMBwoPUB0_private/MockSim/MVTXRecoClusters_HijingMBwoPU0T_MockSim.root";
    TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/HijingAuAuMB_NoPileup_0T_TrackletMinitree/TrackletAna_minitree_layer%d_Evt%dto%d_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s_GapNorth2p9_GapUpper1p45_CentShift0p0.root", layer_, iniEvt, iniEvt + NevtToRun_, randhit_case, clusplit_case, misalignment_num,
                               number_to_string(dRCut).c_str());
    // Simple pions
    // TString EvtVtx_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/Pion_RecoVtx_Optimization/SimplePion_Npart500_VtxZmean0cm_VtxZwidth0cm_pT0p001to5GeV/TrackletAna_RecoClusters_SimplePion_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root";
    // TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/Pion/Ntuple_Pions_Npart500_VtxZmean0cm_VtxZwidth0cm_pT0p001to5GeV.root";
    // TString outfilename = Form("/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer%d_Evt%dto%d_RandhitCase%d_ClusSplitCase%d_MisAlignNum%d_dRcut%s.root", layer_, iniEvt, iniEvt + NevtToRun_, randhit_case, clusplit_case, misalignment_num, number_to_string(dRCut).c_str());

    // Centrality bin (temporary, will move to use the official centrality variable once it's a thing)
    std::fstream is("/sphenix/user/hjheng/TrackletAna/analysis/plot/Centrality_bin.txt", std::ios_base::in);
    vector<float> centbin;
    float number;
    while (is >> number)
    {
        // printf("%f ", number);
        centbin.push_back(number);
    }
    printf("\n");

    TrackletData tkldata = {};
    tkldata.MVTXl1 = l1;
    tkldata.MVTXl2 = l2;

    std::map<int, float> EvtVtx_map = EvtVtx_map_tklcluster(EvtVtx_map_filename.Data());

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("EventTree");
    t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    int event, NTruthVtx;
    float TruthPV_trig_z, TruthPV_mostNpart_z;
    vector<int> *ClusLayer = 0;
    vector<float> *ClusX = 0, *ClusY = 0, *ClusZ = 0;
    vector<int> *UniqueAncG4P_PID = 0;
    vector<float> *UniqueAncG4P_Pt = 0, *UniqueAncG4P_Eta = 0, *UniqueAncG4P_Phi = 0, *UniqueAncG4P_E = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("NTruthVtx", &NTruthVtx);
    t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
    t->SetBranchAddress("TruthPV_mostNpart_z", &TruthPV_mostNpart_z);
    t->SetBranchAddress("ClusLayer", &ClusLayer);
    t->SetBranchAddress("ClusX", &ClusX);
    t->SetBranchAddress("ClusY", &ClusY);
    t->SetBranchAddress("ClusZ", &ClusZ);
    t->SetBranchAddress("UniqueAncG4P_PID", &UniqueAncG4P_PID);
    t->SetBranchAddress("UniqueAncG4P_Pt", &UniqueAncG4P_Pt);
    t->SetBranchAddress("UniqueAncG4P_Eta", &UniqueAncG4P_Eta);
    t->SetBranchAddress("UniqueAncG4P_Phi", &UniqueAncG4P_Phi);
    t->SetBranchAddress("UniqueAncG4P_E", &UniqueAncG4P_E);

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
        tkldata.event = event;
        tkldata.PV_z = PV_z;
        tkldata.TruthPV_trig_z = TruthPV_trig_z;
        tkldata.TruthPV_mostNpart_z = TruthPV_mostNpart_z;

        cout << "[Event info] Event = " << event << endl;

        // A selection to select events -> mimic 0-pileup environment
        tkldata.pu0_sel = (NTruthVtx == 1) ? true : false;
        tkldata.trig = true;
        tkldata.process = 0; // Exclude single diffractive events

        // Event selection
        // tkldata.vz_sel = (fabs(PV_z) <= 10.) ? true : false;
        // tkldata.o_sel = true;    // CMS XeXe: HLT && (nhfn > 2) && (nhfp > 2)
        // tkldata.proc_sel = true;
        // tkldata.evt_sel = tkldata.vz_sel && tkldata.o_sel;
        // tkldata.gen_sel = tkldata.vz_sel && tkldata.proc_sel;

        // Prepare clusters
        for (size_t ihit = 0; ihit < ClusLayer->size(); ihit++)
        {
            // Mis-alignment study
            if (misalignment_num < 40)
            {
                float clusposx = (ClusLayer->at(ihit) == misalign_layer) ? ClusX->at(ihit) + shift_r * TMath::Cos(shift_angle) : ClusX->at(ihit);
                float clusposy = (ClusLayer->at(ihit) == misalign_layer) ? ClusY->at(ihit) + shift_r * TMath::Sin(shift_angle) : ClusY->at(ihit);
                float clusposz = (ClusLayer->at(ihit) == misalign_layer) ? ClusZ->at(ihit) + shift_z : ClusZ->at(ihit);
                Hit *hit = new Hit(clusposx, clusposy, clusposz, 0., 0., PV_z, ClusLayer->at(ihit));
                tkldata.MVTXlayers[ClusLayer->at(ihit)].push_back(hit);
            }
            else
            {
                // Misalignment study for MVTX
                float x0 = ClusX->at(ihit);
                float y0 = ClusY->at(ihit);
                float z0 = ClusZ->at(ihit);
                vector<float> cpos = {x0, y0, z0};
                UpdatePos_GapTwoHalves(cpos, 2.9, 0, 1.45);
                Hit *hit = new Hit(cpos[0], cpos[1], cpos[2], 0., 0., PV_z, ClusLayer->at(ihit));
                tkldata.MVTXlayers[ClusLayer->at(ihit)].push_back(hit);
            }
        }

        // Centrality
        auto bin = std::lower_bound(centbin.begin(), centbin.end(), tkldata.MVTXlayers[0].size());
        tkldata.CentralityBin = (tkldata.MVTXlayers[0].size() == 0) ? 19 : 20 - (bin - centbin.begin());
        cout << "Truth PV z = " << TruthPV_trig_z << "; PV_z = " << PV_z << "; (reco PV z - truth PV z) = " << (PV_z - TruthPV_trig_z) << endl
             << "# of clusters on MVTX 1st layer = " << tkldata.MVTXlayers[0].size() << " -> Centrality = " << tkldata.CentralityBin * 5 << "-" << (tkldata.CentralityBin + 1) * 5 << "% (Centrality bin = " << tkldata.CentralityBin << ")" << endl;
        tkldata.NhitsLayer1 = tkldata.MVTXlayers[0].size();

        // Random hits: for background study
        if (randhit_case != 0)
        {
            int Nrandhit_l1 = (int)randhit_prec[randhit_case] * 0.01 * tkldata.MVTXlayers[tkldata.MVTXl1].size();
            int Nrandhit_l2 = (int)randhit_prec[randhit_case] * 0.01 * tkldata.MVTXlayers[tkldata.MVTXl2].size();

            for (int i = 0; i < Nrandhit_l1; i++)
            {
                tkldata.MVTXlayers[tkldata.MVTXl1].push_back(RandomHit(-3.0, 3.0, -TMath::Pi(), TMath::Pi()));
            }

            for (int i = 0; i < Nrandhit_l2; i++)
            {
                tkldata.MVTXlayers[tkldata.MVTXl2].push_back(RandomHit(-3.0, 3.0, -TMath::Pi(), TMath::Pi()));
            }
        }

        // Clusters spliting
        if (clusplit_case != 0)
        {
            float split = clusplit_prec[clusplit_case] * 0.01;

            for (int i = 0; i < tkldata.MVTXlayers[tkldata.MVTXl1].size(); i++)
            {
                if (gRandom->Rndm() < split)
                {
                    Hit *splithit = new Hit(tkldata.MVTXlayers[tkldata.MVTXl1][i]->Eta() + gRandom->Gaus(0, 0.006), tkldata.MVTXlayers[tkldata.MVTXl1][i]->Phi() + gRandom->Gaus(0, 0.006));
                    tkldata.MVTXlayers[tkldata.MVTXl1].push_back(splithit);
                }
            }

            for (int i = 0; i < tkldata.MVTXlayers[tkldata.MVTXl2].size(); i++)
            {
                if (gRandom->Rndm() < split)
                {
                    Hit *splithit = new Hit(tkldata.MVTXlayers[tkldata.MVTXl2][i]->Eta() + gRandom->Gaus(0, 0.006), tkldata.MVTXlayers[tkldata.MVTXl2][i]->Phi() + gRandom->Gaus(0, 0.006));
                    tkldata.MVTXlayers[tkldata.MVTXl2].push_back(splithit);
                }
            }
        }

        // Drop clusters, for cluster reconstruction inefficiency

        RecoclusToMinitree(tkldata);

        // Tracklet reconstruction: proto-tracklets -> reco-tracklets -> gen-hadron matching
        ProtoTracklets(tkldata, dRCut);
        RecoTracklets(tkldata);

        // Generated charged hadrons
        for (size_t ihad = 0; ihad < UniqueAncG4P_PID->size(); ihad++)
        {
            if (is_chargedHadron(UniqueAncG4P_PID->at(ihad)) == false)
                continue;

            GenHadron *genhadron = new GenHadron(UniqueAncG4P_Pt->at(ihad), UniqueAncG4P_Eta->at(ihad), UniqueAncG4P_Phi->at(ihad), UniqueAncG4P_E->at(ihad));
            tkldata.GenHadrons.push_back(genhadron);
            tkldata.GenHadron_Pt.push_back(UniqueAncG4P_Pt->at(ihad));
            tkldata.GenHadron_eta.push_back(UniqueAncG4P_Eta->at(ihad));
            tkldata.GenHadron_phi.push_back(UniqueAncG4P_Phi->at(ihad));
            tkldata.GenHadron_E.push_back(UniqueAncG4P_E->at(ihad));
        }
        tkldata.NGenHadron = tkldata.GenHadrons.size();

        GenMatch_Recotkl(tkldata);
        tkldata.NRecotkl_GenMatched = tkldata.RecoTkls_GenMatched.size();

        minitree->Fill();
        ResetVec(tkldata);
        cout << "----------" << endl;
    }

    outfile->cd();
    minitree->Write("", TObject::kOverwrite);
    outfile->Close();

    f->Close();

    return 0;
}
