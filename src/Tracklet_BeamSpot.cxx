#include <TFile.h>
#include <TMath.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TTreeIndex.h>
#include <TStyle.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Tracklet.h"
#include "Utilities.h"
#include "Vertex.h"

int main(int argc, char *argv[])
{
    TString EvtVtx_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root";
    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters.root";

    std::map<int, float> EvtVtx_map = EvtVtx_map_tklcluster(EvtVtx_map_filename.Data());

    vector<vector<Hit *>> MVTXlayer = {{}, {}, {}};

    TH2F *hM_VtxX_TruthReco = new TH2F("hM_VtxX_TruthReco", "hM_VtxX_TruthReco", 100, -0.2, 0.2, 100, -0.2, 0.2);
    TH2F *hM_VtxY_TruthReco = new TH2F("hM_VtxY_TruthReco", "hM_VtxY_TruthReco", 100, -0.2, 0.2, 100, -0.2, 0.2);

    TFile *f = new TFile(infilename, "READ");
    TTree *t = (TTree *)f->Get("EventTree");
    t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
    TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
    int event, NTruthVtx;
    float TruthPV_trig_x, TruthPV_trig_y, TruthPV_trig_z;
    vector<int> *ClusLayer = 0;
    vector<float> *ClusX = 0, *ClusY = 0, *ClusZ = 0;
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("NTruthVtx", &NTruthVtx);
    t->SetBranchAddress("TruthPV_trig_x", &TruthPV_trig_x);
    t->SetBranchAddress("TruthPV_trig_y", &TruthPV_trig_y);
    t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
    t->SetBranchAddress("ClusLayer", &ClusLayer);
    t->SetBranchAddress("ClusX", &ClusX);
    t->SetBranchAddress("ClusY", &ClusY);
    t->SetBranchAddress("ClusZ", &ClusZ);

    for (int i = 0; i < index->GetN(); i++)
    {
        Long64_t local = t->LoadTree(index->GetIndex()[i]);
        t->GetEntry(local);

        float PV_z = EvtVtx_map[event];

        cout << "event = " << event << " PV_z = " << PV_z << " Truth trig vtx (x,y,z)=(" << TruthPV_trig_x << "," << TruthPV_trig_y << "," << TruthPV_trig_z
             << ") # of clusters = " << ClusLayer->size() << endl;

        for (size_t i = 0; i < MVTXlayer.size(); i++)
        {
            CleanVec(MVTXlayer[i]);
        }

        for (size_t ihit = 0; ihit < ClusLayer->size(); ihit++)
        {
            Hit *hit = new Hit(ClusX->at(ihit), ClusY->at(ihit), ClusZ->at(ihit), 0., 0., PV_z, ClusLayer->at(ihit));
            MVTXlayer[ClusLayer->at(ihit)].push_back(hit);
        }

        vector<float> vPVx, vPVy;
        CleanVec(vPVx);
        CleanVec(vPVy);

        for (auto hitl1 : MVTXlayer[0])
        {
            for (auto hitl2 : MVTXlayer[1])
            {
                float dEta = hitl2->Eta() - hitl1->Eta();
                float dPhi = deltaPhi(hitl2->Phi(), hitl1->Phi());
                float dR = sqrt(dEta * dEta + dPhi * dPhi);
                if (fabs(deltaPhi(hitl1->Phi(), hitl2->Phi())) > 0.03)
                    continue;
                // float rhoPV = hitl1->rho() - (PV_z - hitl1->posZ()) / (hitl1->posZ() - hitl2->posZ()) * (hitl2->rho() - hitl1->rho());
                bool validcomb = ((PV_z > hitl1->posZ()) && (hitl1->posZ() > hitl2->posZ())) || ((PV_z < hitl1->posZ()) && (hitl1->posZ() < hitl2->posZ()));
                if (!validcomb)
                    continue;

                float extrapolated_x = (PV_z - hitl1->posZ()) / (hitl1->posZ() - hitl2->posZ()) * (hitl1->posX() - hitl2->posX()) + hitl1->posX();
                float extrapolated_y = (PV_z - hitl1->posZ()) / (hitl1->posZ() - hitl2->posZ()) * (hitl1->posY() - hitl2->posY()) + hitl1->posY();
                float extrapolated_rho = sqrt(extrapolated_x * extrapolated_x + extrapolated_y * extrapolated_y);

                if (extrapolated_rho > 2.461)
                    continue;

                vPVx.push_back(extrapolated_x);
                vPVy.push_back(extrapolated_y);
                // if (MVTXlayer[0].size() < 20)
                // {
                //     // cout << "PV_z=" << PV_z << " hit1 (x,y,z)=(" << hitl1->posX() << "," << hitl1->posY() << "," << hitl1->posZ() << ") hit2 (x,y,z)=(" << hitl2->posX() << "," << hitl2->posY() << "," << hitl2->posZ() << ") extrapolated_x=" << extrapolated_x << " extrapolated_y=" << extrapolated_y << endl;
                // }
            }
        }

        float sum_vPVx, avg_vPVx, sum_vPVy, avg_vPVy, median_vPVx, median_vPVy, sq_sum_vPVx, stdev_vPVx, sq_sum_vPVy, stdev_vPVy;
        if (ClusLayer->size() == 0 || vPVx.size() == 0 || vPVy.size() == 0)
        {
            sum_vPVx = 0;
            avg_vPVx = -999.;
            sum_vPVy = 0;
            avg_vPVy = -999.;
            median_vPVx = 0;
            median_vPVy = 0;
            sq_sum_vPVx = 0;
            stdev_vPVx = 0;
            sq_sum_vPVy = 0;
            stdev_vPVy = 0;
        }
        else
        {
            // Calculate the average
            sum_vPVx = accumulate(vPVx.begin(), vPVx.end(), 0.0);
            avg_vPVx = sum_vPVx / vPVx.size();
            sum_vPVy = accumulate(vPVy.begin(), vPVy.end(), 0.0);
            avg_vPVy = sum_vPVy / vPVy.size();

            // Calculate the median
            sort(vPVx.begin(), vPVx.end());
            median_vPVx = (vPVx.size() % 2 == 0) ? (vPVx[vPVx.size() / 2 - 1] + vPVx[vPVx.size() / 2]) / 2 : vPVx[vPVx.size() / 2];
            sort(vPVy.begin(), vPVy.end());
            median_vPVy = (vPVy.size() % 2 == 0) ? (vPVy[vPVy.size() / 2 - 1] + vPVy[vPVy.size() / 2]) / 2 : vPVy[vPVy.size() / 2];

            // Calculate the standard deviation
            sq_sum_vPVx = inner_product(vPVx.begin(), vPVx.end(), vPVx.begin(), 0.0);
            stdev_vPVx = sqrt(sq_sum_vPVx / vPVx.size() - avg_vPVx * avg_vPVx);
            sq_sum_vPVy = inner_product(vPVy.begin(), vPVy.end(), vPVy.begin(), 0.0);
            stdev_vPVy = sqrt(sq_sum_vPVy / vPVy.size() - avg_vPVy * avg_vPVy);
        }

        // hM_VtxX_TruthReco->Fill(TruthPV_trig_x, avg_vPVx);
        // hM_VtxY_TruthReco->Fill(TruthPV_trig_y, avg_vPVy);
        hM_VtxX_TruthReco->Fill(TruthPV_trig_x, median_vPVx);
        hM_VtxY_TruthReco->Fill(TruthPV_trig_y, median_vPVy);

        // cout << "Truth trig vtx (x,y)=(" << TruthPV_trig_x << "," << TruthPV_trig_y << "); Reco PV (x,y)=(" << avg_vPVx << "," << avg_vPVy << ")" << endl;
    }

    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    gPad->SetLeftMargin(0.13);
    hM_VtxX_TruthReco->SetContour(1000);
    hM_VtxX_TruthReco->GetXaxis()->SetTitle("Truth vtx x (cm)");
    hM_VtxX_TruthReco->GetYaxis()->SetTitle("Reco PV x (cm)");
    hM_VtxX_TruthReco->Draw("colz");
    c1->SaveAs("./plot/PV_TruthReco/hM_VtxX_TruthReco_mean.pdf");
    c1->cd();
    hM_VtxY_TruthReco->SetContour(1000);
    hM_VtxY_TruthReco->GetXaxis()->SetTitle("Truth vtx y (cm)");
    hM_VtxY_TruthReco->GetYaxis()->SetTitle("Reco PV y (cm)");
    hM_VtxY_TruthReco->Draw("colz");
    c1->SaveAs("./plot/PV_TruthReco/hM_VtxY_TruthReco_mean.pdf");
}