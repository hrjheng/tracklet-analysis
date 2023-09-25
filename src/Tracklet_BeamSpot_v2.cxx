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

#include "Tracklet.h"
#include "Utilities.h"
#include "Vertex.h"

struct Point
{
    double x;
    double y;
};

void findCircle(Point A, Point B, Point C, double &h, double &k, double &r)
{
    // Find the midpoint of AB and BC
    Point Mab = {(A.x + B.x) / 2, (A.y + B.y) / 2};
    Point Mbc = {(B.x + C.x) / 2, (B.y + C.y) / 2};

    // Find the slope of AB and BC
    double slopeAB = (B.y - A.y) / (B.x - A.x);
    double slopeBC = (C.y - B.y) / (C.x - B.x);

    // Find the equation of the perpendicular bisectors of AB and BC
    double x1 = Mab.x;
    double y1 = Mab.y;
    double m1 = -1 / slopeAB;
    double b1 = y1 - m1 * x1;
    double x2 = Mbc.x;
    double y2 = Mbc.y;
    double m2 = -1 / slopeBC;
    double b2 = y2 - m2 * x2;

    // Solve for the intersection point of the two perpendicular bisectors
    h = (b2 - b1) / (m1 - m2);
    k = m1 * h + b1;

    // Find the radius of the circle
    r = sqrt(pow(h - A.x, 2) + pow(k - A.y, 2));
}

int main(int argc, char *argv[])
{
    TString EvtVtx_map_filename = "/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10.root";
    TString infilename = "/sphenix/user/hjheng/TrackletAna/data/MVTXRecoClusters/NoPileup_Nevt500_ana325private_singleEvtDst/MVTXRecoClusters.root";

    std::map<int, float> EvtVtx_map = EvtVtx_map_tklcluster(EvtVtx_map_filename.Data());

    float dPhiCut = 0.03;

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
    // for (int i = 0; i < 100; i++)
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
                // float dEta21 = hitl2->Eta() - hitl1->Eta();
                float dPhi21 = deltaPhi(hitl2->Phi(), hitl1->Phi());
                // float dR21 = sqrt(dEta21 * dEta21 + dPhi21 * dPhi21);
                if (fabs(dPhi21) > dPhiCut)
                    continue;

                for (auto hitl3 : MVTXlayer[2])
                {
                    // float dEta = hitl3->Eta() - hitl1->Eta();
                    float dPhi32 = deltaPhi(hitl3->Phi(), hitl2->Phi());
                    // Search window for the hit on the 3rd layer
                    float dPhi32_max = dPhi21 * (1 + 0.2);
                    float dPhi32_min = dPhi21 * (1 - 0.2);
                    // float dR = sqrt(dEta * dEta + dPhi * dPhi);
                    // cout << "dPhi21 = " << dPhi21 << " dPhi32_min = " << dPhi32_min << " dPhi32_max = " << dPhi32_max << endl;
                    // cout << "event = " << event << " # of clusters=" << ClusLayer->size() << " hit1 phi=" << hitl1->Phi() << " hit2 phi=" << hitl2->Phi() << " hit3 phi=" << hitl3->Phi()
                    //      << " dPhi21 = " << dPhi21 << " dPhi32 = " << dPhi32 << " dPhi32_min = " << dPhi32_min << " dPhi32_max = " << dPhi32_max << endl;
                    if (dPhi32 > dPhi32_max || dPhi32 < dPhi32_min)
                        continue;
                    // if (fabs(dPhi32) > dPhiCut)
                    // continue;

                    // cout << "event = " << event << "# of clusters=" << ClusLayer->size() << " hit1 phi=" << hitl1->Phi() << " hit2 phi=" << hitl2->Phi() << " hit3 phi=" << hitl3->Phi() << endl;

                    double h, k, r;
                    Point hit1 = {hitl1->posX(), hitl1->posY()};
                    Point hit2 = {hitl2->posX(), hitl2->posY()};
                    Point hit3 = {hitl3->posX(), hitl3->posY()};
                    findCircle(hit1, hit2, hit3, h, k, r);
                    float phi1_circ = std::atan2(hitl1->posY() - k, hitl1->posX() - h);
                    float phi2_circ = std::atan2(hitl2->posY() - k, hitl2->posX() - h);
                    float phi3_circ = std::atan2(hitl3->posY() - k, hitl3->posX() - h);
                    float dPhi21_circ = phi2_circ - phi1_circ;
                    float dPhi32_circ = phi3_circ - phi2_circ;
                    float ratio_hitz321 = (hitl3->posZ() - hitl2->posZ()) / (hitl2->posZ() - hitl1->posZ());
                    float ratio_dphi321 = dPhi32_circ / dPhi21_circ;
                    float ratio_21PV = (hitl2->posZ() - hitl1->posZ()) / (hitl1->posZ() - PV_z);
                    float dPhi1PV_circ = dPhi21_circ * ((hitl1->posZ() - PV_z) / (hitl2->posZ() - hitl1->posZ()));
                    float phiPV_circ = phi1_circ - dPhi1PV_circ;
                    float PVx = h + r * cos(phiPV_circ);
                    float PVy = k + r * sin(phiPV_circ);
                    // cout << "event = " << event << " # of clusters=" << ClusLayer->size() << " hit1 phi=" << hitl1->Phi() << " hit2 phi=" << hitl2->Phi() << " hit3 phi=" << hitl3->Phi()
                    //      << " dPhi21=" << dPhi21 << " dPhi32=" << dPhi32 << " dPhi21_circ=" << dPhi21_circ << " dPhi32_circ=" << dPhi32_circ << " ratio_hitz321=" << ratio_hitz321
                    //      << " ratio_dphi321=" << ratio_dphi321 << " ratio_21PV=" << ratio_21PV << " dPhi1PV_circ=" << dPhi1PV_circ << " phiPV_circ=" << phiPV_circ << " PVx=" << PVx << " PVy=" <<
                    //      PVy
                    //      << endl;
                    // cout << "event = " << event << " # of clusters=" << ClusLayer->size() << " The equation of the circle is: (x - " << h << ")^2 + (y - " << k << ")^2 = " << r * r << endl;
                    vPVx.push_back(PVx);
                    vPVy.push_back(PVy);
                }
            }
        }

        // calculate the median of the PVx and PVy
        float median_vPVx, median_vPVy;
        if (ClusLayer->size() == 0 || vPVx.size() == 0 || vPVy.size() == 0)
        {
            median_vPVx = -999;
            median_vPVy = -999;
        }
        else
        {
            sort(vPVx.begin(), vPVx.end());
            median_vPVx = (vPVx.size() % 2 == 0) ? (vPVx[vPVx.size() / 2 - 1] + vPVx[vPVx.size() / 2]) / 2 : vPVx[vPVx.size() / 2];
            sort(vPVy.begin(), vPVy.end());
            median_vPVy = (vPVy.size() % 2 == 0) ? (vPVy[vPVy.size() / 2 - 1] + vPVy[vPVy.size() / 2]) / 2 : vPVy[vPVy.size() / 2];
        }

        hM_VtxX_TruthReco->Fill(TruthPV_trig_x, median_vPVx);
        hM_VtxY_TruthReco->Fill(TruthPV_trig_y, median_vPVy);

        // print the truth PV x/y and median of the PVx and PVy
        // cout << "event = " << event << " # of clusters=" << ClusLayer->size() << " Truth PV (x,y)=(" << TruthPV_trig_x << "," << TruthPV_trig_y << ") median PV (x,y)=(" << median_vPVx << ","
        //      << median_vPVy << ")" << endl;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    gPad->SetLeftMargin(0.13);
    hM_VtxX_TruthReco->SetContour(1000);
    hM_VtxX_TruthReco->GetXaxis()->SetTitle("Truth vtx x (cm)");
    hM_VtxX_TruthReco->GetYaxis()->SetTitle("Reco PV x (cm)");
    hM_VtxX_TruthReco->Draw("colz");
    c1->SaveAs("./plot/PV_TruthReco/hM_VtxX_TruthReco_v2.pdf");
    c1->cd();
    hM_VtxY_TruthReco->SetContour(1000);
    hM_VtxY_TruthReco->GetXaxis()->SetTitle("Truth vtx y (cm)");
    hM_VtxY_TruthReco->GetYaxis()->SetTitle("Reco PV y (cm)");
    hM_VtxY_TruthReco->Draw("colz");
    c1->SaveAs("./plot/PV_TruthReco/hM_VtxY_TruthReco_v2.pdf");
}