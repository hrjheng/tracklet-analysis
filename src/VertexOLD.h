#ifndef VERTEXOLD_H
#define VERTEXOLD_H

#include <vector>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <stdlib.h>
#include <map>
#include <numeric>

#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>

#include "/sphenix/user/hjheng/TrackletAna/analysis/plot/sPHENIXStyle/sPhenixStyle.C"

using namespace std;

struct VtxData
{
    int event, NhitsLayer1, NTruthVtx;
    float PV_z, TruthPV_trig_z, TruthPV_mostNpart_z;
};

float RMS(vector<float> v)
{
    float sum = std::accumulate(v.begin(), v.end(), 0.0);
    float mean = sum / v.size();

    std::vector<float> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](float x)
                   { return x - mean; });
    float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / v.size());
    return stdev;
}

TCanvas* Canv_RecoVtxZTklZRho(vector<TLine *> v_line, vector<Hit *> hit_l1, vector<Hit *> hit_l2, TH2F *hM_PVz, TH2F *hM_truthPVz)
{
    TH2F *hM = new TH2F("hM", "hM", 2000,-20,20,900,-0.5,4);
    for(size_t i = 0; i < hit_l1.size(); i++)
    {
        hM->Fill(hit_l1[i]->posZ(), hit_l1[i]->rho());
    }
    for(size_t i = 0; i < hit_l2.size(); i++)
    {
        hM->Fill(hit_l2[i]->posZ(), hit_l2[i]->rho());
    }

    TCanvas *c = new TCanvas("c", "c", 800, 700);
    c->cd();
    gPad->SetRightMargin(0.08);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    hM->GetXaxis()->SetTitle("Position Z (cm)");
    hM->GetYaxis()->SetTitle("Radius (cm)");
    hM->GetXaxis()->SetTickSize(0.03);
    hM->GetYaxis()->SetTickSize(0.03);
    hM->GetXaxis()->SetTitleSize(0.05);
    hM->GetYaxis()->SetTitleSize(0.05);
    hM->GetXaxis()->SetLabelSize(0.04);
    hM->GetYaxis()->SetLabelSize(0.04);
    hM->GetXaxis()->SetTitleOffset(1.1);
    hM->GetYaxis()->SetTitleOffset(1.3);
    hM->SetMarkerSize(0.5);
    hM->Draw("SCAT");

    for(size_t i = 0; i < v_line.size(); i++)
    {
        v_line[i]->SetLineColorAlpha(14,0.5);
        v_line[i]->Draw("same");
    }

    hM->Draw("SCATsame");

    // MVTX dimension
    float MVTXLength = 27.1; // cm 
    float MVTXL1_min = 2.461; 
    float MVTXL1_mid = 2.523;
    float MVTXL1_max = 2.793;
    float MVTXL2_min = 3.198;
    float MVTXL2_mid = 3.335;
    float MVTXL2_max = 3.625;
    TLine *beamline = new TLine(-20,0,20,0);
    TLine *MVTXl1_in = new TLine(-(MVTXLength/2),MVTXL1_min,MVTXLength/2,MVTXL1_min);
    TLine *MVTXl1_out = new TLine(-(MVTXLength/2),MVTXL1_max,MVTXLength/2,MVTXL1_max);
    TLine *MVTXl2_in = new TLine(-(MVTXLength/2),MVTXL2_min,MVTXLength/2,MVTXL2_min);
    TLine *MVTXl2_out = new TLine(-(MVTXLength/2),MVTXL2_max,MVTXLength/2,MVTXL2_max);
    MVTXl1_in->SetLineStyle(kDashed);
    MVTXl1_out->SetLineStyle(kDashed);
    MVTXl2_in->SetLineStyle(kDashed);
    MVTXl2_out->SetLineStyle(kDashed);
    beamline->Draw("same");
    MVTXl1_in->Draw("same");
    MVTXl1_out->Draw("same");
    MVTXl2_in->Draw("same");
    MVTXl2_out->Draw("same");

    hM_PVz->SetMarkerSize(0.8);
    hM_PVz->SetMarkerColorAlpha(46,0.9);
    hM_PVz->Draw("SCATsame");

    hM_truthPVz->SetMarkerStyle(34);
    hM_truthPVz->SetMarkerSize(0.8);
    hM_truthPVz->SetMarkerColorAlpha(38,0.9);
    hM_truthPVz->Draw("SCATsame");

    TText *t_l1 = new TText(MVTXLength/2,MVTXL1_min,"Layer 1");
    TText *t_l2 = new TText(MVTXLength/2,MVTXL2_min,"Layer 2");
    t_l1->SetTextSize(0.04);
    t_l2->SetTextSize(0.04);
    t_l1->SetTextAlign(11);
    t_l2->SetTextAlign(11);
    t_l1->Draw("same");
    t_l2->Draw("same");

    TLegend *leg = new TLegend(0.5,0.15,0.88,0.2);
    leg->SetNColumns(3);
    leg->SetTextSize(0.03);
    leg->AddEntry(hM,"Clusters","p");
    leg->AddEntry(hM_PVz,"Reco PV","p");
    leg->AddEntry(hM_truthPVz,"Truth PV","p");
    leg->Draw("same");

    return c;
}

// Old implementation
float TrackletPV_cluster_OBSOLETE(int evt, float TruthPVz, vector<Hit *> layer1, vector<Hit *> layer2, int dPhiCutbin, int dZCutbin, bool verbose)
{
    float dPhi_cut = dPhiCutbin * 0.01;
    float dZ_cut = dZCutbin * 0.01;

    double maxNz = 0;
    double maxTotalZ = 0;
    double maxRMS = 10e10;
    double nRecoZ = 0;

    vector<double> vectorZ;
    vectorZ.clear();
    vector<TLine *> v_tlkprojZ;
    v_tlkprojZ.clear();

    if (verbose)
    {
        int nhitsl1_calc = (layer1.size() > 1000) ? (int)layer1.size() / 5 : layer1.size();
        int nhitsl2_calc = (layer2.size() > 1000) ? (int)layer2.size() / 5 : layer2.size();
        cout << __LINE__ << " [PVfinder-TklCluster Verbose] number of hits to be used in layer1 = " << nhitsl1_calc << "; number of hits to be used in layer2 = " << nhitsl2_calc << endl;
    }

    // int runsf = ((int)layer1.size() / 4000);
    // cout << runsf << endl;

    for (size_t i = 0; i < layer1.size(); i++)
    {
        // if (fabs(layer1[i]->Eta()) > 1.)
        //     continue;
        // If there are more than 1000 first layer hits in the event, only 1/5 of the hits are used for proto-tracklet reconstruction to reduce the reconstruction time.
        if (i % 5 != 0 && layer1.size() > 1000)
            continue;
        // if (i % (runsf + 1) != 0)
        //     continue;
            
        for (size_t j = 0; j < layer2.size(); j++)
        {
            // if (fabs(layer2[j]->Eta()) > 1.)
            //     continue;
            if (j % 5 != 0 && layer1.size() > 1000)
                continue;
            // if (j % (runsf + 1) != 0)
            //     continue;

            if (fabs(deltaPhi(layer1[i]->Phi(), layer2[j]->Phi())) > dPhi_cut)
                continue;

            float z = layer1[i]->posZ() - (layer2[j]->posZ() - layer1[i]->posZ()) / (layer2[j]->rho() - layer1[i]->rho()) * layer1[i]->rho();
            
            if (fabs(z) < 20.)
            {
                nRecoZ++;
                vectorZ.push_back(z);

                // Check
                if (layer1.size() < 101)
                {
                    TLine *l = new TLine(layer2[j]->posZ(), layer2[j]->rho(), z , 0.0);
                    v_tlkprojZ.push_back(l);
                }
            }
        }
    }

    sort(vectorZ.begin(), vectorZ.end());

    if (verbose)
        cout << __LINE__ << " [PVfinder-TklCluster Verbose] size of vectorZ = " << vectorZ.size() << endl;

    TH1F *hM_rmsAll = new TH1F("hM_rmsAll", "hM_rmsAll", 100, 0, 0.02);
    TH1F *hM_vtxclusterZ = new TH1F("hM_vtxclusterZ", "hM_vtxclusterZ", 200,-20,20);
    for (size_t i = 0; i < vectorZ.size(); i++)
    {
        double nz = 0;
        double totalZ = 0.;
        double rms = 0.;
        // vector<float> tmp_v_vtxZ;
        nz++;

        TH1F *h = new TH1F("h", "", 100, vectorZ[i] - dZ_cut, vectorZ[i] + dZ_cut);
        totalZ += vectorZ[i];
        h->Fill(vectorZ[i]);
        hM_vtxclusterZ->Fill(vectorZ[i]);
        // tmp_v_vtxZ.push_back(vectorZ[i]);

        int flag = 0;
        for (size_t j = 0; j < vectorZ.size(); j++)
        {
            if (fabs(vectorZ[j] - vectorZ[i]) < dZ_cut && i != j)
            {
                nz++;
                totalZ += vectorZ[j];
                // tmp_v_vtxZ.push_back(vectorZ[j]);
                h->Fill(vectorZ[j]);
                hM_vtxclusterZ->Fill(vectorZ[j]);
                flag = 1;
            }
            else
            {
                if (flag == 1)
                    continue;
            }
        }

        // vtx_rms = RMS(tmp_v_vtxZ);
        // CleanVec(tmp_v_vtxZ);
        rms = h->GetRMS();
        hM_rmsAll->Fill(rms);
        delete h;

        if (verbose)
            cout << __LINE__ << " [PVfinder-TklCluster Verbose] (nz,totalZ,rms) = (" << nz << "," << totalZ << "," << rms << ")" << endl;

        if (nz > maxNz || (nz == maxNz && rms < maxRMS))
        {
            maxNz = nz;
            maxTotalZ = totalZ;
            maxRMS = rms;
        }
    }

    system(Form("mkdir -p /sphenix/user/hjheng/TrackletAna/analysis/plot/RecoPV_optimization/dPhiCutBin%d_dZCutBin%d/event/ClusZRho_PVClusZ", dPhiCutbin, dZCutbin));
    if (layer1.size() < 101)
    {
        float final_RecoPVz;
        if (maxNz == 0 || nRecoZ == 0)
        {
            final_RecoPVz = -999.;
        }
        else
        {
            final_RecoPVz = maxTotalZ / maxNz;
        }
        TH2F *hM_finalPVz = new TH2F("hM_finalPVz","hM_finalPVz",2000,-20,20,900,-0.5,4);
        hM_finalPVz->Fill(final_RecoPVz, 0.0);
        TH2F *hM_TruthPVz = new TH2F("hM_TruthPVz","hM_TruthPVz",2000,-20,20,900,-0.5,4);
        hM_TruthPVz->Fill(TruthPVz, 0.0);
        cout << "final RecoVtxZ=" << final_RecoPVz << endl;

        SetsPhenixStyle();
        TCanvas *canv_tlkprojZ = Canv_RecoVtxZTklZRho(v_tlkprojZ, layer1, layer2, hM_finalPVz, hM_TruthPVz);
        canv_tlkprojZ->RedrawAxis();
        canv_tlkprojZ->Draw();
        canv_tlkprojZ->SaveAs(Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/RecoPV_optimization/dPhiCutBin%d_dZCutBin%d/event/ClusZRho_PVClusZ/RecoVtxZTklZRho_evt%d.pdf",dPhiCutbin, dZCutbin, evt));
        canv_tlkprojZ->SaveAs(Form("/sphenix/user/hjheng/TrackletAna/analysis/plot/RecoPV_optimization/dPhiCutBin%d_dZCutBin%d/event/ClusZRho_PVClusZ/RecoVtxZTklZRho_evt%d.png",dPhiCutbin, dZCutbin, evt));
        canv_tlkprojZ->Close();
        delete canv_tlkprojZ;
        delete hM_finalPVz;
        delete hM_TruthPVz;
    }
    
    // store the vertex clusters histogram
    system(Form("mkdir -p /sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/hists/dPhiCutbin%d_dZCutbin%d/", dPhiCutbin, dZCutbin));
    TFile *f = new TFile(Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/hists/dPhiCutbin%d_dZCutbin%d/VtxCluster_hist_dPhiCutbin%d_dZCutbin%d_event%d.root", dPhiCutbin, dZCutbin, dPhiCutbin, dZCutbin, evt), "RECREATE");
    hM_rmsAll->Write();
    hM_vtxclusterZ->Write();
    f->Close();
    delete hM_rmsAll;
    delete hM_vtxclusterZ;

    if (verbose)
        cout << __LINE__ << " [PVfinder-TklCluster Verbose] " << maxNz << " " << maxTotalZ << " " << maxRMS << endl;

    CleanVec(vectorZ);
    if (maxNz == 0 || nRecoZ == 0)
    {
        if (verbose)
        {
            cout << __LINE__ << " [PVfinder-TklCluster Verbose] final vtxZ = -999." << endl;
            cout << "---------------------------------------" << endl;
        }
        return -999.;
    }
    else
    {
        if (verbose)
        {
            cout << __LINE__ << " [PVfinder-TklCluster Verbose] final vtxZ = " << (maxTotalZ / maxNz) << endl;
            cout << "---------------------------------------" << endl;
        } 
        return (maxTotalZ / maxNz);
    }
}

#endif