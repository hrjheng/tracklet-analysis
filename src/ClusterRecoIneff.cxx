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

#include "Utilities.h"

using namespace std;

int main(int argc, char *argv[])
{
    // The probability for a pixel hit to not be reconstructed is estimated by studying the fraction of reconstructed trackets from the 1st and 3rd layers that have a corresponding pixel hit in the second layer.
    vector<TString> infilename = {"/sphenix/user/hjheng/TrackletAna/minitree/SimplePion/TrackletAna_minitree_layer13_Evt0to500_RandhitCase0_ClusSplitCase0_MisAlignNum0_dRcut0p5.root"};
    vector<TString> outfilename = {"/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/SimplePions/Hists_ClusterRecoIneff_SimplePions.root"};

    // vector<TString> infilename = {"/sphenix/user/hjheng/TrackletAna/minitree/AuAu_Nominal_NoPileup/TrackletAna_minitree_layer13_Evt0to2000_RandhitCase0_ClusSplitCase0_MisAlignNum0_dRcut0p5.root"};
    // vector<TString> outfilename = {"/sphenix/user/hjheng/TrackletAna/analysis/plot/hists/HijingAuAuNoPileup_ana325private/Hists_ClusterRecoIneff_HijingAuAuNoPileup_ana325private.root"};

    vector<TH1F *> hM_recotkll2_eta_all = {new TH1F("hM_recotkll2_eta_sim", "hM_recotkll2_eta_sim", 70, -3.5, 3.5)};
    vector<TH1F *> hM_recotkll2_eta_nonmatched = {new TH1F("hM_recotkll2_eta_nonmatched_sim", "hM_recotkll2_eta_nonmatched_sim", 70, -3.5, 3.5)};

    for (size_t i = 0; i < infilename.size(); i++)
    {
        TFile *f = new TFile(infilename[i], "READ");
        TTree *t = (TTree *)f->Get("minitree_13");
        t->BuildIndex("event"); // Reference: https://root-forum.cern.ch/t/sort-ttree-entries/13138
        TTreeIndex *index = (TTreeIndex *)t->GetTreeIndex();
        int event;
        vector<float> *recotklraw_eta = 0, *recotklraw_phi = 0;
        vector<float> *recoclus_eta_l2 = 0, *recoclus_phi_l2 = 0;
        t->SetBranchAddress("event", &event);
        t->SetBranchAddress("recotklraw_eta", &recotklraw_eta);
        t->SetBranchAddress("recotklraw_phi", &recotklraw_phi);
        t->SetBranchAddress("recoclus_eta_l2", &recoclus_eta_l2);
        t->SetBranchAddress("recoclus_phi_l2", &recoclus_phi_l2);
        for (int ev = 0; ev < t->GetEntriesFast(); ev++)
        {
            Long64_t local = t->LoadTree(index->GetIndex()[ev]);
            t->GetEntry(local);
            // Match the tracklets from the 1st and 3rd layers with the clusters from the 2nd layer

            int Nrecotklraw = recotklraw_eta->size(), Nrecoclus_l2 = recoclus_eta_l2->size();
            cout << event << " " << Nrecotklraw << " " << Nrecoclus_l2 << endl;
            bool tklmatched[Nrecotklraw] = {false}, clusmatched[Nrecoclus_l2] = {false};

            for (size_t j = 0; j < recotklraw_eta->size(); j++)
            {
                if (tklmatched[j])
                    continue;

                for (size_t k = 0; k < recoclus_eta_l2->size(); k++)
                {
                    if (clusmatched[k])
                        continue;
                    if (deltaR(recotklraw_eta->at(j), recotklraw_phi->at(j), recoclus_eta_l2->at(k), recoclus_phi_l2->at(k)) < 0.01)
                    {
                        // cout << "Matched: (j,k)=(" << j << "," << k << "), (tkl eta, tkl phi, clus eta, clus phi) = (" << recotklraw_eta->at(j) << ", " << recotklraw_phi->at(j) << ", " << recoclus_eta_l2->at(k) << ", " << recoclus_phi_l2->at(k) << ")" << endl;
                        // the 1st and 3rd layers tracklets are matched with the 2nd layer clusters
                        tklmatched[j] = true;
                        clusmatched[k] = true;
                    }
                }
            }

            // Fill the histograms
            for (size_t j = 0; j < recotklraw_eta->size(); j++)
            {
                hM_recotkll2_eta_all[i]->Fill(recotklraw_eta->at(j));
                if (!tklmatched[j])
                    hM_recotkll2_eta_nonmatched[i]->Fill(recotklraw_eta->at(j));
            }
        }

        f->Close();

        TFile *fout = new TFile(outfilename[i], "RECREATE");
        fout->cd();
        hM_recotkll2_eta_all[i]->Write();
        hM_recotkll2_eta_nonmatched[i]->Write();
        fout->Close();
    }

    return 0;
}