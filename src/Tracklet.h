#ifndef TRACKLET_H
#define TRACKLET_H

#include <stdio.h>
#include <iostream>
#include <numeric>

#include <TObject.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TFile.h>

#include "Utilities.h"
#include "Hit.h"
#include "Hists.h"
#include "GenHadron.h"
#include "TrackletData.h"

using namespace std;

class Tracklet : public TObject
{
public:
    Tracklet();
    Tracklet(Hit *hit1, Hit *hit2);
    ~Tracklet();

    void SetHits(Hit *hit1, Hit *hit2);
    Hit *Hit1();
    Hit *Hit2();
    float dEta();
    float dPhi();
    float dR();
    float Eta();
    float Phi();
    float tklVtxZ();
    int LayerComb(); // Layer combination: can only be 12, 23, 13
    void SetMatchedGenHardon();
    bool IsMatchedGenHadron();

private:
    Hit *_hit1;
    Hit *_hit2;
    float _deta;
    float _dphi;
    bool _matched_genhadron;
};

Tracklet::Tracklet()
{
    _hit1 = nullptr;
    _hit2 = nullptr;
    _matched_genhadron = false;
}

Tracklet::Tracklet(Hit *hit1, Hit *hit2)
{
    _hit1 = hit1;
    _hit2 = hit2;
    _deta = _hit1->Eta() - _hit2->Eta();
    _dphi = deltaPhi(_hit1->Phi(), _hit2->Phi());
    _matched_genhadron = false;
}

Tracklet::~Tracklet()
{
}

void Tracklet::SetHits(Hit *hit1, Hit *hit2)
{
    _hit1 = hit1;
    _hit2 = hit2;
    _deta = _hit1->Eta() - _hit2->Eta();
    _dphi = deltaPhi(_hit1->Phi(), _hit2->Phi());
    _matched_genhadron = false;
}

Hit *Tracklet::Hit1()
{
    return (_hit1);
}

Hit *Tracklet::Hit2()
{
    return (_hit2);
}

float Tracklet::dEta()
{
    return (_deta);
}

float Tracklet::dPhi()
{
    return (_dphi);
}

float Tracklet::dR()
{
    return (sqrt(_deta * _deta + _dphi * _dphi));
}

float Tracklet::Eta()
{
    return (_hit1->Eta());
}

float Tracklet::Phi()
{
    return (_hit1->Phi());
}

float Tracklet::tklVtxZ()
{
    return (_hit1->posZ() - _hit1->rho() * (_hit2->posZ() - _hit1->posZ()) / (_hit2->rho() - _hit1->rho()));
}

int Tracklet::LayerComb()
{
    return (_hit1->Layer()) * 10 + (_hit2->Layer());
}

void Tracklet::SetMatchedGenHardon() { _matched_genhadron = true; }

bool Tracklet::IsMatchedGenHadron() { return _matched_genhadron; }

bool compare_dR(Tracklet *a, Tracklet *b) { return a->dR() < b->dR(); }

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

float TrackletPV_cluster(int evt, vector<Hit *> layer1, vector<Hit *> layer2, int dPhiCutbin, int dZCutbin, bool verbose)
{
    float dPhi_cut = dPhiCutbin * 0.01;
    float dZ_cut = dZCutbin * 0.01;

    double maxNz = 0;
    double maxTotalZ = 0;
    double maxRMS = 10e10;
    double nRecoZ = 0;

    vector<double> vectorZ;

    // if (verbose)
    // {
    //     int nhitsl1_calc = (layer1.size() > 1000) ? (int)layer1.size() / 5 : layer1.size();
    //     int nhitsl2_calc = (layer2.size() > 1000) ? (int)layer2.size() / 5 : layer2.size();
    //     cout << __LINE__ << " [PVfinder-TklCluster Verbose] number of hits to be used in layer1 = " << nhitsl1_calc << "; number of hits to be used in layer2 = " << nhitsl2_calc << endl;
    // }

    int runsf = ((int)layer1.size() / 4000);
    // cout << runsf << endl;

    for (size_t i = 0; i < layer1.size(); i++)
    {
        if (fabs(layer1[i]->Eta()) > 1.)
            continue;
        // If there are more than 1000 first layer hits in the event, only 1/5 of the hits are used for proto-tracklet reconstruction to reduce the reconstruction time.
        // if (i % 5 != 0 && layer1.size() > 1000)
        if (i % (runsf + 1) != 0 && layer1.size() > 1000)
            continue;
        for (size_t j = 0; j < layer2.size(); j++)
        {
            if (fabs(layer2[j]->Eta()) > 1.)
                continue;
            // if (j % 5 != 0 && layer1.size() > 1000)
            if (j % (runsf + 1) != 0 && layer1.size() > 1000)
                continue;

            if (fabs(deltaPhi(layer1[i]->Phi(), layer2[j]->Phi())) > dPhi_cut)
                continue;

            float z = layer1[i]->posZ() - (layer2[j]->posZ() - layer1[i]->posZ()) / (layer2[j]->rho() - layer1[i]->rho()) * layer1[i]->rho();

            if (fabs(z) < 20.)
            {
                nRecoZ++;
                vectorZ.push_back(z);
            }
        }
    }

    sort(vectorZ.begin(), vectorZ.end());

    if (verbose)
        cout << __LINE__ << " [PVfinder-TklCluster Verbose] size of vectorZ = " << vectorZ.size() << endl;

    TH1F *hM_rmsAll = new TH1F("hM_rmsAll", "hM_rmsAll", 100, 0, 0.02);
    TH1F *hM_vtxclusterZ = new TH1F("hM_vtxclusterZ", "hM_vtxclusterZ", 200,-10,10);
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

        if (nz > maxNz || (nz == maxNz && rms < maxRMS))
        {
            maxNz = nz;
            maxTotalZ = totalZ;
            maxRMS = rms;
        }
        // cout << "(i,sumZ,N,sumZ/N,vtx_rms)=(" << i << "," << sumZ << "," << N << "," << (sumZ / N) << "," << vtx_rms << ")" << endl;
    }

    // store the vertex clusters histogram
    system(Form("mkdir -p /sphenix/user/hjheng/TrackletAna/minitree/AuAu_NoPileup_RecoVtx_Optimization/hists/dPhiCutbin%d_dZCutbin%d/", dPhiCutbin, dZCutbin));
    TFile *f = new TFile(Form("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_NoPileup_RecoVtx_Optimization/hists/dPhiCutbin%d_dZCutbin%d/VtxCluster_hist_dPhiCutbin%d_dZCutbin%d_event%d.root", dPhiCutbin, dZCutbin, dPhiCutbin, dZCutbin, evt), "RECREATE");
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
            cout << __LINE__ << " [PVfinder-TklCluster Verbose] final vtxZ = -999." << endl;
        return -999.;
    }
    else
    {
        if (verbose)
            cout << __LINE__ << " [PVfinder-TklCluster Verbose] final vtxZ = " << (maxTotalZ / maxNz) << endl;
        return (maxTotalZ / maxNz);
    }
}

vector<Tracklet *> ProtoTracklets(vector<Hit *> layer1, vector<Hit *> layer2, Hists &Hists)
{
    vector<Tracklet *> _ProtoTkls;
    _ProtoTkls.clear();

    // int nhits1 = layer1.size();
    // int nhits2 = layer2.size();

    float Cut_dEta = 0.5;
    float Cut_dR = 0.5;

    // "1st/2nd" means the first/second layer to look at here
    // fprintf(stderr, "Nhits in (1st Layer, 2nd Layer2) = (%d, %d). # of all combinations = %d \n", nhits1, nhits2, nhits1 * nhits2);

    int iComb = 0;
    // Get all combinations
    for (auto ihitl1 : layer1)
    {
        for (auto ihitl2 : layer2)
        {
            // iComb++;
            // if (iComb % 500000 == 0)
            //     fprintf(stderr, "processing %d of %d combinations (%.3f %%)\n", iComb, nhits1 * nhits2, (float)iComb / (nhits1 * nhits2) * 100.);

            float dEta = ihitl2->Eta() - ihitl1->Eta();
            float dPhi = deltaPhi(ihitl1->Phi(), ihitl2->Phi());
            float dR = sqrt(dEta * dEta + dPhi * dPhi);

            Hists.hM_dPhi_dEta_woCuts->Fill(dPhi, dEta);

            if (dR < Cut_dR)
            {
                Tracklet *tmptkl = new Tracklet(ihitl1, ihitl2);
                _ProtoTkls.push_back(tmptkl);

                // cout << ihitl1->Eta() << " " << tmptkl->Eta() << endl;
                Hists.hM_Eta_hit1_proto->Fill(ihitl1->Eta());
                Hists.hM_Phi_hit1_proto->Fill(ihitl1->Phi());
                Hists.hM_Eta_hit2_proto->Fill(ihitl2->Eta());
                Hists.hM_Phi_hit2_proto->Fill(ihitl2->Phi());

                Hists.hM_Eta_proto->Fill(tmptkl->Eta());
                if (ihitl1->vtxZ() < -7)
                    Hists.hM_Eta_proto_PVzRange1->Fill(tmptkl->Eta());
                else if (ihitl1->vtxZ() >= -7 && ihitl1->vtxZ() < -2.5)
                    Hists.hM_Eta_proto_PVzRange2->Fill(tmptkl->Eta());
                else if (ihitl1->vtxZ() >= -2.5 && ihitl1->vtxZ() < 2.5)
                    Hists.hM_Eta_proto_PVzRange3->Fill(tmptkl->Eta());
                else if (ihitl1->vtxZ() >= 2.5 && ihitl1->vtxZ() < 7)
                    Hists.hM_Eta_proto_PVzRange4->Fill(tmptkl->Eta());
                else if (ihitl1->vtxZ() >= 7)
                    Hists.hM_Eta_proto_PVzRange5->Fill(tmptkl->Eta());
                else
                    fprintf(stderr, "Weird! PVz = %f \n", ihitl1->vtxZ());

                Hists.hM_Phi_proto->Fill(tmptkl->Phi());
                Hists.hM_dEta_proto->Fill(tmptkl->dEta());
                Hists.hM_dEta_proto_altrange->Fill(tmptkl->dEta());
                Hists.hM_dEta_proto_altrange2->Fill(tmptkl->dEta());
                Hists.hM_dPhi_proto->Fill(tmptkl->dPhi());
                Hists.hM_dPhi_proto_altrange->Fill(tmptkl->dPhi());
                Hists.hM_dR_proto->Fill(tmptkl->dR());
                Hists.hM_dR_proto_altrange->Fill(tmptkl->dR());
                Hists.hM_dR_proto_LogX->Fill(tmptkl->dR());
                Hists.hM_hit1Eta_hit2Eta_proto->Fill(ihitl1->Eta(), ihitl2->Eta());
                Hists.hM_hit1Phi_hit2Phi_proto->Fill(ihitl1->Phi(), ihitl2->Phi());
                Hists.hM_Eta_vtxZ_proto->Fill(tmptkl->Eta(), ihitl1->vtxZ());
                Hists.hM_Eta_vtxZ_anabin_proto->Fill(tmptkl->Eta(), ihitl1->vtxZ());
                Hists.hM_dPhi_dEta_proto->Fill(tmptkl->dPhi(), tmptkl->dEta());
                Hists.hM_dPhi_dEta_proto_altrange->Fill(tmptkl->dPhi(), tmptkl->dEta());
            }
            else
                continue;
        }
    }

    return _ProtoTkls;
}

vector<Tracklet *> RecoTracklets(vector<Tracklet *> &ProtoTkls, Hists &Hists)
{
    sort(ProtoTkls.begin(), ProtoTkls.end(), compare_dR);

    vector<Tracklet *> _RecoTkls;
    _RecoTkls.clear();

    int itkl = 0;
    for (auto &tkl : ProtoTkls)
    {
        // printf("ith prototkl = %d,\t dR = %f,\t hit1 eta = %f,\t hit1 phi = %f,\t IsMatchedTkl = %s,\t hit2 eta = %f,\t hit2 phi = %f,\t IsMatchedTkl = %s,\t # of matched hits @ L1 = %d\n", itkl, tkl->dR(), tkl->Hit1()->Eta(), tkl->Hit1()->Phi(), tkl->Hit1()->IsMatchedTkl() ? "true" : "false", tkl->Hit2()->Eta(), tkl->Hit2()->Phi(), tkl->Hit2()->IsMatchedTkl() ? "true" : "false", tkl->Hit2()->MatchedL1Hits().size());

        if (tkl->Hit1()->IsMatchedTkl() || tkl->Hit2()->IsMatchedTkl())
        {
            continue;
        }
        else
        {
            _RecoTkls.push_back(tkl);
            tkl->Hit1()->SetMatchedTkl();
            tkl->Hit2()->SetMatchedTkl();

            Hists.hM_Eta_hit1_reco->Fill(tkl->Hit1()->Eta());
            Hists.hM_Phi_hit1_reco->Fill(tkl->Hit1()->Phi());
            Hists.hM_Eta_hit2_reco->Fill(tkl->Hit2()->Eta());
            Hists.hM_Phi_hit2_reco->Fill(tkl->Hit2()->Phi());
            Hists.hM_Eta_reco->Fill(tkl->Eta());
            Hists.hM_Phi_reco->Fill(tkl->Phi());
            Hists.hM_dEta_reco->Fill(tkl->dEta());
            Hists.hM_dEta_reco_altrange->Fill(tkl->dEta());
            Hists.hM_dEta_reco_altrange2->Fill(tkl->dEta());
            Hists.hM_dPhi_reco->Fill(tkl->dPhi());
            Hists.hM_dPhi_reco_altrange->Fill(tkl->dPhi());
            Hists.hM_dR_reco->Fill(tkl->dR());
            Hists.hM_dR_reco_altrange->Fill(tkl->dR());
            Hists.hM_dR_reco_LogX->Fill(tkl->dR());
            Hists.hM_hit1Eta_hit2Eta_reco->Fill(tkl->Hit1()->Eta(), tkl->Hit2()->Eta());
            Hists.hM_hit1Phi_hit2Phi_reco->Fill(tkl->Hit1()->Phi(), tkl->Hit2()->Phi());
            Hists.hM_Eta_vtxZ_reco->Fill(tkl->Eta(), tkl->Hit1()->vtxZ());
            Hists.hM_Eta_vtxZ_anabin_reco->Fill(tkl->Eta(), tkl->Hit1()->vtxZ());
            Hists.hM_dPhi_dEta_reco->Fill(tkl->dPhi(), tkl->dEta());
            Hists.hM_dPhi_dEta_reco_altrange->Fill(tkl->dPhi(), tkl->dEta());
        }
    }

    return _RecoTkls;
}

vector<Tracklet *> GenMatch_Recotkl(vector<Tracklet *> &RecoTkls, vector<GenHadron *> &GenHadrons, Hists &Hists)
{
    vector<Tracklet *> _RecoTkls_GMatched;
    _RecoTkls_GMatched.clear();

    for (auto &tkl : RecoTkls)
    {
        if (tkl->IsMatchedGenHadron())
            continue;
        for (auto &ghadron : GenHadrons)
        {
            if (ghadron->IsMatchedToRecotkl())
                continue;
            // Machting criteria
            if (deltaR(tkl->Eta(), tkl->Phi(), ghadron->Eta(), ghadron->Phi()) > 0.05)
                continue;
            else
            {
                tkl->SetMatchedGenHardon();
                ghadron->SetMatchedToRecotkl();

                _RecoTkls_GMatched.push_back(tkl);

                Hists.hM_Eta_vtxZ_reco_genmatched->Fill(tkl->Eta(), tkl->Hit1()->vtxZ());
                Hists.hM_Eta_vtxZ_anabin_reco_genmatched->Fill(tkl->Eta(), tkl->Hit1()->vtxZ());
            }
        }
    }

    return _RecoTkls_GMatched;
}

#endif