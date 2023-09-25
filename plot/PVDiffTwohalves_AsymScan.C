#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "/sphenix/user/hjheng/TrackletAna/analysis/plot/sPHENIXStyle/sPhenixStyle.C"

float sigmaEff(vector<float> v, float threshold, float &xmin, float &xmax)
{

    std::sort(v.begin(), v.end());

    int total = v.size();
    int max = (int)(threshold * total);

    vector<float> start;
    vector<float> stop;
    vector<float> width;

    unsigned i = 0;
    while (i != v.size() - 1)
    {

        int count = 0;
        unsigned j = i;
        while (j != v.size() - 1 && count < max)
        {

            ++count;
            ++j;
        }

        if (j != v.size() - 1)
        {
            start.push_back(v[i]);
            stop.push_back(v[j]);
            width.push_back(v[j] - v[i]);
        }
        ++i;
    }

    float minwidth = *min_element(width.begin(), width.end());

    unsigned pos = min_element(width.begin(), width.end()) - width.begin();

    xmin = start[pos];
    xmax = stop[pos];
    return minwidth;
}

std::tuple<float, float, float, float, float, float, float, float, float, float, float, float> getPVzDiffSigma(TString fname)
{
    std::vector<float> list_diffPVx, list_diffPVy, list_diffPVz, list_ResPVx, list_ResPVy, list_ResPVz;
    list_diffPVx.clear();
    list_diffPVy.clear();
    list_diffPVz.clear();
    list_ResPVx.clear();
    list_ResPVy.clear();
    list_ResPVz.clear();

    TFile *f = new TFile(fname.Data(), "READ");
    TTree *t = (TTree *)f->Get("minitree");
    int NTruthVtx;
    float PV_x_halves1, PV_x_halves2, PV_y_halves1, PV_y_halves2, PV_z_halves1, PV_z_halves2, TruthPV_mostNpart_x, TruthPV_mostNpart_y, TruthPV_mostNpart_z, TruthPV_trig_x, TruthPV_trig_y, TruthPV_trig_z;
    t->SetBranchAddress("NTruthVtx", &NTruthVtx);
    t->SetBranchAddress("PV_x_halves1", &PV_x_halves1);
    t->SetBranchAddress("PV_x_halves2", &PV_x_halves2);
    t->SetBranchAddress("PV_y_halves1", &PV_y_halves1);
    t->SetBranchAddress("PV_y_halves2", &PV_y_halves2);
    t->SetBranchAddress("PV_z_halves1", &PV_z_halves1);
    t->SetBranchAddress("PV_z_halves2", &PV_z_halves2);
    t->SetBranchAddress("TruthPV_mostNpart_x", &TruthPV_mostNpart_x);
    t->SetBranchAddress("TruthPV_mostNpart_y", &TruthPV_mostNpart_y);
    t->SetBranchAddress("TruthPV_mostNpart_z", &TruthPV_mostNpart_z);
    t->SetBranchAddress("TruthPV_trig_x", &TruthPV_trig_x);
    t->SetBranchAddress("TruthPV_trig_y", &TruthPV_trig_y);
    t->SetBranchAddress("TruthPV_trig_z", &TruthPV_trig_z);
    for (Long64_t idx = 0; idx < t->GetEntries(); ++idx)
    {
        t->GetEntry(idx);
        
        if (NTruthVtx != 1)
            continue;

        list_diffPVx.push_back(PV_x_halves1 - PV_x_halves2);
        list_diffPVy.push_back(PV_y_halves1 - PV_y_halves2);
        list_diffPVz.push_back(PV_z_halves1 - PV_z_halves2);
        float finalPVx = (PV_x_halves1 + PV_x_halves2) / 2.;
        float finalPVy = (PV_y_halves1 + PV_y_halves2) / 2.;
        float finalPVz = (PV_z_halves1 + PV_z_halves2) / 2.;
        // list_ResPVx.push_back(finalPVx - TruthPV_trig_x);
        // list_ResPVy.push_back(finalPVy - TruthPV_trig_y);
        // list_ResPVz.push_back(finalPVz - TruthPV_trig_z);
        list_ResPVx.push_back(finalPVx - TruthPV_mostNpart_x);
        list_ResPVy.push_back(finalPVy - TruthPV_mostNpart_y);
        list_ResPVz.push_back(finalPVz - TruthPV_mostNpart_z);
    }

    f->Close();

    float datax_min = -999.0, datax_max = 999.0, datay_min = -990., detay_max = 999., dataz_min = -999.0, dataz_max = 999.0, dataresx_min = -999.0, dataresx_max = 999.0, dataresy_min = -999.0,
          dataresy_max = 999.0, dataresz_min = -999.0, dataresz_max = 999.0;
    float threshold = 0.685;
    float se_diffPVx = sigmaEff(list_diffPVx, threshold, datax_min, datax_max);
    float se_diffPVy = sigmaEff(list_diffPVy, threshold, datay_min, detay_max);
    float se_diffPVz = sigmaEff(list_diffPVz, threshold, dataz_min, dataz_max);
    float se_resPVx = sigmaEff(list_ResPVx, threshold, dataresx_min, dataresx_max);
    float se_resPVy = sigmaEff(list_ResPVy, threshold, dataresy_min, dataresy_max);
    float se_resPVz = sigmaEff(list_ResPVz, threshold, dataresz_min, dataresz_max);

    // Calculate the median
    std::sort(list_diffPVx.begin(), list_diffPVx.end());
    float median_diffPVx = (list_diffPVx.size() % 2 == 0) ? (list_diffPVx[list_diffPVx.size() / 2 - 1] + list_diffPVx[list_diffPVx.size() / 2]) / 2 : list_diffPVx[list_diffPVx.size() / 2];
    std::sort(list_diffPVy.begin(), list_diffPVy.end());
    float median_diffPVy = (list_diffPVy.size() % 2 == 0) ? (list_diffPVy[list_diffPVy.size() / 2 - 1] + list_diffPVy[list_diffPVy.size() / 2]) / 2 : list_diffPVy[list_diffPVy.size() / 2];
    std::sort(list_diffPVz.begin(), list_diffPVz.end());
    float median_diffPVz = (list_diffPVz.size() % 2 == 0) ? (list_diffPVz[list_diffPVz.size() / 2 - 1] + list_diffPVz[list_diffPVz.size() / 2]) / 2 : list_diffPVz[list_diffPVz.size() / 2];
    std::sort(list_ResPVx.begin(), list_ResPVx.end());
    float median_ResPVx = (list_ResPVx.size() % 2 == 0) ? (list_ResPVx[list_ResPVx.size() / 2 - 1] + list_ResPVx[list_ResPVx.size() / 2]) / 2 : list_ResPVx[list_ResPVx.size() / 2];
    std::sort(list_ResPVy.begin(), list_ResPVy.end());
    float median_ResPVy = (list_ResPVy.size() % 2 == 0) ? (list_ResPVy[list_ResPVy.size() / 2 - 1] + list_ResPVy[list_ResPVy.size() / 2]) / 2 : list_ResPVy[list_ResPVy.size() / 2];
    std::sort(list_ResPVz.begin(), list_ResPVz.end());
    float median_ResPVz = (list_ResPVz.size() % 2 == 0) ? (list_ResPVz[list_ResPVz.size() / 2 - 1] + list_ResPVz[list_ResPVz.size() / 2]) / 2 : list_ResPVz[list_ResPVz.size() / 2];

    return std::make_tuple(median_diffPVx, se_diffPVx, median_diffPVy, se_diffPVy, median_diffPVz, se_diffPVz, median_ResPVx, se_resPVx, median_ResPVy, se_resPVy, median_ResPVz, se_resPVz);
}

void Draw2Dhisttable(TH2F *hist, const TString &XaxisName, const TString &YaxisName, const TString &ZaxisName, const TString &DrawOpt, const TString &outname)
{
    TCanvas *c = new TCanvas("c", "c", 4500, 3500);
    c->cd();
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    // ROOT.gStyle.SetPaintTextFormat('1.5f')
    // hist.SetMarkerSize(0.4)
    hist->GetXaxis()->SetTitle(XaxisName);
    hist->GetYaxis()->SetTitle(YaxisName);
    hist->GetZaxis()->SetTitle(ZaxisName);
    hist->GetXaxis()->SetTickSize(0.03);
    hist->GetYaxis()->SetTickSize(0.03);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetZaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->GetZaxis()->SetTitleOffset(1.5);
    hist->GetZaxis()->SetLabelSize(0.04);
    hist->GetZaxis()->SetRangeUser(hist->GetMinimum(), hist->GetMaximum());
    // hist.LabelsOption("v")
    hist->SetContour(10000);
    hist->Draw(DrawOpt);
    c->RedrawAxis();
    c->Draw();
    c->SaveAs(outname + ".pdf");
    c->SaveAs(outname + ".png");
    if (c)
    {
        c->Close();
        gSystem->ProcessEvents();
        delete c;
        c = nullptr;
    }
}

void PVDiffTwohalves_AsymScan()
{
    std::string plotpath = "./PV_TruthReco/AsymScan_FieldOff";
    mkdir(plotpath.c_str(), 0777);

    float gap_north = 3.5;
    float gap_atClips = 1.5;
    float gap_upper_min = gap_atClips / 2.;
    float gap_upper_max = gap_north - (gap_atClips / 2.);
    int Nbinscan_gap = 100;
    float centshift_min = -0.5;
    float centshift_max = 0.5;
    int Nbinscan_centshift = 50;
    std::vector<float> gap_upper(Nbinscan_gap + 1);
    std::vector<float> cent_shift(Nbinscan_centshift + 1);

    for (int i = 0; i <= Nbinscan_gap; ++i)
    {
        gap_upper[i] = gap_upper_min + (gap_upper_max - gap_upper_min) / Nbinscan_gap * i;
    }

    for (int i = 0; i <= Nbinscan_centshift; ++i)
    {
        cent_shift[i] = centshift_min + (centshift_max - centshift_min) / Nbinscan_centshift * i;
    }

    TH2F *hM_sigmaeff_diffPVx_asym = new TH2F("hM_sigmaeff_diffPVx_asym", "hM_sigmaeff_diffPVx_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_sigmaeff_diffPVy_asym = new TH2F("hM_sigmaeff_diffPVy_asym", "hM_sigmaeff_diffPVy_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_sigmaeff_diffPVz_asym = new TH2F("hM_sigmaeff_diffPVz_asym", "hM_sigmaeff_diffPVz_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_sigmaeff_resPVx_asym = new TH2F("hM_sigmaeff_resPVx_asym", "hM_sigmaeff_resPVx_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_sigmaeff_resPVy_asym = new TH2F("hM_sigmaeff_resPVy_asym", "hM_sigmaeff_resPVy_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_sigmaeff_resPVz_asym = new TH2F("hM_sigmaeff_resPVz_asym", "hM_sigmaeff_resPVz_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_median_diffPVx_asym = new TH2F("hM_median_diffPVx_asym", "hM_median_diffPVx_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_median_diffPVy_asym = new TH2F("hM_median_diffPVy_asym", "hM_median_diffPVy_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_median_diffPVz_asym = new TH2F("hM_median_diffPVz_asym", "hM_median_diffPVz_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_median_resPVx_asym = new TH2F("hM_median_resPVx_asym", "hM_median_resPVx_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_median_resPVy_asym = new TH2F("hM_median_resPVy_asym", "hM_median_resPVy_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);
    TH2F *hM_median_resPVz_asym = new TH2F("hM_median_resPVz_asym", "hM_median_resPVz_asym", Nbinscan_gap, &gap_upper[0], Nbinscan_centshift, &cent_shift[0]);

    for (int i = 0; i < Nbinscan_gap; ++i)
    {
        for (int j = 0; j < Nbinscan_centshift; ++j)
        {
            float gapupper = gap_upper[i];
            float centshift = cent_shift[j];
            TString gapnorth_str = TString::Format("%.1f", gap_north).ReplaceAll(".", "p");
            TString gapupper_str = TString::Format("%.2f", gapupper).ReplaceAll(".", "p");
            TString centshift_str = TString::Format("%.2f", centshift).ReplaceAll(".", "p");

            // TString fname = TString::Format("/sphenix/user/hjheng/TrackletAna/minitree/AuAu_ana325private_NoPileup_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth%s_GapUpper%s_CentShift%s_3DVertex_twohalves.root",
            //                                 gapnorth_str.Data(), gapupper_str.Data(), centshift_str.Data());
            TString fname = TString::Format("/sphenix/user/hjheng/TrackletAna/minitree/HijingAuAuMB_NoPileup_0T_RecoVtx_Optimization/TrackletAna_RecoClusters_RecoVtx_TklCluster_dPhiCutbin3_dZCutbin10_GapNorth%s_GapUpper%s_CentShift%s_3DVertex_twohalves.root",
                                            gapnorth_str.Data(), gapupper_str.Data(), centshift_str.Data());
            auto [mdn_x, se_x, mdn_y, se_y, mdn_z, se_z, mdn_resPVx, se_resPVx, mdn_resPVy, se_resPVy, mdn_resPVz, se_resPVz] = getPVzDiffSigma(fname.Data());
            std::cout << "(gap_north, gap_upper, cent_shift)=(" << gap_north << ", " << gapupper_str << ", " << centshift_str << ")" << std::endl;
            hM_sigmaeff_diffPVx_asym->SetBinContent(i + 1, j + 1, se_x);
            hM_median_diffPVx_asym->SetBinContent(i + 1, j + 1, mdn_x);
            hM_sigmaeff_diffPVy_asym->SetBinContent(i + 1, j + 1, se_y);
            hM_median_diffPVy_asym->SetBinContent(i + 1, j + 1, mdn_y);
            hM_sigmaeff_diffPVz_asym->SetBinContent(i + 1, j + 1, se_z);
            hM_median_diffPVz_asym->SetBinContent(i + 1, j + 1, mdn_z);
            hM_sigmaeff_resPVx_asym->SetBinContent(i + 1, j + 1, se_resPVx);
            hM_median_resPVx_asym->SetBinContent(i + 1, j + 1, mdn_resPVx);
            hM_sigmaeff_resPVy_asym->SetBinContent(i + 1, j + 1, se_resPVy);
            hM_median_resPVy_asym->SetBinContent(i + 1, j + 1, mdn_resPVy);
            hM_sigmaeff_resPVz_asym->SetBinContent(i + 1, j + 1, se_resPVz);
            hM_median_resPVz_asym->SetBinContent(i + 1, j + 1, mdn_resPVz);
        }
    }

    SetsPhenixStyle();
    Draw2Dhisttable(hM_sigmaeff_diffPVx_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "#DeltaX(= PV_{x}^{1} - PV_{x}^{2}) #sigma_{eff} (cm)", "colz",
                    Form("%s/PVxDiff_GapNorth%s_sigmaeff_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_median_diffPVx_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "#DeltaX(= PV_{x}^{1} - PV_{x}^{2}) median (cm)", "colz",
                    Form("%s/PVxDiff_GapNorth%s_median_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_sigmaeff_diffPVy_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "#DeltaY(= PV_{y}^{1} - PV_{y}^{2}) #sigma_{eff} (cm)", "colz",
                    Form("%s/PVyDiff_GapNorth%s_sigmaeff_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_median_diffPVy_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "#DeltaY(= PV_{y}^{1} - PV_{y}^{2}) median (cm)", "colz",
                    Form("%s/PVyDiff_GapNorth%s_median_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_sigmaeff_diffPVz_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "#DeltaZ(= PV_{z}^{1} - PV_{z}^{2}) #sigma_{eff} (cm)", "colz",
                    Form("%s/PVzDiff_GapNorth%s_sigmaeff_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_median_diffPVz_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "#DeltaZ(= PV_{z}^{1} - PV_{z}^{2}) median (cm)", "colz",
                    Form("%s/PVzDiff_GapNorth%s_median_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_sigmaeff_resPVx_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "(PV_{x}^{Reco} - PV_{x}^{Truth}) #sigma_{eff} (cm)", "colz",
                    Form("%s/PVxRes_GapNorth%s_sigmaeff_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_median_resPVx_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "(PV_{x}^{Reco} - PV_{x}^{Truth}) median (cm)", "colz",
                    Form("%s/PVxRes_GapNorth%s_median_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_sigmaeff_resPVy_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "(PV_{y}^{Reco} - PV_{y}^{Truth}) #sigma_{eff} (cm)", "colz",
                    Form("%s/PVyRes_GapNorth%s_sigmaeff_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_median_resPVy_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "(PV_{y}^{Reco} - PV_{y}^{Truth}) median (cm)", "colz",
                    Form("%s/PVyRes_GapNorth%s_median_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_sigmaeff_resPVz_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "(PV_{z}^{Reco} - PV_{z}^{Truth}) #sigma_{eff} (cm)", "colz",
                    Form("%s/PVzRes_GapNorth%s_sigmaeff_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
    Draw2Dhisttable(hM_median_resPVz_asym, "East-half gap (mm)", "Shift in MVTX center (mm)", "(PV_{z}^{Reco} - PV_{z}^{Truth}) median (cm)", "colz",
                    Form("%s/PVzRes_GapNorth%s_median_asFuncOfGapUpperCentShift", plotpath.c_str(), TString::Format("%.1f", gap_north).ReplaceAll(".", "p").Data()));
}