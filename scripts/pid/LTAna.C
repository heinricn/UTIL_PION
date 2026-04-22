/****************************
*
*   Code originally writen by N. Hamdi to produce pid plots for pionLT and KaonLT data
*
*****************************/
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TText.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TFormula.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TLine.h"
#include "TArrow.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

TH1F *randomSub(TChain *ch, TString var,TString cut,int nbins,double xmin,double xmax)
{
    static unsigned int ncall = 0;
    ncall++;

    int nWindows = 6;
    double prompt_min = -0.5;
    double prompt_max = 1.5;
    double random_min1 = -10.5;
    double random_max1 = -4.5;
    double random_min2 = 5.5;
    double random_max2 = 11.5;

    TString promptCut = Form("((CTime.eKCoinTime_ROC1 > %f) && (CTime.eKCoinTime_ROC1 < %f))",prompt_min, prompt_max);
    TString randomCut = Form("((CTime.eKCoinTime_ROC1 > %f) && (CTime.eKCoinTime_ROC1 < %f)) | ((CTime.eKCoinTime_ROC1 > %f) && (CTime.eKCoinTime_ROC1 < %f))",random_min1,random_max1,random_min2,random_max2);

    TH1F *h_mmk_prompt = new TH1F("h_mmk_prompt",Form(";%s;Counts",var.Data()),nbins,xmin,xmax);
    TH1F *h_mmk_random = new TH1F("h_mmk_random",Form(";%s;Counts",var.Data()),nbins,xmin,xmax);
    TH1F *h_mmk_sub = new TH1F(Form("h_mmk_sub_%i",ncall),Form(";%s;Counts",var.Data()),nbins,xmin,xmax);
    ch->Project("h_mmk_prompt",var.Data(),Form("%s&&%s",cut.Data(),promptCut.Data()));
    ch->Project("h_mmk_random",var.Data(),Form("%s&&%s",cut.Data(),randomCut.Data()));
    h_mmk_random->Scale(1.0/nWindows);
    h_mmk_sub->Add(h_mmk_prompt, h_mmk_random, 1, -1);

    delete h_mmk_prompt;
    delete h_mmk_random;

    return h_mmk_sub;
}

void LTAna()
{
    // Open the file containing the tree.
    TFile *fin = NULL;
    TString pathIn = "/volatile/hallc/c-pionlt/heinricn/ROOTfiles/Analysis/PionLT/";
    TString pathOut = "/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/UTIL_PION/scripts/pid";
    // TString pathIn = "/home/nacer/hallc_uofr/kaonLT/input";
    // if(dataORmc == "data") fin = TFile::Open(Form("%s/input/20220604/Kaon_coin_replay_production_6890_-1.root", fpath.Data()));

    ifstream fruns("/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/UTIL_BATCH/InputRunLists/PionLT_2021_2022/Q5p00_W2p95_mideps_center");
    if(!fruns.is_open())
    {
        cout << "Could not open run list file" << endl;
        return;
    }

    vector<int> runs;
    int run;  
    while(fruns>>run)  
        runs.push_back(run);

    TChain *ch = new TChain("T");
    for(size_t i = 0; i < runs.size(); i++)
    {
        ch->AddFile(Form("%s/PionLT_ProdCoin_replay_production_%i_-1.root", pathIn.Data(), runs[i]));
    }

    // TTree *t = (TTree *)fin->Get("T");
    if(ch == nullptr)
        cout<<"something went wrong! abort!"<<endl;
    
    ch->Print();

    // ================ CUTS ================= end
    
    //Save histograms to output ROOT file
    TFile *fout = new TFile(Form("%s/kaonLTAna.root",pathOut.Data()), "RECREATE");

    gStyle->SetPalette(kRainBow); // choose different Palette color scheme
    // gStyle->SetOptStat("e"); // remove stat box
    gStyle->SetOptStat(0); // remove stat box 
    gStyle->SetOptFit(1);
    gStyle->SetTitleSize(.05, "xyz");
    gStyle->SetTitleOffset(1.4, "xy");
    gStyle->SetLabelSize(.05, "xyz");
    gStyle->SetPadTopMargin(0.12);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.17); // 2d:0.17 , 1d:0.05
    gStyle->SetPadLeftMargin(0.16);  // 2d:0.16, 1d:0.15
    gStyle->SetLineWidth(2);
    // gStyle->SetTextSize(0.03);
    // gStyle->SetPadGridX(true);
    // gStyle->SetPadGridY(true);

    // ================ CUTS ================= start

    // ===== Acceptance ===== start

    double H_del_low = -8.0;
    double H_del_high = 8.0;
    double P_del_low = -10.0;
    double P_del_high = 20.0;
    double H_xpfp_low = -0.08;
    double H_xpfp_high = 0.08;
    double H_ypfp_low = -0.045;
    double H_ypfp_high = 0.045;
    double P_xpfp_low = -0.06;
    double P_xpfp_high = 0.06;
    double P_ypfp_low = -0.04;
    double P_ypfp_high = 0.04;

    TString deltaCut = Form("((H.gtr.dp > %f) && (H.gtr.dp < %f)) && ((P.gtr.dp > %f) && (P.gtr.dp < %f))",H_del_low,H_del_high,P_del_low,P_del_high);

    TString h_pfpCut = Form("((H.gtr.th > %f) && (H.gtr.th < %f)) && ((H.gtr.ph > %f) && (H.gtr.ph < %f))",H_xpfp_low,H_xpfp_high,H_ypfp_low,H_ypfp_high);

    TString p_pfpCut = Form("((P.gtr.th > %f) && (P.gtr.th < %f)) && ((P.gtr.ph > %f) && (P.gtr.ph < %f))",P_xpfp_low,P_xpfp_high,P_ypfp_low,P_ypfp_high);

    TString acceptCuts = Form("%s && %s && %s",deltaCut.Data(),h_pfpCut.Data(),p_pfpCut.Data());

    // ===== Acceptance ===== end

    // ===== PID ===== Start

    double P_kcut_P_beta = 0.3;
    double P_kcut_P_hgcer = 1.5;
    double P_kcut_P_aero = 3;
    double H_ecut_H_cal = 0.5;
    double P_kcut_H_cer = 2.0;
    double KRFHigh = 0.95;
    double KRFLow = 0.15;

    // --- corrected RF time
    TString RFtimeDist_corr = "RFTime.SHMS_RFtimeDist+(3.05*P.gtr.th)";
    TString pidCuts = Form("(H.cal.etottracknorm > %f) && (H.cer.npeSum > %f) && ((%s > %f) | (%s < %f))",H_ecut_H_cal,P_kcut_H_cer,RFtimeDist_corr.Data(),KRFHigh,RFtimeDist_corr.Data(),KRFLow);

    double prompt_min = -0.5;
    double prompt_max = 1.5;
    TString promptCut = Form("((CTime.eKCoinTime_ROC1 > %f) && (CTime.eKCoinTime_ROC1 < %f))",prompt_min, prompt_max);

    TString totalCuts = Form("%s && %s",acceptCuts.Data(),pidCuts.Data());

    // ===== PID ===== end

    // ================ Histograms ================

    // ================= H.cal.etottracknorm
    auto c_HcalEtottracknorm = new TCanvas("c_HcalEtottracknorm","c_HcalEtottracknorm",900,600);
    c_HcalEtottracknorm->cd();
    c_HcalEtottracknorm->SetLogy();
    auto h_HcalEtottracknorm = new TH1F("h_HcalEtottracknorm",";H.cal.etottracknorm;Counts", 300,  0.0, 2.0);
    ch->Project("h_HcalEtottracknorm", "H.cal.etottracknorm", Form("%s && (H.cer.npeSum > %f) && ((%s > %f) | (%s < %f)) && %s",acceptCuts.Data(),P_kcut_H_cer,RFtimeDist_corr.Data(),KRFHigh,RFtimeDist_corr.Data(),KRFLow,promptCut.Data()));
    h_HcalEtottracknorm->Draw();
    h_HcalEtottracknorm->Write();
    auto l1_HcalEtottracknorm = new TLine(H_ecut_H_cal,0,H_ecut_H_cal,h_HcalEtottracknorm->GetMaximum());
    l1_HcalEtottracknorm->SetLineColor(2);
    l1_HcalEtottracknorm->Draw("same");
    auto ar_HcalEtottracknorm = new TArrow(H_ecut_H_cal,h_HcalEtottracknorm->GetMaximum(),H_ecut_H_cal+0.1,h_HcalEtottracknorm->GetMaximum(),0.03,"|>");
    ar_HcalEtottracknorm->SetFillColor(2);
    ar_HcalEtottracknorm->Draw("same");
    c_HcalEtottracknorm->Print(Form("%s/c_HcalEtottracknorm.root",pathOut.Data()),"root");
    c_HcalEtottracknorm->Print(Form("%s/c_HcalEtottracknorm.pdf",pathOut.Data()),"pdf");

    // ================= H.cer.npeSum

    auto c_HcerNpeSum = new TCanvas("c_HcerNpeSum",";H.cer.npeSum;Counts",900,600);
    c_HcerNpeSum->cd();
    c_HcerNpeSum->SetLogy();
    auto h_HcerNpeSum = new TH1F("h_HcerNpeSum",";H.cer.npeSum;Counts", 300, 0.1, 10.0);
    ch->Project("h_HcerNpeSum", "H.cer.npeSum", Form("%s && (H.cal.etottracknorm > %f) && (H.cer.npeSum > %f) && ((%s > %f) | (%s < %f)) && %s",acceptCuts.Data(),H_ecut_H_cal,P_kcut_H_cer,RFtimeDist_corr.Data(),KRFHigh,RFtimeDist_corr.Data(),KRFLow,promptCut.Data()));
    h_HcerNpeSum->Draw();
    h_HcerNpeSum->Write();
    auto l1_HcerNpeSum = TLine(P_kcut_H_cer,0, P_kcut_H_cer, h_HcerNpeSum->GetMaximum());
    l1_HcerNpeSum.SetLineColor(2);
    l1_HcerNpeSum.Draw("same");
    auto ar_HcerNpeSum = TArrow(P_kcut_H_cer,h_HcerNpeSum->GetMaximum(),P_kcut_H_cer+1.5, h_HcerNpeSum->GetMaximum(),0.03,"|>");
    ar_HcerNpeSum.SetFillColor(2);
    ar_HcerNpeSum.Draw("same");
    c_HcerNpeSum->Print(Form("%s/c_HcerNpeSum.root",pathOut.Data()),"root");
    c_HcerNpeSum->Print(Form("%s/c_HcerNpeSum.pdf",pathOut.Data()),"pdf");

    // ================= Missing Mass ================

    auto h_mmk_sub = randomSub(ch,"P.kin.secondary.MMK",totalCuts,300,0.85,1.4);

    TString pionCut = Form("%s && (H.cal.etottracknorm > %f) && (H.cer.npeSum > %f) && ((%s > 0.3) && (%s < 0.8))",acceptCuts.Data(),H_ecut_H_cal,P_kcut_H_cer,RFtimeDist_corr.Data(),RFtimeDist_corr.Data());

    auto h_mmk_sub_pion = randomSub(ch,"P.kin.secondary.MMK",pionCut,300,0.85,1.4);
    auto c_mmk = new TCanvas("c_mmk", "c_mmk", 600, 400);
    c_mmk->cd();
    h_mmk_sub->SetLineColor(1);
    h_mmk_sub_pion->SetLineColor(2);
    h_mmk_sub->Draw();
    h_mmk_sub_pion->Draw("same");
    h_mmk_sub->Write("h_mmk_sub",TObject::kWriteDelete);
    c_mmk->Print(Form("%s/c_mmk.root",pathOut.Data()),"root");
    c_mmk->Print(Form("%s/c_mmk.pdf",pathOut.Data()),"pdf");

    // ================= Cut Efficiency ================ start

    // ===== RF
    TString cut_noRF = Form("%s && (H.cal.etottracknorm > %f) && (H.cer.npeSum > %f) && (P.kin.secondary.MMK < 0.95) && ((%s > 0.0) && (%s < 2.0))",acceptCuts.Data(),H_ecut_H_cal,P_kcut_H_cer,RFtimeDist_corr.Data(),RFtimeDist_corr.Data());
    TString cut_RF = Form("%s && (H.cal.etottracknorm > %f) && (H.cer.npeSum > %f) && (P.kin.secondary.MMK < 0.95) && ((%s > 0.0) && (%s < 1.0))",acceptCuts.Data(),H_ecut_H_cal,P_kcut_H_cer,RFtimeDist_corr.Data(),RFtimeDist_corr.Data());

    auto h_mmk_noRF = randomSub(ch,"P.kin.secondary.MMK",cut_noRF,300,0.85,0.95);
    auto h_mmk_RF = randomSub(ch,"P.kin.secondary.MMK",cut_RF,300,0.85,0.95);
    
    int xbin_min = 1;//h_mmk->FindFixBin(1.08);
    int xbin_max = 300;//h_mmk->FindFixBin(1.15);
    double dN_noRF = 0.0;
    double dN_RF = 0.0;
    double N_noRF = h_mmk_noRF->IntegralAndError(xbin_min,xbin_max,dN_noRF);
    double N_RF = h_mmk_RF->IntegralAndError(xbin_min,xbin_max,dN_RF);
    double eff_noRF = N_RF/N_noRF;

    auto c_mmk_eff_rf = new TCanvas("c_mmk_eff_rf", "c_mmk_eff_rf", 600, 400);
    c_mmk_eff_rf->cd();
    h_mmk_noRF->SetLineColor(1);
    h_mmk_RF->SetLineColor(6);
    h_mmk_noRF->Draw("hist");
    h_mmk_RF->Draw("hist same");
    h_mmk_noRF->Write("h_mmk_noRF",TObject::kWriteDelete);
    h_mmk_RF->Write("h_mmk_RF",TObject::kWriteDelete);
    auto lg_eff_rf = TLegend(0.55,0.6,0.85,0.8);
    lg_eff_rf.AddEntry(h_mmk_noRF,"0 - 2 ns","l");
    lg_eff_rf.AddEntry(h_mmk_RF,"0 - 1 ns","l");
    lg_eff_rf.Draw("same");
    auto lat_eff_rf = TLatex();
    lat_eff_rf.SetTextSize(0.05);
    lat_eff_rf.SetTextAlign(13);
    lat_eff_rf.SetNDC();
    lat_eff_rf.DrawLatex(0.2,0.85,Form("#color[6]{efficiency = %0.3f}",eff_noRF));
    c_mmk_eff_rf->Print(Form("%s/c_mmk_eff_rf.root",pathOut.Data()),"root");
    c_mmk_eff_rf->Print(Form("%s/c_mmk_eff_rf.pdf",pathOut.Data()),"pdf");

    // printf("N_RF = %0.3f, N_noRF = %0.3f\n",N_RF,N_noRF);

    // ===== HCal

    auto h_HcalEtottracknorm_d = new TH1F("h_HcalEtottracknorm_d",";H.cal.etottracknorm;Counts",300,0.0,2.0);
    auto h_HcalEtottracknorm_s = new TH1F("h_HcalEtottracknorm_s",";H.cal.etottracknorm;Counts",300,0.0,2.0);
    ch->Project("h_HcalEtottracknorm_d", "H.cal.etottracknorm", Form("%s && (H.cal.etottracknorm > %f) && (H.cer.npeSum > %f)",acceptCuts.Data(),H_ecut_H_cal,P_kcut_H_cer));
    ch->Project("h_HcalEtottracknorm_s", "H.cal.etottracknorm", Form("%s && (H.cer.npeSum > %f)",acceptCuts.Data(),P_kcut_H_cer));

    double dNs_Hcal = 0.0;
    double dNd_Hcal = 0.0;
    double Nd_Hcal = h_HcalEtottracknorm_d->IntegralAndError(xbin_min,xbin_max,dNd_Hcal);
    double Ns_Hcal = h_HcalEtottracknorm_s->IntegralAndError(xbin_min,xbin_max,dNs_Hcal);
    double eff_Hcal = Nd_Hcal/Ns_Hcal;

    auto c_HcalEtottracknorm_eff = new TCanvas("c_HcalEtottracknorm_eff","c_HcalEtottracknorm_eff",900,600);
    c_HcalEtottracknorm_eff->cd();
    h_HcalEtottracknorm_d->SetLineColor(2);
    h_HcalEtottracknorm_s->SetLineColor(1);
    c_HcalEtottracknorm_eff->SetLogy();
    h_HcalEtottracknorm_s->Draw();
    h_HcalEtottracknorm_d->Draw("same");
    auto leg_HcalEtottracknorm_eff = TLegend(0.55,0.6,0.75,0.8);
    leg_HcalEtottracknorm_eff.AddEntry(h_HcalEtottracknorm_s,Form("Hcer > %0.1f + Hcal > %0.1f",P_kcut_H_cer,H_ecut_H_cal),"l");
    leg_HcalEtottracknorm_eff.AddEntry(h_HcalEtottracknorm_d,Form("Hcer > %0.1f",P_kcut_H_cer),"l");
    leg_HcalEtottracknorm_eff.Draw("same");
    auto lat_HcalEtottracknorm_eff = TLatex();
    lat_HcalEtottracknorm_eff.SetTextSize(0.05);
    lat_HcalEtottracknorm_eff.SetTextAlign(13);
    lat_HcalEtottracknorm_eff.SetNDC();
    lat_HcalEtottracknorm_eff.DrawLatex(0.2,0.85,Form("#color[6]{efficiency = %0.3f}",eff_Hcal));
    c_HcalEtottracknorm_eff->Print(Form("%s/c_HcalEtottracknorm_eff.root",pathOut.Data()),"root");
    c_HcalEtottracknorm_eff->Print(Form("%s/c_HcalEtottracknorm_eff.pdf",pathOut.Data()),"pdf");

    // ===== HCer

    auto h_HcerNpeSum_d = new TH1F("h_HcerNpeSum_d",";H.cer.npeSum;Counts", 300, 0.1, 10.0);
    auto h_HcerNpeSum_s = new TH1F("h_HcerNpeSum_s",";H.cer.npeSum;Counts", 300, 0.1, 10.0);
    ch->Project("h_HcerNpeSum_d", "H.cer.npeSum", Form("%s && (H.cal.etottracknorm > %f) && (H.cer.npeSum > %f)",acceptCuts.Data(),H_ecut_H_cal,P_kcut_H_cer));
    ch->Project("h_HcerNpeSum_s", "H.cer.npeSum", Form("%s && (H.cal.etottracknorm > %f)",acceptCuts.Data(),H_ecut_H_cal));

    double dNs_Hcer = 0.0;
    double dNd_Hcer = 0.0;
    double Nd_Hcer = h_HcerNpeSum_d->IntegralAndError(xbin_min,xbin_max,dNd_Hcer);
    double Ns_Hcer = h_HcerNpeSum_s->IntegralAndError(xbin_min,xbin_max,dNs_Hcer);
    double eff_Hcer = Nd_Hcer/Ns_Hcer;

    auto c_HcerNpeSum_eff = new TCanvas("c_HcerNpeSum_eff","c_HcerNpeSum_eff",900,600);
    c_HcerNpeSum_eff->cd();
    h_HcerNpeSum_d->SetLineColor(2);
    h_HcerNpeSum_s->SetLineColor(1);
    c_HcerNpeSum_eff->SetLogy();
    h_HcerNpeSum_s->Draw();
    h_HcerNpeSum_d->Draw("same");
    auto leg_HcerNpeSum_eff = TLegend(0.55,0.6,0.75,0.8);
    leg_HcerNpeSum_eff.AddEntry(h_HcerNpeSum_s,Form("Hcer > %0.1f + Hcal > %0.1f",P_kcut_H_cer,H_ecut_H_cal),"l");
    leg_HcerNpeSum_eff.AddEntry(h_HcerNpeSum_d,Form("Hcal > %0.1f",H_ecut_H_cal),"l");
    leg_HcerNpeSum_eff.Draw("same");
    auto lat_HcerNpeSum_eff = TLatex();
    lat_HcerNpeSum_eff.SetTextSize(0.05);
    lat_HcerNpeSum_eff.SetTextAlign(13);
    lat_HcerNpeSum_eff.SetNDC();
    lat_HcerNpeSum_eff.DrawLatex(0.2,0.85,Form("#color[6]{efficiency = %0.3f}",eff_Hcer));
    c_HcerNpeSum_eff->Print(Form("%s/c_HcerNpeSum_eff.root",pathOut.Data()),"root");
    c_HcerNpeSum_eff->Print(Form("%s/c_HcerNpeSum_eff.pdf",pathOut.Data()),"pdf");

    // ================= Cut Efficiency ================ end

}