### 27/02/2025 - Abdennacer Hamdi, University of Regina
### This code is used to study the event selection cuts for the kaonLT experiment
################################################################################################################################################

# Import relevant packages
import sys, os, subprocess
import array
import csv
from matplotlib.pylab import f
import uproot as up
import numpy as np
import math as m
import root_numpy as rnp
import pandas as pd
import ROOT as r
from ROOT import TCanvas,TH1F,TH2D,RDataFrame,TChain,TFile,TString,TLegend,TStyle,TText,TLine,TArrow,TF1,TGraph,TGraphErrors,TFitResultPtr,TLatex,gStyle,TCutG
from ctypes import c_double
import uncertainties as u

np.bool = bool
np.float = float

gStyle.SetLineWidth(2)
### Set ROOT to batch mode explicitly, does not splash anything to screen
r.gROOT.SetBatch(r.kTRUE) 
# Removes stat box
r.gStyle.SetOptStat(0)

r.EnableImplicitMT() ### Enable ROOT's implicit multi-threading

dataSet = sys.argv[1]

SIMCPATH = "/volatile/hallc/c-pionlt/heinricn/OUTPUT/Analysis/SIMC"
RUNLISTPATH = "/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/UTIL_BATCH/InputRunLists/PionLT_2021_2022/"
CSVPATH = "/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/UTIL_PION/efficiencies/"
REPLAYPATH = "/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/ROOTfiles/Analysis/PionLT"
OUTPATH = "/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/UTIL_PION/scripts/pid/OUTPUT"
DATA_RUN_LIST = f'{dataSet}'
DUMMY_RUN_LIST = f'{dataSet}_dummy'
CSV_FILE = "PionLT_coin_production_Prod_efficiency_data_2025_11_21"
SIMC_Suffix_lambda = f'Prod_Coin_Q0p5W2p40{dataSet}_lambda'
SIMC_Suffix_sigma = f'Prod_Coin_Q0p5W2p40{dataSet}_sigma'
data_run_list = "%s/%s" % (RUNLISTPATH, DATA_RUN_LIST)
dummy_run_list = "%s/%s" % (RUNLISTPATH, DUMMY_RUN_LIST)
csv_file = "%s/%s.csv" % (CSVPATH, CSV_FILE)

print("Running over replay path assumed as %s" % (REPLAYPATH))

# Read CSV File to use later for total charge in mC calculation
# Input runlists and csv files
data_run_list_file = (data_run_list)
dummy_run_list_file = (dummy_run_list)
csv_file_name = (csv_file)

# Read run numbers from the run list file
with open(data_run_list_file, 'r') as run_list_file_1:
    data_runs = [line.strip() for line in run_list_file_1 if line.strip()]
with open(dummy_run_list_file, 'r') as run_list_file_2:
    dummy_runs = [line.strip() for line in run_list_file_2 if line.strip()]

# Read CSV file using pandas
df = pd.read_csv(csv_file_name)

# Filter DataFrame to include only rows with run numbers in the run list
filtered_data_df = df[df['Run_Number'].astype(str).str.replace('.0', '', regex=False).isin(data_runs)]
filtered_dummy_df = df[df['Run_Number'].astype(str).str.replace('.0', '', regex=False).isin(dummy_runs)]

print ('\nBeam Energy = ',filtered_data_df['Beam_Energy'].values[0], '\n')
print (f'\nData Setting : {dataSet} \n')
print("-"*40)

### ====== add all the ROOT trees to the TChain for Data    
ch_data =  TChain("T")
for run in data_runs:
    ch_data.AddFile(f'{REPLAYPATH}/PionLT_ProdCoin_replay_production_{run}_-1.root')
    print(f'{run}')

if ch_data is None:
    print("something went wrong! abort!")

### Create ROOT dataframe from a TChain
rdf = RDataFrame(ch_data)

### save the output histograms to ROOT file
output_file = TFile(f'{OUTPATH}/{dataSet}/evtSelect_{dataSet}.root', 'RECREATE')

### ======== koan Cuts ========= start
### ===== Acceptance ===== start
H_del_low = -8.0
H_del_high = 8.0
H_xpfp_low = -0.08
H_xpfp_high = 0.08
H_ypfp_low = -0.045
H_ypfp_high = 0.045
P_del_low = -10.0
P_del_high = 20.0
P_xpfp_low = -0.06
P_xpfp_high = 0.06
P_ypfp_low = -0.04
P_ypfp_high = 0.04

### === Data
HMS_Acceptance_data = f'(H.gtr.dp>={H_del_low}) && (H.gtr.dp<={H_del_high}) && (H.gtr.th>={H_xpfp_low}) && (H.gtr.th<={H_xpfp_high}) && (H.gtr.ph>={H_ypfp_low}) && (H.gtr.ph<={H_ypfp_high})'
SHMS_Acceptance_data = f'(P.gtr.dp>={P_del_low}) && (P.gtr.dp<={P_del_high}) && (P.gtr.th>={P_xpfp_low}) && (P.gtr.th<={P_xpfp_high}) && (P.gtr.ph>={P_ypfp_low}) && (P.gtr.ph<={P_ypfp_high})'
acceptCuts_data = f'({HMS_Acceptance_data} && {SHMS_Acceptance_data})'

### ===== Acceptance ===== end

### ===== PID ===== Start
P_kcut_P_beta = 0.3
P_kcut_P_hgcer = 1.5
P_kcut_P_aero = 1.5
H_ecut_H_cal = 0.7
P_kcut_H_cer = 1.5
KRFHigh = 1.5
KRFLow = 0.1

### === corrected RF time
RFtimeDist_corr = f'RFTime.SHMS_RFtimeDist+(3.05*P.gtr.th)'

pidCuts = f'(H.cal.etottracknorm > {H_ecut_H_cal}) && (H.cer.npeSum > {P_kcut_H_cer}) && (({RFtimeDist_corr} < {KRFHigh}) && ({RFtimeDist_corr} > {KRFLow}))'

### ===== PID ===== end

### === total cut
cuts = f'(({acceptCuts_data}) && ({pidCuts}))'

### ========= random events subtraction time window for kaons
nWindows = 6
prompt_min = -2.004
prompt_max = 2.004
random_min1 = -18.036
random_max1 = -6.012
random_min2 = 6.012
random_max2 = 18.036

promptCut = f'(CTime.ePiCoinTime_ROC1 > {prompt_min}) && (CTime.ePiCoinTime_ROC1 < {prompt_max})'
randomCut = f'((CTime.ePiCoinTime_ROC1 > {random_min1}) && (CTime.ePiCoinTime_ROC1 < {random_max1})) | ((CTime.ePiCoinTime_ROC1 > {random_min2}) && (CTime.ePiCoinTime_ROC1 < {random_max2}))'

### ================ pion cuts to subtract from kaon missing mass later
## -- pion pid cuts
pion_cut = f'({acceptCuts_data} && H.cal.etottracknorm > {H_ecut_H_cal}) && (({RFtimeDist_corr} < {KRFHigh}) && ({RFtimeDist_corr} > {KRFLow}))'

## --- Coin Time for pions
prompt_pi_min = -0.5
prompt_pi_max = 1.5
random_pi_min1 = -12.5
random_pi_max1 = -6.5
random_pi_min2 = 7.5
random_pi_max2 = 13.5
# prompt_pi_min = -1.0
# prompt_pi_max = 1.0
# random_pi_min1 = -13.0
# random_pi_max1 = -7.0
# random_pi_min2 = 7.0
# random_pi_max2 = 13.0

promptCut_pi = f'(CTime.ePiCoinTime_ROC1 > {prompt_pi_min}) && (CTime.ePiCoinTime_ROC1 < {prompt_pi_max})'
randomCut_pi = f'((CTime.ePiCoinTime_ROC1 > {random_pi_min1}) && (CTime.ePiCoinTime_ROC1 < {random_pi_max1})) | ((CTime.ePiCoinTime_ROC1 > {random_pi_min2}) && (CTime.ePiCoinTime_ROC1 < {random_pi_max2}))'

### ============== Event Selection ================= start

## --- electron selection using electron total track normalized in HMS Calorimeter
cuts_noHcal = f'(({acceptCuts_data}) && H.cer.npeSum > {P_kcut_H_cer}) && (({RFtimeDist_corr} < {KRFHigh}) && ({RFtimeDist_corr} > {KRFLow})) && {promptCut}'
c_HcalEtottracknorm = TCanvas("c_HcalEtottracknorm","c_HcalEtottracknorm",900,600)
c_HcalEtottracknorm.cd()
c_HcalEtottracknorm.SetLogy()
h_HcalEtottracknorm = rdf.Filter(cuts_noHcal).Histo1D(("h_HcalEtottracknorm", ";H.cal.etottracknorm;Counts", 300, 0.0, 2.0), "H.cal.etottracknorm")
h_HcalEtottracknorm.Draw()
l1_HcalEtottracknorm = TLine(H_ecut_H_cal, 0, H_ecut_H_cal, h_HcalEtottracknorm.GetMaximum())
l1_HcalEtottracknorm.SetLineColor(2)
l1_HcalEtottracknorm.Draw("same")
ar_HcalEtottracknorm = TArrow(H_ecut_H_cal,h_HcalEtottracknorm.GetMaximum(),H_ecut_H_cal+0.1,h_HcalEtottracknorm.GetMaximum(),0.03,"|>")
ar_HcalEtottracknorm.SetFillColor(2)
ar_HcalEtottracknorm.Draw("same")
h_HcalEtottracknorm.Write("h_HcalEtottracknorm",r.TObject.kWriteDelete)
c_HcalEtottracknorm.SaveAs(f'{OUTPATH}/{dataSet}/c_HcalEtottracknorm_{dataSet}.root')
c_HcalEtottracknorm.SaveAs(f'{OUTPATH}/{dataSet}/c_HcalEtottracknorm_{dataSet}.pdf')

## --- electron selection using electron aerogel cherenkov
cuts_noHcer = f'(({acceptCuts_data}) && H.cal.etottracknorm > {H_ecut_H_cal}) && (({RFtimeDist_corr} < {KRFHigh}) && ({RFtimeDist_corr} > {KRFLow})) && {promptCut}'
c_HcerNpeSum = TCanvas("c_HcerNpeSum","c_HcerNpeSum",900,600)
c_HcerNpeSum.cd()
c_HcerNpeSum.SetLogy()
h_HcerNpeSum = rdf.Filter(cuts_noHcer).Histo1D(("h_HcerNpeSum", ";H.cer.npeSum;Counts", 300, 0.1, 10.0), "H.cer.npeSum")
h_HcerNpeSum.Draw()
l1_HcerNpeSum = TLine(P_kcut_H_cer, 0, P_kcut_H_cer, h_HcerNpeSum.GetMaximum())
l1_HcerNpeSum.SetLineColor(2)
l1_HcerNpeSum.Draw("same")
ar_HcerNpeSum = TArrow(P_kcut_H_cer,h_HcerNpeSum.GetMaximum(),P_kcut_H_cer+0.1,h_HcerNpeSum.GetMaximum(),0.03,"|>")
ar_HcerNpeSum.SetFillColor(2)
ar_HcerNpeSum.Draw("same")
h_HcerNpeSum.Write("h_HcerNpeSum",r.TObject.kWriteDelete)
c_HcerNpeSum.SaveAs(f'{OUTPATH}/{dataSet}/c_HcerNpeSum_{dataSet}.root')
c_HcerNpeSum.SaveAs(f'{OUTPATH}/{dataSet}/c_HcerNpeSum_{dataSet}.pdf')

## --- RF time cut to select kaons
## RFTime.SHMS_RFtimeDist = (RFTime − StartTime + RF_Offset)%Bunch_Spacing
cuts_noRF = f'(({acceptCuts_data}) && (H.cal.etottracknorm > {H_ecut_H_cal}) && (H.cer.npeSum > {P_kcut_H_cer}) && {promptCut})'
c_RFtimeDistVsMmk = TCanvas("c_RFtimeDistVsMmk","c_RFtimeDistVsMmk",900,600)
c_RFtimeDistVsMmk.cd()
h_RFtimeDistVsMmk = rdf.Define("RFtimeDist_corr", f'{RFtimeDist_corr}').Filter(cuts_noRF).Histo2D(("h_RFtimeDistVsMmk", ";P.kin.secondary.MMpi;RFTime.SHMS_RFtimeDist_corr", 600, 0.85, 1.42, 600, 0.0, 4.008), "P.kin.secondary.MMpi","RFtimeDist_corr")
h_RFtimeDistVsMmk.Draw("colz")
l_RFLow = TLine(0.85,KRFLow,1.42,KRFLow)
l_RFHigh = TLine(0.85,KRFHigh,1.42,KRFHigh)
l_RFLow.SetLineColor(2)
l_RFHigh.SetLineColor(2)
l_RFLow.Draw("same")
l_RFHigh.Draw("same")
ar_RFLow = TArrow(1.35,KRFLow,1.35,KRFLow-0.1,0.02,"|>")
ar_RFHigh = TArrow(1.35,KRFHigh,1.35,KRFHigh+0.1,0.02,"|>")
ar_RFLow.SetFillColor(2)
ar_RFHigh.SetFillColor(2)
ar_RFLow.Draw("same")
ar_RFHigh.Draw("same")
h_RFtimeDistVsMmk.Write("h_RFtimeDistVsMmk",r.TObject.kWriteDelete)
c_RFtimeDistVsMmk.SaveAs(f'{OUTPATH}/{dataSet}/c_RFtimeDistVsMmk_{dataSet}.root')
c_RFtimeDistVsMmk.SaveAs(f'{OUTPATH}/{dataSet}/c_RFtimeDistVsMmk_{dataSet}.pdf')

## --- ekCoinTime Vs. mmpi
c_ekCoinTimeVsmmpi = TCanvas("c_ekCoinTimeVsmmpi","c_ekCoinTimeVsmmpi", 900, 600)
c_ekCoinTimeVsmmpi.cd()
h_ekCoinTimeVsmmpi = rdf.Filter(cuts).Histo2D(("h_ekCoinTimeVsmmpi", ";P.kin.secondary.MMpi;CTime.ePiCoinTime_ROC1", 600, 0.85, 1.42, 600, -20, 20), "P.kin.secondary.MMpi","CTime.ePiCoinTime_ROC1")
h_ekCoinTimeVsmmpi.Draw("colz")
h_ekCoinTimeVsmmpi.Write("h_ekCoinTimeVsmmpi",r.TObject.kWriteDelete)
c_ekCoinTimeVsmmpi.SaveAs(f'{OUTPATH}/{dataSet}/c_ekCoinTimeVsmmpi_{dataSet}.root')
c_ekCoinTimeVsmmpi.SaveAs(f'{OUTPATH}/{dataSet}/c_ekCoinTimeVsmmpi_{dataSet}.pdf')

## --- Coin Time for lambda region
cuts_lambda = f'{cuts} && (P.kin.secondary.MMpi>1.10) && (P.kin.secondary.MMpi<1.14)'
c_ekCoinTime_lambda = TCanvas("c_ekCoinTime_lambda","c_ekCoinTime_lambda",900,600)
c_ekCoinTime_lambda.cd()
h_ekCoinTime_lambda = rdf.Filter(cuts_lambda).Histo1D(("h_CoinTime_lambda", ";CTime.ePiCoinTime_ROC1;Counts", 600, -20, 20), "CTime.ePiCoinTime_ROC1")
h_ekCoinTime_lambda.SetLineColor(1)
h_ekCoinTime_lambda.SetMinimum(0.0)
h_ekCoinTime_lambda.Draw()
l1_prompt_min = TLine(prompt_min, 0, prompt_min, h_ekCoinTime_lambda.GetMaximum())
l1_prompt_max = TLine(prompt_max, 0, prompt_max, h_ekCoinTime_lambda.GetMaximum())
l1_random_min1 = TLine(random_min1, 0, random_min1, h_ekCoinTime_lambda.GetMaximum())
l1_random_max1 = TLine(random_max1, 0, random_max1, h_ekCoinTime_lambda.GetMaximum())
l1_random_min2 = TLine(random_min2, 0, random_min2, h_ekCoinTime_lambda.GetMaximum())
l1_random_max2 = TLine(random_max2, 0, random_max2, h_ekCoinTime_lambda.GetMaximum())
l1_prompt_min.SetLineColor(2)
l1_prompt_max.SetLineColor(2)
l1_random_min1.SetLineColor(4)
l1_random_max1.SetLineColor(4)
l1_random_min2.SetLineColor(4)
l1_random_max2.SetLineColor(4)
l1_prompt_min.Draw("same")
l1_prompt_max.Draw("same")
l1_random_min1.Draw("same")
l1_random_max1.Draw("same")
l1_random_min2.Draw("same")
l1_random_max2.Draw("same")
h_ekCoinTime_lambda.Write("h_CoinTime_lambda",r.TObject.kWriteDelete)
c_ekCoinTime_lambda.SaveAs(f'{OUTPATH}/{dataSet}/c_ekCoinTime_lambda_{dataSet}.root')
c_ekCoinTime_lambda.SaveAs(f'{OUTPATH}/{dataSet}/c_ekCoinTime_lambda_{dataSet}.pdf')

## --- pions selection for later subtraction 

##----------- random events subtraction time window for pions

c_ePiCoinTime = TCanvas("c_ePiCoinTime","c_ePiCoinTime", 900, 600)
c_ePiCoinTime.cd()

# ePiCTime_offset = 0.0
# if dataSet == 'right_highe' or dataSet == 'center_highe' or dataSet == 'left_highe':
#     ePiCTime_offset = -0.5
# h_ePiCoinTime = rdf.Define("ePiCTime_shift", f'CTime.ePiCoinTime_ROC1+{ePiCTime_offset}').Filter(pion_cut).Histo1D(("h_ePiCoinTime", ";CTime.ePiCoinTime_ROC1;Counts", 600, -20, 20), "ePiCTime_shift")

h_ePiCoinTime = rdf.Filter(pion_cut).Histo1D(("h_ePiCoinTime", ";CTime.ePiCoinTime_ROC1;Counts", 600, -20, 20), "CTime.ePiCoinTime_ROC1")
h_ePiCoinTime.SetLineColor(1)
h_ePiCoinTime.SetMinimum(0.0)
h_ePiCoinTime.Draw()
l1_prompt_pi_min = TLine(prompt_pi_min, 0, prompt_pi_min, h_ePiCoinTime.GetMaximum())
l1_prompt_pi_max = TLine(prompt_pi_max, 0, prompt_pi_max, h_ePiCoinTime.GetMaximum())
l1_random_pi_min1 = TLine(random_pi_min1, 0, random_pi_min1, h_ePiCoinTime.GetMaximum())
l1_random_pi_max1 = TLine(random_pi_max1, 0, random_pi_max1, h_ePiCoinTime.GetMaximum())
l1_random_pi_min2 = TLine(random_pi_min2, 0, random_pi_min2, h_ePiCoinTime.GetMaximum())
l1_random_pi_max2 = TLine(random_pi_max2, 0, random_pi_max2, h_ePiCoinTime.GetMaximum())
l1_prompt_pi_min.SetLineColor(2)
l1_prompt_pi_max.SetLineColor(2)
l1_random_pi_min1.SetLineColor(4)
l1_random_pi_max1.SetLineColor(4)
l1_random_pi_min2.SetLineColor(4)
l1_random_pi_max2.SetLineColor(4)
l1_prompt_pi_min.Draw("same")
l1_prompt_pi_max.Draw("same")
l1_random_pi_min1.Draw("same")
l1_random_pi_max1.Draw("same")
l1_random_pi_min2.Draw("same")
l1_random_pi_max2.Draw("same")
h_ePiCoinTime.Write("h_ePiCoinTime",r.TObject.kWriteDelete)
c_ePiCoinTime.SaveAs(f'{OUTPATH}/{dataSet}/c_ePiCoinTime_{dataSet}.root')
c_ePiCoinTime.SaveAs(f'{OUTPATH}/{dataSet}/c_ePiCoinTime_{dataSet}.pdf')

### ----- pion missing mass
c_mmpi = TCanvas("c_mmpi", "c_mmpi", 600, 400)
c_mmpi.cd()
h_mmpi_precut = rdf.Histo1D(("h_mmpi_precut", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi = rdf.Filter(f'{pion_cut} && {promptCut_pi}').Histo1D(("h_mmpi", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_random = rdf.Filter(f'{pion_cut} && {randomCut_pi}').Histo1D(("h_mmpi_random", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_random.Scale(1.0/nWindows)
h_mmpi.Add(h_mmpi_random.GetPtr(), -1)
h_mmpi_precut.SetLineColor(r.kBlack)
h_mmpi.SetLineColor(r.kRed)
lg_mmpi = TLegend(0.4,0.7,0.6,0.9)
lg_mmpi.AddEntry(h_mmpi_precut.GetPtr(),"pre-cuts","l")
lg_mmpi.AddEntry(h_mmpi.GetPtr(),"post-cuts","l")
h_mmpi_precut.SetMinimum(0.0)
h_mmpi_precut.Draw("hist")
h_mmpi.Draw("hist same")
lg_mmpi.Draw()
h_mmpi.Write("h_mmpi_prompt",r.TObject.kWriteDelete)
c_mmpi.SaveAs(f'{OUTPATH}/{dataSet}/c_mmpi_{dataSet}.root')
c_mmpi.SaveAs(f'{OUTPATH}/{dataSet}/c_mmpi_{dataSet}.pdf')

# sys.exit("stopped process !")

### ----- kaon missing mass

c_mmpi = TCanvas("c_mmpi", "c_mmpi", 1000, 400)
c_mmpi.Divide(2,1)
c_mmpi.cd(1)
# c_mmpi.SetLogy()
# mmpi_offset = 0.1
# h_mmpi_precut = rdf.Define("MMpi_shift", f'P.kin.secondary.MMpi+{mmpi_offset}').Histo1D(("h_mmpi_precut", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "MMpi_shift")
h_mmpi_precut = rdf.Histo1D(("h_mmpi_precut", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_hCalCut = rdf.Filter(f'(({acceptCuts_data}) && H.cal.etottracknorm > {H_ecut_H_cal})').Histo1D(("h_mmpi_hCalCut", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_hCerCut = rdf.Filter(f'(({acceptCuts_data}) && H.cal.etottracknorm > {H_ecut_H_cal} && H.cer.npeSum > {P_kcut_H_cer})').Histo1D(("h_mmpi_hCerCut", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_RFtimeCut = rdf.Filter(f'(({acceptCuts_data}) && H.cal.etottracknorm > {H_ecut_H_cal} && H.cer.npeSum > {P_kcut_H_cer} && (({RFtimeDist_corr} < {KRFHigh}) && ({RFtimeDist_corr} > {KRFLow})))').Histo1D(("h_mmpi_RFtimeCut", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi = rdf.Filter(f'({cuts} && {promptCut})').Histo1D(("h_mmpi", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_random = rdf.Filter(f'({cuts} && {randomCut})').Histo1D(("h_mmpi_random", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_random.Scale(1.0/nWindows)
h_mmpi.Add(h_mmpi_random.GetPtr(), -1)

### ----- Normalize neutron peak with pion selection to kaon missing mass
h_mmpi_pionPID = rdf.Filter(f'({pion_cut} && {promptCut_pi})').Histo1D(("h_mmpi_pionPID", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_pionPID_random = rdf.Filter(f'({pion_cut} && {randomCut_pi})').Histo1D(("h_mmpi_pionPID_random", ";P.kin.secondary.MMpi;Counts", 300,0.85,1.4), "P.kin.secondary.MMpi")
h_mmpi_pionPID_random.Scale(1.0/nWindows)
h_mmpi_pionPID.Add(h_mmpi_pionPID_random.GetPtr(), -1)
normfac_pion = h_mmpi.Integral(h_mmpi.GetXaxis().FindBin(0.89),h_mmpi.GetXaxis().FindBin(0.92))/h_mmpi_pionPID.Integral(h_mmpi_pionPID.GetXaxis().FindBin(0.89),h_mmpi_pionPID.GetXaxis().FindBin(0.92))
h_mmpi_pionPID.Scale(normfac_pion)
h_mmpi_pionSub = TH1F("h_mmpi_pionSub", ";P.kin.secondary.MMpi; Counts",300,0.85,1.4)
h_mmpi_pionSub.Add(h_mmpi.GetPtr(),h_mmpi_pionPID.GetPtr(),1,-1)

h_mmpi_precut.SetLineColor(r.kBlack)
h_mmpi_hCalCut.SetLineColor(r.kRed)
h_mmpi_hCerCut.SetLineColor(r.kGreen)
h_mmpi_RFtimeCut.SetLineColor(r.kMagenta)
h_mmpi.SetLineColor(28)
h_mmpi_pionPID.SetLineColor(r.kRed)
h_mmpi_pionSub.SetLineColor(r.kBlue)
lg_mmpi = TLegend(0.4,0.7,0.6,0.9)
lg_mmpi.AddEntry(h_mmpi_precut.GetPtr(),"pre-cuts","l")
lg_mmpi.AddEntry(h_mmpi_hCalCut.GetPtr(),"H.cal.etottracknorm cut","l")
lg_mmpi.AddEntry(h_mmpi_hCerCut.GetPtr(),"H.cer.npeSum cut","l")
lg_mmpi.AddEntry(h_mmpi_RFtimeCut.GetPtr(),"RF time cut","l")
lg_mmpi.AddEntry(h_mmpi.GetPtr(),"all cuts","l")
lg_mmpi.AddEntry(h_mmpi_pionPID.GetPtr(),"MM_{#pi}","l")
lg_mmpi.AddEntry(h_mmpi_pionSub,"#pi subtracted","l")
h_mmpi_precut.SetMinimum(0.0)
h_mmpi_precut.Draw()
h_mmpi_hCalCut.Draw("same")
h_mmpi_hCerCut.Draw("same")
h_mmpi_RFtimeCut.Draw("same")
lg_mmpi.Draw()
c_mmpi.cd(2)
h_mmpi.SetMinimum(-50.0)
h_mmpi.Draw("hist")
h_mmpi_pionPID.Draw("hist same")
h_mmpi_pionSub.Draw("hist same")
h_mmpi.Write("h_mmpi",r.TObject.kWriteDelete)
h_mmpi_pionSub.Write("h_mmpi_pionSub",r.TObject.kWriteDelete)
c_mmpi.SaveAs(f'{OUTPATH}/{dataSet}/c_mmpi_{dataSet}.root')
c_mmpi.SaveAs(f'{OUTPATH}/{dataSet}/c_mmpi_{dataSet}.pdf')

### ============== Event Selection ================= end

### ================= Cut Efficiency ================ start

### ===== RF

cut_noRF = f'({acceptCuts_data} && (H.cal.etottracknorm > {H_ecut_H_cal}) && (H.cer.npeSum > {P_kcut_H_cer}) && (P.kin.secondary.MMpi < 0.95) && (({RFtimeDist_corr} > 0.0) && ({RFtimeDist_corr} < 2.0)))'
cut_RF = f'({acceptCuts_data} && (H.cal.etottracknorm > {H_ecut_H_cal}) && (H.cer.npeSum > {P_kcut_H_cer}) && (P.kin.secondary.MMpi < 0.95) && (({RFtimeDist_corr} > 0.0) && ({RFtimeDist_corr} < 1.0)))'

h_mmpi_noRF = rdf.Filter(f'({cut_noRF} && {promptCut})').Histo1D(("h_mmpi_noRF", ";P.kin.secondary.MMpi;Counts",300,0.85,0.95), "P.kin.secondary.MMpi")
h_mmpi_noRF_random = rdf.Filter(f'({cut_noRF} && {randomCut})').Histo1D(("h_mmpi_noRF_random", ";P.kin.secondary.MMpi;Counts",300,0.85,0.95), "P.kin.secondary.MMpi")
h_mmpi_noRF_random.Scale(1.0/nWindows)
h_mmpi_noRF.Add(h_mmpi_noRF_random.GetPtr(), -1)

h_mmpi_RF = rdf.Filter(f'({cut_RF} && {promptCut})').Histo1D(("h_mmpi_RF", ";P.kin.secondary.MMpi;Counts",300,0.85,0.95), "P.kin.secondary.MMpi")
h_mmpi_RF_random = rdf.Filter(f'({cut_RF} && {randomCut})').Histo1D(("h_mmpi_RF_random", ";P.kin.secondary.MMpi;Counts",300,0.85,0.95), "P.kin.secondary.MMpi")
h_mmpi_RF_random.Scale(1.0/nWindows)
h_mmpi_RF.Add(h_mmpi_RF_random.GetPtr(), -1)

xbin_min = 1
xbin_max = 300
dN_noRF = array.array('d', [0.0])
dN_RF = array.array('d', [0.0])
N_noRF = h_mmpi_noRF.IntegralAndError(xbin_min,xbin_max,dN_noRF,"")
N_RF = h_mmpi_RF.IntegralAndError(xbin_min,xbin_max,dN_RF,"")
eff_noRF = N_RF/N_noRF
deff_noRF = eff_noRF*m.sqrt((dN_RF[0]/N_RF)**2 + (dN_noRF[0]/N_noRF)**2)
c_mmpi_eff_rf = TCanvas("c_mmpi_eff_rf", "c_mmpi_eff_rf", 600, 400)
c_mmpi_eff_rf.cd()
h_mmpi_noRF.SetLineColor(1)
h_mmpi_RF.SetLineColor(6)
h_mmpi_noRF.Draw("hist")
h_mmpi_RF.Draw("hist same")
h_mmpi_noRF.Write("h_mmpi_noRF",r.TObject.kWriteDelete)
h_mmpi_RF.Write("h_mmpi_RF",r.TObject.kWriteDelete)
lg_eff_rf = TLegend(0.2,0.6,0.4,0.7)
lg_eff_rf.AddEntry(h_mmpi_noRF.GetPtr(),"]0 - 2[ ns","l")
lg_eff_rf.AddEntry(h_mmpi_RF.GetPtr(),"]0 - 1[ ns","l")
lg_eff_rf.Draw("same")
lat_eff_rf = TLatex()
lat_eff_rf.SetTextSize(0.05)
lat_eff_rf.SetTextAlign(13)
lat_eff_rf.SetNDC()
lat_eff_rf.DrawLatex(0.15,0.8,f'efficiency = %.3f #pm %.3f' % (eff_noRF,deff_noRF))
c_mmpi_eff_rf.Print(f'{OUTPATH}/{dataSet}/c_mmpi_eff_{dataSet}.root',"root")
c_mmpi_eff_rf.Print(f'{OUTPATH}/{dataSet}/c_mmpi_eff_{dataSet}.pdf',"pdf")
print("-"*40)
print(f'N_RF = {N_RF}, N_noRF = {N_noRF}, eff_RF = {eff_noRF}+/-{deff_noRF}\n')
print("-"*40)



