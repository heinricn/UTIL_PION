#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-03-13 01:28:15 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

# Import relevant packages
import uproot
import uproot as up
import numpy as np

np.bool = bool
np.float = float

import root_numpy as rnp
import pandas as pd
import root_pandas as rpd
import ROOT
import scipy
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import sys, math, os, subprocess, csv
import awkward as ak

##################################################################################################################################################

# Check the number of arguments provided to the script
if len(sys.argv)-1!=3:
    print("!!!!! ERROR !!!!!\n Expected 3 arguments\n Usage is with - ROOTfilePrefix RunNumber MaxEvents \n!!!!! ERROR !!!!!")
    sys.exit(1)

##################################################################################################################################################

# Input params - run number and max number of events
ROOTPrefix = sys.argv[1]
runNum = sys.argv[2]
MaxEvent = sys.argv[3]

################################################################################################################################################

# Import package for cuts
from ltsep import Root

#Define and set up cuts
cut_f = '/DB/CUTS/run_type/coin_prod.cuts'

# defining Cuts
cuts = ["coin_epi_cut_accpt", "coin_epi_cut_all_RF", "coin_epi_cut_prompt_RF", "coin_epi_cut_rand_RF"]
#cuts = ["coin_epi_cut_accpt", "coin_epi_cut_all_noRF", "coin_epi_cut_prompt_noRF", "coin_epi_cut_rand_noRF"]
lt=Root(os.path.realpath(__file__),"Prod",ROOTPrefix,runNum,MaxEvent,cut_f,cuts)

USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
HOLEPATH="/u/group/c-pionlt/USERS/heinricn/hallc_replay_lt/UTIL_PION/scripts/pid/src/holeBoundry_60.csv" 

proc_root = lt.setup_ana()
c = proc_root[0] # Cut object
tree = proc_root[1] # Dictionary of branches
strDict = proc_root[2] # Dictionary of cuts as strings

###############################################################################################################################################

def loadHole(HoleCsv): #simple function, pass string with path to file and returns a fully constructed Root Polygon for doing HGC hole cut. - NH
    x = []
    y = []
    print("Attempting to open hole cut cevs at: ", HoleCsv)
    with open (HoleCsv) as f:
        for line in f:
            x.append(float(line[:line.find(",")]))
            y.append(float(line[line.find(",")+1:]))

    return [x,y]

###############################################################################################################################################

def coin_pions():
    # Define the list of relevant branches for COIN Pion analysis
    COIN_Pion_Data_Header = [
        "H_gtr_beta", "H_gtr_xp", "H_gtr_yp", "H_gtr_dp", "H_gtr_p", "H_dc_x_fp", "H_dc_y_fp", "H_dc_xp_fp", "H_dc_yp_fp",
        "H_hod_goodscinhit", "H_hod_goodstarttime", "H_cal_etotnorm", "H_cal_etottracknorm", "H_cer_npeSum",
        "CTime_ePiCoinTime_ROC1", "P_gtr_beta", "P_gtr_xp", "P_gtr_yp", "P_gtr_p", "P_gtr_dp", "P_dc_x_fp", "P_dc_y_fp",
        "P_dc_xp_fp", "P_dc_yp_fp", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_cal_etotnorm", "P_cal_etottracknorm",
        "P_aero_npeSum", "P_aero_xAtAero", "P_aero_yAtAero", "P_hgcer_npeSum", "P_hgcer_xAtCer", "P_hgcer_yAtCer",
        "P_ngcer_npeSum", "P_ngcer_xAtCer", "P_ngcer_yAtCer", "MMpi", "H_RF_Dist", "P_RF_Dist", "Q2", "W", "epsilon",
        "ph_q", "th_q", "MandelT", "pmiss", "pmiss_x", "pmiss_y", "pmiss_z", "emiss", "Erecoil", "Mrecoil",
    ]

    # Extract the data arrays for each branch
    NoCut_COIN_Pions = [tree[branch] for branch in COIN_Pion_Data_Header]

    # Create cut arrays for different scenarios
    Cut_COIN_Pions_accpt_tmp = [c.add_cut(arr, "coin_epi_cut_accpt") for arr in NoCut_COIN_Pions]
    Cut_COIN_Pions_all_tmp = [c.add_cut(arr, "coin_epi_cut_all_RF") for arr in NoCut_COIN_Pions]
    Cut_COIN_Pions_prompt_tmp = [c.add_cut(arr, "coin_epi_cut_prompt_RF") for arr in NoCut_COIN_Pions]
    Cut_COIN_Pions_rand_tmp = [c.add_cut(arr, "coin_epi_cut_rand_RF") for arr in NoCut_COIN_Pions]
    
    iXname = 0
    iYname = 0 
    for i in range(len(COIN_Pion_Data_Header)):
        if COIN_Pion_Data_Header[i] == "P_hgcer_xAtCer":
            iXname = i
        if COIN_Pion_Data_Header[i] == "P_hgcer_yAtCer":
            iYname = i
    #print("pre hole cut:")
    #print(Cut_COIN_Pions_prompt_tmp)      
    HGChole = loadHole(HOLEPATH)
    hgcHoleX = np.asarray(HGChole[0], dtype=float) #force the data type to be double*
    hgcHoleY = np.asarray(HGChole[1], dtype=float)
    HoleLen = len(hgcHoleX)
    hgcHolePoly = ROOT.TGeoPolygon(HoleLen) #constructer
    print("initializing hole. Number of verticies in Hole: ", HoleLen)
    hgcHolePoly.SetXY(hgcHoleX, hgcHoleY) #define verticies
    hgcHolePoly.FinishPolygon() #bake the polygon
    #temp = [[]]*len(COIN_Pion_Data_Header)
    temp = []
    for event in range(0,len(Cut_COIN_Pions_prompt_tmp[iXname])):
        px = Cut_COIN_Pions_prompt_tmp[iXname][event]
        py = Cut_COIN_Pions_prompt_tmp[iYname][event]
        hgcPoint = np.array([px,py],dtype=float)
        #hgcPoint = np.asarray(hgcPoint, dtype=float)
        #if(abs(px) < 2 and abs(py) < 2):
        #    print("potential event in hole: ", px, ", ", py)
        #    print("hgcPoint: ", hgcPoint)
        #    print(hgcHolePoly.Contains(hgcPoint))
        if hgcHolePoly.Contains(hgcPoint):
            temp.append(event)
    #print("length before removal: ", len(Cut_COIN_Pions_prompt_tmp[iXname]))
    print("number entries to remove: ", len(temp))
    
    if len(temp) != 0: 
        for Iarr in range(len(Cut_COIN_Pions_prompt_tmp)):
            ti = Cut_COIN_Pions_prompt_tmp[Iarr][:temp[0]]
            for i in range(1,len(temp)):
                tj = Cut_COIN_Pions_prompt_tmp[Iarr][temp[i-1]+1:temp[i]]
                ti = ak.concatenate([ti, tj],0)
            tj = Cut_COIN_Pions_prompt_tmp[Iarr][temp[len(temp)-1]+1:]
            Cut_COIN_Pions_prompt_tmp[Iarr] = ak.concatenate([ti, tj],0)

    #print("new length of array: ", len(Cut_COIN_Pions_prompt_tmp[iXname]))

    temp = []
    for event in range(0,len(Cut_COIN_Pions_rand_tmp[iXname])):
        px = Cut_COIN_Pions_rand_tmp[iXname][event]
        py = Cut_COIN_Pions_rand_tmp[iYname][event]
        hgcPoint = [px,py]
        hgcPoint = np.asarray(hgcPoint, dtype=float)
        #if(abs(px) < 2 and abs(py) < 2):
        #    print("potential event in hole: ", px, ", ", py)
        #    print("hgcPoint: ", hgcPoint)
        if hgcHolePoly.Contains(hgcPoint):
        #if (px < 2) and (px > -2):
            temp.append(event)
    #print("length before removal: ", len(Cut_COIN_Pions_rand_tmp[iXname]))
    print("number entries to remove: ", len(temp))
    
    if len(temp) != 0: 
        for Iarr in range(len(Cut_COIN_Pions_rand_tmp)):
            ti = Cut_COIN_Pions_rand_tmp[Iarr][:temp[0]]
            for i in range(1,len(temp)):
                tj = Cut_COIN_Pions_rand_tmp[Iarr][temp[i-1]+1:temp[i]]
                ti = ak.concatenate([ti, tj],0)
            tj = Cut_COIN_Pions_rand_tmp[Iarr][temp[len(temp)-1]+1:]
            Cut_COIN_Pions_rand_tmp[Iarr] = ak.concatenate([ti, tj],0)
    
    # Transpose to create event-wise records
    Uncut_COIN_Pions = list(zip(*NoCut_COIN_Pions))
    Cut_COIN_Pions_accpt = list(zip(*Cut_COIN_Pions_accpt_tmp))
    Cut_COIN_Pions_all = list(zip(*Cut_COIN_Pions_all_tmp))
    Cut_COIN_Pions_prompt = list(zip(*Cut_COIN_Pions_prompt_tmp))
    Cut_COIN_Pions_random = list(zip(*Cut_COIN_Pions_rand_tmp))

    COIN_Pions = {
        "Uncut_Pion_Events": Uncut_COIN_Pions,
        "Cut_Pion_Events_Accpt": Cut_COIN_Pions_accpt,
        "Cut_Pion_Events_All": Cut_COIN_Pions_all,
        "Cut_Pion_Events_Prompt": Cut_COIN_Pions_prompt,
        "Cut_Pion_Events_Random": Cut_COIN_Pions_random,
#        "Cut_Pion_Events_Random_Subtracted": Cut_COIN_Pions_random_subtracted
    }

    return COIN_Pions

##################################################################################################################################################################

def main():
    COIN_Pion_Data = coin_pions()

    COIN_Pion_Data_Header = [
        "H_gtr_beta", "H_gtr_xp", "H_gtr_yp", "H_gtr_dp", "H_gtr_p", "H_dc_x_fp", "H_dc_y_fp", "H_dc_xp_fp", "H_dc_yp_fp", "H_hod_goodscinhit", "H_hod_goodstarttime", "H_cal_etotnorm", "H_cal_etottracknorm", "H_cer_npeSum", "CTime_ePiCoinTime_ROC1", "P_gtr_beta", "P_gtr_xp", "P_gtr_yp", "P_gtr_p", "P_gtr_dp", "P_dc_x_fp", "P_dc_y_fp", "P_dc_xp_fp", "P_dc_yp_fp", "P_hod_goodscinhit", "P_hod_goodstarttime", "P_cal_etotnorm", "P_cal_etottracknorm", "P_aero_npeSum", "P_aero_xAtAero", "P_aero_yAtAero", "P_hgcer_npeSum", "P_hgcer_xAtCer", "P_hgcer_yAtCer", "P_ngcer_npeSum", "P_ngcer_xAtCer", "P_ngcer_yAtCer", "MMpi", "H_RF_Dist", "P_RF_Dist", "Q2", "W", "epsilon", "ph_q", "th_q", "MandelT", "pmiss", "pmiss_x", "pmiss_y", "pmiss_z", "emiss", "Erecoil", "Mrecoil"
    ]

    data = {}
    data.update(COIN_Pion_Data)
    data_keys = list(data.keys())

    for i in range(len(data_keys)):
        if "Pion" in data_keys[i]:
            DFHeader = list(COIN_Pion_Data_Header)
        else:
            continue

        df = pd.DataFrame(data.get(data_keys[i]), columns=DFHeader, index=None)
        # Convert all columns to numeric to avoid dtype('O') issues
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        if i == 0:
            df.to_root(
                "%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key="%s" % data_keys[i])
        elif i != 0:
            df.to_root(
                "%s/%s_%s_Analysed_Data.root" % (OUTPATH, runNum, MaxEvent), key="%s" % data_keys[i], mode='a')

if __name__ == '__main__':
    main()
print("Processing Complete")
