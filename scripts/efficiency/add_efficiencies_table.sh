#! /bin/bash

#
# Description:
# ================================================================
# Time-stamp: "2024-02-09 12:05:00 junaid"
# ================================================================
#
# Copied: Richard L. Trotta III <trotta@cua.edu>
# Created: Muhammad junaid  <mjo147@uregina.ca>
# Copyright (c) trottar & junaid
#

while getopts 'hprs' flag; do
    case "${flag}" in
        h) 
        echo "-------------------------------------------------------------------"
        echo "./add_efficiencies_table.sh -{flags} {variable arguments, see help}"
        echo "-------------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
	echo "    If no flags called arguments are..."
	echo "        coin -> RunType=arg1"
	echo "        sing -> RunType=arg1 SPEC=arg2 (requires -s flag)"		
        echo "    -h, help"
        echo "    -p, plot efficiencies (code exits upon completion, table required)"
	echo "        coin -> RunType=arg1 DATE=arg2"
	echo "        sing -> RunType=arg1 DATE=arg2 SPEC=arg3 (requires -s flag)"
	echo "    -r, run hgcer root analysis"
	echo "        coin -> RunType=arg1"
	echo "        sing -> RunType=arg1 SPEC=arg2 (requires -s flag)"	
	echo "    -s, single arm"
        exit 0
        ;;
        p) p_flag='true' ;;
	r) r_flag='true' ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

RunList=$2

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ "${HOSTNAME}" = *"farm"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi
#echo $PWD
#echo
#echo $PATHFILE_INFO
#echo
# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
#echo $ANATYPE

cd "${SCRIPTPATH}/efficiency/src/"

if [[ $p_flag = "true" ]]; then
    RunType=$3
    DATE=$4
    if [[ $RunType = "HeePCoin" ]]; then
	ROOTPREFIX=PionLT_replay_HeeP_coin
#	HGCERPREFIX=${ANATYPE}_coin_replay_production
    elif [[ $RunType = "LumiCoin" ]]; then
	ROOTPREFIX=PionLT_replay_luminosity_coin
#	HGCERPREFIX=${ANATYPE}_coin_replay_production	
    elif [[ $RunType = "pTRIG6" ]]; then
        ROOTPREFIX=PionLT_replay_coin_production_pTRIG6
#       HGCERPREFIX=${ANATYPE}_coin_replay_production
    elif [[ $RunType = "Prod" ]]; then
        ROOTPREFIX=PionLT_replay_coin_production
#       HGCERPREFIX=${ANATYPE}_coin_replay_production 
    else
	ROOTPREFIX=PionLT_replay_coin_production
#       HGCERPREFIX=${ANATYPE}_coin_replay_production 
    fi
    python3 plot_efficiency.py ${ROOTPREFIX} ${RunType} ${DATE}
#    python3 plot/plot_efficiency_beam.py ${ROOTPREFIX} ${RunType} ${DATE}
    cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
    convert *.png "${RunType}_${DATE}.pdf"
    evince "${RunType}_${DATE}.pdf"
    rm -f *.png
    exit 1

elif [[ $p_flag = "true" && $s_flag = "true" ]]; then
    RunType=$3
    DATE=$4
    spec=$(echo "$5" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    if [[ $RunType = "HeePSing" ]]; then
	ROOTPREFIX=PionLT_replay_${spec}_HeePSing
#	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RunType = "LumiSing" ]]; then
        ROOTPREFIX=PionLT_replay_luminosity
#       HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
        #inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    else
	ROOTPREFIX=PionLt_replay_${spec}_production
#	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    fi
    python3 plot_efficiency.py ${ROOTPREFIX} ${RunType} ${DATE}
#    python3 plot/plot_efficiency_beam.py ${ROOTPREFIX} ${RunType} ${DATE}
    cd "${SCRIPTPATH}/efficiency/OUTPUTS/plots"
    convert *.png "${RunType}_${DATE}.pdf"
    evince "${RunType}_${DATE}.pdf"
    exit 1

elif [[ $s_flag = "true" ]]; then
    RunType=$3
    spec=$(echo "$4" | tr '[:upper:]' '[:lower:]')
    SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
    if [[ $RunType = "HeePSing" ]]; then
	ROOTPREFIX=PionLT_replay_${SPEC}_HeePSing
#	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RunType = "LumiSing" ]]; then
        ROOTPREFIX=PionLT_replay_luminosity
#       HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
        #inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    else
	ROOTPREFIX=PionLT_replay_${spec}_production
#	HGCERPREFIX=${ANATYPE}_${SPEC}_replay_production	
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    fi    

else
    RunType=$3
    if [[ $RunType = "HeePCoin" ]]; then
        ROOTPREFIX=PionLT_replay_HeeP_coin
#        HGCERPREFIX=${ANATYPE}_coin_replay_production
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
        #inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RunType = "LumiCoin" ]]; then
	ROOTPREFIX=PionLT_replay_luminosity_coin
#	HGCERPREFIX=${ANATYPE}_coin_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RunType = "pTRIG6" ]]; then
        ROOTPREFIX=PionLT_replay_coin_production_pTRIG6
#       HGCERPREFIX=${ANATYPE}_coin_replay_production
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
        #inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    elif [[ $RunType = "Prod" ]]; then
        ROOTPREFIX=PionLT_replay_coin_production
#       HGCERPREFIX=${ANATYPE}_coin_replay_production
        inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
        #inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    else
	ROOTPREFIX=PionLT_replay_coin_production
#	HGCERPREFIX=${ANATYPE}_coin_replay_production
	inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
	#inputFile="${REPLAYPATH}/UTIL_BATCH/InputRunLists/eff_runlist/${RunList}"
    fi    
fi

if [[ $r_flag = "true" ]]; then
    while true; do
	read -p "Do you wish to analyse hgcer efficiencies with run list ${inputFile}? (Please answer yes or no) " yn
	case $yn in
	    [Yy]* )
		i=-1
		(
		##Reads in input file##
		while IFS='' read -r line || [[ -n "$line" ]]; do
		    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		    echo "Run number read from file: $line"
		    echo ""
		    cd "${SCRIPTPATH}/efficiency/src/hgcer"
		    python3 hgcer.py ${HGCERPREFIX} $line -1
		done < "$inputFile"
		)
		break;;
	    [Nn]* ) 
		exit;;
	    * ) echo "Please answer yes or no.";;
	esac
    done
else
    while true; do
	read -p "Do you wish to create efficiency table for run list ${inputFile}? (Please answer yes or no) " yn
	case $yn in
	    [Yy]* )
		i=-1
		(
		##Reads in input file##
		while IFS='' read -r line || [[ -n "$line" ]]; do
		    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		    echo "Run number read from file: $line"
		    echo ""
		    python3 efficiency_main.py ${ROOTPREFIX} $RunType $line -1
		done < "$inputFile"
		)
		break;;
	    [Nn]* ) 
		exit;;
	    * ) echo "Please answer yes or no.";;
	esac
    done
fi
