#! /bin/bash

eval `scramv1 runtime -sh`

das_client.py --limit=0 --query="dataset dataset=/*Mu*/*16Dec2015*/AOD* "  &> data_all.txt
das_client.py --limit=0 --query="dataset dataset=/*Cosm*/*/* release=CMSSW_7_6_*"  &> cosmics_all.txt
das_client.py --limit=0 --query="dataset=/*/*DR76*/AOD*" &> mc_all.txt

echo "[$0] Status of SingleMuon data reprocessing:" 
cat data_all.txt | grep Single
echo
echo

echo "[$0] Status of DoubleMuon data reprocessing:" 
cat data_all.txt | grep Double
echo
echo

echo "[$0] Status of Cosmics reprocessing:" 
cat cosmics_all.txt | grep -v RelVal
echo
echo

echo "[$0] Status of ZToMuMu MC reprocessing:" 
cat mc_all.txt | grep -i ZToMuMu
echo
echo

echo "[$0] Status of DY MC reprocessing:" 
cat mc_all.txt | grep -i DY
echo
echo

echo "[$0] Status of JPsi reprocessing:" 
cat mc_all.txt | grep -i JPsi
echo
echo

echo "[$0] Status of MuEnriched reprocessing:" 
cat mc_all.txt | grep -i MuEnriched
echo
echo





 
