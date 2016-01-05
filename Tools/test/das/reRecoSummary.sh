#! /bin/bash

eval `scramv1 runtime -sh`

das_client.py --limit=0 --query="dataset=/*Mu*/*16Dec2015*/AOD*" &> data_all.txt
das_client.py --limit=0 --query="dataset=/*/*DR76*/AOD*" &> mc_all.txt

echo "[$0] Status of SingleMuon data reprocessing:" 
cat data_all.txt | grep Single
echo
echo

echo "[$0] Status of DoubleMuon data reprocessing:" 
cat data_all.txt | grep Double
echo
echo

echo "[$0] Status of ZMuMu MC reprocessing:" 
cat mc_all.txt | grep -i ZMuMu
echo
echo

echo "[$0] Status of DY MC reprocessing:" 
cat mc_all.txt | grep -i DY
echo
echo

echo "[$0] Status of ZMuMu JPsi reprocessing:" 
cat mc_all.txt | grep -i JPsi
echo
echo

echo "[$0] Status of ZMuMu MuEnriched reprocessing:" 
cat mc_all.txt | grep -i MuEnriched
echo
echo





 
