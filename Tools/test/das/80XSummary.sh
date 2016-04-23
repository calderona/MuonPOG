#! /bin/bash

eval `scramv1 runtime -sh`
LAST_UPDATE="Last update: "`date`

#echo $LAST_UPDATE &> data_all.txt
#echo $LAST_UPDATE &> cosmics_all.txt
echo $LAST_UPDATE &> mc_all.txt
echo $LAST_UPDATE &> summary.txt

#echo >> data_all.txt
#echo >> cosmics_all.txt
echo >> mc_all.txt
echo >> summary.txt

#das_client.py --limit=0 --query="dataset dataset=/*Mu*/*16Dec2015*/AOD* status=*"  >> data_all.txt
#das_client.py --limit=0 --query="dataset dataset=/*Cosm*/*/* release=CMSSW_8_0_* status=*"  >> cosmics_all.txt
das_client.py --limit=0 --query="dataset dataset=/*/*DR80*/AOD* status=*" >> mc_all.txt

{
#    echo "[$0] Status of SingleMuon data reprocessing:" 
#    cat data_all.txt | grep Single
#    echo
#    echo

#    echo "[$0] Status of DoubleMuon data reprocessing:" 
#    cat data_all.txt | grep Double
#    echo
#    echo

#    echo "[$0] Status of Cosmics reprocessing:" 
#    cat cosmics_all.txt | grep -v RelVal | grep -v update 
#    echo
#    echo

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
} >> summary.txt

mkdir -p ~/www/80X_Reco/
cp *txt ~/www/80X_Reco/



 
