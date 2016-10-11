# MuonPOG
Collection of MuonPOG related tools

## Installation instructions

```bash
cmsrel CMSSW_8_0_20 # Just an example release, works in CMSSW >= 74X at present 
cd CMSSW_8_0_20/src/

git clone https://github.com/calderona/MuonPOG/ -b 80X-23SepReReco

cmsenv
scramv1 b -j 5
```

## Event dumper:

A duper, printing event information relevant for muons is available in :  MuonPOG/Tools/plugins/EventDumper.cc

It prints HLT, GEN level, beam spot, vertex and muon information. It works both in AOD and miniAOD.

To dump information from a given event :

```bash
cd MuonPOG/Tools/test/
cmsRun muonPrinter_cfg.py # modify files according to your needs
```

## Ntuples

The interface of muon Ntuples is defined in : MuonPOG/Tools/src/MuonPogTree.h

The code filling ntuples is available in : MuonPOG/Tools/plugins/MuonPogTreeProducer.cc

It fills HLT, GEN level, beam spot, vertex and muon information. It works both in AOD and miniAOD (NOTE: trigger information not filled when running in miniAOD).


To create some ntuples :

```bash
cd MuonPOG/Tools/test/
python muonPogNtuples_cfg.py --print # this will give you the default input parameters of the filler. 
                                     # As the ntuple cfg is based on VarParsing you can customise the
                                     # ntuple production via command line [1] or in a crab cfg [2] 

[1] 
cmsRun muonPogNtuples_cfg.py globalTag=80X_mcRun2_asymptotic_v5 \\
  eosInputFolder=/store/relval/CMSSW_8_0_0_patch2/RelValZMM_13/GEN-SIM-RECO/PU25ns_80X_mcRun2_asymptotic_v5_refGT-v1/10000

[2] 
https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters (find pyCfgParams)
```


## Running on CRAB

For running trhough crab you can go to: 

    MuonPOG/Tools/test/crab/

The CRAB client can be sourced using the command below after cmsenv.

    source /cvmfs/cms.cern.ch/crab3/crab.sh
  
Check if you have writing permissions in the common area if you already asked for that. 

    crab checkwrite --site=T2_CH_CERN --lfn=/store/group/phys_muon/

For running on the most recent re-reco data you can open file example: 

    crab_SingleMuon_Run2016B-23SepReReco_Run273158.py    

whit the configuration pset parameters: 

    config.JobType.pyCfgParams = ['globalTag=80X_dataRun2_2016SeptRepro_v3',
                                  'ntupleName=muonPOGNtuple_SingleMuonRun2016B_23SepReReco.root',
                                  'nEvents=-1',
                                  'runOnMC=False',
                                  'hltPathFilter=all',
                                  'minMuPt=10.0',
                                  'minNMu=2' ]


The ntuple producer gets loaded by :

```python
from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple
appendMuonPogNtuple(process,False,"HLT","MY_NTUPLE_NAME.root")
```

Where arguments are :

1. The CMS configuration process
2. A bool to say whether you are running on MC
3. The label of the process giving HLT results
4. The name of the output ntuple

## List of macros running on ntuples

The presently committed list of macros are:

| *Macro*        | *Description*  |
| -------------- | -------------- |
| MuonPOG/Tools/invariant_mass/  | makes dimuon invariant mas plots for different resonances  |
| MuonPOG/Tools/variables_comparison/  | performs a cut'n'count tnp study of commissioning and isolation variables, muon IDs and muon scale and resolution using muons from Z  |

Both macros use ini configuration files to allow freedom to chose parameters at run time, the are stored under the config/ directory of each macro package.

The syntax to run the macros is:

```bash
cd MuonPOG/Tools/invariant_mass/
./invariantMassPlots PATH_TO_INPUT_FILE PAT_TO_CONFIG_FILE(s)
```

```bash
MuonPOG/Tools/variables_comparison/
./variableComparisonPlots PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR
```

