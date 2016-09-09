# MuonPOG
Collection of MuonPOG related tools

## Installation instructions

```bash
cmsrel CMSSW_7_6_3 # Just an example release, works in CMSSW >= 74X at present 
cd CMSSW_7_6_3/src/

git clone https://github.com/battibass/MuonPOG/

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

The CRAB client can be sourced using the command below after cmsenv.

    source /cvmfs/cms.cern.ch/crab3/crab.sh
  
Check if you have writing permissions in the common area.

    crab checkwrite --site=T2_CH_CERN --lfn=/store/group/phys_muon/Ntuples/2016/
    
Submit jobs.

    python multicrab.py samples/samples_spring15_miniaodv2_25ns.py
    python multicrab.py samples/samples_dataD_05Oct2015_25ns.py

Resubmit jobs.

    python multicrab.py crab_projects_21October resubmit

Check status.

    python multicrab.py crab_projects_21October status



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

