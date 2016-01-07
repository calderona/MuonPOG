# MuonPOG
Collection of MuonPOG related tools

## Installation instructions

cmsrel CMSSW_7_6_3 # Just an example release, works in CMSSW >= 74X at present 

cd CMSSW_7_6_3/src/

git clone git@github.com:battibass/MuonPOG.git

cmsenv

scramv1 b -j 5

## Event dumper:

A duper, printing event information relevant for muons is available in :  MuonPOG/Tools/plugins/EventDumper.cc

It prints HLT, GEN level, beam spot, vertex and muon information. It works both in AOD and miniAOD.

To dump information from a given event :

cd MuonPOG/Tools/test/

cmsRun muonPrinter_cfg.py # modify files according to your needs

## Ntuples

The interface of muon Ntuples is defined in : MuonPOG/Tools/src/MuonPogTree.h

The code filling ntuples is available in : MuonPOG/Tools/plugins/MuonPogTreeProducer.cc

It fills HLT, GEN level, beam spot, vertex and muon information. It works both in AOD and miniAOD.

To create some ntuples :

cd MuonPOG/Tools/test/

cmsRun muonPogNtuples_cfg.py # modify files according to your needs

The ntuple producer gets loaded by :

from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple

appendMuonPogNtuple(process,False,"HLT","MY_NTUPLE_NAME.root")

Where arguments are :

1. The CMS configuration process
2. A bool to say whether you are running on MC
3. The label of the process giving HLT results
4. The name of the output ntuple

## List of macros running on ntuples

The presently committed list of macros are:

MuonPOG/Tools/invariant_mass/ : makes dimuon invariant mas plots for different resonances

MuonPOG/Tools/variables_comparison/ : performs a cut'n'count tnp study of commissioning and isolation variables, muon IDs and muon scale and resolution using muons from Z

Both macros use ini configuration files to allow freedom to chose some configuration parameters, the are stored under the config/ directory in each macro package.

The syntax to run the macros is:

cd MuonPOG/Tools/invariant_mass/

./invariantMassPlots PATH_TO_INPUT_FILE PAT_TO_CONFIG_FILE(s)

MuonPOG/Tools/variables_comparison/

./variableComparisonPlots PATH_TO_CONFIG_FILE PATH_TO_OUTPUT_DIR

