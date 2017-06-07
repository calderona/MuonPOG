# MuonPOG
Collection of MuonPOG related tools

## Installation instructions

```bash
cmsrel CMSSW_9_2_1 # Just an example release, works in CMSSW >= 74X at present 
cd CMSSW_9_2_1/src/

git clone https://github.com/calderona/MuonPOG/ -b 92X

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

