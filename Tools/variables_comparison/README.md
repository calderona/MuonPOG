# What does the varables_comparison macro do?
It compares some ID and isolation variables in data and MC, using a cut and count tag and probe method.
It selects a tag muon matched to the trigger and a set of probes that fall in the Z mass peak range, then computes histogram for a set of variables of such probes.
It runs on DATA and on a programmable set of MC samles and overlays the DATA histograms to SHStacks from the sum of themontecarlo sample.

## How do I run it?
Simply by something like:

./variableComparisonPlots config_z/config.ini myResult

The code will compile, pharse the cfg, run and produce the results in a rootfile in the output directory.

## How do I configure it?
Using an INI file like the one in config_z/config.ini .
The cfg is rather self explanatory, it consist in different parts:
1.
A TagAndProbe section, defining the cuts on the tag, the probe eta binning and the cuts on the probe for isolation studies, as well as the Z mass window used for the study.
2. A sample section where ones has to specify the sample name (in the name of the section), where the ntuple of such sample is located, and the MC process corss section.
One can add as many samples as needed, the one with name [Data] is of course recognised and used differently, there the cross section value exist but is ignored.

## How do I add a variable to be monitored?
To add a variable to be monitored you should:

1. Define the variable in the Muon object of the tree in: ../src/MuonPogTree.h
2. Fill it properly in in the fillMuons method in: ../plugins/MuonPogTreeProducer.cc
3. ReRun the ntuple production using the cfg in: ../test/muonPogNtuples_cfg.py

Please think to all the variables you need before running so you avoid to run many times!

4. Book a plot in the Plotter::book method in: ./variableComparisonPlots.C
5. Fill the plot in the Plotter::fill method in: ./variableComparisonPlots.C
6. Ask the code to make a comparison plot by adding a muon_pog::comparisonPlot() call in the main function in: ./variableComparisonPlots.C

Please note that 3 types of plots exist:
a. Control plots, they are booked once for each plotter and not in multiple eta bins
b. Id variables plots, they are booked for each eta range the plotter gets in the configuration and filled for all probe muons that are tracker or global
c. Isolation variables plots, they are booked for each eta range the plotter gets in the configuration and filled for all probe muons after an id selection (programmable)

## What are the caveat, missing parts?

1.
The plots are normalised by the integral of the DATA plot and each MC contributes to tothe stack with its cross section.
No luminosity information is used, and plots filled before the invariant mass cut (for example the invariant mass cut plot) are known to be badly normalized.
You can change the normalization logic in the muon_pog::comparisonPlot() function of:  ./variableComparisonPlots.C

2.
The ratios below the plots are missing

3.
No PU reweighting is applied to the plots!
To change the weights look into the main function of ./variableComparisonPlots.C

4.
No attempt to plot variables for background (e.g. outside the Z peak) is made in the macro.

5.
The muon_pog::comparisonPlot() has colors setup only for 5 MC samples, for more MC you have to add them in the colorMap array.





