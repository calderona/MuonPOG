#! /usr/bin/python

import json, sys
from pprint import pprint
from datetime import datetime, timedelta

secPerLumi = 23.3

if len(sys.argv) != 2 :
    print sys.argv[0], "running with arglist:", str(sys.argv), \
        "\ncorrect syntax :", sys.argv[0], "INPUT_JSON_FILE"
    sys.exit(100)

fileData=open(sys.argv[1]).read()
json = json.loads(fileData)

nRuns = 0.;
nLumis = 0.;

for run, lumiRanges in json.iteritems() :
    nRuns += 1
    for lumiRange in lumiRanges :
        nLumis += lumiRange[1] - lumiRange[0] + 1

nSeconds = timedelta(seconds=nLumis*secPerLumi)
time = datetime(1,1,1) + nSeconds

print

print("Total number of runs JSON : %i" % nRuns)
print("For a total number of GOOD lumisections:seconds : %i:%i\n" % (nLumis, nLumis*secPerLumi))
print("Which corresponds to a total running time of :\nMONTHS:DAYS:HOURS:MIN:SEC: %d:%d:%d:%d:%d" % (time.month-1, time.day-1, time.hour, time.minute, time.second))

print
