#! /usr/bin/env bash

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC10b/LHC10b/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC10c/LHC10c/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC10d/LHC10d/ AnalysisResults.root | grep AnalysisResults.root`; do

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC10e/LHC10e/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC14j4b/LHC14j4b/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC14j4c/LHC14j4c/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC14j4d/LHC14j4d/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done

for I in `alien_find /alice/cern.ch/user/n/nrubini/LHC14j4e/LHC14j4e/ AnalysisResults.root | grep AnalysisResults.root`; do
    
    ((ifile+=1))

while [ $(ps -ef | grep alien_cp | wc -l) -ge 20 ]; do
    sleep 1s;
done

    dir=./${I%AnalysisResults.root}
    mkdir -p $dir
    alien_cp $I file://$dir/AnalysisResults.root &
    
done
