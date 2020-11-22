#! /bin/bash

./bash/GeneratorMC_Linux.sh
mkdir Pythia8_genev || exit 1

strun=0
nruns=10000
njobs=100
nevents=100000
nTotal=1000000000

echo "[INFO] Starting production in Pythia8 with $nTotal events in $nruns runs"

for run in $(seq $strun $(($strun + $nruns - 1))); do

    ### wait if there are too many jobs running
    while true; do
        bkgjobs=$(jobs | grep GeneratorMC | wc -l | xargs)
        if [ $bkgjobs -lt $njobs ]; then
            break
        fi
        echo "[INFO] Sleep while waiting for a free job slot"
        sleep 10
    done
    
    runid=$(printf "%05d" $run)
    seed=$((123456789 + $run * 2))
    
    echo "[INFO] Starting run: $runid"

    ./exe/Anls_MonteCarloGenerator Pythia8_genev/outGeneratorMC_$runid $nevents $seed >& /dev/null & #&& \
	#./dropbox_uploader.sh upload result/outGeneratorMC_$runid.root Rubini/. && \
	#rm -rf result/outGeneratorMC_$runid.root &

    sleep 1s

done
exit 0
