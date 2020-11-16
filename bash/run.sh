#! /bin/bash

./bash/GeneratorMC.sh
mkdir Pythia8_genev || exit 1

strun=0
nruns=1000
njobs=16
nevents=1000

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

    ./exe/Anls_MonteCarloGenerator Pythia8_genev/outGeneratorMC_$runid $nevents >& /dev/null & #&& \
	#./dropbox_uploader.sh upload result/outGeneratorMC_$runid.root Rubini/. && \
	#rm -rf result/outGeneratorMC_$runid.root &

    sleep 1s

done
exit 0
