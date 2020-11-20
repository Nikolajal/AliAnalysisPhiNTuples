mkdir -p ./exe

    g++ ./src/Anls_MonteCarloGenerator.C \
	-o ./exe/GeneratorMC \
	-std=c++11 \
	-lpythia8 \
	-L$PYTHIA8/lib/ \
	-I$PYTHIA8/include/ \
	-I$ROOTSYS/include/ \
	-L$ROOTSYS/lib/ \
	-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lm -ldl
