mkdir -p ./exe

    g++ ./src/GeneratorMC.C \
	-o ./exe/GeneratorMC \
	-std=c++11 \
	-lpythia8 \
	-L$PYTHIA_ROOT/lib/ \
	-I$PYTHIA_ROOT/include/ \
	-I$ROOTSYS/include/ \
	-L$ROOTSYS/lib/ \
	-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lm -ldl
