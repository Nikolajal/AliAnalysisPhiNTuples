mkdir -p ./exe

if [ `hostname` == "bownalice07.bo.infn.it" ]; then

    g++ ./root/GeneratorMC.C \
	-o ./exe/GeneratorMC \
	-std=c++11 \
	-lpythia8 \
	-L$PYTHIA_ROOT/lib/ \
	-I$PYTHIA_ROOT/include/ \
	-I$ROOTSYS/include/ \
	-L$ROOTSYS/lib/ \
	-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lm -ldl
    
else
  
    g++ ./root/GeneratorMC.C \
	-o ./exe/GeneratorMC \
	-std=c++11 \
	-lpythia8 \
	-L/Applications/pythia8243/lib/ \
	-I/Applications/pythia8243/include/ \
	-I/Applications/root-build/include/ \
	-L/Applications/root-build/lib/ \
	-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -stdlib=libc++ -lm -ldl
  
fi