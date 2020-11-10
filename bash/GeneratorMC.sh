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
  
    g++ -g -Wall `root-config --cflags --libs` \
    ./root/Anls_MonteCarloGenerator.C \
    -lpythia8 \
    -L/Applications/pythia8303/lib\
    -I/Applications/pythia8303/include\
	-o ./exe/Anls_MonteCarloGenerator -std=c++11 -ldl
 
 install_name_tool -change @rpath/libpythia8.dylib /Applications/pythia8303/lib/libpythia8.dylib /Users/nikolajal/alidock/AliAnalysisPhiPair/exe/Anls_MonteCarloGenerator 
fi
