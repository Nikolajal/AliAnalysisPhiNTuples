g++ ./root/AnalysisCount.C \
-o ./exe/AnalysisCount \
-std=c++11 \
-I$ROOTSYS/include/ \
-L$ROOTSYS/lib/ \
-lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lm -ldl 
