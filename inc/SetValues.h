// Global Values and constants file
// !TODO: All set!

#ifndef SETVALUES_H
#define SETVALUES_H

// C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>

// ROOT // Trees
#include "TTree.h"
#include "TLorentzVector.h"

// ROOT // Histograms
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

// ROOT // Functions
#include "TF1.h"

// ROOT // Graphs
#include "TGraphErrors.h"

// ROOT // I/O
#include "TFile.h"
#include "TRandom.h"

using namespace std;

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//

// - // File Names

auto const  fInvMasHist         =   "./result/InvariantMassHistograms.root";
auto const  fEfficiHist         =   "./result/Efficiencies_MCTruth.root";
auto const  fFitResHist         =   "./result/InvariantMassFitResultsPlots.root";
auto const  fFitResults         =   "./result/InvariantMassFitResults.root";
auto const  fAnlResHist         =   "./result/AnalysisResultsPlots.root";
auto const  fAnlResults         =   "./result/AnalysisResults.root";
auto const  fSystError_         =   "./result/Syst_SigExt.root";

// - // Efficiencies & Errors

auto const  kRapidityInterval   =   1.000;

// Uniform Bin by Bin

// Signal Extraction
auto const  kSyst_SigExtr1D     =   0.030;
auto const  kSyst_SigExtr2D     =   0.100;

// Branching Ratio
auto const  kBranchingRatio__   =   0.489;
auto const  kSyst_BRatio1D      =   0.005;
auto const  kSyst_BRatio2D      =   0.010;

// Vertex
auto const  kVertexEfficiency   =   0.990;

// Tracking
auto const  kSyst_TrackEff1D    =   0.080;
auto const  kSyst_TrackEff2D    =   0.160;

// PID
auto const  kSyst_PID1D         =   0.015;
auto const  kSyst_PID2D         =   0.030;

// On the final result

// Event
auto const  kEventEfficiency_   =   0.850;
auto const  kEventEfficienERP   =   0.062;
auto const  kEventEfficienERM   =   0.030;

// Total
auto const  kSystematical1D_   =   sqrt(kSyst_SigExtr1D*kSyst_SigExtr1D + kSyst_BRatio1D*kSyst_BRatio1D + kSyst_TrackEff1D*kSyst_TrackEff1D + kSyst_PID1D*kSyst_PID1D);
auto const  kSystematical2D_   =   sqrt(kSyst_SigExtr2D*kSyst_SigExtr2D + kSyst_BRatio2D*kSyst_BRatio2D + kSyst_TrackEff2D*kSyst_TrackEff2D + kSyst_PID2D*kSyst_PID2D);

// Pythia8
auto const  kPythia1DEfficien   =   0.970;
auto const  kPythia2DEfficien   =   0.943;

auto const  fFileDTAnal     = "./result/DTFinalAnalysisResults.root";
auto const  fFileZZ1D1D     = "./result/ZZFitResults1D-1D.root";
auto const  fFileZZ2D1D     = "./result/ZZFitResults1D-2D.root";
auto const  fFileZZ2D2D     = "./result/ZZFitResults2D-2D.root";
auto const  fFileZZ1DPT     = "./result/ZZFitResults1D-PT.root";
auto const  fFileZZ2DPT     = "./result/ZZFitResults2D-PT.root";
auto const  fSystError2     = "./result/bSystematicErrors.root";
auto const  fFileDT_10b     = "./result/LHC10b.root";
auto const  fFileDT_10c     = "./result/LHC10c.root";
auto const  fFileDT_10d     = "./result/LHC10d.root";
auto const  fFileDT_10e     = "./result/LHC10e.root";
auto const  fFileMC_10b     = "./result/LHC14j4b.root";
auto const  fFileMC_10c     = "./result/LHC14j4c.root";
auto const  fFileMC_10d     = "./result/LHC14j4d.root";
auto const  fFileMC_10e     = "./result/LHC14j4e.root";
auto const  bPythiaTest     = false;
auto const  bSave           = true;

// - // Physics Constants


// - // DATA
auto const  fFileDTForm     = "./result/DTZFetched.root";
auto const  fFileDTPreP     = "./result/DTPreProcessed.root";
auto const  fFileDTAn1D     = "./result/DTAnalysis1D.root";
auto const  fFileDTAn2D     = "./result/DTAnalysis2D.root";
auto const  fTreeSigName    = "SIG_Kaon_Tree";

// - // MONTECARLO
auto const  fFileMCForm     = "/Volumes/[HD][Nikolajal]_Toshiba/Rubini/outAAA.root";
//auto const  fFileMCForm     = "./result/MCZGenerated.root";
auto const  fFileMCPreP     = "./result/MCPreProcessed.root";
auto const  fTreeBkgName    = "BKG_Kaon_Tree";

// - // EFFICIENCY
auto const  fFileEFPreP     = "./result/EFPreProcessed.root";
auto const  fTreeTruName    = "TRU_Phi__Tree";


#endif
