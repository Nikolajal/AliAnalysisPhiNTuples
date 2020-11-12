// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H

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

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//

// Tree Names

auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaon_Tree              =   "KaonCandidate";
auto const  fKaonEff_Tree           =   "KaonEfficiency";

#endif
