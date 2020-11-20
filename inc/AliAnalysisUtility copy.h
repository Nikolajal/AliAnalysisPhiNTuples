// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISUTILITY_H
#define ALIANALYSISUTILITY_H

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

// Random Generator
TRandom *   fRandomGen      =   new TRandom();

// Title and Name for histograms
auto        hName           =   "Name";
auto        hTitle          =   "Title";

//---------------------------------------//
//- Physics is by defualt in GeV, we    -//
//- define some constants to evaluate   -//
//- results in other units              -//
//---------------------------------------//

auto const  KeV             =   1e6;
auto const  MeV             =   1e3;
auto const  GeV             =   1;
auto const  TeV             =   1e-3;

#endif
