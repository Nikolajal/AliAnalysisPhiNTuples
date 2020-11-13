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

// Data Structures

typedef struct
{
    UChar_t nPhi,           iKaon[1024],    jKaon[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_PhiCandidate;

typedef struct
{
    UChar_t nKaon,          Charge[1024];
    Char_t  SigmaTOF[1024], SigmaTPC[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
} Struct_KaonCandidate;

typedef struct
{
    UChar_t nPhi,           Selection[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  fTru,           fGen,           fRec;
} Struct_PhiEfficiency;

typedef struct
{
    UChar_t nKaon,          Charge[1024],   Selection[1024];
    Float_t Multiplicity,   Px[1024],       Py[1024],       Pz[1024],   InvMass[1024];
    Bool_t  ftru;
} Struct_KaonEfficiency;

#endif
