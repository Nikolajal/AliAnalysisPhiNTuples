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

// TRandom
TRandom *   fRandomGen      =   new TRandom();

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
auto        hName           = "Name";
auto        hTitle          = "Title";
auto const  bPythiaTest     = false;
auto const  bSave           = true;

// - // Physics Constants

//---------------------------------------//
//- Physics is by defualt in GeV, we    -//
//- define some constants to evaluate   -//
//- results in other units              -//
//---------------------------------------//

auto const  kPMas           = 1.019455; // 1.019455 ± 0.000020 GeV
auto const  kPWid           = 0.004260; // 0.004260 ± 0.000040 GeV
auto const  KeV             = 1e6;
auto const  MeV             = 1e3;
auto const  GeV             = 1;
auto const  TeV             = 1e-3;

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

// INVARIANT MASS
//-// InvMass range Pythia MC
const   Float_t   fMinIMMC  = 0.75;
const   Float_t   fMaxIMMC  = 1.25;

//-// InvMass range 1D
const   Int_t     nBinIM1D  = 200;
const   Float_t   fMinIM1D  = 0.99;
const   Float_t   fMaxIM1D  = 1.05;
        Float_t * fArrIM1D  = new Float_t [nBinIM1D+1];

//-// InvMass range 2D
const   Int_t     nBinIM2D  = 100;
const   Float_t   fMinIM2D  = 0.99;
const   Float_t   fMaxIM2D  = 1.05;
        Float_t * fArrIM2D  = new Float_t [nBinIM2D+1];

// MOMENTUM
// pT cuts 1D
        Int_t     nBinPT1D  = 32;
const   Float_t   fMinPT1D  = 0.0;
const   Float_t   fMaxPT1D  = 10.0;
        Float_t * fArrPT1D  = new Float_t [nBinPT1D+1];

// pT cuts 2D
        Int_t     nBinPT2D  = 12;
const   Float_t   fMinPT2D  = 0.0;
const   Float_t   fMaxPT2D  = 10.0;
        Float_t * fArrPT2D  = new Float_t [nBinPT2D+1];

// DATA Structures
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

//------------------------------//
//    VARIABLES UTILITIES       //
//------------------------------//

void    fSetBinIM1D ()
{
    for (int i = 0; i <= nBinIM1D; i++ )
    {
        fArrIM1D[i] = fMinIM1D+(i)*(fMaxIM1D - fMinIM1D)/(static_cast<Float_t>(nBinIM1D));
    }
}

Int_t   fGetBinIM1D (Float_t input_value )
{
    if ( input_value > fMaxIM1D ) return -1;
    for ( Int_t iBin = 0; iBin <= nBinIM1D; iBin++ )
    {
        if ( input_value <= fArrIM1D[iBin] )
        {
            return iBin -1;
        }
    }
    return nBinIM1D-1;
}

void    fSetBinIM2D ()
{
    for (int i = 0; i <= nBinIM2D; i++ )
    {
        fArrIM2D[i] = fMinIM2D+(i)*(fMaxIM2D - fMinIM2D)/(static_cast<Float_t>(nBinIM2D));
    }
}

Int_t   fGetBinIM2D (Float_t input_value )
{
    if ( input_value > fMaxIM2D ) return -1;
    for ( Int_t iBin = 0; iBin <= nBinIM2D; iBin++ )
    {
        if ( input_value <= fArrIM2D[iBin] )
        {
            return iBin -1;
        }
    }
    return nBinIM2D-1;
}

void    fSetBinPT1D ()
{
    fArrPT1D[0] =   0.0;
    fArrPT1D[1] =   0.1;
    fArrPT1D[2] =   0.2;
    fArrPT1D[3] =   0.3;
    fArrPT1D[4] =   0.4;
    fArrPT1D[5] =   0.5;
    fArrPT1D[6] =   0.6;
    fArrPT1D[7] =   0.7;
    fArrPT1D[8] =   0.8;
    fArrPT1D[9] =   0.9;
    fArrPT1D[10]=   1.0;
    fArrPT1D[11]=   1.1;
    fArrPT1D[12]=   1.2;
    fArrPT1D[13]=   1.3;
    fArrPT1D[14]=   1.4;
    fArrPT1D[15]=   1.5;
    fArrPT1D[16]=   1.6;
    fArrPT1D[17]=   1.7;
    fArrPT1D[18]=   1.8;
    fArrPT1D[19]=   1.9;
    fArrPT1D[20]=   2.0;
    fArrPT1D[21]=   2.2;
    fArrPT1D[22]=   2.4;
    fArrPT1D[23]=   2.6;
    fArrPT1D[24]=   2.8;
    fArrPT1D[25]=   3.0;
    fArrPT1D[26]=   3.5;
    fArrPT1D[27]=   4.0;
    fArrPT1D[28]=   4.5;
    fArrPT1D[29]=   5.0;
    fArrPT1D[30]=   6.0;
    fArrPT1D[31]=   8.0;
    fArrPT1D[32]=   10.;
}

void    fSetBinPT2D ()
{
    fArrPT2D[0] =   0.0;
    fArrPT2D[1] =   0.2;
    fArrPT2D[2] =   0.40;
    fArrPT2D[3] =   0.68;
    fArrPT2D[4] =   0.82;
    fArrPT2D[5] =   0.95;
    fArrPT2D[6] =   1.1;
    fArrPT2D[7] =   1.3;
    fArrPT2D[8] =   1.6;
    fArrPT2D[9] =   2.3;
    fArrPT2D[10] =  3.0;
    fArrPT2D[11] =  5.0;
    fArrPT2D[12] =  10.;
 }

Int_t   fGetBinPT1D (Float_t input_value )
{
    for ( Int_t iBin = 0; iBin <= nBinPT1D; iBin++ )
    {
        if ( input_value <= fArrPT1D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}

Int_t   fGetBinPT2D (Float_t input_value )
{
    for ( Int_t iBin = 0; iBin <= nBinPT2D; iBin++ )
    {
        if ( input_value <= fArrPT2D[iBin] )
        {
            return iBin -1;
        }
    }
    return -1;
}

#endif
