// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H

// Analysis Utility
#include "AliAnalysisUtility.h"

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

// Analysis Values

//-// File Names
auto const  fInvMasHist         =   "./result/InvariantMassHistograms.root";

//-// Tree Names
auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaonCandidate_Tree     =   "KaonCandidate";
auto const  fKaonCandidateEff_Tree  =   "KaonEfficiency";

//-// InvMass range Pythia MC
const   Float_t   fMinIMMC  =   0.75;
const   Float_t   fMaxIMMC  =   1.25;

//-// InvMass range 1D
const   Int_t     nBinIM1D  =   200;
const   Float_t   fMinIM1D  =   0.99;
const   Float_t   fMaxIM1D  =   1.05;
        Float_t  *fArrIM1D  =   new Float_t [nBinIM1D+1];

//-// InvMass range 2D
const   Int_t     nBinIM2D  =   100;
const   Float_t   fMinIM2D  =   0.99;
const   Float_t   fMaxIM2D  =   1.05;
        Float_t  *fArrIM2D  =   new Float_t [nBinIM2D+1];

//-// pT cuts 1D
        Int_t     nBinPT1D  =   32;
const   Float_t   fMinPT1D  =   0.0;
const   Float_t   fMaxPT1D  =   10.0;
        Float_t  *fArrPT1D  =   new Float_t [nBinPT1D+1];

//-// pT cuts 2D
        Int_t     nBinPT2D  =   12;
const   Float_t   fMinPT2D  =   0.0;
const   Float_t   fMaxPT2D  =   10.0;
        Float_t  *fArrPT2D  =   new Float_t [nBinPT2D+1];

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

void    fSetBinIM2D ()
{
    for (int i = 0; i <= nBinIM2D; i++ )
    {
        fArrIM2D[i] = fMinIM2D+(i)*(fMaxIM2D - fMinIM2D)/(static_cast<Float_t>(nBinIM2D));
    }
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

//------------------------------//
//    HISTOGRAM UTILITIES       //
//------------------------------//

template < class Tclass >
void    SetAxis             ( Tclass * aTarget, string aOption = "" )
{
    if ( aOption.find("IM") != -1 )
    {
        if ( aOption.find("2D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("m_{K^{+}K^{-}} candidate #phi_{1} (GeV/c^{2})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("m_{K^{+}K^{-}} candidate #phi_{2} (GeV/c^{2})");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("1D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
        }
    }
    else if ( aOption.find("PT") != -1 )
    {
        if ( aOption.find("DD") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi#phi}}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("2D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
                
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
                
            // Z-Axis formatting
            aTarget->GetZaxis()->SetTitle("#frac{d^{3}N_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
            aTarget->GetZaxis()->SetTitleOffset(1.15);
        }
        else if ( aOption.find("1D") != -1 )
        {
            // X-Axis formatting
            aTarget->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
            aTarget->GetXaxis()->SetTitleOffset(1.15);
            
            // Y-Axis formatting
            aTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
            aTarget->GetYaxis()->SetTitleOffset(1.15);
        }
    }
}

bool    fRapidityCut        ( Double_t  dRapidity )
{
    if ( fabs(dRapidity) <= 0.5 ) return true;
    return false;
}

bool    fTransverseMomCut   ( Double_t  dTransverseMom )
{
    if ( dTransverseMom < fMinPT1D ) return false;
    if ( dTransverseMom > fMaxPT1D ) return false;
    return true;
}

#endif
