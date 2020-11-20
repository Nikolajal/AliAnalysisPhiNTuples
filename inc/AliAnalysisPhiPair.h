// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISPHIPAIR_H
#define ALIANALYSISPHIPAIR_H

// Analysis Utility
#include "AliAnalysisUtility.h"

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//

// Analysis Values

auto const  bPythiaTest             =   kTRUE;

//-// File Names
auto const  fInvMasHist             =   "./result/InvariantMassHistograms.root";
auto const  fEfficiHist             =   "./result/Efficiencies_MCTruth.root";
auto const  fFitResHist         =   "./result/InvariantMassFitResultsPlots.root";
auto const  fFitResults         =   "./result/InvariantMassFitResults.root";
auto const  fAnlResHist         =   "./result/AnalysisResultsPlots.root";
auto const  fAnlResults         =   "./result/AnalysisResults.root";
auto const  fSystError_         =   "./result/Syst_SigExt.root";

//-// Tree Names
auto const  fPhiCandidate_Tree      =   "PhiCandidate";
auto const  fPhiCandidateEff_Tree   =   "PhiEfficiency";
auto const  fKaonCandidate_Tree     =   "KaonCandidate";
auto const  fKaonCandidateEff_Tree  =   "KaonEfficiency";

auto const  kParticleMass_          =   1.019455;   //  1.019455    +- 0.000020
auto const  kParticleWidth          =   0.00426;    //  0.00426     +- 0.00004
auto const  kDetectorSlope          =   1.;

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

//------------------------------//
//    ANALYSISI SEPCIFIC Fncs   //
//------------------------------//

Int_t kColor[] = {38,kBlue,kBlue+3,46,38};
Int_t kStyle[] = {26,9,10,25,22};
Int_t kWidth[] = {1,3,3,1,1};

int             fLegendSelect                   ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 2;
    else return -1;
}

void            fLegendMaker                    ( RooPlot * fRooPlot, const char * fSelect, TLegend * fLegend )
{
    switch (fLegendSelect(fSelect))
    {
        case 1:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg)",            "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        case 2:
            fLegend                     ->SetFillColor(kWhite);
            fLegend                     ->SetLineColor(kWhite);
            fLegend                     ->AddEntry(fRooPlot->findObject("RooData"), "Data",                 "EP");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSS"),   "Fit (Sig #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBS"),   "Fit (Bkg #times Sig)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooSB"),   "Fit (Sig #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooBB"),   "Fit (Bkg #times Bkg)", "L");
            fLegend                     ->AddEntry(fRooPlot->findObject("RooMod"),  "Fit (Model)",          "L");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

int             fAxisSelect                     ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 3;
    else return -1;
}
 
void            fAxisMaker                      ( RooPlot * fRooPlot, const char * fSelect )
{
    switch (fAxisSelect(fSelect))
    {
        case 1:
            fRooPlot                    ->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
            break;
        case 2:
            fRooPlot                    ->GetXaxis()->SetTitle("m^{x}_{K^{+}K^{-}} (GeV/c^{2})");
            break;
        case 3:
            fRooPlot                    ->GetXaxis()->SetTitle("m^{y}_{K^{+}K^{-}} (GeV/c^{2})");
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

int             fPlotterSelect                  ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 2;
    else return -1;
}

void            fRooPlotPlotter                 ( RooPlot * fRooPlot, const char * fSelect, RooAddPdf fModel , RooDataHist * fData )
{
    switch (fPlotterSelect(fSelect))
    {
        case 1:
            fData                           ->plotOn(fRooPlot,      MarkerColor(38),                MarkerStyle(26),    Name("RooData"));
            fModel                          .plotOn (fRooPlot,      LineColor(4),                   LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      Components("fBkg"),             LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                          .plotOn (fRooPlot,      Components("fSig"),             LineColor(2),       Name("RooSS"));
            break;
        case 2:
            fData                           ->plotOn(fRooPlot,      CutRange("fDrawRange"),         MarkerColor(38),    MarkerStyle(26) ,   Name("RooData"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  LineColor(4),       LineStyle(kDashed), Name("RooMod"));
            fModel                          .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkg"), LineStyle(kDashed), LineColor(38),      Name("RooBB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigSig"),LineColor(2),       Name("RooSS"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fSigBkg"),LineStyle(kDashed), LineColor(33),    Name("RooSB"));
            fModel                      .plotOn (fRooPlot,      ProjectionRange("fDrawRange"),  Components("fBkgSig"),LineStyle(kDashed), LineColor(36),    Name("RooBS"));
            break;
        default:
            cout << "Improper option, no changes made" << endl;
            break;
    }
}

void            fRooPlotMaker                   ( RooPlot * fRooPlot, TLegend * fLegend, RooAddPdf fModel , RooDataHist * fData, const char * fSelect )
{
    fRooPlotPlotter(fRooPlot,fSelect,fModel,fData);
    fLegendMaker(fRooPlot,fSelect,fLegend);
    fAxisMaker(fRooPlot,fSelect);
}

void            fCoreFitModelSetBoundaries      ( string fOption, Double_t &aValMin, Double_t &aValMax )
{
    aValMin = 0.99;
    aValMax = 1.05;
    if ( fOption.find("RA") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.040;
    }
    if ( fOption.find("RB") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.060;
    }
    if ( fOption.find("RC") != -1 )
    {
        aValMin =   0.990;
        aValMax =   1.070;
    }
    if ( fOption.find("RD") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.040;
    }
    if ( fOption.find("RE") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.050;
    }
    if ( fOption.find("RF") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.060;
    }
    if ( fOption.find("RG") != -1 )
    {
        aValMin =   0.995;
        aValMax =   1.070;
    }
    if ( fOption.find("RH") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.040;
    }
    if ( fOption.find("RI") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.050;
    }
    if ( fOption.find("RJ") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.060;
    }
    if ( fOption.find("RK") != -1 )
    {
        aValMin =   1.000;
        aValMax =   1.070;
    }
}

Double_t*       fCoreFitModelOptionSelect       ( string fOption )
{
    Double_t   *fResult     =   new Double_t [7];
    fResult = { 0., 0., 1., 0., 0., 0., 1. };
    fCoreFitModelSetBoundaries ( fOption, fResult[0], fResult[1] );
    if  ( fOption.find("CH3") != -1 ) { fResult[2] = 0.; };
    if  ( fOption.find("CH5") != -1 ) { fResult[3] = 1.; };
    if  ( fOption.find("Wdt") != -1 ) { fResult[4] = kParticleMass_*0.1; };
    if  ( fOption.find("Mss") != -1 ) { fResult[5] = kParticleWidth*0.1; };
    if  ( bPythiaTest               ) { fResult[6] = 0.; };
}

RooFitResult*   fCoreFitModel                   ( RooDataHist *fDataHist, RooRealVar fVariable, RooAddPdf &fModel_, Double_t kCh4Limits, Double_t kCh5Limits, Double_t kVmsLimits, Double_t kVwdLimits, Double_t kVslLimits )
{
    //------ Define what your model is made of -----//
    
    Int_t       nEntries        =   fDataHist->sumEntries();
    
    // Background PDF Coefficients
    RooRealVar  ChebychevPar0   =   RooRealVar      ("ch0","ch0",   0., -1, 1);
    RooRealVar  ChebychevPar1   =   RooRealVar      ("ch1","ch1",   0., -1, 1);
    RooRealVar  ChebychevPar2   =   RooRealVar      ("ch2","ch2",   0., -1, 1);
    RooRealVar  ChebychevPar3   =   RooRealVar      ("ch3","ch3",   0., -1, 1);
    RooRealVar  ChebychevPar4   =   RooRealVar      ("ch4","ch4",   0., 0. - kCh4Limits, 0. + kCh4Limits);
    RooRealVar  ChebychevPar5   =   RooRealVar      ("ch5","ch5",   0., 0. - kCh5Limits, 0. + kCh5Limits);

    //Signal
    RooRealVar  VoigtianMass_   =   RooRealVar      ("Vms","Vms",   kParticleMass_, kParticleMass_ - kVmsLimits, kParticleMass_ + kVmsLimits);
    RooRealVar  VoigtianWidth   =   RooRealVar      ("Vwd","Vwd",   kParticleWidth, kParticleWidth - kVwdLimits, kParticleWidth + kVwdLimits);
    RooRealVar  VoigtianSlope   =   RooRealVar      ("Vsl","Vsl",   kDetectorSlope, kDetectorSlope - kVslLimits, kDetectorSlope + kVslLimits);
    
    // Normalisation coefficients
    RooRealVar  SignalMagnit    =   RooRealVar      ("1SS","1SS",   0.5*nEntries, 0., nEntries);
    RooRealVar  BackgroundMg    =   RooRealVar      ("1BB","1BB",   0.5*nEntries, 0., nEntries);
    
    // PDFs
    RooVoigtian     fSignal     =   RooVoigtian     ("fSig","fSig", fVariable,  VoigtianMass_,  VoigtianWidth,  VoigtianSlope);
    RooChebychev    fBkgrnd     =   RooChebychev    ("fBkg","fBkg", fVariable,  RooArgSet( ChebychevPar0, ChebychevPar1, ChebychevPar2, ChebychevPar3, ChebychevPar4, ChebychevPar5 ));
                    fModel_     =   RooAddPdf       ("fMod","fMod", RooArgList( fSignal, fBkgrnd ),RooArgList( SignalMagnit, BackgroundMg ));
    
    return fModel_.fitTo(*fDataHist,Extended(kTRUE),SumW2Error(kTRUE),Save());
}

RooFitResult*   FitModel                        ( TH1F * _h_Data, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
   
    // Defining the Fit Options
    Double_t       *fFitOptions =   fCoreFitModelOptionSelect( fOption );
    
    //------ General Information on the Data Histogram -----//
    
    // Defining the RooFit Variable
    RooRealVar      fVariable   =       RooRealVar  ("fVariable","fVariable",fFitOptions[0],fFitOptions[1]);
    RooDataHist    *fDataHist   =   new RooDataHist ("","",fVariable,Import(*_h_Data));
    RooAddPdf       fModel_;
    
    // Fitting histogram and
    RooFitResult   *fResult     =   fCoreFitModel( fDataHist, fVariable, fModel_, fFitOptions[2], fFitOptions[3], fFitOptions[4], fFitOptions[5], fFitOptions[6] );
    
    // Saving to canvas on file if requested
    if ( fOption.find("S") != -1 )
    {
        Float_t fPTMax      =   fMaxPT1D;
        Float_t fPTMin      =   fMaxPT2D;
        Int_t   fPTindex    =   0;
        if ( fOption.find("PT=")    != -1 ) { fPTindex  =   10*(fOption.at(fOption.find("PT=")+3)-'0')+(fOption.at(fOption.find("PT=")+4)-'0'); }
        if ( fOption.find("1D")     != -1 ) { fPTMin    =   fArrPT1D[fPTindex];  fPTMax    =   fArrPT1D[fPTindex+1]; };
        if ( fOption.find("12D")    != -1 ) { fPTMin    =   fArrPT2D[fPTindex];  fPTMax    =   fArrPT2D[fPTindex+1]; };
        
        hName           =   "DT";
        hTitle          =   Form( "Invariant Mass of Kaons in pT %.2f-%.2f GeV/c", fPTMin, fPTMax );
        if ( bPythiaTest )  { hTitle +=  " (MC)";   hName = "MC"; };
        
        // Canvas to plot
        TCanvas * fSaveToCanvas     =   new TCanvas(Form( "PT_%.1f_%.1f_1D_%s", fPTMin, fPTMax, hName ),
                                                    Form( "PT_%.1f_%.1f_1D_%s", fPTMin, fPTMax, hName ) );
        
        RooPlot * fSaveToFrame      =   fVariable.frame(Name(hName),Title(hTitle));
        TLegend * fLegend           =   new TLegend   (0.12,0.60,0.30,0.85);
        
        fRooPlotMaker(fSaveToFrame,fLegend,fModel_,fDataHist,"InvMass1D");
        
        fSaveToFrame                ->Draw("same");
        fLegend                     ->Draw("same");
        fSaveToCanvas               ->Write ();
        delete fSaveToCanvas;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    return fResult;
}

/*
RooFitResult*   FitModel        ( TH2F * THdata, RooFitResult * fFitShapeX, RooFitResult * fFitShapeY, string fHistName = "", Bool_t fSaveToFile = false, Int_t PTindex = -1, Int_t PTjndex = -1, string fOption = "" )
{
    // Silencing TCanvas Pop-Up
    gROOT->SetBatch();
    
    Bool_t  fBackg  =   false;
    if  ( fOption.find("BK") != -1 )     fBackg  =   true;
    
    Double_t fInvMassValMax, fInvMassValMin;
    SetBoundaries(fOption,fInvMassValMin,fInvMassValMax);
    
    // Global Variables
    Int_t nEntries      = THdata->GetEntries();
    RooRealVar varx     = RooRealVar        ("xInvMass2D","xInvMass2D",fInvMassValMin,fInvMassValMax);
    RooRealVar vary     = RooRealVar        ("yInvMass2D","yInvMass2D",fInvMassValMin,fInvMassValMax);
    RooDataHist* data   = new RooDataHist   (fHistName.c_str(),fHistName.c_str(),RooArgList(varx,vary),Import(*THdata));
    Int_t kNCycle       = 5;
    
    RooArgSet  *utilityx    =   new RooArgSet(fFitShapeX->floatParsFinal(),fFitShapeX->constPars());
    RooArgSet  *utilityy    =   new RooArgSet(fFitShapeY->floatParsFinal(),fFitShapeY->constPars());
    
    
    // Background
    RooRealVar ch0x, ch1x, ch2x, ch3x, ch4x, ch5x, ch0y, ch1y, ch2y, ch3y, ch4y, ch5y;
    
    if ( fBackg )
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0),utilityx->getRealValue("ch0",0)*0.75,utilityx->getRealValue("ch0",0)*1.25);
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0),utilityx->getRealValue("ch1",0)*0.75,utilityx->getRealValue("ch1",0)*1.25);
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0),utilityx->getRealValue("ch2",0)*0.75,utilityx->getRealValue("ch2",0)*1.25);
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0),utilityx->getRealValue("ch3",0)*0.75,utilityx->getRealValue("ch3",0)*1.25);
        ch4x     = RooRealVar ("ch4x","ch4x"     ,utilityx->getRealValue("ch4",0),utilityx->getRealValue("ch4",0)*0.75,utilityx->getRealValue("ch4",0)*1.25);
        ch5x     = RooRealVar ("ch5x","ch5x"     ,utilityx->getRealValue("ch5",0),utilityx->getRealValue("ch5",0)*0.75,utilityx->getRealValue("ch5",0)*1.25);
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0),utilityy->getRealValue("ch0",0)*0.75,utilityy->getRealValue("ch0",0)*1.25);
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0),utilityy->getRealValue("ch1",0)*0.75,utilityy->getRealValue("ch1",0)*1.25);
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0),utilityy->getRealValue("ch2",0)*0.75,utilityy->getRealValue("ch2",0)*1.25);
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0),utilityy->getRealValue("ch3",0)*0.75,utilityy->getRealValue("ch3",0)*1.25);
        ch4y     = RooRealVar ("ch4y","ch4y"     ,utilityy->getRealValue("ch4",0),utilityy->getRealValue("ch4",0)*0.75,utilityy->getRealValue("ch4",0)*1.25);
        ch5y     = RooRealVar ("ch5y","ch5y"     ,utilityy->getRealValue("ch5",0),utilityy->getRealValue("ch5",0)*0.75,utilityy->getRealValue("ch5",0)*1.25);
    }
    else
    {
        ch0x     = RooRealVar ("ch0x","ch0x"     ,utilityx->getRealValue("ch0",0));
        ch1x     = RooRealVar ("ch1x","ch1x"     ,utilityx->getRealValue("ch1",0));
        ch2x     = RooRealVar ("ch2x","ch2x"     ,utilityx->getRealValue("ch2",0));
        ch3x     = RooRealVar ("ch3x","ch3x"     ,utilityx->getRealValue("ch3",0));
        ch4x     = RooRealVar ("ch4x","ch4x"     ,utilityx->getRealValue("ch4",0));
        ch5x     = RooRealVar ("ch5x","ch5x"     ,utilityx->getRealValue("ch5",0));
        ch0y     = RooRealVar ("ch0y","ch0y"     ,utilityy->getRealValue("ch0",0));
        ch1y     = RooRealVar ("ch1y","ch1y"     ,utilityy->getRealValue("ch1",0));
        ch2y     = RooRealVar ("ch2y","ch2y"     ,utilityy->getRealValue("ch2",0));
        ch3y     = RooRealVar ("ch3y","ch3y"     ,utilityy->getRealValue("ch3",0));
        ch4y     = RooRealVar ("ch4y","ch4y"     ,utilityy->getRealValue("ch4",0));
        ch5y     = RooRealVar ("ch5y","ch5y"     ,utilityy->getRealValue("ch5",0));
    }
    
    //Signal
    RooRealVar pMassx   = RooRealVar ("pMassx","pMassx" ,utilityx->getRealValue("sMass",0));
    RooRealVar pWidthx  = RooRealVar ("pWidtx","pWidtx" ,utilityx->getRealValue("sWidt",0));
    RooRealVar pSlopex  = RooRealVar ("pSlopx","pSlopx" ,utilityx->getRealValue("sSlop",0));
    RooRealVar pMassy   = RooRealVar ("pMassy","pMassy" ,utilityy->getRealValue("sMass",0));
    RooRealVar pWidthy  = RooRealVar ("pWidty","pWidty" ,utilityy->getRealValue("sWidt",0));
    RooRealVar pSlopey  = RooRealVar ("pSlopy","pSlopy" ,utilityy->getRealValue("sSlop",0));
    
    // Coefficients
    auto fS__x  =   utilityx->getRealValue("1nSS",0);
    auto fB__x  =   utilityx->getRealValue("1nBB",0);
    auto fS__y  =   utilityy->getRealValue("1nSS",0);
    auto fB__y  =   utilityy->getRealValue("1nBB",0);
    auto fTot_  =   fS__x*fS__y+fB__x*fB__y+fS__x*fB__y+fB__x*fS__y;
    RooRealVar n0       = RooRealVar ("anSS2D","anSS2D" ,nEntries*(fS__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n1       = RooRealVar ("anBB2D","anBB2D" ,nEntries*(fB__x*fB__y)/fTot_,0.,nEntries);
    RooRealVar n2       = RooRealVar ("anBS2D","anBS2D" ,nEntries*(fB__x*fS__y)/fTot_,0.,nEntries);
    RooRealVar n3       = RooRealVar ("anSB2D","anSB2D" ,nEntries*(fS__x*fB__y)/fTot_,0.,nEntries);
    
    // PDFs
    RooChebychev        fBkgx ("fBkgx","fBkgx"          ,varx,RooArgSet(ch0x,ch1x,ch2x,ch3x,ch4x,ch5x));
    RooVoigtian         fSigx ("fSigx","fSigx"          ,varx,pMassx,pWidthx,pSlopex);
    RooChebychev        fBkgy ("fBkgy","fBkgy"          ,vary,RooArgSet(ch0y,ch1y,ch2y,ch3y,ch4y,ch5x));
    RooVoigtian         fSigy ("fSigy","fSigy"          ,vary,pMassy,pWidthy,pSlopey);
    RooProdPdf          fBB   ("fBkg","fBkg"            ,fBkgx,fBkgy);
    RooProdPdf          fSB   ("fSigBkg","fSBWBkg"      ,fSigx,fBkgy);
    RooProdPdf          fBS   ("fBkgSig","fBkgSig"      ,fBkgx,fSigy);
    RooProdPdf          fSS   ("fSigSig","fSigSig"      ,fSigx,fSigy);
    RooAddPdf           fMod  ("fMod2D","fMod2D"        ,RooArgList(fBB,fSS,fSB,fBS),RooArgList(n1,n0,n3,n2));
    

    RooFitResult* FitResults;
    for ( Int_t iCycle = 0; iCycle < kNCycle; iCycle++ )
    {
       FitResults = fMod.fitTo(*data,Extended(kTRUE),SumW2Error(kTRUE),Save());
    }
    
    // Save to file
    if ( fSaveToFile )
    {
        int         nBinsPrint      =   3;
        double      dIncrement      =   (fMaxIM2D-fMinIM2D)/nBinsPrint;
        TLatex*     latext          =   new TLatex();
        TCanvas*    cTotal          =   new TCanvas("","",0,45,1440,855);
                    cTotal          ->  SetTitle(Form("Slices of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]));
                    cTotal          ->  SetName(Form("PT_%.1f_%.1f__%.1f_%.1f_%s",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()));
                    cTotal          ->  Divide(2,nBinsPrint);
        
                            varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                            vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()),
                                                    Form("xInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame  =   vary.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);

                            varx.setRange("fDrawRange",fMinIM2D+i*dIncrement,fMinIM2D+(i+1)*dIncrement);
                            vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);

            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"yInvMass2D");
            
            cTotal->cd( i+1 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.3f < m^{x}_{K^{+}K^{-}} < %.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.3f < m^{x}_{K^{+}K^{-}} < %.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                                        vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
        for (int i = 0; i < nBinsPrint; i++)
        {
            hName                       = "Slice of 2D Invariant Mass of Kaons";
            hTitle                      = "Slice of 2D Invariant Mass of Kaons";
            if ( PTindex != -1 && !bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            if ( PTindex != -1 &&  bPythiaTest ) hTitle = Form("Slice of 2D Invariant Mass of Kaons in pT %.1f-%.1f GeV, %.1f-%.1f GeV for MC",fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1]);
            
            TCanvas * fSaveToCanvas =   new TCanvas(
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f_%s",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1],fHistName.c_str()),
                                                    Form("yInvMass_%.3f_%.3f_PTx_%.3f_%.3f_PTy_%.3f_%.3f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1),fArrPT2D[PTindex],fArrPT2D[PTindex+1],fArrPT2D[PTjndex],fArrPT2D[PTjndex+1])
                                                    );
            
            RooPlot * fSaveToFrame      =   varx.frame(Name(hName),Title(hTitle));
            TLegend * fLegend           = new TLegend   (0.12,0.60,0.30,0.85);
            
                                        varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                                        vary.setRange("fDrawRange",fMinIM2D+i*dIncrement,fMinIM2D+(i+1)*dIncrement);
                                                                            
            fRooPlotMaker(fSaveToFrame,fLegend,fMod,data,"xInvMass2D");
            
            cTotal->cd( i+1 +3 );
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.2f < m^{y}_{K^{+}K^{-}} < %.2f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            
            fSaveToCanvas->cd();
            fSaveToFrame                ->Draw("same");
            fLegend                     ->Draw("same");
            latext                      ->DrawLatexNDC(0.6, 0.85, Form("%.2f < m^{y}_{K^{+}K^{-}} < %.2f",fMinIM2D+dIncrement*i,fMinIM2D+dIncrement*(i+1)));
            fSaveToCanvas               ->Write();
            delete fSaveToCanvas;
        }
                                        varx.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
                                        vary.setRange("fDrawRange",fMinIM2D,fMaxIM2D);
        cTotal ->Write();
        delete cTotal;
    }
    
    // Un-Silencing TCanvas Pop-Up
    gROOT->SetBatch(false);
    
    // Fit
    return FitResults;
}
*/

#endif
