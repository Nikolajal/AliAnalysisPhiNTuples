// Global Functions and data structures file
// !TODO: General housekeeping and new manage rules
// !TODO: Headers clean-up
// !TODO: 2D plot options
// !TODO: Transfer to RooFit the Fit PT Count
// !TODO: Legend as in ...Count.C the slice plots

#ifndef SETFUNCTIONS_H
#define SETFUNCTIONS_H
#include "SetValues.h"
#include "SpectraUtils.C"
#include "Util_Functions.cpp"
#include "Util_Histograms.cpp"

// C++
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <chrono>

// ROOT
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"

// RooFit
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"

// RooFitFunction
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooUniform.h"

using namespace std;
using namespace RooFit;

//--------------------------------------------//
//  Utilities for histogram creations         //
//--------------------------------------------//


//--------------------------------------------//
//  Signal Extraction from Invariant Mass hst //
//--------------------------------------------//

enum            fitresults1D
{
    Background, Signal
};

enum            fitresults2D
{

    BackgBackg, BackgSignl, SignlBackg, SignlSignl
};

//--------------------------------------------//
//  Final Yield Extraction and extrapolation  //
//--------------------------------------------//

 // Utility

void                SetLevyTsalPT1D ( )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,6.0,7.4);
    fLevyFit1D  ->  SetParameter(1,6.7); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.25,.30);
    fLevyFit1D  ->  SetParameter(2,.272); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,0.028,0.036);
    fLevyFit1D  ->  SetParameter(3,0.032);
}
 
void                SetLevyTsalPT2D ( )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,3.5,7.5);
    fLevyFit1D  ->  SetParameter(1,4.5); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.18,.45);
    fLevyFit1D  ->  SetParameter(2,.272); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,0.5e-6,1.e-3);
    fLevyFit1D  ->  SetParameter(3,1.e-6);
}

void                SetLevyTsalis   ( )
{
    // - // Setting up Fit parameters
    
    // Mass
    fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
    fLevyFit1D  ->  SetParameter(0,kPMas);
    
    // n-Parameter
    fLevyFit1D  ->  SetParLimits(1,2.1,7.5);
    fLevyFit1D  ->  SetParameter(1,6.); // 6.7
    
    // T-Parameter
    fLevyFit1D  ->  SetParLimits(2,.21,.750);
    fLevyFit1D  ->  SetParameter(2,.272); // .272
    
    // dN/dy
    fLevyFit1D  ->  SetParLimits(3,1.e-7,1.e-1);
    fLevyFit1D  ->  SetParameter(3,0.032);
    
    if ( bPythiaTest )
    {
        // Mass
        fLevyFit1D  ->  SetParLimits(0,kPMas,kPMas);
        fLevyFit1D  ->  SetParameter(0,kPMas);
           
        // n-Parameter
        fLevyFit1D  ->  SetParLimits(1,2.1,7.5);
        fLevyFit1D  ->  SetParameter(1,4.); // 6.7
           
        // T-Parameter
        fLevyFit1D  ->  SetParLimits(2,.01,.750);
        fLevyFit1D  ->  SetParameter(2,.21); // .272
           
        // dN/dy
        fLevyFit1D  ->  SetParLimits(3,1.e-9,1.);
        fLevyFit1D  ->  SetParameter(3,.04);
    }
}

 // Measurements

Double_t *          EvalStErInteg   ( TH1F * aTarget,   string aName__ = "",  Bool_t aSaveFit = false, Int_t aIntBin = -1  )
{
    Int_t       fPoints =   1.e3;
    Double_t *  fReturn = new Double_t [2];
    Double_t *  fUtilPt = new Double_t [fPoints];
    Double_t *  fUtilMn = new Double_t [fPoints];
    Double_t *  fUtilP1 = new Double_t [fPoints];
    Double_t *  fUtilP2 = new Double_t [fPoints];
    Double_t *  fUtilP3 = new Double_t [fPoints];
    Double_t *  fUtilW0 = new Double_t [fPoints];
    
    // Speeding multiple fits
    gROOT->SetBatch();
    
    TH1F *  fRandPt,    *fTotal_;
    for ( Int_t iTer = 0; iTer < fPoints; iTer++ )
    {
        fRandPt =   SetRandPoints(aTarget);
        fTotal_ =   SetSystErrorsh(fRandPt);
        SetLevyTsalis();
        fTotal_->Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        if ( bPythiaTest ) fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,1.6);
        fUtilPt[iTer]   =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W x^1");
        fUtilMn[iTer]   =   fLevyFit1D->Integral(0.,0.4);
        fUtilP1[iTer]   =   fLevyFit1D->GetParameter(1);
        fUtilP2[iTer]   =   fLevyFit1D->GetParameter(2);
        fUtilP3[iTer]   =   fLevyFit1D->GetParameter(3);
        fUtilW0[iTer]   =   1;
        if ( aSaveFit )
        {
            TCanvas * c1    =   new TCanvas();
            gPad->SetLogy();
            fTotal_     ->Draw();
            fLevyFit1D  ->Draw("same");
            c1          ->Write();
            delete  c1;
        }
    }

    // Speeding multiple fits
    gROOT->SetBatch(false);

    TVectorD fUtilVP    (fPoints,fUtilPt);
    TVectorD fUtilVM    (fPoints,fUtilMn);
    TVectorD fUtilV1    (fPoints,fUtilP1);
    TVectorD fUtilV2    (fPoints,fUtilP2);
    TVectorD fUtilV3    (fPoints,fUtilP3);
    
    TH1F *  fUtilHP     =   new TH1F    (Form("gp_%s",aName__.c_str()),"",25,fUtilVP.Min(),fUtilVP.Max());
    fUtilHP             ->  FillN       (fPoints,fUtilPt,fUtilW0);
    fUtilHP             ->  Write();
    
    TH1F *  fUtilHM     =   new TH1F    (Form("gm_%s",aName__.c_str()),"",25,fUtilVM.Min(),fUtilVM.Max());
    fUtilHM             ->  FillN       (fPoints,fUtilMn,fUtilW0);
    fUtilHM             ->  Write();
    
    TH1F *  fUtilH1     =   new TH1F    (Form("P1_%s",aName__.c_str()),"",25,fUtilV1.Min(),fUtilV1.Max());
    fUtilH1             ->  FillN       (fPoints,fUtilP1,fUtilW0);
    fUtilH1             ->  Write();
    
    TH1F *  fUtilH2     =   new TH1F    (Form("P2_%s",aName__.c_str()),"",25,fUtilV2.Min(),fUtilV2.Max());
    fUtilH2             ->  FillN       (fPoints,fUtilP2,fUtilW0);
    fUtilH2             ->  Write();
    
    TH1F *  fUtilH3     =   new TH1F    (Form("P3_%s",aName__.c_str()),"",25,fUtilV3.Min(),fUtilV3.Max());
    fUtilH3             ->  FillN       (fPoints,fUtilP3,fUtilW0);
    fUtilH3             ->  Write();
    
    fUtilHM             ->  Fit ("gaus","IMREQ0S");
    fUtilHP             ->  Fit ("gaus","IMREQ0S");
    fReturn[0]  =   fUtilHM      ->  GetRMS();//GetFunction ("gaus")->GetParameter(2);
    fReturn[1]  =   fUtilHP      ->  GetRMS();//GetFunction ("gaus")->GetParameter(2);
    return fReturn;
}

Double_t *          EvalStErInteg   ( TH1D * aTarget,   string aName__ = "",  Bool_t aSaveFit = false, Int_t aIntBin = -1  )
{
    Int_t       fPoints =   1.e2;
    Double_t *  fReturn = new Double_t [2];
    Double_t *  fUtilPt = new Double_t [fPoints];
    Double_t *  fUtilMn = new Double_t [fPoints];
    Double_t *  fUtilP1 = new Double_t [fPoints];
    Double_t *  fUtilP2 = new Double_t [fPoints];
    Double_t *  fUtilP3 = new Double_t [fPoints];
    Double_t *  fUtilW0 = new Double_t [fPoints];
    
    Double_t kMinInt, kMaxInt;
    if ( aIntBin == 12 )
    {
        kMinInt = 0.4;
        kMaxInt = 0.68;
    }
    
    // Speeding multiple fits
    gROOT->SetBatch();
    
    TH1D *  fRandPt,    *fTotal_;
    for ( Int_t iTer = 0; iTer < fPoints; iTer++ )
    {
        fRandPt =   SetRandPoints(aTarget);
        fTotal_ =   SetSystErrorsh(fRandPt);
        SetLevyTsalis();
        fTotal_->Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        if ( bPythiaTest ) fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,1.6);
        fUtilPt[iTer]   =   fLevyFit1D->Moment(1,0.0,0.4) + fLevyFit1D->Moment(1,kMinInt,kMaxInt);
        fUtilMn[iTer]   =   fLevyFit1D->Integral(0.0,0.4) + fLevyFit1D->Integral(kMinInt,kMaxInt);
        fUtilP1[iTer]   =   fLevyFit1D->GetParameter(1);
        fUtilP2[iTer]   =   fLevyFit1D->GetParameter(2);
        fUtilP3[iTer]   =   fLevyFit1D->GetParameter(3);
        fUtilW0[iTer]   =   1;
        
        if ( aSaveFit )
        {
            TCanvas * c1    =   new TCanvas();
            gPad->SetLogy();
            fTotal_     ->Draw();
            fLevyFit1D  ->Draw("same");
            c1          ->Write();
            delete  c1;
        }
    }

    // Speeding multiple fits
    gROOT->SetBatch(false);

    TVectorD fUtilVP    (fPoints,fUtilPt);
    TVectorD fUtilVM    (fPoints,fUtilMn);
    TVectorD fUtilV1    (fPoints,fUtilP1);
    TVectorD fUtilV2    (fPoints,fUtilP2);
    TVectorD fUtilV3    (fPoints,fUtilP3);
    TVectorD fUtilV0    (fPoints,fUtilW0);
    
    TH1F *  fUtilHP     =   new TH1F    (Form("gp_%s",aName__.c_str()),"",25,fUtilVP.Min(),fUtilVP.Max());
    fUtilHP             ->  FillN       (fPoints,fUtilPt,fUtilW0);
    fUtilHP             ->  Write();
    
    TH1F *  fUtilHM     =   new TH1F    (Form("gm_%s",aName__.c_str()),"",25,fUtilVM.Min(),fUtilVM.Max());
    fUtilHM             ->  FillN       (fPoints,fUtilMn,fUtilW0);
    fUtilHM             ->  Write();
    
    TH1F *  fUtilH1     =   new TH1F    (Form("P1_%s",aName__.c_str()),"",25,fUtilV1.Min(),fUtilV1.Max());
    fUtilH1             ->  FillN       (fPoints,fUtilP1,fUtilW0);
    fUtilH1             ->  Write();
    
    TH1F *  fUtilH2     =   new TH1F    (Form("P2_%s",aName__.c_str()),"",25,fUtilV2.Min(),fUtilV2.Max());
    fUtilH2             ->  FillN       (fPoints,fUtilP2,fUtilW0);
    fUtilH2             ->  Write();
    
    TH1F *  fUtilH3     =   new TH1F    (Form("P3_%s",aName__.c_str()),"",25,fUtilV3.Min(),fUtilV3.Max());
    fUtilH3             ->  FillN       (fPoints,fUtilP3,fUtilW0);
    fUtilH3             ->  Write();
    
    fUtilHM             ->  Fit ("gaus","IMREQ0S");
    fUtilHP             ->  Fit ("gaus","IMREQ0S");
    fReturn[0]  =   fUtilHM      ->  GetFunction ("gaus")->GetParameter(2);
    fReturn[1]  =   fUtilHP      ->  GetFunction ("gaus")->GetParameter(2);
    return fReturn;
}

Double_t *          EvalStErInteg   ( TH1D * aTarget,  TH1D * aTarge2,   string aName__ = "",  Bool_t aSaveFit = false, Int_t aIntBin = -1  )
{
    Int_t       fPoints =   1.e2;
    Double_t *  fReturn = new Double_t [2];
    Double_t *  fUtilPt = new Double_t [fPoints];
    Double_t *  fUtilMn = new Double_t [fPoints];
    Double_t *  fUtilP1 = new Double_t [fPoints];
    Double_t *  fUtilP2 = new Double_t [fPoints];
    Double_t *  fUtilP3 = new Double_t [fPoints];
    Double_t *  fUtilW0 = new Double_t [fPoints];
    
    Double_t kMinInt, kMaxInt;
    if ( aIntBin == 12 )
    {
        kMinInt = 0.4;
        kMaxInt = 0.68;
    }
    
    // Speeding multiple fits
    gROOT->SetBatch();
    
    TH1D *  fRandPt,    *fTotal_;
    for ( Int_t iTer = 0; iTer < fPoints; iTer++ )
    {
        fRandPt =   SetRandPoints(aTarget);
        fTotal_ =   SetSystErrorsh(fRandPt,aTarge2);
        SetLevyTsalis();
        fTotal_->Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        if ( bPythiaTest ) fRandPt->Fit(fLevyFit1D,"IMREQ0S","",0.4,1.6);
        fUtilPt[iTer]   =   fLevyFit1D->Moment(1,0.0,0.4) + fLevyFit1D->Moment(1,kMinInt,kMaxInt);
        fUtilMn[iTer]   =   fLevyFit1D->Integral(0.0,0.4) + fLevyFit1D->Integral(kMinInt,kMaxInt);
        fUtilP1[iTer]   =   fLevyFit1D->GetParameter(1);
        fUtilP2[iTer]   =   fLevyFit1D->GetParameter(2);
        fUtilP3[iTer]   =   fLevyFit1D->GetParameter(3);
        fUtilW0[iTer]   =   1;
        
        if ( aSaveFit )
        {
            TCanvas * c1    =   new TCanvas();
            gPad->SetLogy();
            fTotal_     ->Draw();
            fLevyFit1D  ->Draw("same");
            c1          ->Write();
            delete  c1;
        }
    }

    // Speeding multiple fits
    gROOT->SetBatch(false);

    TVectorD fUtilVP    (fPoints,fUtilPt);
    TVectorD fUtilVM    (fPoints,fUtilMn);
    TVectorD fUtilV1    (fPoints,fUtilP1);
    TVectorD fUtilV2    (fPoints,fUtilP2);
    TVectorD fUtilV3    (fPoints,fUtilP3);
    TVectorD fUtilV0    (fPoints,fUtilW0);
    
    TH1F *  fUtilHP     =   new TH1F    (Form("gp_%s",aName__.c_str()),"",25,fUtilVP.Min(),fUtilVP.Max());
    fUtilHP             ->  FillN       (fPoints,fUtilPt,fUtilW0);
    fUtilHP             ->  Write();
    
    TH1F *  fUtilHM     =   new TH1F    (Form("gm_%s",aName__.c_str()),"",25,fUtilVM.Min(),fUtilVM.Max());
    fUtilHM             ->  FillN       (fPoints,fUtilMn,fUtilW0);
    fUtilHM             ->  Write();
    
    TH1F *  fUtilH1     =   new TH1F    (Form("P1_%s",aName__.c_str()),"",25,fUtilV1.Min(),fUtilV1.Max());
    fUtilH1             ->  FillN       (fPoints,fUtilP1,fUtilW0);
    fUtilH1             ->  Write();
    
    TH1F *  fUtilH2     =   new TH1F    (Form("P2_%s",aName__.c_str()),"",25,fUtilV2.Min(),fUtilV2.Max());
    fUtilH2             ->  FillN       (fPoints,fUtilP2,fUtilW0);
    fUtilH2             ->  Write();
    
    TH1F *  fUtilH3     =   new TH1F    (Form("P3_%s",aName__.c_str()),"",25,fUtilV3.Min(),fUtilV3.Max());
    fUtilH3             ->  FillN       (fPoints,fUtilP3,fUtilW0);
    fUtilH3             ->  Write();
    
    fUtilHM             ->  Fit ("gaus","IMREQ0S");
    fUtilHP             ->  Fit ("gaus","IMREQ0S");
    fReturn[0]  =   fUtilHM      ->  GetFunction ("gaus")->GetParameter(2);
    fReturn[1]  =   fUtilHP      ->  GetFunction ("gaus")->GetParameter(2);
    return fReturn;
}

Double_t *          ExtrapolateVl   ( TH1F * aTarget,   string aName__ = "", Bool_t aSaveFit = false,  Int_t aIntBin = -1  )
{
    gROOT->SetBatch(true);
    
    Double_t *  fReturn =   new Double_t    [6];    // Result of the Process
                                                    // 1. Mean Value    2. Stat Err     3. Syst Err     4. Mean PT
    // Prepping the FIT
    SetLevyTsalis();
    TH1F   *hTargetSyst = new TH1F (*SetSystErrorsh(aTarget));
    TH1F   *aTarge2 = new TH1F (*SetSystErrors2(aTarget));
    hTargetSyst                     ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,10.);
    if ( bPythiaTest ) hTargetSyst  ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,1.6);
    if ( true )
    {
        aTarget     ->Write();
        fLevyFit1D  ->Write();
        TCanvas * fSaveToCanvas = new TCanvas(aName__.c_str());
        gPad->SetLogy();
        gStyle->SetOptStat(0);
        TLegend * lLegend1 = new TLegend(0.6,0.75,0.9,0.9);
        hTargetSyst ->SetMarkerStyle(22);
        hTargetSyst ->SetMarkerColor(38);
        hTargetSyst->SetTitle("");
        hTargetSyst->GetYaxis()->SetTitle("#frac{d^{2}N#phi}{dydp_{T}#phi}(GeV/c)^{-1}");
        hTargetSyst->GetXaxis()->SetTitle("P_{T}#phi (GeV/c)");
        hTargetSyst     ->Draw();
        if ( bPythiaTest )  fLevyFit1D->SetRange(0.,2.);
        fLevyFit1D  ->Draw("SAME");
        lLegend1->AddEntry(hTargetSyst,"Data","EP");
        lLegend1->AddEntry(fLevyFit1D,"Fit","L");
        lLegend1->Draw("same");
        fSaveToCanvas->Write();
        fSaveToCanvas->SaveAs(Form("./graphs/FITLEVY_%s.pdf",aName__.c_str()));
        fSaveToCanvas->SaveAs(Form("./graphs/FITLEVY_%s.png",aName__.c_str()));
        delete fSaveToCanvas;
    }
    Double_t    fHIntWidth, fHIntWidtE, fHIntWidSh, fHIntWidSE;
    fHIntWidth  =   aTarget->IntegralAndError(-1,100,fHIntWidtE,"width");
    fHIntWidSh  =   aTarge2->IntegralAndError(-1,100,fHIntWidSE,"width");
    Double_t    fHIntWidx1      =   HistIntegrals( aTarget,     "W x^1" );
    Double_t    fHIntWidxE      =   HistIntegrals( aTarget,     "WE x^1");
    Double_t    fFIntWidth      =   fLevyFit1D->Integral(0.,0.4); //FuncIntegrals( fLevyFit1D,  0.,0.4,"W");
    Double_t    fFIntWidtE      =   fLevyFit1D->IntegralError(0.,0.4);
    Double_t    fFIntWidx1      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W x^1");
    Double_t *  fStatErr        =   EvalStErInteg(aTarget,aName__,false,aIntBin);
    
    // Mean Value
    fReturn[0]  =   fHIntWidth + fFIntWidth;
    
    //Stat Error
    fReturn[1]  =   sqrt( fStatErr[0]*fStatErr[0] + fHIntWidtE*fHIntWidtE);
    
    //Syst Error +
    fReturn[2]  =   fReturn[0]*kSystematical1D_;
    
    // Mean PT
    fReturn[3]  =   (fHIntWidx1+fFIntWidx1)/(fHIntWidth+fFIntWidth);
    
    //Stat Error
    auto fValError1 =   fHIntWidx1*fHIntWidx1*fStatErr[1]*fStatErr[1];
    auto fValError2 =   fFIntWidx1*fFIntWidx1*fHIntWidxE*fHIntWidxE;
    auto fValError3 =   fReturn[3]*fReturn[3]*fStatErr[0]*fStatErr[0];
    auto fValError4 =   fReturn[3]*fReturn[3]*fHIntWidtE*fHIntWidtE;
    fReturn[4]  =   (1/(fHIntWidth+fFIntWidth))*sqrt(fValError1+fValError2+fValError3+fValError4);
    
    //Syst Error
    fReturn[5]  =   0;
    
    gROOT->SetBatch(false);
    return fReturn;
    
    gROOT->SetBatch(false);
    return fReturn;
}

Double_t *          ExtrapolateVl   ( TH1D * aTarget,   string aName__ = "", Bool_t aSaveFit = false, Int_t aIntBin = -1  )
{
    gROOT->SetBatch(true);
    
    Double_t *  fReturn =   new Double_t    [6];    // Result of the Process
                                                    // 1. Mean Value    2. Stat Err     3. Syst Err     4. Mean PT
    // Prepping the FIT
    SetLevyTsalis();
    TH1D   *hTargetSyst = new TH1D (*SetSystErrorsh(aTarget));
    TH1D   *aTarge2 = new TH1D (*SetSystErrors2(aTarget));
    hTargetSyst                     ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,10.);
    if ( bPythiaTest ) hTargetSyst  ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,1.6);
    if ( true )
    {
        aTarget     ->Write();
        fLevyFit1D  ->Write();
        TCanvas * fSaveToCanvas = new TCanvas(aName__.c_str());
        gPad->SetLogy();
        TLegend * lLegend1 = new TLegend(0.6,0.75,0.9,0.9);
        hTargetSyst ->SetMarkerStyle(22);
        hTargetSyst ->SetMarkerColor(38);
        hTargetSyst->SetTitle("");
        hTargetSyst->GetYaxis()->SetTitle("#frac{d^{3}N#phi#phi}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
        hTargetSyst->GetXaxis()->SetTitle("P_{T}#phi_{1} (GeV/c)");
        hTargetSyst     ->Draw();
        if ( bPythiaTest )  fLevyFit1D->SetRange(0.,2.);
        fLevyFit1D  ->Draw("SAME");
        lLegend1->AddEntry(hTargetSyst,"Data","EP");
        lLegend1->AddEntry(fLevyFit1D,"Fit","L");
        lLegend1->Draw("same");
        fSaveToCanvas->Write();
        fSaveToCanvas->SaveAs(Form("./graphs/FITLEVY_%s.pdf",aName__.c_str()));
        fSaveToCanvas->SaveAs(Form("./graphs/FITLEVY_%s.png",aName__.c_str()));
        delete fSaveToCanvas;
    }
    Double_t    fHIntWidth, fHIntWidtE, fHIntWidSh, fHIntWidSE;
    fHIntWidth  =   aTarget->IntegralAndError(-1,100,fHIntWidtE,"width");
    fHIntWidSh  =   aTarge2->IntegralAndError(-1,100,fHIntWidSE,"width");
    Double_t    fHIntWidx1      =   HistIntegrals( aTarget,     "W x^1" );
    Double_t    fHIntWidxE      =   HistIntegrals( aTarget,     "WE x^1");
    Double_t    fFIntWidth      =   fLevyFit1D->Integral(0.,0.4); //FuncIntegrals( fLevyFit1D,  0.,0.4,"W");
    Double_t    fFIntWidtE      =   fLevyFit1D->IntegralError(0.,0.4);
    Double_t    fFIntWidx1      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W x^1");
    if ( aIntBin == 12 && !bPythiaTest )
    {
        fFIntWidth      =   fLevyFit1D->Integral(0.,0.68);
        fFIntWidtE      =   fLevyFit1D->IntegralError(0.,0.68);
        fFIntWidx1      =   FuncIntegrals( fLevyFit1D,  0.,0.68,"W x^1");
    }
    if ( aIntBin == 3 && !bPythiaTest )
    {
        fFIntWidth     +=   fLevyFit1D->Integral(5.,10.);
        fFIntWidtE     +=   fLevyFit1D->IntegralError(5.,10.);
        fFIntWidx1     +=   FuncIntegrals( fLevyFit1D,  5.,10.,"W x^1");
    }
    Double_t *  fStatErr        =   EvalStErInteg(aTarget,aName__,false,aIntBin);
    
    // Mean Value
    fReturn[0]  =   fHIntWidth + fFIntWidth;
    
    //Stat Error
    fReturn[1]  =   sqrt( fStatErr[0]*fStatErr[0] + fHIntWidtE*fHIntWidtE);
    
    //Syst Error +
    fReturn[2]  =   fReturn[0]*kSystematical2D_;
    
    // Mean PT
    fReturn[3]  =   (fHIntWidx1+fFIntWidx1)/(fHIntWidth+fFIntWidth);
    
    //Stat Error
    auto fValError1 =   fHIntWidx1*fHIntWidx1*fStatErr[1]*fStatErr[1];
    auto fValError2 =   fFIntWidx1*fFIntWidx1*fHIntWidxE*fHIntWidxE;
    auto fValError3 =   fReturn[3]*fReturn[3]*fStatErr[0]*fStatErr[0];
    auto fValError4 =   fReturn[3]*fReturn[3]*fHIntWidtE*fHIntWidtE;
    fReturn[4]  =   (1/(fHIntWidth+fFIntWidth))*sqrt(fValError1+fValError2+fValError3+fValError4);
    
    //Syst Error
    fReturn[5]  =   0;
    
    gROOT->SetBatch(false);
    return fReturn;
}

Double_t *          ExtrapolateVl   ( TH1D * aTarget,   TH1D * aTarge2, string aName__ = "", Bool_t aSaveFit = false, Int_t aIntBin = -1  )
{
    gROOT->SetBatch(true);
    
    Double_t *  fReturn =   new Double_t    [6];    // Result of the Process
                                                    // 1. Mean Value    2. Stat Err     3. Syst Err     4. Mean PT
    // Prepping the FIT
    SetLevyTsalis();
    TH1D   *hTargetSyst = new TH1D ("Total","Total",nBinPT2D,fArrPT2D);
    for ( Int_t i = 3; i <= nBinPT1D; i++ )
    {
        hTargetSyst->SetBinContent(i,aTarget->GetBinContent(i));
        hTargetSyst->SetBinError(i,sqrt(aTarget->GetBinError(i)*aTarget->GetBinError(i)+aTarge2->GetBinError(i)*aTarge2->GetBinError(i)));
    }
    hTargetSyst                     ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,10.);
    if ( bPythiaTest ) hTargetSyst  ->  Fit         (fLevyFit1D,"IMREQ0S","",0.4,1.6);
    if ( true )
    {
        aTarget     ->Write();
        fLevyFit1D  ->Write();
        TCanvas * fSaveToCanvas = new TCanvas(aName__.c_str());
        gPad->SetLogy();
        TLegend * lLegend1 = new TLegend(0.6,0.75,0.9,0.9);
        hTargetSyst ->SetMarkerStyle(22);
        hTargetSyst ->SetMarkerColor(38);
        hTargetSyst->SetTitle("");
        hTargetSyst->GetYaxis()->SetTitle("#frac{d^{2}N#phi#phi}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
        hTargetSyst->GetXaxis()->SetTitle("P_{T}#phi_{2} (GeV/c)");
        hTargetSyst     ->Draw();
        if ( bPythiaTest )  fLevyFit1D->SetRange(0.,2.);
        fLevyFit1D  ->Draw("SAME");
        lLegend1->AddEntry(hTargetSyst,"Data","EP");
        lLegend1->AddEntry(fLevyFit1D,"Fit","L");
        lLegend1->Draw("same");
        fSaveToCanvas->Write();
        fSaveToCanvas->SaveAs(Form("./graphs/FITLEVY_%s.pdf",aName__.c_str()));
        fSaveToCanvas->SaveAs(Form("./graphs/FITLEVY_%s.png",aName__.c_str()));
        delete fSaveToCanvas;
    }
    Double_t    fHIntWidth, fHIntWidtE, fHIntWidSh, fHIntWidSE;
    fHIntWidth  =   aTarget->IntegralAndError(-1,100,fHIntWidtE,"width");
    fHIntWidSh  =   aTarge2->IntegralAndError(-1,100,fHIntWidSE,"width");
    Double_t    fHIntWidx1      =   HistIntegrals( aTarget,     "W x^1" );
    Double_t    fHIntWidxE      =   HistIntegrals( aTarget,     "WE x^1");
    Double_t    fFIntWidth      =   fLevyFit1D->Integral(0.,0.4); //FuncIntegrals( fLevyFit1D,  0.,0.4,"W");
    Double_t    fFIntWidtE      =   fLevyFit1D->IntegralError(0.,0.4);
    Double_t    fFIntWidx1      =   FuncIntegrals( fLevyFit1D,  0.,0.4,"W x^1");
    Double_t *  fStatErr        =   EvalStErInteg(aTarget,aTarge2,aName__,false,aIntBin);
    
    // Mean Value
    fReturn[0]  =   fHIntWidth + fFIntWidth;
    
    //Stat Error
    fReturn[1]  =   sqrt( fStatErr[0]*fStatErr[0] + fHIntWidtE*fHIntWidtE);
    
    //Syst Error +
    fReturn[2]  =   fHIntWidSE + fFIntWidtE;
    
    // Mean PT
    fReturn[3]  =   (fHIntWidx1+fFIntWidx1)/(fHIntWidth+fFIntWidth);
    
    //Stat Error
    auto fValError1 =   fHIntWidx1*fHIntWidx1*fStatErr[1]*fStatErr[1];
    auto fValError2 =   fFIntWidx1*fFIntWidx1*fHIntWidxE*fHIntWidxE;
    auto fValError3 =   fReturn[3]*fReturn[3]*fStatErr[0]*fStatErr[0];
    auto fValError4 =   fReturn[3]*fReturn[3]*fHIntWidtE*fHIntWidtE;
    fReturn[4]  =   (1/(fHIntWidth+fFIntWidth))*sqrt(fValError1+fValError2+fValError3+fValError4);
    
    //Syst Error
    fReturn[5]  =   0;
    
    gROOT->SetBatch(false);
    return fReturn;
}

Double_t ***        ExtrapolateVl   ( TH2F * aTarget )
{
    gROOT->SetBatch();
    Double_t***     fReturn =   new Double_t**  [nBinPT2D+1];
    TF1  **         fLevyMm =   new TF1 *       [nBinPT2D+1];
    for ( Int_t iTer = 0; iTer < nBinPT2D+1; iTer++ )
    {
        fReturn[iTer]       =   new Double_t*   [2];
        fReturn[iTer][0]    =   new Double_t    [6];
        fReturn[iTer][1]    =   new Double_t    [6];
        fLevyMm[iTer]       =   new TF1         [2];
    }
    
    TH1D *      hSliceFX    =   new TH1D("hSliceFX_","hSliceFX_",nBinPT2D,fArrPT2D);
    TH1D *      hSliceFY    =   new TH1D("hSliceFY_","hSliceFY_",nBinPT2D,fArrPT2D);
    TH1D *      hSlice2X    =   new TH1D("hSlice2X_","hSlice2X_",nBinPT2D,fArrPT2D);
    TH1D *      hSlice2Y    =   new TH1D("hSlice2Y_","hSlice2Y_",nBinPT2D,fArrPT2D);
    
    TCanvas *   fDrawAllX   =   new TCanvas("cDrawAllX","cDrawAllX");
    fDrawAllX               ->  Divide(2,5);
    TCanvas *   fDrawAllY   =   new TCanvas("cDrawAllY","cDrawAllY");
    fDrawAllY               ->  Divide(2,5);
    TCanvas *   fDrawFull   =   new TCanvas("fDrawFull","fDrawFull");
    fDrawFull               ->  Divide(2,1);
    
    Int_t iHisto = 0;
    for ( Int_t iFit = 1; iFit <= nBinPT2D-2; iFit++ )
    {
        cout << "X-" << iFit << endl;
        // X-Projection
        fDrawAllX           ->  cd(iFit);
        gStyle              ->  SetOptStat(0);
        gPad                ->  SetLogy();
        
        Int_t aIntBin = -1;
        if ( iFit == 1 )
        {
            aIntBin = 3;
        }
        if ( iFit == 10 )
        {
            aIntBin = 12;
        }
        
        cout << "XXX" << endl;
        hName = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit+1],fArrPT2D[iFit+2]);
        TH1D *  h1D_ResX    =   new TH1D    (*(aTarget->ProjectionX(hName,iFit+2,iFit+2)));
        h1D_ResX            ->  SetTitle(Form("Slice in #phi_{1} P_{T} from %.1f to %.1f GeV",fArrPT2D[iFit+1],fArrPT2D[iFit+2]));
        h1D_ResX            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
        h1D_ResX            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
        h1D_ResX            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        fLevyMm[iFit][0]    =   *fLevyFit1D;
        h1D_ResX            ->  Draw();
        fLevyMm[iFit][0]    .   Draw("same");
        
        cout << "XXX" << endl;
        fReturn[iFit][0]    =   ExtrapolateVl   (h1D_ResX,Form("X_%d",iFit),true,aIntBin);
        hSliceFY            ->  SetBinContent   (iFit+2,fReturn[iFit][0][0]);
        hSliceFY            ->  SetBinError     (iFit+2,fReturn[iFit][0][1]);
        hSlice2Y            ->  SetBinContent   (iFit+2,fReturn[iFit][0][0]);
        hSlice2Y            ->  SetBinError     (iFit+2,fReturn[iFit][0][2]);
        
        cout << "Y-" << iFit << endl;
        // Y-Projection
        fDrawAllY           ->  cd(iFit);
        gStyle              ->  SetOptStat(0);
        gPad                ->  SetLogy();
        hName = Form("XProjection_PT_%.1f_%.1f",fArrPT2D[iFit+1],fArrPT2D[iFit+2]);
        TH1D *  h1D_ResY    =   new TH1D    (*(aTarget->ProjectionY(hName,iFit+2,iFit+2)));
        h1D_ResY            ->  SetTitle(Form("Slice in #phi_{1} P_{T} from %.1f to %.1f GeV",fArrPT2D[iFit+1],fArrPT2D[iFit+2]));
        h1D_ResY            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
        h1D_ResY            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
        h1D_ResY            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
        fLevyMm[iFit][1]    =   *fLevyFit1D;
        h1D_ResY            ->  Draw();
        fLevyMm[iFit][1]    .   Draw("same");
        
        fReturn[iFit][1]    =   ExtrapolateVl   (h1D_ResY,Form("Y_%d",iFit),true,aIntBin);
        hSliceFX            ->  SetBinContent   (iFit+2,fReturn[iFit][1][0]);
        hSliceFX            ->  SetBinError     (iFit+2,fReturn[iFit][1][1]);
        hSlice2X            ->  SetBinContent   (iFit+2,fReturn[iFit][1][0]);
        hSlice2X            ->  SetBinError     (iFit+2,fReturn[iFit][1][2]);
    }
    
    fDrawFull           ->  cd(0);
    gStyle              ->  SetOptStat(0);
    gPad                ->  SetLogy();
    hSliceFY            ->  SetTitle("Slice in #phi_{1} P_{T} Y");
    hSliceFY            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
    hSliceFY            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
    hSliceFY            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
    fLevyMm[0][0]       =   *fLevyFit1D;
    hSliceFY            ->  Draw();
    fLevyMm[0][0]       .   Draw("same");
    
    fReturn[0][0]       =   ExtrapolateVl   (hSliceFY,"Y_Full",true,-1);
    
    fDrawFull           ->  cd(1);
    gStyle              ->  SetOptStat(0);
    gPad                ->  SetLogy();
    hSliceFX            ->  SetTitle("Slice in #phi_{1} P_{T} X");
    hSliceFX            ->  GetXaxis()  ->  SetTitle("#phi_{2} P_{T} GeV");
    hSliceFX            ->  GetYaxis()  ->  SetTitle("#frac{d^{3}N #phi_{1} }{dydp_{T}d#phi_{2}}(GeV/c)^{-1}");
    hSliceFX            ->  Fit(fLevyFit1D,"IMREQ0S","",0.4,10.);
    fLevyMm[0][1]       =   *fLevyFit1D;
    hSliceFX            ->  Draw();
    fLevyMm[0][1]       .   Draw("same");
    
    fReturn[0][1]       =   ExtrapolateVl   (hSliceFX,"X_Full",true,-1);
    
    fDrawFull->Write();
    fDrawAllX->Write();
    fDrawAllY->Write();
    delete fDrawAllX;
    delete fDrawAllY;
    
    gROOT->SetBatch(false);
    return fReturn;
}
 
TGraphErrors **     BuildTGraphEr   ( Double_t * aExtrapolateVl1D, Double_t *** aExtrapolateVl2D )
{
    TGraphErrors ** fResult =   new TGraphErrors* [16];
    for ( Int_t iPnt = 0; iPnt < 16; iPnt++ )
    {
        fResult[iPnt]   =   new TGraphErrors();
    }
    
    // Final Results yields 1D
    fResult[0]      ->  SetNameTitle    ("1D_dN_Stat","1D_dN_Stat");
    fResult[0]      ->  SetPoint        (0,1,aExtrapolateVl1D[0]);
    fResult[0]      ->  SetPointError   (0,0,aExtrapolateVl1D[1]);
    fResult[1]      ->  SetNameTitle    ("1D_dN_Syst","1D_dN_Syst");
    fResult[1]      ->  SetPoint        (0,1,aExtrapolateVl1D[0]);
    fResult[1]      ->  SetPointError   (0,0,aExtrapolateVl1D[2]);
    
    // Final Results yields 2D
    fResult[2]      ->  SetNameTitle    ("2D_dN_Stat","2D_dN_Stat");
    fResult[2]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[2]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][1]);
    fResult[3]      ->  SetNameTitle    ("2D_dN_Syst","2D_dN_Syst");
    fResult[3]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[3]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][2]);
    /*
    fResult[1]      ->  SetNameTitle    ("2D_dN_Stat","2D_dN_Stat");
    fResult[2]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[2]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][1][1]);
    fResult[1]      ->  SetNameTitle    ("2D_dN_Syst","2D_dN_Syst");
    fResult[3]      ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][0]);
    fResult[3]      ->  SetPointError   (0,0,aExtrapolateVl2D[0][1][2]);
    */
     
    // Final Results yields 2D - XY Proj
    for ( Int_t iPnt = 1; iPnt <= 10; iPnt++ )
    {
        fResult[4]   ->  SetNameTitle   ("2D_dN_Stat_X","2D_dN_Stat_X");
        fResult[4]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][0]);
        fResult[5]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][0]);
        fResult[5]   ->  SetNameTitle   ("2D_dN_Stat_Y","2D_dN_Stat_Y");
        fResult[4]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][1]);
        fResult[5]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][1]);
        fResult[6]   ->  SetNameTitle   ("2D_dN_Syst_X","2D_dN_Syst_X");
        fResult[6]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][0]);
        fResult[7]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][0]);
        fResult[7]   ->  SetNameTitle   ("2D_dN_Syst_Y","2D_dN_Syst_Y");
        fResult[6]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][2]);
        fResult[7]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][2]);
    }
    
    // Final Results PT 1D
    fResult[8]      ->  SetNameTitle    ("1D_PT_Stat","1D_PT_Stat");
    fResult[8]      ->  SetPoint        (0,1,aExtrapolateVl1D[3]);
    fResult[8]      ->  SetPointError   (0,0,aExtrapolateVl1D[4]);
    fResult[9]      ->  SetNameTitle    ("1D_PT_Syst","1D_PT_Syst");
    fResult[9]      ->  SetPoint        (0,1,aExtrapolateVl1D[3]);
    fResult[9]      ->  SetPointError   (0,0,aExtrapolateVl1D[5]);
    
    // Final Results PT 2D
    fResult[10]     ->  SetNameTitle    ("2D_PT_Stat","2D_PT_Stat");
    fResult[10]     ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][3]);
    fResult[10]     ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][4]);
    fResult[11]     ->  SetNameTitle    ("2D_PT_Syst","2D_PT_Syst");
    fResult[11]     ->  SetPoint        (0,2,aExtrapolateVl2D[0][0][3]);
    fResult[11]     ->  SetPointError   (0,0,aExtrapolateVl2D[0][0][5]);
    
    // Final Results PT 2D - XY Proj
    for ( Int_t iPnt = 1; iPnt <= 10; iPnt++ )
    {
        fResult[12]   ->  SetNameTitle   ("2D_PT_Stat_X","2D_PT_Stat_X");
        fResult[12]   ->  SetPoint       (iPnt-1,  .5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][3]);
        fResult[13]   ->  SetPoint       (iPnt-1,  .5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][3]);
        fResult[13]   ->  SetNameTitle   ("2D_PT_Stat_Y","2D_PT_Stat_Y");
        fResult[12]   ->  SetPointError  (iPnt-1,  .5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][4]);
        fResult[13]   ->  SetPointError  (iPnt-1,  .5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][4]);
        fResult[14]   ->  SetNameTitle   ("2D_PT_Syst_X","2D_PT_Syst_X");
        fResult[14]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][3]);
        fResult[15]   ->  SetPoint       (iPnt-1,.5*fArrPT2D[iPnt+1]+.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][3]);
        fResult[15]   ->  SetNameTitle   ("2D_PT_Syst_Y","2D_PT_Syst_Y");
        fResult[14]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][0][5]);
        fResult[15]   ->  SetPointError  (iPnt-1,.5*fArrPT2D[iPnt+1]-.5*fArrPT2D[iPnt+2],aExtrapolateVl2D[iPnt][1][5]);
    }
    
    return fResult;
}

#endif
