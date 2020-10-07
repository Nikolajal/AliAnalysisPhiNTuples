//
//  Util_Histograms.cpp
//  Utilities for Histograms and Graphs general managment in Root and Roofit
//  Nicola Rubini (R&D)
//
// TODO: Set options for systematic add and square
// TODO: Generic Ratio plot
// TODO: Clean-up the Mosaic, slice check

#include "SetValues.h"
#include "SetFunctions.h"

using namespace std;
using namespace RooFit;

// Error Management

template <class Tclass>                     // You should have an header file SetValues.h where the systematic uncertainties are defined.
Tclass*         AddSystematics              ( Tclass* aTarget, string aOption = "" )
{
    Tclass*     fReturn     =   new Tclass  (*aTarget);
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fReturn ->  SetBinContent           (iBin,   aTarget ->GetBinContent (iBin) );
        fReturn ->  SetBinError             (iBin,   SquareSystematics( aTarget ->GetBinError   (iBin) ) );
    }
    return fReturn;
}

template <class Tclass>    
Tclass*         AddSystematics              ( Tclass* aTarget, Tclass* aTarge2 )
{
    Tclass*     fReturn     =   new Tclass  (*aTarget);
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fReturn ->  SetBinContent           (iBin,   aTarget ->GetBinContent (iBin) );
        fReturn ->  SetBinError             (iBin,   sqrt( pow( aTarget ->GetBinError   (iBin), 2 ) + pow( aTarge2 ->GetBinError   (iBin), 2 ) )  );
    }
    return fReturn;
}
                            
Double_t        SquareSystematics           ( Double_t aTarget )
{
    return 1.;
}
                            
TH1F *          SetSystErrorsh              ( TH1F * aTarget ) // Togliere l'h
{
    TH1F *  fReturn =   new TH1F (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   sqrt(fError*fError+kSystematical1D_*fValue*kSystematical1D_*fValue));
    }
    return fReturn;
}
                            
TH1D *          SetSystErrorsh              ( TH1D * aTarget ) // Stat + systematics
{
    TH1D *  fReturn =   new TH1D (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   sqrt(fError*fError+kSystematical2D_*fValue*kSystematical2D_*fValue));
    }
    return fReturn;
}

TH1D *          SetSystErrorsh              ( TH1D * aTarget, TH1D * aTarge2 ) // Stat + systematics
{
    TH1D *  fReturn =   new TH1D (*aTarget);
    Double_t    fError,    fValue,  fErro2;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fErro2  =   aTarge2 ->GetBinError   (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   sqrt(fError*fError+fErro2*fErro2));
    }
    return fReturn;
}
                                    
TH1F *          SetSystErrors2              ( TH1F * aTarget ) // Solo systematics
{
    TH1F *  fReturn =   new TH1F (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   sqrt(kSystematical1D_*fValue*kSystematical1D_*fValue));
    }
    return fReturn;
}
                            
TH1D *          SetSystErrors2              ( TH1D * aTarget ) // Togliere l'h
{
    TH1D *  fReturn =   new TH1D (*aTarget);
    Double_t    fError,    fValue;
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fError  =   aTarget ->GetBinError   (iBin);
        fValue  =   aTarget ->GetBinContent (iBin);
        fReturn ->  SetBinContent           (iBin,   fValue );
        fReturn ->  SetBinError             (iBin,   sqrt(kSystematical2D_*fValue*kSystematical2D_*fValue));
    }
    return fReturn;
}

// General Histogram management
                            
Double_t        FuncIntegrals               ( TF1 * aFunction, Double_t aLowBound, Double_t aHigBound, string aOption = "", Double_t aEpsilon = 1.e-4 )
{
    Double_t    fResult = 0;
    Int_t       fiTer,  fnPowr;
    Bool_t      fWidth, fPower;
    
    if ( aOption.find("W") != -1 )  fWidth  =   true;
    else                                fWidth  =   false;
    
    if ( aOption.find("x^") != -1 )   { fPower  =   true;   fnPowr  =   aOption.at(aOption.find("x^")+2)-'0'; }
    else                                fPower  =   false;
    
    // Starting iterations at 0
    fiTer = 0;
    while ( true )
    {
        // Determining the pT point to calculate
        auto fPoint = aLowBound + aEpsilon*fiTer;
        
        // Exiting the loop if out of bound
        if ( fPoint >= aHigBound ) break;
        
        Double_t        fCycleAdd   =   ( aFunction->Eval(fPoint) );
        
        if ( fWidth )   fCycleAdd  *=   ( aEpsilon );
        if ( fPower )   fCycleAdd  *=   pow( ( fPoint ) , fnPowr );
        
        fResult += fCycleAdd;
        fiTer++;
    }
    
    return fResult;
}

template <class Tclass>
Double_t        HistIntegrals               ( Tclass* aTarget, string aOption = "", Double_t aLowBound = 0, Double_t aHigBound = 0 )
{
    
    Double_t    fResult = 0;
    Int_t       fnPowr;
    Bool_t      fWidth, fPower, fError;
    
    if ( aOption.find("W") != -1 )      fWidth  =   true;
    else                                fWidth  =   false;
    
    if ( aOption.find("x^") != -1 )   { fPower  =   true;   fnPowr  =   aOption.at(aOption.find("x^")+2)-'0'; }
    else                                fPower  =   false;
    
    if ( aOption.find("E") != -1 )      fError  =   true;
    else                                fError  =   false;
    
    // Starting iterations at 0
    for ( Int_t fiTer = 1; fiTer <= aTarget->GetNbinsX(); fiTer++ )
    {
        // Determining the pT point to calculate
        auto fPoint =   aTarget->GetBinCenter(fiTer);
        
        // Exiting the loop if out of bound
        if ( aTarget->GetBinLowEdge(fiTer)   <= aLowBound && aLowBound != aHigBound ) continue;
        if ( aTarget->GetBinLowEdge(fiTer+1) >= aHigBound && aLowBound != aHigBound ) break;
        
        Double_t        fCycleAdd   =   ( aTarget->GetBinContent(fiTer) );
        if ( fError )   fCycleAdd   =   ( aTarget->GetBinError(fiTer) );
        if ( fWidth )   fCycleAdd  *=   ( aTarget->GetBinWidth(fiTer) );
        if ( fPower )   fCycleAdd  *=   pow( ( aTarget->GetBinCenter(fiTer) ) , fnPowr );
        if ( fError )   fCycleAdd  *=   fCycleAdd;
        fResult += fCycleAdd;
    }
    
    if ( fError )   return sqrt(fResult);
    return fResult;
}

template <class Tclass>
Tclass*         SetRandPoints               ( Tclass* aTarget, string aOption = "" )
{
    Tclass*     fReturn     =   new Tclass  (*aTarget);
    for ( Int_t iBin = 1; iBin <= aTarget->GetNbinsX(); iBin++ )
    {
        fReturn ->  SetBinContent           (iBin,   aTarget ->GetBinContent (iBin) + fRandomGen->Gaus(0.,aTarget ->GetBinError   (iBin)) );
        fReturn ->  SetBinError             (iBin,   aTarget ->GetBinError   (iBin) );
    }
    return fReturn;
}

void ratioplot ( )  // To Be Done
{
    /*
    // Setting Systematics
    TGraphErrors *hDataSyst = new TGraphErrors (SetSystErrors2(hData));
    hDataSyst->RemovePoint(0);
    hDataSyst->RemovePoint(0);
    if ( bPythiaTest ) aName += "P8";
    
    // Setting Aspect
    setLine                 (hPythia6,1);
    setLine                 (hPythia8,2);
    
    hData       ->SetMarkerStyle(22);
    hData       ->SetMarkerColor(38);
    hData       ->SetFillColorAlpha(33,0.66);
    hData       ->SetLineColorAlpha(kRed,1.0);
    hData       ->SetTitle("");
    hData       ->GetYaxis()->SetTitleSize(20);
    hData       ->GetYaxis()->SetTitleFont(43);
    hData       ->GetYaxis()->SetTitle("#frac{d^{3}N_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hData       ->GetYaxis()->SetTitle("#frac{d^{3}N^{RAW}_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
    hData       ->GetYaxis()->SetTitleOffset(1.5);
    
    hDataSyst       ->SetMarkerStyle(22);
    hDataSyst       ->SetMarkerColor(38);
    hDataSyst       ->SetLineColor(33);

    
    // Define the Canvas
    TCanvas *cRatioPlot = new TCanvas( aName.c_str(), aName.c_str(), 800, 800);
    
    TLegend *Legend     = new TLegend(0.7,0.75,0.9,0.9);
    Legend      ->AddEntry(hData, "Data", "P");
    Legend      ->AddEntry(hData, "Stat. Err.", "E");
    if ( !bPythiaTest ) Legend      ->AddEntry(hDataSyst, "Syst. Err.", "F");
    if ( !bPythiaTest ) Legend      ->AddEntry(hPythia6, "Pythia6", "L");
    Legend      ->AddEntry(hPythia8, "Pythia8", "L");
    
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1, 1.0);
    pad1    ->SetLogy();
    pad1    ->SetBottomMargin(0); // Upper and lower plot are joined
    pad1    ->Draw();             // Draw the upper pad: pad1
    pad1    ->cd();               // pad1 becomes the current pad
    
    hData   ->SetStats(0);        // No statistics on upper plot
    hData   ->Draw("EP");
    if ( !bPythiaTest ) hDataSyst->Draw("SAME P 5");// Draw h1
    hData   ->Draw("SAME EP");
    if ( !bPythiaTest ) hPythia6->Draw("SAME HIST L");
    hPythia8->Draw("SAME HIST L");       // Draw h2 on top of h1
    Legend  ->Draw("same");
    
    //hData->GetYaxis()->SetLabelSize(0.);
    
    TGaxis *axis = new TGaxis( 0., 20, 10., 220, 20,220,510,"");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();

    // lower plot will be in pad
    cRatioPlot->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0., 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    TH1F *hPythia6Ratio = (TH1F*)hPythia6->Clone("hPythia6Ratio");
    hPythia6Ratio->Sumw2();     // No statistics on lower plot
    hPythia6Ratio->Divide(hData);
    hPythia6Ratio->SetLineStyle(1);
    TH1F *hPythia8Ratio = (TH1F*)hPythia8->Clone("hPythia8Ratio");
    hPythia8Ratio->Sumw2();
    hPythia8Ratio->Divide(hData);
    hPythia8Ratio->SetLineStyle(1);
    TH1F *hDataRatio = new TH1F(Form("hDataRatio_%s",aName.c_str()),"hDataRatio",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValue =   hData  ->GetBinContent(iTer);
        auto fError =   hData  ->GetBinError(iTer);
        hDataRatio->SetBinContent   (iTer,1);
        hDataRatio->SetBinError     (iTer,0.062);
        if ( fValue == 0 ) hDataRatio->SetBinError     (iTer,0);
    }
    TH1F * hDatS = SetSystErrorsh(hDataRatio);
    if ( bPythiaTest ) hDatS = hDataRatio;
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValu6 =   hPythia6Ratio  ->GetBinContent(iTer);
        auto fValu8 =   hPythia8Ratio  ->GetBinContent(iTer);
        auto fError =   hDatS  ->GetBinError(iTer);
        hPythia6Ratio->SetBinError     (iTer,fError*fValu6);
        hPythia8Ratio->SetBinError     (iTer,fError*fValu8);
    }
    
    hDataRatio->SetStats(0);  // Define Y ..
    hDataRatio->SetTitle("");
    hDataRatio->SetLineColorAlpha(kBlack,1.00);
    hDataRatio->SetFillColorAlpha(33,0.66);
    hDataRatio->Draw("HIST L E0 SAME");
    
    // Setting Axes
    // - // XAxis
    hDataRatio->GetXaxis()->SetTitle("P_{T}#phi_{1} (GeV/c)");
    hDataRatio->GetXaxis()->SetTitleSize(20);
    hDataRatio->GetXaxis()->SetTitleFont(43);
    hDataRatio->GetXaxis()->SetTitleOffset(3.);
    hDataRatio->GetXaxis()->SetLabelFont(43);
    hDataRatio->GetXaxis()->SetLabelSize(15);
    hDataRatio->SetMinimum(0.);                     // Range set
    hDataRatio->SetMaximum(2.9);
    // - // YAxis
    hDataRatio->GetYaxis()->SetTitle("MODEL / DATA");
    if ( bPythiaTest ) hDataRatio->GetYaxis()->SetTitle("Tru / Res");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hDataRatio->GetYaxis()->SetTitle("Rec / Raw");
    hDataRatio->GetYaxis()->SetNdivisions(505);
    hDataRatio->GetYaxis()->SetTitleSize(20);
    hDataRatio->GetYaxis()->SetTitleFont(43);
    hDataRatio->GetYaxis()->SetTitleOffset(2.);
    hDataRatio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hDataRatio->GetYaxis()->SetLabelSize(15);
    
    // Define the ratio plot
    hPythia6Ratio->SetMarkerStyle(0);
    hPythia6Ratio->SetStats(0);// No statistics on lower plot
    if ( !bPythiaTest ) hPythia6Ratio->Draw("SAME EP");       // Draw the ratio plot
    
    hPythia8Ratio->SetMarkerStyle(0);
    hPythia8Ratio->SetStats(0);
    hPythia8Ratio->Draw("SAME EP");       // Draw the ratio plot

    // Ratio plot (h3) settings
    hPythia6Ratio->SetTitle(""); // Remove the first ratio title
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    hPythia8Ratio->SetTitle(""); // Remove the second ratio presence
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    cRatioPlot->Write();
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.pdf", aName.c_str()));
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.png", aName.c_str()));
    delete cRatioPlot;
    */
}

void fSliceCheck ( TH1F * hCheck1, TH1F * hCheck2, const char* hName1 = "1", const char* hName2 = "2", const char* fName = "-1", const char* fTitle = "" )
{
    TLatex * latext             =   new TLatex();
    
    TH1F hUtility1 = *hCheck1;
    TH1F hUtility2 = *hCheck2;
    hUtility1.SetMarkerStyle(29);
    hUtility1.SetMarkerColor(2);
    hUtility1.SetName("1");
    hUtility2.SetMarkerStyle(33);
    hUtility2.SetMarkerColor(4);
    hUtility2.SetName("2");
    
    TCanvas * fSaveToCanvas;
    fSaveToCanvas               =   new TCanvas(fName,fTitle);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hUtility1                   .Draw("same");
    hUtility2                   .Draw("same");
    TLegend * fLegend           = new TLegend   ();
    fLegend                     ->SetFillColor(kWhite);
    fLegend                     ->SetLineColor(kWhite);
    fLegend                     ->AddEntry("1",    hName1,      "EP");
    fLegend                     ->AddEntry("2",    hName2,      "EP");
    fLegend                     ->Draw("same");
    fSaveToCanvas->Write();
    
    delete fSaveToCanvas;
}

void fSliceCheck ( TH1D * hCheck1, TH1D * hCheck2, const char* hName1 = "1", const char* hName2 = "2", const char* fName = "-1", const char* fTitle = "")
{
    TLatex * latext             =   new TLatex();
    
    TH1D hUtility1 = *hCheck1;
    TH1D hUtility2 = *hCheck2;
    TH1D hUtility3 = *hCheck1;
    hUtility1.SetName("1");
    hUtility2.SetName("2");
    hUtility1.SetMarkerStyle(29);
    hUtility1.SetMarkerColor(2);
    hUtility3.SetMarkerStyle(29);
    hUtility3.SetMarkerColor(2);
    hUtility2.SetMarkerStyle(33);
    hUtility2.SetMarkerColor(4);
    
    TCanvas * fSaveToCanvas;
    fSaveToCanvas               =   new TCanvas(fName,fTitle);
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hUtility1                   .Draw("same");
    hUtility2                   .Draw("same");
    TLegend * fLegend           = new TLegend   ();
    fLegend                     ->SetFillColor(kWhite);
    fLegend                     ->SetLineColor(kWhite);
    fLegend                     ->AddEntry("1",    hName1,      "EP");
    fLegend                     ->AddEntry("2",    hName2,      "EP");
    fLegend                     ->Draw("same");
    fSaveToCanvas->Write();
    
    delete fSaveToCanvas;
}

void fSliceCheckPT2D ( TH2F * hCheck1, TH2F * hCheck2, const char* hName1 = "1", const char* hName2 = "2" )
{
    TH1D * hUtility1X;
    TH1D * hUtility2X;
    TH1D * hUtility1Y;
    TH1D * hUtility2Y;
    for ( int iChk = 0; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        hUtility1X  = (hCheck1->ProjectionX(Form("%s_X_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility2X  = (hCheck2->ProjectionX(Form("%s_X_Profile2_PT_%.2f_%.2f",hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1X->SetTitle(Form("X_%s",hName1));
        hUtility2X->SetTitle(Form("X_%s",hName2));
        fSliceCheck(hUtility1X,hUtility2X,hName1,hName2,Form("%s_%s_X_Profile_PT_%.2f_%.2f",hName1,hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),Form("X Profile in pT %.2f to %.2f GeV",fArrPT2D[iChk],fArrPT2D[iChk+1]));
                                         
        // Y-Projection
        hUtility1Y  = (hCheck1->ProjectionY(Form("%s_Y_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility2Y  = (hCheck2->ProjectionY(Form("%s_Y_Profile2_PT_%.2f_%.2f",hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1Y->SetTitle(Form("Y_%s",hName1));
        hUtility2Y->SetTitle(Form("Y_%s",hName2));
        fSliceCheck(hUtility1Y,hUtility2Y,hName1,hName2,Form("%s_%s_Y_Profile_PT_%.2f_%.2f",hName1,hName2,fArrPT2D[iChk],fArrPT2D[iChk+1]),Form("Y Profile in pT %.2f to %.2f GeV",fArrPT2D[iChk],fArrPT2D[iChk+1]));
    }
}

void fSlicePT2D ( TH2F * hCheck1, const char* hName1 = "1")
{
    TH1D * hUtility1X;
    TH1D * hUtility1Y;
    for ( int iChk = 0; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        hUtility1X  = (hCheck1->ProjectionX(Form("%s_X_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1X->SetTitle(Form("X_%s",hName1));
        hUtility1X->SetMaximum(1.5);
        hUtility1X->SetMinimum(0.5);
        hUtility1X->Write();
                                         
        // Y-Projection
        hUtility1Y  = (hCheck1->ProjectionY(Form("%s_Y_Profile1_PT_%.2f_%.2f",hName1,fArrPT2D[iChk],fArrPT2D[iChk+1]),iChk+1,iChk+1));
        hUtility1Y->SetTitle(Form("Y_%s",hName1));
        hUtility1Y->SetMaximum(1.5);
        hUtility1Y->SetMinimum(0.5);
        hUtility1Y->Write();
    }
}

void fMosaicCanvas ( TH1F ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TF1 ** hFunc = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
        
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
        
        // Function Draw
        if ( !hFunc ) continue;
        if ( !hFunc[iHisto-1] ) continue;
        hFunc[iHisto-1]->SetLineColor(2);
        hFunc[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH1D ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TF1 ** hFunc = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
        
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
        
        // Function Draw
        if ( !hFunc ) continue;
        if ( !hFunc[iHisto-1] ) continue;
        hFunc[iHisto-1]->SetLineColor(2);
        hFunc[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH1F ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TH1F ** hDat2 = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
        
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
        
        // Data 2 Draw
        if ( !hDat2 ) continue;
        if ( !hDat2[iHisto-1] ) continue;
        hDat2[iHisto-1]->SetMarkerColor(2);
        hDat2[iHisto-1]->SetMarkerStyle(33);
        hDat2[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH1D ** hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, TH1D ** hDat2 = nullptr, const char*  fFTOpt = "" )
{
    TCanvas * cResult = new TCanvas (fName,fName);
    cResult->Divide(xDiv,yDiv);
    for ( Int_t iHisto = 1; iHisto <= xDiv*yDiv; iHisto++ )
    {
        cResult->cd(iHisto);
        if ( fLogy ) gPad->SetLogy();
               
        // Data Draw
        if ( !hData ) continue;
        if ( !hData[iHisto-1] ) continue;
        hData[iHisto-1]->SetMarkerColor(1);
        hData[iHisto-1]->SetMarkerStyle(33);
        hData[iHisto-1]->Draw(fDTOpt);
               
        // Data 2 Draw
        if ( !hDat2 ) continue;
        if ( !hDat2[iHisto-1] ) continue;
        hDat2[iHisto-1]->SetMarkerColor(2);
        hDat2[iHisto-1]->SetMarkerStyle(33);
        hDat2[iHisto-1]->Draw(Form("same %s",fFTOpt));
    }
    if ( fSaveToFile )
    {
        cResult->Write();
        cResult->SaveAs(Form("%s.pdf",fName));
        cResult->SaveAs(Form("%s.png",fName));
        delete cResult;
    }
}

void fMosaicCanvas ( TH2F * hData, const char*  fDTOpt = "", const char*  fName = "", Int_t xDiv = 1, Int_t yDiv = 1, Bool_t fSaveToFile = false, Bool_t fLogy = false, Bool_t fPT2D = false, Double_t fMax = -1., Double_t fMin = -1. )
{
    // Building utility histograms
    TH1D ** hSlcX;
    if ( fPT2D )    hSlcX = new TH1D * [nBinPT2D];
    else            hSlcX = new TH1D * [xDiv*yDiv];
    TH1D ** hSlcY;
    if ( fPT2D )    hSlcY = new TH1D * [nBinPT2D];
    else            hSlcY = new TH1D * [xDiv*yDiv];
    TH1D ** fAmbg   = nullptr;
    
    auto iHist2 = 0;
    auto nTotal = 0;
    if ( fPT2D )    nTotal = nBinPT2D;
    else            nTotal = xDiv*yDiv;
    for ( Int_t iHisto = 0; iHisto < nTotal; iHisto++ )
    {
        if ( fPT2D && fArrPT2D[iHisto+1] <= 0.41 ) continue;
        
        // X-Projection
        if ( fPT2D )    hName = Form("XProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        else            hName = Form("XProjection_Bin_%i_%i_%s",iHisto,iHisto+1,fName);
        if ( fPT2D )    hSlcX[iHist2] = static_cast<TH1D*>(hData->ProjectionX("",iHisto+1,iHisto+1));
        else            hSlcX[iHisto] = static_cast<TH1D*>(hData->ProjectionX("",iHisto+1,iHisto+1));
        if ( fPT2D )    hSlcX[iHist2]->SetNameTitle(hName,hName);
        else            hSlcX[iHisto]->SetNameTitle(hName,hName);
        if ( fPT2D && fMax != -1. )  hSlcX[iHist2]->SetMaximum(fMax);
        if ( fPT2D && fMin != -1. )  hSlcX[iHist2]->SetMinimum(fMin);
        if ( !fPT2D&& fMax != -1. )  hSlcX[iHisto]->SetMaximum(fMax);
        if ( !fPT2D&& fMin != -1. )  hSlcX[iHisto]->SetMinimum(fMin);
        
        // Y Projection
        if ( fPT2D )    hName = Form("YProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        else            hName = Form("YProjection_Bin_%i_%i_%s",iHisto,iHisto+1,fName);
        if ( fPT2D )    hSlcY[iHist2] = static_cast<TH1D*>(hData->ProjectionY("",iHisto+1,iHisto+1));
        else            hSlcY[iHisto] = static_cast<TH1D*>(hData->ProjectionY("",iHisto+1,iHisto+1));
        if ( fPT2D )    hSlcY[iHist2]->SetNameTitle(hName,hName);
        else            hSlcY[iHisto]->SetNameTitle(hName,hName);
        if ( fPT2D && fMax != -1. )  hSlcY[iHist2]->SetMaximum(fMax);
        if ( fPT2D && fMin != -1. )  hSlcY[iHist2]->SetMinimum(fMin);
        if ( !fPT2D&& fMax != -1. )  hSlcY[iHisto]->SetMaximum(fMax);
        if ( !fPT2D&& fMin != -1. )  hSlcY[iHisto]->SetMinimum(fMin);
        
        iHist2++;
    }
    
    fMosaicCanvas(hSlcX,"",Form("XProjection_%s",fName),xDiv,yDiv,fSaveToFile,fLogy,fAmbg);
    fMosaicCanvas(hSlcY,"",Form("YProjection_%s",fName),xDiv,yDiv,fSaveToFile,fLogy,fAmbg);
    return;
}

// Analysis-Targeted functions

// --- // Setters for graphical styles

template < class Tclass >
void            SetGraphicStyle             ( Tclass * aTarget, string aOption = "" )
{
    if ( aOption.find("DT") != -1 )
    {
        aTarget->SetMarkerStyle(22);
        aTarget->SetMarkerColor(38);
        aTarget->SetLineColor(kRed);
    }
    if ( aOption.find("P8") != -1 )
    {
        aTarget->SetLineStyle(9);
        aTarget->SetLineColor(kBlue+3);
        aTarget->SetLineWidth(3);
    }
    if ( aOption.find("P6") != -1 )
    {
        aTarget->SetLineStyle(10);
        aTarget->SetLineColor(kBlue);
        aTarget->SetLineWidth(3);
    }
}

template < class Tclass >
void            SetDescription              ( Tclass * aTarget, string aOption = "" )
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


// --- // Generate histogram for results
                            
TH1D*           GenResultHistogram          ( Double_t***   aInput, string aOption = "" )
{
    return nullptr;
}

// --------- //


// --------- //

// Color Style and Width, Fill Marker Line
Int_t kColor[] = {38,kBlue,kBlue+3,46,38};
Int_t kStyle[] = {26,9,10,25,22};
Int_t kWidth[] = {1,3,3,1,1};

int             fLegendSelect       ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 2;
    else return -1;
}

void            fLegendMaker        ( RooPlot * fRooPlot, const char * fSelect, TLegend * fLegend )
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

int             fAxisSelect         ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 3;
    else return -1;
}
 
void            fAxisMaker          ( RooPlot * fRooPlot, const char * fSelect )
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

int             fPlotterSelect      ( string fOption )
{
    if ( !fOption.compare("InvMass1D") )   return 1;
    if ( !fOption.compare("xInvMass2D") )  return 2;
    if ( !fOption.compare("yInvMass2D") )  return 2;
    else return -1;
}

void            fRooPlotPlotter     ( RooPlot * fRooPlot, const char * fSelect, RooAddPdf fModel , RooDataHist * fData )
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

void            fRooPlotMaker       ( RooPlot * fRooPlot, TLegend * fLegend, RooAddPdf fModel , RooDataHist * fData, const char * fSelect )
{
    fRooPlotPlotter(fRooPlot,fSelect,fModel,fData);
    fLegendMaker(fRooPlot,fSelect,fLegend);
    fAxisMaker(fRooPlot,fSelect);
}

void            setMarker           ( TH1F * hTarget, Int_t iType = 0 )
{
    hTarget->SetMarkerColor(kColor[iType]);
    hTarget->SetMarkerStyle(kStyle[iType]);
    hTarget->SetMarkerSize(kWidth[iType]);
}

void            setLine             ( TH1F * hTarget, Int_t iType = 0 )
{
    hTarget->SetLineColor(kColor[iType]);
    hTarget->SetLineStyle(kStyle[iType]);
    hTarget->SetLineWidth(kWidth[iType]);
}

void            setMarker           ( TH1D * hTarget, Int_t iType = 0 )
{
    hTarget->SetMarkerColor(kColor[iType]);
    hTarget->SetMarkerStyle(kStyle[iType]);
    hTarget->SetMarkerSize(kWidth[iType]);
}

void            setMarker           ( TGraphAsymmErrors * hTarget, Int_t iType = 0 )
{
    hTarget->SetMarkerColor(kColor[iType]);
    hTarget->SetMarkerStyle(kStyle[iType]);
    hTarget->SetMarkerSize(kWidth[iType]);
}

void            setLine             ( TH1D * hTarget, Int_t iType = 0 )
{
    hTarget->SetLineColor(kColor[iType]);
    hTarget->SetLineStyle(kStyle[iType]);
    hTarget->SetLineWidth(kWidth[iType]);
}

void            TH1_SetDescription ( TH1F * hTarget )
{
    hTarget->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
    hTarget->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
    hTarget->GetXaxis()->SetTitleOffset(1.15);
    hTarget->GetYaxis()->SetTitleOffset(1.15);
}

void            TH2_SetDescription ( TH2F * hTarget )
{
    hTarget->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    hTarget->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
    hTarget->GetXaxis()->SetTitleOffset(1.5);
    hTarget->GetYaxis()->SetTitleOffset(1.5);
}

void                TGrCompare1D    ( TGraphAsymmErrors * gDataStat, TGraphAsymmErrors * gDataSyst, TH1F * h1D_Tru_P6, TH1F * h1D_Tru_P8, TGraphAsymmErrors * gCompare = nullptr )
{
    // Setting Canvas
    TCanvas *       cCompare1D  =   new TCanvas ("cCompare1D","cCompare1D");
    gPad->SetLogy();
    
    // Setting Standard Colors
    if ( !!gCompare ) setMarker               (gCompare,3);
    setLine                 (h1D_Tru_P6,1);
    setLine                 (h1D_Tru_P8,2);
    
    // Setting Statistical error
    gDataStat->SetMarkerStyle(22);
    gDataStat->SetMarkerColor(38);
    gDataStat->SetFillColorAlpha(33,0.33);
    gDataStat->SetLineColorAlpha(33,0.66);
    gDataSyst->SetMarkerStyle(22);
    gDataSyst->SetMarkerColor(38);
    gDataSyst->SetLineColor(kBlack);
    
    // Setting Legend
    TLegend * lLegend1          =   new TLegend(0.65,0.65,0.85,0.85);
    lLegend1                    ->SetFillColor(kWhite);
    lLegend1                    ->SetLineColor(kWhite);
    lLegend1                    ->AddEntry(gDataStat,   "Data",             "P");
    lLegend1                    ->AddEntry(gDataStat,   "Stat. Err.",       "EF");
    lLegend1                    ->AddEntry(gDataSyst,   "Syst. Err.",       "EP[]");
    if ( !!gCompare ) lLegend1                    ->AddEntry(gCompare,    "Previous Paper",   "EP");
    lLegend1                    ->AddEntry(h1D_Tru_P6,  "Pythia 6",         "EP");
    lLegend1                    ->AddEntry(h1D_Tru_P8,  "Pythia 8",         "EP");

    // Setting Multigraph
    TMultiGraph *   mCompare1D  =   new TMultiGraph();
    mCompare1D      ->Add   (gDataStat, "AP35");
    mCompare1D      ->Add   (gDataSyst, "EP[]");
    if ( !!gCompare )    mCompare1D      ->Add   (gCompare,  "P");
    
    mCompare1D      ->GetXaxis()->SetLimits(0.,10.);
    //mCompare1D      ->GetYaxis()->SetLimits(3.e-6,1.e-1);
    mCompare1D      ->SetMaximum(5.e-2);
    mCompare1D      ->SetMinimum(3.e-6);
    mCompare1D      ->GetXaxis()->SetTitle("P_{T}#phi (GeV/c)");
    mCompare1D      ->GetYaxis()->SetTitle("#frac{d^{2}N#phi}{dydp_{T}#phi}(GeV/c)^{-1}");
    
    // Drawing
    mCompare1D      ->Draw  ("AP");
    h1D_Tru_P6      ->Draw  ("SAME HIST L");
    h1D_Tru_P8      ->Draw  ("SAME HIST L");
    lLegend1        ->Draw  ("SAME");
    
    cCompare1D->Write();
    cCompare1D->SaveAs("./graphs/cCompare1D.png");
    cCompare1D->SaveAs("./graphs/cCompare1D.pdf");
    gPad->SetLogx();
    cCompare1D->Write();
    cCompare1D->SaveAs("./graphs/cCompare1D_logx.png");
    cCompare1D->SaveAs("./graphs/cCompare1D_logx.pdf");
    gPad->SetLogx(false);
    delete cCompare1D;
}

void                TGrCompare2D    ( TH2F* h2D_Res, TH2F * h2D_Tru_P6, TH2F * h2D_Tru_P8)
{
    gROOT->SetBatch(true);
    TCanvas*cCompare2D  =   new TCanvas("cCompare2D","cCompare2D");
    cCompare2D->Divide(3,3);
    TH1D  **hProfileP8  =   new TH1D*   [nBinPT2D];
    TH1D  **hProfileP6  =   new TH1D*   [nBinPT2D];
    TH1D  **hProfileDT  =   new TH1D*   [nBinPT2D];
    TGraphAsymmErrors  **gProfileDT_Stat  =   new TGraphAsymmErrors*   [nBinPT2D];
    TGraphAsymmErrors  **gProfileDT_Syst  =   new TGraphAsymmErrors*   [nBinPT2D];
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        cCompare2D->cd(iTer-1);
        hProfileP8[iTer]        =   new TH1D(*(h2D_Tru_P8  ->ProjectionX(Form("P8_%i",iTer+1),iTer+1,iTer+1)));
        hProfileP6[iTer]        =   new TH1D(*(h2D_Tru_P6  ->ProjectionX(Form("P6_%i",iTer+1),iTer+1,iTer+1)));
        hProfileDT[iTer]        =   new TH1D(*(h2D_Res     ->ProjectionX(Form("DT_%i",iTer+1),iTer+1,iTer+1)));
        gProfileDT_Stat[iTer]   =   new TGraphAsymmErrors(hProfileDT[iTer]);
        gProfileDT_Syst[iTer]   =   new TGraphAsymmErrors(SetSystErrorsh(hProfileDT[iTer]));
        
        for ( Int_t jTer = 0; jTer < 2; jTer++ )
        {
            gProfileDT_Stat[iTer]->RemovePoint(0);
            gProfileDT_Syst[iTer]->RemovePoint(0);
        }
        
        // Setting Canvas
        TCanvas *       cCompare2D_bin  =   new TCanvas (Form("cCompare2D_%i",iTer),Form("cCompare2D_%i",iTer));
        gPad->SetLogy();
        
        // Setting Standard Colors
        setLine                 (hProfileP6[iTer],1);
        setLine                 (hProfileP8[iTer],2);
        
        // Setting Statistical error
        gProfileDT_Stat[iTer]->SetMarkerStyle(22);
        gProfileDT_Stat[iTer]->SetMarkerColor(38);
        gProfileDT_Stat[iTer]->SetFillColorAlpha(33,0.33);
        gProfileDT_Stat[iTer]->SetLineColorAlpha(33,0.66);
        gProfileDT_Syst[iTer]->SetMarkerStyle(22);
        gProfileDT_Syst[iTer]->SetMarkerColor(38);
        gProfileDT_Syst[iTer]->SetLineColor(kBlack);
        
        // Setting Legend
        TLegend * lLegend1          =   new TLegend(0.65,0.65,0.85,0.85);
        lLegend1                    ->SetFillColor(kWhite);
        lLegend1                    ->SetLineColor(kWhite);
        lLegend1                    ->AddEntry(gProfileDT_Stat[iTer],   "Data",             "P");
        lLegend1                    ->AddEntry(gProfileDT_Stat[iTer],   "Stat. Err.",       "EF");
        lLegend1                    ->AddEntry(gProfileDT_Syst[iTer],   "Syst. Err.",       "EP[]");
        lLegend1                    ->AddEntry(hProfileP6[iTer],  "Pythia 6",         "EP");
        lLegend1                    ->AddEntry(hProfileP8[iTer],  "Pythia 8",         "EP");

        // Setting Multigraph
        TMultiGraph *   mCompare1D  =   new TMultiGraph();
        mCompare1D      ->Add   (gProfileDT_Stat[iTer], "AP35");
        mCompare1D      ->Add   (gProfileDT_Syst[iTer], "EP[]");
        
        // Drawing
        mCompare1D      ->Draw  ("AP");
        hProfileP6[iTer]      ->Draw  ("same HIST L");
        hProfileP8[iTer]      ->Draw  ("same HIST L");
        lLegend1        ->Draw  ("same");
        
        cCompare2D_bin->Write();
        cCompare2D_bin->SaveAs(Form("./graphs/cCompare2D_%i.png",iTer));
        cCompare2D_bin->SaveAs(Form("./graphs/cCompare2D_%i.pdf",iTer));
        gPad->SetLogx();
        cCompare2D_bin->Write();
        cCompare2D_bin->SaveAs(Form("./graphs/cCompare2D_%i_logx.png",iTer));
        cCompare2D_bin->SaveAs(Form("./graphs/cCompare2D_%i_logx.pdf",iTer));
        gPad->SetLogx(false);
        
        cCompare2D->cd(iTer-2);
        mCompare1D      ->Draw  ("AP");
        hProfileP6[iTer]      ->Draw  ("same HIST L");
        hProfileP8[iTer]      ->Draw  ("same HIST L");
        lLegend1        ->Draw  ("same");
    }
    cCompare2D->Write();
    cCompare2D->SaveAs(Form("./graphs/cCompare2D.png"));
    cCompare2D->SaveAs(Form("./graphs/cCompare2D.pdf"));
    gROOT->SetBatch(false);
}

void  TGraphAEGeneratorPT   ( Double_t***   aInput, TH2F * h2D_Tru_P6, TH2F * h2D_Tru_P8)
{
    TCanvas           * c1      =   new TCanvas             ("c1_e","c1_e");
    TLegend           * Legend  =   new TLegend             (0.6,0.15,0.9,0.3);
    TMultiGraph       * fPTTotal=   new TMultiGraph         ("fPTTotal","fPTTotal");
    TGraphAsymmErrors * fPTX    =   new TGraphAsymmErrors   ();
    TGraphAsymmErrors * fPTY    =   new TGraphAsymmErrors   ();
    TGraphAsymmErrors * fPTXSyst=   new TGraphAsymmErrors   ();
    TGraphAsymmErrors * fPTYSyst=   new TGraphAsymmErrors   ();
    TGraphAsymmErrors * fPTXP6  =   new TGraphAsymmErrors   ();
    TGraphAsymmErrors * fPTYP8  =   new TGraphAsymmErrors   ();
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        fPTX->SetPoint       (iTer-2,  (fArrPT2D[iTer]+fArrPT2D[iTer+1])*0.5,aInput[iTer-1][0][3]);
        fPTX->SetPointError  (iTer-2, -(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,-(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,aInput[iTer-1][0][4],aInput[iTer-1][0][4]);
        fPTY->SetPoint       (iTer-2,  (fArrPT2D[iTer]+fArrPT2D[iTer+1])*0.5,aInput[iTer-1][1][3]);
        fPTY->SetPointError  (iTer-2, -(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,-(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,aInput[iTer-1][1][4],aInput[iTer-1][1][4]);
        fPTXSyst->SetPoint       (iTer-2,  (fArrPT2D[iTer]+fArrPT2D[iTer+1])*0.5,aInput[iTer-1][0][3]);
        fPTXSyst->SetPointError  (iTer-2, -(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,-(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,aInput[iTer-1][0][4],aInput[iTer-1][0][4]*1.2);
        fPTYSyst->SetPoint       (iTer-2,  (fArrPT2D[iTer]+fArrPT2D[iTer+1])*0.5,aInput[iTer-1][1][3]);
        fPTYSyst->SetPointError  (iTer-2, -(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,-(fArrPT2D[iTer]-fArrPT2D[iTer+1])*0.5,aInput[iTer-1][1][4],aInput[iTer-1][1][4]*1.2);
        auto xP6 = h2D_Tru_P6->ProjectionX("dd",iTer+1,iTer+1);
        auto xP8 = h2D_Tru_P8->ProjectionX("xx",iTer+1,iTer+1);
        fPTXP6->SetPoint(iTer-2,(fArrPT2D[iTer]+fArrPT2D[iTer+1])*0.5,(HistIntegrals( xP6,"W x^1"))/(xP6->Integral("width")));
        fPTYP8->SetPoint(iTer-2,(fArrPT2D[iTer]+fArrPT2D[iTer+1])*0.5,(HistIntegrals( xP8,"W x^1"))/(xP8->Integral("width")));
    }
    fPTX->SetMarkerStyle(22);
    fPTX->SetMarkerColor(38);
    fPTX->SetFillColorAlpha(33,0.33);
    fPTX->SetLineColorAlpha(kRed,1.);
    fPTXP6->SetLineColor(kBlue);
    fPTXP6->SetLineWidth(3);
    fPTXP6->SetLineStyle(10);
    fPTYP8->SetLineColor(kBlue+3);
    fPTYP8->SetLineWidth(3);
    fPTYP8->SetLineStyle(9);
    TLegend * PTLegend = new TLegend (0.7,0.7,0.9,0.9);
    fPTTotal->SetTitle("#LT p_{T}#phi_{2} #GT as a function of p_{T} #phi_{1}");
    fPTTotal->SetTitle("");
    fPTTotal->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    fPTTotal->GetYaxis()->SetTitle("#LT p_{T}#phi_{2} #GT (GeV/c)");
    //Legend      ->AddEntry(fPTX,"Stat. Err.","EP");
    fPTTotal    ->Add(fPTX,"EP");
    fPTTotal    ->Add(fPTXP6,"L");
    fPTTotal    ->Add(fPTYP8,"L");
    fPTTotal    ->Draw("AP");
    PTLegend->AddEntry(fPTX,"Data","EP");
    PTLegend->AddEntry(fPTXP6,"Pythia 6","L");
    PTLegend->AddEntry(fPTYP8,"Pythia 8","L");
    PTLegend->Draw("same");
    //Legend      ->Draw("same");
    c1->Write();
    c1->SaveAs("./graphs/PTTOTAL.pdf");
    c1->SaveAs("./graphs/PTTOTAL.png");
}

TH1D * TH1DGeneratorYield ( Double_t***   aInput, Bool_t aSystBool )
{
    TH1D * fResult = new TH1D ("","",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        fResult->SetBinContent  (iTer+1, aInput[iTer-1][0][0] );
        fResult->SetBinError    (iTer+1, aInput[iTer-1][0][1] );
        if ( aSystBool ) fResult->SetBinError    (iTer+1, aInput[iTer-1][0][2] );
    }
    return fResult;
}

void ratioplot ( TH1D * hData, TH1D * hPythia6, TH1D * hPythia8, string aName = "" )
{
    // Setting Systematics
    TGraphErrors *hDataSyst = new TGraphErrors (SetSystErrors2(hData));
    hDataSyst->RemovePoint(0);
    hDataSyst->RemovePoint(0);
    if ( bPythiaTest ) aName += "P8";
    
    // Setting Aspect
    setLine                 (hPythia6,1);
    setLine                 (hPythia8,2);
    
    hData       ->SetMarkerStyle(22);
    hData       ->SetMarkerColor(38);
    hData       ->SetFillColorAlpha(33,0.66);
    hData       ->SetLineColorAlpha(kRed,1.0);
    hData       ->SetTitle("");
    hData       ->GetYaxis()->SetTitleSize(20);
    hData       ->GetYaxis()->SetTitleFont(43);
    hData       ->GetYaxis()->SetTitle("#frac{d^{3}N_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hData       ->GetYaxis()->SetTitle("#frac{d^{3}N^{RAW}_{#phi#phi}}{dydp_{T}#phi_{1}dp_{T}#phi_{2}}(GeV/c)^{-1}");
    hData       ->GetYaxis()->SetTitleOffset(1.5);
    
    hDataSyst       ->SetMarkerStyle(22);
    hDataSyst       ->SetMarkerColor(38);
    hDataSyst       ->SetLineColor(33);

    
    // Define the Canvas
    TCanvas *cRatioPlot = new TCanvas( aName.c_str(), aName.c_str(), 800, 800);
    
    TLegend *Legend     = new TLegend(0.7,0.75,0.9,0.9);
    Legend      ->AddEntry(hData, "Data", "P");
    Legend      ->AddEntry(hData, "Stat. Err.", "E");
    if ( !bPythiaTest ) Legend      ->AddEntry(hDataSyst, "Syst. Err.", "F");
    if ( !bPythiaTest ) Legend      ->AddEntry(hPythia6, "Pythia6", "L");
    Legend      ->AddEntry(hPythia8, "Pythia8", "L");
    
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1, 1.0);
    pad1    ->SetLogy();
    pad1    ->SetBottomMargin(0); // Upper and lower plot are joined
    pad1    ->Draw();             // Draw the upper pad: pad1
    pad1    ->cd();               // pad1 becomes the current pad
    
    hData   ->SetStats(0);        // No statistics on upper plot
    hData   ->Draw("EP");
    if ( !bPythiaTest ) hDataSyst->Draw("SAME P 5");// Draw h1
    hData   ->Draw("SAME EP");
    if ( !bPythiaTest ) hPythia6->Draw("SAME HIST L");
    hPythia8->Draw("SAME HIST L");       // Draw h2 on top of h1
    Legend  ->Draw("same");
    
    //hData->GetYaxis()->SetLabelSize(0.);
    
    TGaxis *axis = new TGaxis( 0., 20, 10., 220, 20,220,510,"");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();

    // lower plot will be in pad
    cRatioPlot->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0., 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    TH1F *hPythia6Ratio = (TH1F*)hPythia6->Clone("hPythia6Ratio");
    hPythia6Ratio->Sumw2();     // No statistics on lower plot
    hPythia6Ratio->Divide(hData);
    hPythia6Ratio->SetLineStyle(1);
    TH1F *hPythia8Ratio = (TH1F*)hPythia8->Clone("hPythia8Ratio");
    hPythia8Ratio->Sumw2();
    hPythia8Ratio->Divide(hData);
    hPythia8Ratio->SetLineStyle(1);
    TH1F *hDataRatio = new TH1F(Form("hDataRatio_%s",aName.c_str()),"hDataRatio",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValue =   hData  ->GetBinContent(iTer);
        auto fError =   hData  ->GetBinError(iTer);
        hDataRatio->SetBinContent   (iTer,1);
        hDataRatio->SetBinError     (iTer,0.062);
        if ( fValue == 0 ) hDataRatio->SetBinError     (iTer,0);
    }
    TH1F * hDatS = SetSystErrorsh(hDataRatio);
    if ( bPythiaTest ) hDatS = hDataRatio;
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValu6 =   hPythia6Ratio  ->GetBinContent(iTer);
        auto fValu8 =   hPythia8Ratio  ->GetBinContent(iTer);
        auto fError =   hDatS  ->GetBinError(iTer);
        hPythia6Ratio->SetBinError     (iTer,fError*fValu6);
        hPythia8Ratio->SetBinError     (iTer,fError*fValu8);
    }
    
    hDataRatio->SetStats(0);  // Define Y ..
    hDataRatio->SetTitle("");
    hDataRatio->SetLineColorAlpha(kBlack,1.00);
    hDataRatio->SetFillColorAlpha(33,0.66);
    hDataRatio->Draw("HIST L E0 SAME");
    
    // Setting Axes
    // - // XAxis
    hDataRatio->GetXaxis()->SetTitle("P_{T}#phi_{1} (GeV/c)");
    hDataRatio->GetXaxis()->SetTitleSize(20);
    hDataRatio->GetXaxis()->SetTitleFont(43);
    hDataRatio->GetXaxis()->SetTitleOffset(3.);
    hDataRatio->GetXaxis()->SetLabelFont(43);
    hDataRatio->GetXaxis()->SetLabelSize(15);
    hDataRatio->SetMinimum(0.);                     // Range set
    hDataRatio->SetMaximum(2.9);
    // - // YAxis
    hDataRatio->GetYaxis()->SetTitle("MODEL / DATA");
    if ( bPythiaTest ) hDataRatio->GetYaxis()->SetTitle("Tru / Res");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hDataRatio->GetYaxis()->SetTitle("Rec / Raw");
    hDataRatio->GetYaxis()->SetNdivisions(505);
    hDataRatio->GetYaxis()->SetTitleSize(20);
    hDataRatio->GetYaxis()->SetTitleFont(43);
    hDataRatio->GetYaxis()->SetTitleOffset(2.);
    hDataRatio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hDataRatio->GetYaxis()->SetLabelSize(15);
    
    // Define the ratio plot
    hPythia6Ratio->SetMarkerStyle(0);
    hPythia6Ratio->SetStats(0);// No statistics on lower plot
    if ( !bPythiaTest ) hPythia6Ratio->Draw("SAME EP");       // Draw the ratio plot
    
    hPythia8Ratio->SetMarkerStyle(0);
    hPythia8Ratio->SetStats(0);
    hPythia8Ratio->Draw("SAME EP");       // Draw the ratio plot

    // Ratio plot (h3) settings
    hPythia6Ratio->SetTitle(""); // Remove the first ratio title
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    hPythia8Ratio->SetTitle(""); // Remove the second ratio presence
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    cRatioPlot->Write();
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.pdf", aName.c_str()));
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.png", aName.c_str()));
    delete cRatioPlot;
}

void ratioplot ( TH1F * hData, TH1F * hPythia6, TH1F * hPythia8, string aName = "" )
{
    // Setting Systematics
    TGraphErrors *hDataSyst = new TGraphErrors (SetSystErrors2(hData));
    hDataSyst->RemovePoint(0);
    hDataSyst->RemovePoint(0);
    if ( bPythiaTest ) aName += "P8";
    
    // Setting Aspect
    setLine                 (hPythia6,1);
    setLine                 (hPythia8,2);
    
    hData       ->SetMarkerStyle(22);
    hData       ->SetMarkerColor(38);
    hData       ->SetFillColorAlpha(33,0.66);
    hData       ->SetLineColorAlpha(kRed,1.0);
    hData       ->SetTitle("");
    hData       ->GetYaxis()->SetTitleSize(20);
    hData       ->GetYaxis()->SetTitleFont(43);
    hData       ->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}#phi}(GeV/c)^{-1}");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hData       ->GetYaxis()->SetTitle("#frac{d^{2}N^{RAWS}_{#phi}}{dydp_{T}#phi}(GeV/c)^{-1}");
    hData       ->GetYaxis()->SetTitleOffset(1.5);
    
    hDataSyst       ->SetMarkerStyle(22);
    hDataSyst       ->SetMarkerColor(38);
    hDataSyst       ->SetLineColor(33);

    
    // Define the Canvas
    TCanvas *cRatioPlot = new TCanvas( aName.c_str(), aName.c_str(), 800, 800);
    
    TLegend *Legend     = new TLegend(0.7,0.75,0.9,0.9);
    Legend      ->AddEntry(hData, "Data", "P");
    Legend      ->AddEntry(hData, "Stat. Err.", "E");
    if ( !bPythiaTest ) Legend      ->AddEntry(hDataSyst, "Syst. Err.", "F");
    if ( !bPythiaTest ) Legend      ->AddEntry(hPythia6, "Pythia6", "L");
    Legend      ->AddEntry(hPythia8, "Pythia8", "L");
    
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1, 1.0);
    pad1    ->SetLogy();
    pad1    ->SetBottomMargin(0); // Upper and lower plot are joined
    pad1    ->Draw();             // Draw the upper pad: pad1
    pad1    ->cd();               // pad1 becomes the current pad
    
    hData   ->SetStats(0);        // No statistics on upper plot
    hData   ->Draw("EP");
    if ( !bPythiaTest ) hDataSyst->Draw("SAME P 5");// Draw h1
    hData   ->Draw("SAME EP");
    if ( !bPythiaTest ) hPythia6->Draw("SAME HIST L");
    hPythia8->Draw("SAME HIST L");       // Draw h2 on top of h1
    Legend  ->Draw("same");
    
    //hData->GetYaxis()->SetLabelSize(0.);
    
    TGaxis *axis = new TGaxis( 0., 20, 10., 220, 20,220,510,"");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();

    // lower plot will be in pad
    cRatioPlot->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0., 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    TH1F *hPythia6Ratio = (TH1F*)hPythia6->Clone("hPythia6Ratio");
    hPythia6Ratio->Sumw2();     // No statistics on lower plot
    hPythia6Ratio->Divide(hData);
    hPythia6Ratio->SetLineStyle(1);
    TH1F *hPythia8Ratio = (TH1F*)hPythia8->Clone("hPythia8Ratio");
    hPythia8Ratio->Sumw2();
    hPythia8Ratio->Divide(hData);
    hPythia8Ratio->SetLineStyle(1);
    TH1F *hDataRatio = new TH1F(Form("hDataRatio_%s",aName.c_str()),"hDataRatio",nBinPT1D,fArrPT1D);
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValue =   hData  ->GetBinContent(iTer);
        auto fError =   hData  ->GetBinError(iTer);
        hDataRatio->SetBinContent   (iTer,1);
        hDataRatio->SetBinError     (iTer,0.062);
        if ( fValue == 0 ) hDataRatio->SetBinError     (iTer,0);
    }
    TH1F * hDatS = SetSystErrorsh(hDataRatio);
    if ( bPythiaTest ) hDatS = hDataRatio;
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValu6 =   hPythia6Ratio  ->GetBinContent(iTer);
        auto fValu8 =   hPythia8Ratio  ->GetBinContent(iTer);
        auto fError =   hDatS  ->GetBinError(iTer);
        hPythia6Ratio->SetBinError     (iTer,fError*fValu6);
        hPythia8Ratio->SetBinError     (iTer,fError*fValu8);
    }
    
    hDataRatio->SetStats(0);  // Define Y ..
    hDataRatio->SetTitle("");
    hDataRatio->SetLineColorAlpha(kBlack,1.00);
    hDataRatio->SetFillColorAlpha(33,0.66);
    hDataRatio->Draw("HIST L E0 SAME");
    
    // Setting Axes
    // - // XAxis
    hDataRatio->GetXaxis()->SetTitle("P_{T}#phi (GeV/c)");
    hDataRatio->GetXaxis()->SetTitleSize(20);
    hDataRatio->GetXaxis()->SetTitleFont(43);
    hDataRatio->GetXaxis()->SetTitleOffset(3.);
    hDataRatio->GetXaxis()->SetLabelFont(43);
    hDataRatio->GetXaxis()->SetLabelSize(15);
    hDataRatio->SetMinimum(0.);                     // Range set
    hDataRatio->SetMaximum(2.9);
    // - // YAxis
    hDataRatio->GetYaxis()->SetTitle("MODEL / DATA");
    if ( bPythiaTest ) hDataRatio->GetYaxis()->SetTitle("Tru / Res");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hDataRatio->GetYaxis()->SetTitle("Rec / Raw");
    hDataRatio->GetYaxis()->SetNdivisions(505);
    hDataRatio->GetYaxis()->SetTitleSize(20);
    hDataRatio->GetYaxis()->SetTitleFont(43);
    hDataRatio->GetYaxis()->SetTitleOffset(2.);
    hDataRatio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hDataRatio->GetYaxis()->SetLabelSize(15);
    
    // Define the ratio plot
    hPythia6Ratio->SetMarkerStyle(0);
    hPythia6Ratio->SetStats(0);// No statistics on lower plot
    if ( !bPythiaTest ) hPythia6Ratio->Draw("SAME EP");       // Draw the ratio plot
    
    hPythia8Ratio->SetMarkerStyle(0);
    hPythia8Ratio->SetStats(0);
    hPythia8Ratio->Draw("SAME EP");       // Draw the ratio plot

    // Ratio plot (h3) settings
    hPythia6Ratio->SetTitle(""); // Remove the first ratio title
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    hPythia8Ratio->SetTitle(""); // Remove the second ratio presence
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    cRatioPlot->Write();
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.pdf", aName.c_str()));
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.png", aName.c_str()));
    delete cRatioPlot;
}

void ratioplot ( Double_t***   aInput, TH1D * hPythia6, TH1D * hPythia8, string aName = "" )
{
    // Setting Systematics
    TH1D * hData            = new TH1D (*TH1DGeneratorYield(aInput,false));
    TGraphErrors *hDataSyst = new TGraphErrors (TH1DGeneratorYield(aInput,true));
    hDataSyst->RemovePoint(0);
    hDataSyst->RemovePoint(0);
    if ( bPythiaTest ) aName += "P8";
    
    // Setting Aspect
    setLine                 (hPythia6,1);
    setLine                 (hPythia8,2);
    
    hData       ->SetMarkerStyle(22);
    hData       ->SetMarkerColor(38);
    hData       ->SetFillColorAlpha(33,0.66);
    hData       ->SetLineColorAlpha(kRed,1.0);
    hData       ->SetTitle("");
    hData       ->GetYaxis()->SetTitleSize(20);
    hData       ->GetYaxis()->SetTitleFont(43);
    hData       ->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi#phi}}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
    hData       ->GetYaxis()->SetTitleOffset(1.5);
    
    hDataSyst       ->SetMarkerStyle(22);
    hDataSyst       ->SetMarkerColor(38);
    hDataSyst       ->SetLineColor(33);

    
    // Define the Canvas
    TCanvas *cRatioPlot = new TCanvas( aName.c_str(), aName.c_str(), 800, 800);
    
    TLegend *Legend     = new TLegend(0.7,0.75,0.9,0.9);
    Legend      ->AddEntry(hData, "Data", "P");
    Legend      ->AddEntry(hData, "Stat. Err.", "E");
    if ( !bPythiaTest ) Legend      ->AddEntry(hDataSyst, "Syst. Err.", "F");
    if ( !bPythiaTest ) Legend      ->AddEntry(hPythia6, "Pythia6", "L");
    Legend      ->AddEntry(hPythia8, "Pythia8", "L");
    
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0., 0.3, 1, 1.0);
    pad1    ->SetLogy();
    pad1    ->SetBottomMargin(0); // Upper and lower plot are joined
    pad1    ->Draw();             // Draw the upper pad: pad1
    pad1    ->cd();               // pad1 becomes the current pad
    
    hData   ->SetStats(0);        // No statistics on upper plot
    hData   ->Draw("EP");
   if ( !bPythiaTest )  hDataSyst->Draw("SAME P 5");// Draw h1
    hData   ->Draw("SAME EP");
    if ( !bPythiaTest ) hPythia6->Draw("SAME HIST L");
    hPythia8->Draw("SAME HIST L");       // Draw h2 on top of h1
    Legend  ->Draw("same");
    
    //hData->GetYaxis()->SetLabelSize(0.);
    
    TGaxis *axis = new TGaxis( 0., 20, 10., 220, 20,220,510,"");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();

    // lower plot will be in pad
    cRatioPlot->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0., 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    TH1F *hPythia6Ratio = (TH1F*)hPythia6->Clone("hPythia6Ratio");
    hPythia6Ratio->Sumw2();     // No statistics on lower plot
    hPythia6Ratio->Divide(hData);
    hPythia6Ratio->SetLineStyle(1);
    TH1F *hPythia8Ratio = (TH1F*)hPythia8->Clone("hPythia8Ratio");
    hPythia8Ratio->Sumw2();
    hPythia8Ratio->Divide(hData);
    hPythia8Ratio->SetLineStyle(1);
    TH1F *hDataRatio = new TH1F(Form("hDataRatio_%s",aName.c_str()),"hDataRatio",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValue =   hData  ->GetBinContent(iTer);
        auto fError =   hData  ->GetBinError(iTer);
        hDataRatio->SetBinContent   (iTer,1);
        hDataRatio->SetBinError     (iTer,0.062);
        if ( fValue == 0 ) hDataRatio->SetBinError     (iTer,0);
    }
    TH1F * hDatS = SetSystErrorsh(hDataRatio);
    if ( bPythiaTest ) hDatS = hDataRatio;
    for ( Int_t iTer = 1; iTer <= nBinPT1D; iTer++ )
    {
        auto fValu6 =   hPythia6Ratio  ->GetBinContent(iTer);
        auto fValu8 =   hPythia8Ratio  ->GetBinContent(iTer);
        auto fError =   hDatS  ->GetBinError(iTer);
        hPythia6Ratio->SetBinError     (iTer,fError*fValu6);
        hPythia8Ratio->SetBinError     (iTer,fError*fValu8);
    }
    
    hDataRatio->SetStats(0);  // Define Y ..
    hDataRatio->SetTitle("");
    hDataRatio->SetLineColorAlpha(kBlack,1.00);
    hDataRatio->SetFillColorAlpha(33,0.66);
    hDataRatio->Draw("HIST L E0 SAME");
    
    // Setting Axes
    // - // XAxis
    hDataRatio->GetXaxis()->SetTitle("P_{T}#phi_{2} (GeV/c)");
    hDataRatio->GetXaxis()->SetTitleSize(20);
    hDataRatio->GetXaxis()->SetTitleFont(43);
    hDataRatio->GetXaxis()->SetTitleOffset(3.);
    hDataRatio->GetXaxis()->SetLabelFont(43);
    hDataRatio->GetXaxis()->SetLabelSize(15);
    hDataRatio->SetMinimum(0.);                     // Range set
    hDataRatio->SetMaximum(2.9);
    // - // YAxis
    hDataRatio->GetYaxis()->SetTitle("MODEL / DATA");
    if ( bPythiaTest ) hDataRatio->GetYaxis()->SetTitle("Tru / Res");
    if ( bPythiaTest && aName.find("Raw") != -1 ) hDataRatio->GetYaxis()->SetTitle("Rec / Raw");
    hDataRatio->GetYaxis()->SetNdivisions(505);
    hDataRatio->GetYaxis()->SetTitleSize(20);
    hDataRatio->GetYaxis()->SetTitleFont(43);
    hDataRatio->GetYaxis()->SetTitleOffset(2.);
    hDataRatio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hDataRatio->GetYaxis()->SetLabelSize(15);
    
    // Define the ratio plot
    hPythia6Ratio->SetMarkerStyle(0);
    hPythia6Ratio->SetStats(0);// No statistics on lower plot
    if ( !bPythiaTest ) hPythia6Ratio->Draw("SAME EP");       // Draw the ratio plot
    
    hPythia8Ratio->SetMarkerStyle(0);
    hPythia8Ratio->SetStats(0);
    hPythia8Ratio->Draw("SAME EP");       // Draw the ratio plot

    // Ratio plot (h3) settings
    hPythia6Ratio->SetTitle(""); // Remove the first ratio title
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    hPythia8Ratio->SetTitle(""); // Remove the second ratio presence
    hPythia8Ratio->GetXaxis()->SetTitleSize(0);
    hPythia8Ratio->GetYaxis()->SetTitleSize(0);
    
    cRatioPlot->Write();
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.pdf", aName.c_str()));
    cRatioPlot->SaveAs(Form("./graphs/Ratio_%s.png", aName.c_str()));
    delete cRatioPlot;
}
