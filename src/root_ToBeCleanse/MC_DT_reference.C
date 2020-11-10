// File for 2-Dimensional Analysis:
// !TODO: All Set!

#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

void MCreference (Char_t * MCFileIn, Char_t * MCFileOut)
{
    //Retrieving MC data
    TFile *insFileMC    =   new TFile   (MCFileIn);
    
    // Retrieving Histograms
    hName   = "hNP_1D_Tru_PT_S";
    TH1F *  MC1D_Curve_TRU  =   (TH1F*)(insFileMC->Get(hName));
    hName   = "hNP_2D_Tru_PT_S";
    TH2F *  MC2D_Curve_TRU  =   (TH2F*)(insFileMC->Get(hName));
    
    // Line Attibutes
    Int_t kColor    = kBlue+3;
    Int_t kStyle    = 10;
    Int_t kWidth    = 5;
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    // Creating the histograms-------------------------------------------------------------------------------
    
    TH1D *  MC1D_Util0_TRU  = new TH1D();
    TH1D *  MC1D_Util1_TRU  = new TH1D("","",nBinPT2D,fArrPT2D);
    TH1D *  MC1D_Util2_TRU  = new TH1D("","",nBinPT2D,fArrPT2D);
    
    // Output file
    TFile *outFileMC    =   new TFile   (MCFileOut,"recreate");
    
    gStyle->SetOptStat(0);
    
    TCanvas * c1 = new TCanvas();
    gPad->SetLogy();
    MC1D_Curve_TRU->SetName(Form("XProjection"));
    MC1D_Curve_TRU->SetTitle(Form("XProjection"));
    MC1D_Curve_TRU->GetXaxis()->SetTitle("");
    MC1D_Curve_TRU->GetXaxis()->SetTitle("");
    MC1D_Curve_TRU->SetLineColor(kColor);
    MC1D_Curve_TRU->SetLineStyle(kStyle);
    MC1D_Curve_TRU->SetLineWidth(kWidth);
    MC1D_Curve_TRU->Draw("HIST L");
    c1->Write();
    c1->Clear();
    
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        c1->SetName(Form("X_%i",iHisto));
        MC1D_Util0_TRU = MC2D_Curve_TRU->ProjectionX("",iHisto,iHisto+1);
        MC1D_Util0_TRU->SetName(Form("XProjection"));
        MC1D_Util0_TRU->SetTitle(Form("XProjection"));
        MC1D_Util0_TRU->GetXaxis()->SetTitle("");
        MC1D_Util0_TRU->GetXaxis()->SetTitle("");
        MC1D_Util0_TRU->SetLineColor(kColor);
        MC1D_Util0_TRU->SetLineStyle(kStyle);
        MC1D_Util0_TRU->SetLineWidth(kWidth);
        MC1D_Util0_TRU->Draw("HIST L");
        c1->Write();
        c1->Clear();
    
        MC1D_Util1_TRU->SetBinContent(iHisto+1,MC1D_Util0_TRU->Integral());
        
        c1->SetName(Form("Y_%i",iHisto));
        MC1D_Util0_TRU = MC2D_Curve_TRU->ProjectionX("",iHisto,iHisto+1);
        MC1D_Util0_TRU->SetName(Form("XProjection"));
        MC1D_Util0_TRU->SetTitle(Form("XProjection"));
        MC1D_Util0_TRU->GetXaxis()->SetTitle("");
        MC1D_Util0_TRU->GetXaxis()->SetTitle("");
        MC1D_Util0_TRU->SetLineColor(kColor);
        MC1D_Util0_TRU->SetLineStyle(kStyle);
        MC1D_Util0_TRU->SetLineWidth(kWidth);
        MC1D_Util0_TRU->Draw("HIST L");
        c1->Write();
        c1->Clear();
        
        MC1D_Util2_TRU->SetBinContent(iHisto+1,MC1D_Util0_TRU->Integral());
    }
    
    c1->SetName("FinalX");
    MC1D_Util1_TRU->SetName(Form("XProjection"));
    MC1D_Util1_TRU->SetTitle(Form("XProjection"));
    MC1D_Util1_TRU->GetXaxis()->SetTitle("");
    MC1D_Util1_TRU->GetXaxis()->SetTitle("");
    MC1D_Util1_TRU->SetLineColor(kColor);
    MC1D_Util1_TRU->SetLineStyle(kStyle);
    MC1D_Util1_TRU->SetLineWidth(kWidth);
    MC1D_Util1_TRU->Draw("HIST L");
    c1->Write();
    c1->Clear();
    
    c1->SetName("FinalY");
    MC1D_Util2_TRU->SetName(Form("XProjection"));
    MC1D_Util2_TRU->SetTitle(Form("XProjection"));
    MC1D_Util2_TRU->GetXaxis()->SetTitle("");
    MC1D_Util2_TRU->GetXaxis()->SetTitle("");
    MC1D_Util2_TRU->SetLineColor(kColor);
    MC1D_Util2_TRU->SetLineStyle(kStyle);
    MC1D_Util2_TRU->SetLineWidth(kWidth);
    MC1D_Util2_TRU->Draw("HIST L");
    c1->Write();
    c1->Clear();
    
    delete c1;
    
    return;
}
