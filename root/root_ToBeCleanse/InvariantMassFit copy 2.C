// File for 1-Dimensional Analysis:
// !TODO: All Set!


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"
#include "RooMsgService.h"

void InvariantMassFit ( bool fSilent = false )
{
    // Silencing warnings for smoother 
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    
    //Retrieving PreProcessed Data
    TFile*  insFile_DT              =   new TFile   (fInvMasHist);
    
    // Fit Variables for Roofit
    RooRealVar  xInvMass2D ("xInvMass2D","xInvMass2D",fMinIM2D,fMaxIM2D);
    RooRealVar  yInvMass2D ("yInvMass2D","yInvMass2D",fMinIM2D,fMaxIM2D);
    
    //Recovering histograms
    TH1F ** hIM_1D_Rec_PT_S  = new TH1F * [nBinPT1D];
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        auto hName                  =   Form("hIM_1D_Rec_PT_B_S_%i",iHisto);
        hIM_1D_Rec_PT_S[iHisto]     =   (TH1F*)(insFile_DT->Get(hName));
    }
    
    TH1F **  hdM_dpT_Tot_Rec     = new TH1F *  [nBinPT2D];
    RooDataHist *** hdM_dpT_Tot_Rec2D   = new RooDataHist ** [nBinPT2D];
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        // Importing 1D Invariant Mass Data Histograms for 2D Fit
        hName = Form("hIM_2D_Rec_PT_B_S_%i",iHisto);
        hdM_dpT_Tot_Rec[iHisto] = (TH1F*)(insFile_DT->Get(hName));
        
        // 2D set-up
        hdM_dpT_Tot_Rec2D[iHisto]   = new RooDataHist * [nBinPT2D];
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            // Importing 2D Invariant Mass Data Histograms
            hName = Form("hIM_2D_Rec_PT_PT_BB_BS_SB_SS_%i_%i",iHisto,jHisto);
            hdM_dpT_Tot_Rec2D[iHisto][jHisto]   = new RooDataHist(hName,hName,RooArgList(xInvMass2D,yInvMass2D),Import(*(TH2F*)(insFile_DT->Get(hName))));
        }
    }
    
    TH1F * hUtlEntry    = (TH1F*)(insFile_EF->Get("Entry_MC_DT"));
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Creating the histograms-------------------------------------------------------------------------------
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    //--- 1D Histograms ---//
    TH1F*   hNP_1D_Eff_PT_S     =   (TH1F*)(insFile_EF->Get("hNP_1D_Eff_PT_S"));
    
    TH1F*   hNP_1D_Tru_PT_S     =   (TH1F*)(insFile_EF->Get("hNP_1D_Tru_PT_S"));
    
    TH1F*   hNP_1D_Rec_PT_S     =   (TH1F*)(insFile_EF->Get("hNP_1D_Rec_PT_S"));
    
    hName   = "hNP_1D_Raw_PT_S";
    hTitle  = "Number of #phi found in Fit per |y|<5 in K_{+}K_{-} Decay mode, 1D analysis";
    TH1F*   hNP_1D_Raw_PT_S     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    hNP_1D_Raw_PT_S->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
    hNP_1D_Raw_PT_S->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
    hNP_1D_Raw_PT_S->GetXaxis()->SetTitleOffset(1.15);
    hNP_1D_Raw_PT_S->GetYaxis()->SetTitleOffset(1.15);
    
    hName   = "hNP_1D_Res_PT_S";
    hTitle  = "Number of #phi reconstructed from Fit per |y|<5, 1D analysis";
    TH1F*   hNP_1D_Res_PT_S     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    hNP_1D_Res_PT_S->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
    hNP_1D_Res_PT_S->GetYaxis()->SetTitle("#frac{d^{2}N_{#phi}}{dydp_{T}}(GeV/c)^{-1}");
    hNP_1D_Res_PT_S->GetXaxis()->SetTitleOffset(1.15);
    hNP_1D_Res_PT_S->GetYaxis()->SetTitleOffset(1.15);
    
    hName   = "hNP_1D_Chk_RT_S";
    hTitle  = "Ratio N^{Res}/N^{Tru}, 1D analysis";
    TH1F*   hNP_1D_Chk_RT_S     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    hNP_1D_Chk_RT_S->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
    hNP_1D_Chk_RT_S->GetYaxis()->SetTitle("#frac{N^{Res}}{N^{Tru}}");
    hNP_1D_Chk_RT_S->GetXaxis()->SetTitleOffset(1.15);
    hNP_1D_Chk_RT_S->GetYaxis()->SetTitleOffset(1.15);
    hNP_1D_Chk_RT_S->SetTitleSize(0.75,"t");
    
    hName   = "hNP_1D_Chk_RR_S";
    hTitle  = "Ratio N^{Raw}/N^{Rec}, 1D analysis";
    TH1F*   hNP_1D_Chk_RR_S     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    hNP_1D_Chk_RR_S->GetXaxis()->SetTitle("p_{T} #phi (GeV/c)");
    hNP_1D_Chk_RR_S->GetYaxis()->SetTitle("#frac{N^{Raw}}{N^{Rec}}");
    hNP_1D_Chk_RR_S->GetXaxis()->SetTitleOffset(1.15);
    hNP_1D_Chk_RR_S->GetYaxis()->SetTitleOffset(1.15);
    hNP_1D_Chk_RR_S->SetTitleSize(0.75,"t");
    
    //--- 2D Histograms ---//
 
    TH2F * hNP_2D_Eff_PT_S  = (TH2F*)(insFile_EF->Get("hNP_2D_Eff_PT_S"));
    
    TH2F * hNP_2D_Tru_PT_S  = (TH2F*)(insFile_EF->Get("hNP_2D_Tru_PT_S"));
    
    TH2F * hNP_2D_Rec_PT_S  = (TH2F*)(insFile_EF->Get("hNP_2D_Rec_PT_S"));
    
    hName   = "hNP_2D_Raw_PT_S";
    hTitle  = "Number of #phi found in Fit per |y|<5 in K_{+}K_{-} Decay mode, 2D analysis";
    TH2F*   hNP_2D_Raw_PT_S     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hNP_2D_Raw_PT_S->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    hNP_2D_Raw_PT_S->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
    hNP_2D_Raw_PT_S->GetXaxis()->SetTitleOffset(1.5);
    hNP_2D_Raw_PT_S->GetYaxis()->SetTitleOffset(1.5);
    
    hName   = "hNP_2D_Res_PT_S";
    hTitle  = "Number of #phi reconstructed from Fit per |y|<5, 2D analysis";
    TH2F*   hNP_2D_Res_PT_S     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hNP_2D_Res_PT_S->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    hNP_2D_Res_PT_S->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
    hNP_2D_Res_PT_S->GetXaxis()->SetTitleOffset(1.5);
    hNP_2D_Res_PT_S->GetYaxis()->SetTitleOffset(1.5);
    
    hName   = "hNP_2D_Chk_RT_S";
    hTitle  = "Ratio N^{Res}/N^{Tru}, 2D analysis";
    TH2F*   hNP_2D_Chk_RT_S     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hNP_2D_Chk_RT_S->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    hNP_2D_Chk_RT_S->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
    hNP_2D_Chk_RT_S->GetXaxis()->SetTitleOffset(1.5);
    hNP_2D_Chk_RT_S->GetYaxis()->SetTitleOffset(1.5);
    
    hName   = "hNP_2D_Chk_RR_S";
    hTitle  = "Ratio  N^{Raw}/N^{Rec}, 2D analysis";
    TH2F*   hNP_2D_Chk_RR_S     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    hNP_2D_Chk_RR_S->GetXaxis()->SetTitle("p_{T} #phi_{1} (GeV/c)");
    hNP_2D_Chk_RR_S->GetYaxis()->SetTitle("p_{T} #phi_{2} (GeV/c)");
    hNP_2D_Chk_RR_S->GetXaxis()->SetTitleOffset(1.5);
    hNP_2D_Chk_RR_S->GetYaxis()->SetTitleOffset(1.5);
    
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    // Output File for Fit Check
    TFile*  outFileFit  =   new TFile(fFitResults,"recreate");
    
    // Fit Results and PlotOn object
    RooFitResult *** Results = new RooFitResult **  [nBinPT2D];
    RooFitResult **  utility = new RooFitResult *   [nBinPT2D];
    
    // Fit  to  data, N_raw
    for (Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
    {
        // Not considering pT < 0.4 GeV
        if ( fArrPT1D[iFit+1] <= 0.41 ) continue;
        auto fResults = FitModel(hIM_1D_Rec_PT_S[iFit],"",bSave,iFit,1);
        
        // Building N_Raw histogram
        auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal____));
        hNP_1D_Raw_PT_S->SetBinContent      (iFit+1,N_Raw->getVal());
        hNP_1D_Raw_PT_S->SetBinError        (iFit+1,N_Raw->getError());
    }
    
    // Fit  to  data, N_raw
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        Results[iFit]   = new RooFitResult * [nBinPT2D];
        utility[iFit]   = FitModel(hdM_dpT_Tot_Rec[iFit],"",bSave,iFit,2);
    }
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        //break;
        for (int jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
            if ( fArrPT2D[jFit+1] <= 0.41 ) continue;
            Float_t nTOT = hdM_dpT_Tot_Rec2D[iFit][jFit]->sumEntries();
            Results[iFit][jFit] = FitModel(utility[iFit],utility[jFit],hdM_dpT_Tot_Rec2D[iFit][jFit],xInvMass2D,yInvMass2D,bSave,iFit,jFit);
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(Results[iFit][jFit]->floatParsFinal().at(SignlSignl));
            hNP_2D_Raw_PT_S->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            hNP_2D_Raw_PT_S->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
    
    // // Necessary scaling
    // Scaling in pT
    hNP_1D_Raw_PT_S->Scale(1.,"width");
    hNP_2D_Raw_PT_S->Scale(1.,"width");
    
    // Divided by events
    hNP_1D_Raw_PT_S->Scale(1./(hUtlEntry->GetBinContent(2)));
    hNP_2D_Raw_PT_S->Scale(1./(hUtlEntry->GetBinContent(2)));
    
    // Producing of N_res by applying corrections to N_raw
    hNP_1D_Res_PT_S->Divide(hNP_1D_Raw_PT_S,hNP_1D_Eff_PT_S,1.,kS1DEf*kBRKK,"");
    hNP_2D_Res_PT_S->Divide(hNP_2D_Raw_PT_S,hNP_2D_Eff_PT_S,1.,kS2DEf*kBRKK*kBRKK,"");
    
    // N_Chk
    hNP_1D_Chk_RT_S->Divide(hNP_1D_Res_PT_S,hNP_1D_Tru_PT_S,1.,1.,"");
    hNP_1D_Chk_RR_S->Divide(hNP_1D_Raw_PT_S,hNP_1D_Rec_PT_S,1.,1.,"");
    hNP_2D_Chk_RT_S->Divide(hNP_2D_Res_PT_S,hNP_2D_Tru_PT_S,1.,1.,"");
    hNP_2D_Chk_RR_S->Divide(hNP_2D_Raw_PT_S,hNP_2D_Rec_PT_S,1.,1.,"");
    
    // Output File for Final Results
    TFile*  outFile_RS  =   new TFile(fFitResHist,"recreate");
    
    hNP_1D_Eff_PT_S->Write();
    hNP_1D_Tru_PT_S->Write();
    hNP_1D_Raw_PT_S->Write();
    hNP_1D_Res_PT_S->Write();
    hNP_1D_Chk_RT_S->Write();
    hNP_1D_Chk_RR_S->Write();
    hNP_2D_Eff_PT_S->Write();
    hNP_2D_Tru_PT_S->Write();
    hNP_2D_Raw_PT_S->Write();
    hNP_2D_Res_PT_S->Write();
    hNP_2D_Chk_RT_S->Write();
    hNP_2D_Chk_RR_S->Write();
    
    
    insFile_DT->Close();
    insFile_EF->Close();
    outFile_RS->Close();
    outFileFit->Close();
}
