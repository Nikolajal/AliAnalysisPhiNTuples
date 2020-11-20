// File for 1,2-Dimensional Fit:
// !TODO: All Set!

#include "../inc/AliAnalysisPhiPair.h"

void Anls_InvariantMassFit ( bool fSilent = false )
{
    //-// OPTIONS
    
    // Silencing warnings for smoother 
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(kTRUE);
    }
    
    //---------------------//
    //  Setting up input   //-------------------------------------------------------------------------------
    //---------------------//
    
    // Opening Data File
    TFile      *insFile_DT              =   new TFile   (fInvMasHist);
    
    //-// Recovering Data histograms
    // 1D
    TH1F       *hREC_1D                 =   (TH1F*)(insFile_DT->Get("hREC_1D"));
    TH1F      **hREC_1D_in_PT           =   new TH1F     *[nBinPT1D];
    TH1F       *Entry_DT                =   (TH1F*)(insFile_DT->Get("Entry_DT"));

    // 2D
    TH2F       *hREC_2D                 =   (TH1F*)(insFile_DT->Get("hREC_2D"));
    TH1F      **hREC_1D_in_PT_2D_bin    =   new TH1F     *[nBinPT2D];
    TH2F     ***hREC_2D_in_PT           =   new TH2F    **[nBinPT2D];
    
    //------ 1D Histograms Recovery ------//
    
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hName                           =   Form("hREC_1D_in_PT_%i",iHisto);    // Name of the histogram in the preprocessed file
        hREC_1D_in_PT[iHisto]           =   (TH1F*)(insFile_DT->Get(hName));
    }
    
    //------ 2D Histograms Recovery ------//
    
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        // Importing 1D Invariant Mass Data Histograms for 2D Fit
        hName                           =   Form("hREC_1D_in_PT_2D_bin_%i",iHisto);   // Name of the histogram in the preprocessed file
        hREC_1D_in_PT_2D_bin[iHisto]    =   (TH1F*)(insFile_DT->Get(hName));
        
        // 2D set-up
        hREC_2D_in_PT[iHisto]           =   new TH2F * [nBinPT2D];
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            // Importing 2D Invariant Mass Data Histograms
            hName                           =   Form("hREC_2D_in_PT_%i_%i",iHisto,jHisto);  // Name of the histogram in the preprocessed file
            hREC_2D_in_PT[iHisto][jHisto]   =   (TH2F*)(insFile_DT->Get(hName));
        }
    }
    
    //---------------------//
    //  Setting up output  //-------------------------------------------------------------------------------
    //---------------------//
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    //---------- 1D Histograms ----------//
    
    hName   = "h1D_Raw";
    hTitle  = "Number of #phi found in Invariant Mass Fit per |y|<5 in K_{+}K_{-} Decay mode";
    TH1F*   h1D_Raw     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    TH1_SetDescription(h1D_Raw);
    
    //---------- 2D Histograms ----------//
    
    hName   = "h2D_Raw";
    hTitle  = "Number of #phi found in Invariant Mass Fit per |y|<5 in K_{+}K_{-} Decay mode";
    TH2F*   h2D_Raw     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2_SetDescription(h2D_Raw);
    
    //---------------------//
    //      Analysis       //-------------------------------------------------------------------------------
    //---------------------//
    
    // Output File for Fit Check
    TFile*  outFile_FT          =   new TFile(fFitResHist,"recreate");
    
    // Fit Results and PlotOn object
    RooFitResult  ***Results    =   new RooFitResult  **[nBinPT2D];
    RooFitResult   **utility    =   new RooFitResult   *[nBinPT2D];
    
    //------ 1D Histogram of N_raw ------//
    
    for (Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
    {
        // Not considering pT < 0.4 GeV
        if ( fArrPT1D[iFit+1] <= 0.41 ) continue;
        
        // Fit
        auto fResults = FitModel(hREC_1D_in_PT[iFit],Form("PT=%.2i S 1D",iFit));        // FitModel with Save enabled will write every fit performed on a Canvas
        
        // Building N_Raw histogram
        auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(Signal));
        h1D_Raw->SetBinContent      (iFit+1,N_Raw->getVal());
        h1D_Raw->SetBinError        (iFit+1,N_Raw->getError());
    }
    
    //------ 2D Histogram of N_raw ------//
    
    // Preparing 1D fit for shape in pT bins
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        Results[iFit]   = new RooFitResult * [nBinPT2D];
        utility[iFit]   = FitModel(hREC_1D_in_PT_2D_bin[iFit],Form("PT=%.2i S 12D",iFit));
    }
    
    /*
    // 2D Fits
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
        for (int jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            // Not considering pT < 0.4 GeV
            if ( fArrPT2D[jFit+1] <= 0.41 ) continue;
            if ( iFit == 11 && jFit == 2  ) continue;
            if ( iFit == 2  && jFit == 11 ) continue;
            
            // Fit
            Results[iFit][jFit] = FitModel(hdM_dpT_Tot_Rec2D[iFit][jFit],utility[iFit],utility[jFit],"STD",bSave,iFit,jFit);
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(Results[iFit][jFit]->floatParsFinal().at(SignlSignl));
            h2D_Raw->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            h2D_Raw->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
    */
    //---------------------//
    // Output and wrap up  //-------------------------------------------------------------------------------
    //---------------------//
    
    // Output File for Results
    TFile*  outFile_RS  =   new TFile(fFitResults,"recreate");
    
    // Writing Results to file
    h1D_Raw->Write();
    h2D_Raw->Write();
    Entry_DT->Write();
    
    // Closing opened files
    insFile_DT->Close();
    outFile_RS->Close();
    outFile_FT->Close();
    outFile_RS->Close();
}
