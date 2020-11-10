// File for 1,2-Dimensional Fit:
// !TODO: All Set!

#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

void Anls_InvariantMassFit ( bool fSilent = false )
{
    //---------------------//
    //  Setting up input   //-------------------------------------------------------------------------------
    //---------------------//
    
    //-// OPTIONS
    
    // Silencing warnings for smoother 
    if ( fSilent )
    {
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(fSilent);
    }
    
    // Opening Data File
    TFile*  insFile_DT              =   new TFile   (fInvMasHist);
    
    //-// Recovering Data histograms
    
    // Generating local Histograms to work on
    TH1F *          Entry_DT            = new TH1F          ();
    TH1F **         hIM_1D_Rec_PT_S     = new TH1F *        [nBinPT1D];
    TH1F **         hdM_dpT_Tot_Rec     = new TH1F *        [nBinPT2D];
    TH2F ***        hdM_dpT_Tot_Rec2D   = new TH2F **       [nBinPT2D];
    
    //------ 1D Histograms Recovery ------//
    
    hName = "Entry_DT";
    Entry_DT    =   (TH1F*)(insFile_DT->Get(hName));
    
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hName                       =   Form("hIM_1D_Rec_PT_B_S_%i",iHisto);    // Name of the histogram in the preprocessed file
        hIM_1D_Rec_PT_S[iHisto]     =   (TH1F*)(insFile_DT->Get(hName));
    }
    
    //------ 2D Histograms Recovery ------//
    
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        // Importing 1D Invariant Mass Data Histograms for 2D Fit
        hName                       =   Form("hIM_2D_Rec_PT_B_S_%i",iHisto);    // Name of the histogram in the preprocessed file
        hdM_dpT_Tot_Rec[iHisto]     =   (TH1F*)(insFile_DT->Get(hName));
        
        // 2D set-up
        hdM_dpT_Tot_Rec2D[iHisto]   =   new TH2F * [nBinPT2D];
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            // Importing 2D Invariant Mass Data Histograms
            hName                               =   Form("hIM_2D_Rec_PT_PT_BB_BS_SB_SS_%i_%i",iHisto,jHisto);   // Name of the histogram in the preprocessed file
            hdM_dpT_Tot_Rec2D[iHisto][jHisto]   =   (TH2F*)(insFile_DT->Get(hName));
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
    TFile*  outFile_FT  =   new TFile(fFitResHist,"recreate");
    
    // Fit Results and PlotOn object
    RooFitResult *** Results = new RooFitResult **  [nBinPT2D];
    RooFitResult **  utility = new RooFitResult *   [nBinPT2D];
    
    //------ 1D Histogram of N_raw ------//
    
    for (Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
    {
        // Not considering pT < 0.4 GeV
        if ( fArrPT1D[iFit+1] <= 0.41 ) continue;
        
        // Fit
        auto fResults = FitModel(hIM_1D_Rec_PT_S[iFit],"",bSave,iFit,1);        // FitModel with Save enabled will write every fit performed on a Canvas
        
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
        utility[iFit]   = FitModel(hdM_dpT_Tot_Rec[iFit],"STD",bSave,iFit,2);
    }
    
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
