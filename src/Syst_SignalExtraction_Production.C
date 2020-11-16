
// File for 1-Dimensional Analysis:
// !TODO: All Set!


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"
#include "RooMsgService.h"

void Syst_SignalExtraction_Production ( bool fSilent = false)
{
    if ( fSilent )
    {
        RooMsgService::instance().setSilentMode(fSilent);
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    }
    
    //Retrieving PreProcessed Data
    TFile*  insFile_DT              =   new TFile   (fInvMasHist);
    
    //Target variables
    RooRealVar  xInvMass2D ("xInvMass2D","xInvMass2D",fMinIM2D,fMaxIM2D);
    RooRealVar  yInvMass2D ("yInvMass2D","yInvMass2D",fMinIM2D,fMaxIM2D);
    
    //Recovering histograms
    TH1F ** hIM_1D_Rec_PT_S  = new TH1F * [nBinPT1D];
    for ( Int_t iHisto = 4; iHisto < nBinPT1D; iHisto++)
    {
        hName                  =   Form("hIM_1D_Rec_PT_B_S_%i",iHisto);
        hIM_1D_Rec_PT_S[iHisto]     =   (TH1F*)(insFile_DT->Get(hName));
    }
    
    TH1F **  hdM_dpT_Tot_Rec     = new TH1F *  [nBinPT2D];
    TH2F *** hdM_dpT_Tot_Rec2D   = new TH2F ** [nBinPT2D];
    for (int iHisto = 2; iHisto < nBinPT2D; iHisto++)
    {
        // Importing 1D Histograms
        hName = Form("hIM_2D_Rec_PT_B_S_%i",iHisto);
        hdM_dpT_Tot_Rec[iHisto] = (TH1F*)(insFile_DT->Get(hName));
        
        // 2D set-up
        hdM_dpT_Tot_Rec2D[iHisto]   = new TH2F * [nBinPT2D];
        for (int jHisto = 2; jHisto < nBinPT2D; jHisto++)
        {
            // Importing 2D Histograms
            hName = Form("hIM_2D_Rec_PT_PT_BB_BS_SB_SS_%i_%i",iHisto,jHisto);
            hdM_dpT_Tot_Rec2D[iHisto][jHisto]   = (TH2F*)(insFile_DT->Get(hName));
        }
    }
    hName                           =   "Entry_DT";
    TH1F *  h1D_EDT                 =   (TH1F*)(insFile_DT->Get(hName));

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
    
    hName   = "h1D_Raw";
    hTitle  = "Number of #phi found in Invariant Mass Fit per |y|<5 in K_{+}K_{-} Decay mode";
    TH1F*   h1D_Raw     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    TH1_SetDescription(h1D_Raw);
    
    const   Bool_t  fStandar    =   true;
    const   Bool_t  f1DOptio    =   true;
    const   Bool_t  f2DOptio    =   true;
    const   Int_t   nOptions    =   0; //15
    const   string  sOptions[]  =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ","RK","M","W","CH3","CH5"};
    const   Int_t   nOption2    =   0; //11
    const   string  sOption2[]  =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ","RK"};
    
    TH1F**  h1D_Syst    =   new TH1F*   [nOptions];
    for ( Int_t iTer = 0; iTer < nOptions; iTer++ )
    {
        h1D_Syst[iTer]  =   new TH1F (Form("1D_%s",sOptions[iTer].c_str()),Form("1D_%s",sOptions[iTer].c_str()),nBinPT1D,fArrPT1D);
    }
    
    TH2F**  h2D_Syst    =   new TH2F*   [nOptions+nOption2];
    for ( Int_t iTer = 0; iTer < nOptions; iTer++ )
    {
        h2D_Syst[iTer]  =   new TH2F (Form("1D-2D_%s",sOptions[iTer].c_str()),Form("1D-2D_%s",sOptions[iTer].c_str()),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    }
    for ( Int_t iTer = nOptions; iTer < nOption2+nOptions; iTer++ )
    {
        h2D_Syst[iTer]  =   new TH2F (Form("2D_%s",sOption2[iTer-nOptions].c_str()),Form("2D_%s",sOption2[iTer-nOptions].c_str()),nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    }
    
    //---------- 2D Histograms ----------//
    
    hName   = "h2D_Raw";
    hTitle  = "Number of #phi found in Invariant Mass Fit per |y|<5 in K_{+}K_{-} Decay mode";
    TH2F*   h2D_Raw     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2_SetDescription(h2D_Raw);
    
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    //------ 1D Histogram of N_raw ------//
    
    cout << "Starting evaluating 1D STD" << endl;
    
    TFile*  outFileR1D;
    if ( fStandar ) outFileR1D  =   new TFile(Form("./result/N_Raw_1D.root"),"recreate");
    for (Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
    {
        if ( !fStandar ) break;
        // Not considering pT < 0.4 GeV
        if ( fArrPT1D[iFit+1] <= 0.41 ) continue;
        
        auto fFitResult = FitModel(hIM_1D_Rec_PT_S[iFit],"1D_STD",bSave,iFit,1,"");
        
        // Building N_Raw histogram
        auto N_Raw      = static_cast<RooRealVar*>(fFitResult->floatParsFinal().at(Signal));
        h1D_Raw->SetBinContent      (iFit+1,N_Raw->getVal());
        h1D_Raw->SetBinError        (iFit+1,N_Raw->getError());
    }
    if ( fStandar )
    {
        h1D_Raw->Scale(1./(h1D_EDT->GetBinContent(1)),"width");
        h1D_Raw->Write();
        outFileR1D->Close();
    }
    
    for ( Int_t iSys = 0; iSys < nOptions; iSys++ )
    {
        if ( !f1DOptio )    break;
        TFile*  outFileFit  =   new TFile(Form("./result/Systematics_1D_%s.root",sOptions[iSys].c_str()),"recreate");
        cout << "Starting evaluating 1D " << sOptions[iSys].c_str() << endl;
        for (Int_t iFit = 0; iFit < nBinPT1D; iFit++ )
        {
            // Not considering pT < 0.4 GeV
            if ( fArrPT1D[iFit+1] <= 0.41 ) continue;
            
            auto fFitResult = FitModel(hIM_1D_Rec_PT_S[iFit],Form("1D_%s",sOptions[iSys].c_str()),bSave,iFit,1,sOptions[iSys].c_str());
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fFitResult->floatParsFinal().at(Signal));
            h1D_Syst[iSys]->SetBinContent      (iFit+1,N_Raw->getVal());
            h2D_Syst[iSys]->SetBinError        (iFit+1,N_Raw->getError());
            
        }
        h1D_Syst[iSys]->Scale(1./(h1D_EDT->GetBinContent(1)),"width");
        h1D_Syst[iSys]->Write();
        outFileFit->Close();
    }
    
    
    //------ 2D Histogram of N_raw ------//
    
    // Fit Results and PlotOn object
    RooFitResult  **utility = new RooFitResult *    [nBinPT2D];
    
    TFile*  outFileR2D;
    if ( fStandar ) outFileR2D  =   new TFile(Form("./result/N_Raw_2D.root"),"recreate");
    cout << "Starting evaluating 1D-2D STD" << endl;
    // Preparing 1D fit for shape in pT bins
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        if ( !fStandar && !(nOption2 != 0 && f2DOptio) ) break;
        if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
        utility[iFit]  = FitModel(hdM_dpT_Tot_Rec[iFit],"2D_STD",bSave,iFit,2);
    }
    
    cout << "Starting evaluating 2D STD" << endl;
    // 2D Fits
    for (int iFit = 0; iFit < nBinPT2D; iFit++ )
    {
        if ( !fStandar ) break;
        if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
        for (int jFit = 0; jFit < nBinPT2D; jFit++ )
        {
            // Not considering pT < 0.4 GeV
            if ( fArrPT2D[jFit+1] <= 0.41 ) continue;
            
            // Fit
            auto fFitResult = FitModel(hdM_dpT_Tot_Rec2D[iFit][jFit],utility[iFit],utility[jFit],"STD",bSave,iFit,jFit);
            
            // Building N_Raw histogram
            auto N_Raw      = static_cast<RooRealVar*>(fFitResult->floatParsFinal().at(SignlSignl));
            
            h2D_Raw->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
            h2D_Raw->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
        }
    }
    if ( fStandar )
    {
        h2D_Raw->Scale(1./(h1D_EDT->GetBinContent(1)),"width");
        h2D_Raw->Write();
        outFileR2D->Close();
    }
    
    for ( Int_t iSys = nOptions; iSys < nOption2+nOptions; iSys++ )
    {
        TFile*  outFileFit  =   new TFile(Form("./result/Systematics_2D_%s.root",sOption2[iSys-nOptions].c_str()),"recreate");
        cout << "Starting evaluating 2D " << sOption2[iSys-nOptions].c_str() << endl;
        // 2D Fits
        for (int iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
            for (int jFit = 0; jFit < nBinPT2D; jFit++ )
            {
                // Not considering pT < 0.4 GeV
                if ( fArrPT2D[jFit+1] <= 0.41 ) continue;
                
                // Fit
                auto fResults = FitModel(hdM_dpT_Tot_Rec2D[iFit][jFit],utility[iFit],utility[jFit],Form("2D_%s",sOption2[iSys].c_str()),bSave,iFit,jFit,sOption2[iSys-nOptions].c_str());
                   
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(SignlSignl));
                h2D_Syst[iSys]->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                h2D_Syst[iSys]->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
            }
        }
        h2D_Syst[iSys]->Scale(1./(h1D_EDT->GetBinContent(1)),"width");
        h2D_Syst[iSys]->Write();
        outFileFit->Close();
    }
    
    for ( Int_t iSys = 0; iSys < nOptions; iSys++ )
    {
        if ( !f2DOptio )    break;
        TFile*  outFileFit  =   new TFile(Form("./result/Systematics_2D_1D_%s.root",sOptions[iSys].c_str()),"recreate");
        cout << "Starting evaluating 1D-2D " << sOptions[iSys].c_str() << endl;
        // Preparing 1D fit for shape in pT bins
        for (int iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
            utility[iFit]   = FitModel(hdM_dpT_Tot_Rec[iFit],Form("2D_%s",sOptions[iSys].c_str()),bSave,iFit,2,sOptions[iSys].c_str());
        }
        // 2D Fits
        for (int iFit = 0; iFit < nBinPT2D; iFit++ )
        {
            if ( fArrPT2D[iFit+1] <= 0.41 ) continue;
            for (int jFit = 0; jFit < nBinPT2D; jFit++ )
            {
                // Not considering pT < 0.4 GeV
                if ( fArrPT2D[jFit+1] <= 0.41 ) continue;

                // Fit
                auto fResults = FitModel(hdM_dpT_Tot_Rec2D[iFit][jFit],utility[iFit],utility[jFit],Form("1D-2D_%s",sOptions[iSys].c_str()),bSave,iFit,jFit,sOptions[iSys].c_str());
                
                // Building N_Raw histogram
                auto N_Raw      = static_cast<RooRealVar*>(fResults->floatParsFinal().at(SignlSignl));
                h2D_Syst[iSys]->SetBinContent      (iFit+1,jFit+1,N_Raw->getVal());
                h2D_Syst[iSys]->SetBinError        (iFit+1,jFit+1,N_Raw->getError());
            }
        }
        h2D_Syst[iSys]->Scale(1./(h1D_EDT->GetBinContent(1)),"width");
        h2D_Syst[iSys]->Write();
        outFileFit->Close();
    }
    
    // Output File for Fit Check
    TFile*  outFileFin  =   new TFile(fSystError_,"recreate");
    
    // 1D Raw
    h1D_Raw->Write();
    
    // 1D Systematics
    for ( Int_t iSys = 0; iSys < nOptions; iSys++ )
    {
        h1D_Syst[iSys]->Write();
    }
    
    // 2D Raw
    h2D_Raw->Write();
    
    // 2D Systematics
    for ( Int_t iSys = 0; iSys < nOptions+nOption2; iSys++ )
    {
        h2D_Syst[iSys]->Write();
    }
    
    insFile_DT->Close();
    outFileFin->Close();
}
