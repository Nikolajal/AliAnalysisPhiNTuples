// File for 1-Dimensional Analysis:
// !TODO: All Set!


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"
#include "RooMsgService.h"

void Syst_SignalExtraction_Analysis ( bool fSilent = false)
{
    const   Int_t   nOptions    =   15;
    const   string  sOptions[]  =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ","RK","M","W","CH3","CH5"};
    const   Int_t   nOption2    =   0;
    const   string  sOption2[]  =   {"RA","RB","RC","RD","RE","RF","RG","RH","RI","RJ","RK"};
    
    TFile **insFileHst  =   new TFile*  [2*nOptions+nOption2+2];
    
    //Recovering histograms
    TH1F**  h1D_Syst    =   new TH1F*   [nOptions];
    for ( Int_t iTer = 0; iTer < nOptions; iTer++ )
    {
        insFileHst[iTer]=   new TFile   (Form("./result/Systematics_1D_%s.root",sOptions[iTer].c_str()));
        hName           =   Form("1D_%s",sOptions[iTer].c_str());
        h1D_Syst[iTer]  =   (TH1F*)(insFileHst[iTer]->Get(hName));
    }
    TH2F**  h2D_Syst    =   new TH2F*   [nOptions+nOption2];
    for ( Int_t iTer = 0; iTer < nOptions; iTer++ )
    {
        insFileHst[iTer+nOptions]   =   new TFile   (Form("./result/Systematics_2D_1D_%s.root",sOptions[iTer].c_str()));
        hName           =   Form("1D-2D_%s",sOptions[iTer].c_str());
        h2D_Syst[iTer]  =   (TH2F*)(insFileHst[iTer+nOptions]->Get(hName));
    }
    for ( Int_t iTer = nOptions; iTer < nOption2+nOptions; iTer++ )
    {
        insFileHst[iTer+nOptions]   =   new TFile   (Form("./result/Systematics_2D_%s.root",sOption2[iTer-nOptions].c_str()));
        hName           =   Form("2D_%s",sOption2[iTer-nOptions].c_str());
        h2D_Syst[iTer]  =   (TH2F*)(insFileHst[iTer+nOptions]->Get(hName));
    }
    insFileHst[nOptions*2+nOption2]       =   new TFile   ("./result/N_Raw_1D.root");
    hName                           =   "h1D_Raw";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Raw                 =   (TH1F*)(insFileHst[nOptions*2+nOption2]->Get(hName));
    insFileHst[nOptions*2+nOption2+1]       =   new TFile   ("./result/N_Raw_2D.root");
    hName                           =   "h2D_Raw";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Raw                 =   (TH2F*)(insFileHst[nOptions*2+nOption2+1]->Get(hName));
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Creating the histograms-------------------------------------------------------------------------------
    TH1F  **fh1D_Systematics_Bin    =   new TH1F   *[nBinPT1D];
    for ( Int_t iAll = 0; iAll < nBinPT1D; iAll++ )
    {
        fh1D_Systematics_Bin[iAll]  =   new TH1F    (Form("1D_BIN_%i",iAll+1),Form("1D_BIN_%i",iAll+1),180,0.1,1.9);
        fh1D_Systematics_Bin[iAll]  ->SetTitle(Form("Percentage variation of raw yield for bin %i",iAll+1));
        fh1D_Systematics_Bin[iAll]  ->GetXaxis()->SetTitle("Percentage variation");
    }
    TH1F ***fh2D_Systematics_Bin    =   new TH1F  **[nBinPT2D];
    for ( Int_t iAll = 0; iAll < nBinPT2D; iAll++ )
    {
        fh2D_Systematics_Bin[iAll]  =   new TH1F   *[nBinPT2D];
        for ( Int_t jAll = 0; jAll < nBinPT2D; jAll++ )
        {
            fh2D_Systematics_Bin[iAll][jAll]    =   new TH1F    (Form("2D_BIN_%i_%i",iAll+1,jAll+1),Form("2D_BIN_%i_%i",iAll+1,jAll+1),180,0.1,1.9);
            fh2D_Systematics_Bin[iAll][jAll]    ->SetTitle(Form("Percentage variation of raw yield for bin %i-%i",iAll+1,jAll+1));
            fh2D_Systematics_Bin[iAll][jAll]    ->GetXaxis()->SetTitle("Percentage variation");
        }
    }
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
        
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    // Output File for Fit Check
    TFile*  outFileFit  =   new TFile("./result/Syst.root","recreate");
    
    //------ 1D Histograms ------//
    
    TH1F   *h1D_Systematics_Util    =   new TH1F("1D","1D",nBinPT1D,fArrPT1D);
    for ( Int_t iFll = 0; iFll < nOptions; iFll++ )
    {
        h1D_Systematics_Util    ->  Divide(h1D_Syst[iFll],h1D_Raw);
        for ( Int_t jFll = 0; jFll < nBinPT1D; jFll++ )
        {
            fh1D_Systematics_Bin[jFll]  ->  Fill(h1D_Systematics_Util->GetBinContent(jFll+1));
        }
        h1D_Systematics_Util  ->  SetName(Form("Frac_1D_%s",sOptions[iFll].c_str()));
        h1D_Systematics_Util  ->  Write();
    }
    
    TH1F   *Histogram_1D_Systematics_Error              =   new TH1F ("Systematic_1D",      "Systematic_1D",        nBinPT1D-4, 4,  nBinPT1D);
    TH1F   *Histogram_1D_Systematics_Error_Mean         =   new TH1F ("Systematic_1D_Mean", "Systematic_1D_Mean",   nBinPT1D-4, 4,  nBinPT1D);
    TH1F   *Histogram_1D_Systematics_Error_Widt         =   new TH1F ("Systematic_1D_Widt", "Systematic_1D_Widt",   nBinPT1D-4, 4,  nBinPT1D);
    TH1F   *Histogram_1D_Systematics_Error_Normalised   =   new TH1F ("Systematic_1D_Norm", "Systematic_1D_Norm",   nBinPT1D-4, 4,  nBinPT1D);
    for ( Int_t iTer = 4; iTer < nBinPT1D; iTer++ )
    {
        fh1D_Systematics_Bin[iTer]->GetXaxis()->SetRangeUser(0.9,1.1);
        Histogram_1D_Systematics_Error              ->  SetBinContent   (iTer-3,  fabs(1- fh1D_Systematics_Bin[iTer]->GetMean()) + fabs(fh1D_Systematics_Bin[iTer]->GetRMS()) );
        Histogram_1D_Systematics_Error_Mean         ->  SetBinContent   (iTer-3,  fabs(1- fh1D_Systematics_Bin[iTer]->GetMean()) );
        Histogram_1D_Systematics_Error_Widt         ->  SetBinContent   (iTer-3,  fabs(fh1D_Systematics_Bin[iTer]->GetRMS()) );
        Histogram_1D_Systematics_Error_Normalised   ->  SetBinContent   (iTer-3,  0);
        Histogram_1D_Systematics_Error_Normalised   ->  SetBinError     (iTer-3,  fabs(1- fh1D_Systematics_Bin[iTer]->GetMean()) + fabs(fh1D_Systematics_Bin[iTer]->GetRMS()) );
        fh1D_Systematics_Bin[iTer]->GetXaxis()->SetRangeUser(0.1,1.9);
        fh1D_Systematics_Bin[iTer]->Write();
        if ( iTer == 18 )
        {
            TCanvas *c1 = new TCanvas();
            gStyle->SetOptStat(0);
            fh1D_Systematics_Bin[iTer]->GetXaxis()->SetRangeUser(0.8,1.2);
            fh1D_Systematics_Bin[iTer]->Draw();
            c1->SaveAs("./graphs/SysBin1.pdf");
            delete c1;
        }
    }
    
    Histogram_1D_Systematics_Error              ->Scale(100.);
    Histogram_1D_Systematics_Error_Mean         ->Scale(100.);
    Histogram_1D_Systematics_Error_Widt         ->Scale(100.);
    Histogram_1D_Systematics_Error_Normalised   ->Scale(100.);
    
    TH1F   *fh1D_FnlTT_   = new   TH1F("Total1D",           "Total1D",          100,0.,.2);
    TH1F   *fh1D_FnlTTM   = new   TH1F("Total1D_mean",      "Total1D_mean",     100,0.,.2);
    TH1F   *fh1D_FnlTTW   = new   TH1F("Total1D_widt",      "Total1D_widt",     100,0.,.2);
    for ( Int_t iTer = 4; iTer < nBinPT1D; iTer++ )
    {
        fh1D_FnlTT_->Fill(Histogram_1D_Systematics_Error->GetBinContent(iTer-3));
        fh1D_FnlTTM->Fill(Histogram_1D_Systematics_Error_Mean->GetBinContent(iTer-3));
        fh1D_FnlTTW->Fill(Histogram_1D_Systematics_Error_Widt->GetBinContent(iTer-3));
    }
    
    //------ 2D Histogram of N_raw ------//
    
    TH2F   *h2D_Systematics_Util    =   new TH2F("2D","2D",nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    for ( Int_t iFll = 0; iFll < nOptions; iFll++ )
    {
        h2D_Systematics_Util    ->  Divide(h2D_Syst[iFll],h2D_Raw);
        for ( Int_t jFll = 0; jFll < nBinPT2D; jFll++ )
        {
            for ( Int_t kFll = 0; kFll < nBinPT2D; kFll++ )
            {
                if ( jFll == 0 && kFll == 9 ) continue;
                if ( jFll == 9 && kFll == 0  ) continue;
                fh2D_Systematics_Bin[jFll][kFll]    ->  Fill(h2D_Systematics_Util->GetBinContent(jFll+1,kFll+1));
            }
        }
        h2D_Systematics_Util  ->  SetName(Form("Frac_1D_2D_%s",sOptions[iFll].c_str()));
        h2D_Systematics_Util  ->  Write();
    }
    for ( Int_t iFll = nOptions; iFll < nOptions+nOption2; iFll++ )
    {
        h2D_Systematics_Util    ->  Divide(h2D_Syst[iFll],h2D_Raw);
        for ( Int_t jFll = 0; jFll < nBinPT2D; jFll++ )
        {
            for ( Int_t kFll = 0; kFll < nBinPT2D; kFll++ )
            {
                fh2D_Systematics_Bin[jFll][kFll]    ->  Fill(h2D_Systematics_Util->GetBinContent(jFll,kFll));
            }
        }
        h2D_Systematics_Util  ->  SetName(Form("Frac_2D_%s",sOption2[iFll-nOptions].c_str()));
        h2D_Systematics_Util  ->  Write();
    }
    
    TH2F   *fh2D_FnlBin  = new   TH2F ("Final_Sys2",        "Final_Sys2",       nBinPT2D-2,2,nBinPT2D,nBinPT2D-2,2,nBinPT2D);
    TH2F   *fh2D_FnlBiM  = new   TH2F ("Final_Syst_Mea2",   "Final_Syst_Mea2",  nBinPT2D-2,2,nBinPT2D,nBinPT2D-2,2,nBinPT2D);
    TH2F   *fh2D_FnlBiW  = new   TH2F ("Final_Syst_Wid2",   "Final_Syst_Wid2",  nBinPT2D-2,2,nBinPT2D,nBinPT2D-2,2,nBinPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        for ( Int_t jTer = 2; jTer < nBinPT2D; jTer++ )
        {
            fh2D_Systematics_Bin[iTer][jTer]->GetXaxis()->SetRangeUser(0.4,1.6);
            fh2D_FnlBiM     ->  SetBinContent   (iTer-1,  jTer-1,   fabs(1-fh2D_Systematics_Bin[iTer][jTer]->GetMean()) );
            fh2D_FnlBiW     ->  SetBinContent   (iTer-1,  jTer-1,   fabs(fh2D_Systematics_Bin[iTer][jTer]->GetRMS()) );
            fh2D_FnlBin     ->  SetBinContent   (iTer-1,  jTer-1,   fabs(1-fh2D_Systematics_Bin[iTer][jTer]->GetMean())+ fabs(fh2D_Systematics_Bin[iTer][jTer]->GetRMS()));
            fh2D_Systematics_Bin[iTer][jTer]->GetXaxis()->SetRangeUser(0.1,1.9);
            fh2D_Systematics_Bin[iTer][jTer]->Write();
            if ( iTer == 6 && jTer == 4 )
            {
                TCanvas *c1 = new TCanvas();
                gStyle->SetOptStat(0);
                fh2D_Systematics_Bin[iTer][jTer]->GetXaxis()->SetRangeUser(0.4,1.6);
                fh2D_Systematics_Bin[iTer][jTer]->Draw();
                c1->SaveAs("./graphs/SysBin2.pdf");
                delete c1;
            }
        }
    }
    TH1F   *fh2D_FnlTT_   = new   TH1F("Total2D",           "Total2D",              100,0.,1.);
    TH1F   *fh2D_FnlTTM   = new   TH1F("Total2D_mean",      "Total2D_mean",            100,0.,1.);
    TH1F   *fh2D_FnlTTW   = new   TH1F("Total2D_widt",      "Total2D_widt",            100,0.,1.);
    TH1F   *fh2D_FnlTT1   = new   TH1F("Total2D_1D",        "Total2D_1D",              nBinPT2D,0,nBinPT2D);
    TH1F   *fh2D_FnlTM1   = new   TH1F("Total2D_mean_1D",   "Total2D_mean_1D",    nBinPT2D,0,nBinPT2D);
    TH1F   *fh2D_FnlTW1   = new   TH1F("Total2D_widt_1D",   "Total2D_widt_1D",    nBinPT2D,0,nBinPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        for ( Int_t jTer = 2; jTer < nBinPT2D; jTer++ )
        {
            fh2D_FnlTT_->Fill(fh2D_FnlBin->GetBinContent(iTer-1,jTer-1));
            fh2D_FnlTTM->Fill(fh2D_FnlBiM->GetBinContent(iTer-1,jTer-1));
            fh2D_FnlTTW->Fill(fh2D_FnlBiW->GetBinContent(iTer-1,jTer-1));
        }
    }
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        fh2D_FnlTT1->Fill(iTer,fabs(fh2D_FnlBin->GetBinContent(iTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        fh2D_FnlTM1->Fill(iTer,fabs(fh2D_FnlBiM->GetBinContent(iTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        fh2D_FnlTW1->Fill(iTer,fabs(fh2D_FnlBiW->GetBinContent(iTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        for ( Int_t jTer = iTer+1; jTer < nBinPT2D; jTer++ )
        {
            fh2D_FnlTT1->Fill(iTer,fabs(fh2D_FnlBin->GetBinContent(iTer-1,jTer-1))/(2*(nBinPT2D-iTer)-1));
            fh2D_FnlTM1->Fill(iTer,fabs(fh2D_FnlBiM->GetBinContent(iTer-1,jTer-1))/(2*(nBinPT2D-iTer)-1));
            fh2D_FnlTW1->Fill(iTer,fabs(fh2D_FnlBiW->GetBinContent(iTer-1,jTer-1))/(2*(nBinPT2D-iTer)-1));
            fh2D_FnlTT1->Fill(iTer,fabs(fh2D_FnlBin->GetBinContent(jTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
            fh2D_FnlTM1->Fill(iTer,fabs(fh2D_FnlBiM->GetBinContent(jTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
            fh2D_FnlTW1->Fill(iTer,fabs(fh2D_FnlBiW->GetBinContent(jTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        }
    }
    
    TH2F   *Histogram_2D_Systematics_Error              =   new TH2F ("Systematic_2D",      "Systematic_2D",        nBinPT2D-2, 2,  nBinPT2D,   nBinPT2D-2, 2,  nBinPT2D);
    TH2F   *Histogram_2D_Systematics_Error_Normalised   =   new TH2F ("Systematic_2D_Norm", "Systematic_2D_Norm",   nBinPT2D-2, 2,  nBinPT2D,   nBinPT2D-2, 2,  nBinPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        for ( Int_t jTer = 2; jTer < nBinPT2D; jTer++ )
        {
            fh2D_Systematics_Bin[iTer][jTer]->GetXaxis()->SetRangeUser(0.6,1.4);
            Histogram_2D_Systematics_Error              ->  SetBinContent   (iTer-1,  jTer-1,   fabs(1-fh2D_Systematics_Bin[iTer][jTer]->GetMean())+ fabs(fh2D_Systematics_Bin[iTer][jTer]->GetRMS()));
            Histogram_2D_Systematics_Error_Normalised   ->  SetBinContent   (iTer-1,  jTer-1,   0.);
            Histogram_2D_Systematics_Error_Normalised   ->  SetBinError     (iTer-1,  jTer-1,   fabs(1-fh2D_Systematics_Bin[iTer][jTer]->GetMean())+ fabs(fh2D_Systematics_Bin[iTer][jTer]->GetRMS()));
            fh2D_Systematics_Bin[iTer][jTer]->GetXaxis()->SetRangeUser(0.1,1.9);
        }
    }
    
    Histogram_2D_Systematics_Error              ->Scale(100.);
    Histogram_2D_Systematics_Error_Normalised   ->Scale(100.);
    
    TCanvas *c2 = new TCanvas();
    Histogram_2D_Systematics_Error->Draw("colz text");
    c2->Write();
    c2->SaveAs("./graphs/Sys2D.pdf");
    delete c2;
    
    TH1F   *Histogram_2D_Systematics_Error_Projection               = new   TH1F("Systematics_2D_1D",       "Systematics_2D_1D",        nBinPT2D-2, 2,  nBinPT2D);
    TH1F   *Histogram_2D_Systematics_Error_Projection_Normalised    = new   TH1F("Systematics_2D_1D_Norm",  "Systematics_2D_1D_Norm",   nBinPT2D-2, 2,  nBinPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        for ( Int_t jTer = iTer+1; jTer < nBinPT2D; jTer++ )
        {
            Histogram_2D_Systematics_Error_Projection->Fill(iTer,fabs(Histogram_2D_Systematics_Error->GetBinContent(iTer-1,jTer-1))/(2*(nBinPT2D-iTer)-1));
            Histogram_2D_Systematics_Error_Projection->Fill(iTer,fabs(Histogram_2D_Systematics_Error->GetBinContent(jTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        }
        if ( iTer == 2  ) continue;
        Histogram_2D_Systematics_Error_Projection->Fill(iTer,fabs(Histogram_2D_Systematics_Error->GetBinContent(iTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
    }
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )
    {
        Histogram_2D_Systematics_Error_Projection_Normalised->SetBinError(iTer,Histogram_2D_Systematics_Error_Projection->GetBinContent(iTer));
        Histogram_2D_Systematics_Error_Projection_Normalised->SetBinContent(iTer,0.);
    }
    
    // Statistical Error
    
    TH1F   *Histogram_1D_Statistical_Error              =   new TH1F("Statistical_1D",      "Statistical_1D",       nBinPT1D-4, 4,  nBinPT1D);
    TH1F   *Histogram_1D_Statistical_Error_Normalised   =   new TH1F("Statistical_1D_Norm", "Statistical_1D_Norm",  nBinPT1D-4, 4,  nBinPT1D);
    for ( Int_t iTer = 4; iTer < nBinPT1D; iTer++ )
    {
        auto fVal = h1D_Raw                         ->GetBinContent (iTer+1);
        auto fErr = h1D_Raw                         ->GetBinError   (iTer+1);
        Histogram_1D_Statistical_Error              ->SetBinContent (iTer-3,fErr/fVal);
        Histogram_1D_Statistical_Error_Normalised   ->SetBinContent (iTer-3,0.);
        Histogram_1D_Statistical_Error_Normalised   ->SetBinError   (iTer-3,fErr/fVal);
    }
    Histogram_1D_Statistical_Error                      ->Scale(100.);
    Histogram_1D_Statistical_Error_Normalised           ->Scale(100.);
    
    
    TH2F   *Histogram_2D_Statistical_Error              =   new TH2F("Statistical_2D",      "Statistical_2D",       nBinPT2D-2, 2,  nBinPT2D,   nBinPT2D-2, 2,  nBinPT2D);
    TH2F   *Histogram_2D_Statistical_Error_Normalised   =   new TH2F("Statistical_2D_Norm", "Statistical_2D_Norm",  nBinPT2D-2, 2,  nBinPT2D,   nBinPT2D-2, 2,  nBinPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        for ( Int_t jTer = 2; jTer < nBinPT2D; jTer++ )
        {
            auto fVal = h2D_Raw                         ->GetBinContent(iTer+1,jTer+1);
            auto fErr = h2D_Raw                         ->GetBinError(iTer+1,jTer+1);
            Histogram_2D_Statistical_Error              ->SetBinContent(iTer-1,jTer-1,fErr/fVal);
            Histogram_2D_Statistical_Error_Normalised   ->SetBinContent(iTer-1,jTer-1,0.);
            Histogram_2D_Statistical_Error_Normalised   ->SetBinError(iTer-1,jTer-1,fErr/fVal);
        }
    }
    Histogram_2D_Statistical_Error                      ->Scale(100.);
    Histogram_2D_Statistical_Error_Normalised           ->Scale(100.);
    
    TCanvas *c1 = new TCanvas();
    Histogram_2D_Statistical_Error->Draw("colz text");
    c1->Write();
    c1->SaveAs("./graphs/Stat2D.pdf");
    delete c1;
    
    TH1F   *Histogram_2D_Statistical_Error_Projection               = new   TH1F("Statistical_2D_1D",       "Statistical_2D_1D",        nBinPT2D-2, 2,  nBinPT2D);
    TH1F   *Histogram_2D_Statistical_Error_Projection_Normalised    = new   TH1F("Statistical_2D_1D_Norm",  "Statistical_2D_1D_Norm",   nBinPT2D-2, 2,  nBinPT2D);
    for ( Int_t iTer = 2; iTer < nBinPT2D; iTer++ )
    {
        Histogram_2D_Statistical_Error_Projection->Fill(iTer,fabs(Histogram_2D_Statistical_Error->GetBinContent(iTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        for ( Int_t jTer = iTer+1; jTer < nBinPT2D; jTer++ )
        {
            Histogram_2D_Statistical_Error_Projection->Fill(iTer,fabs(Histogram_2D_Statistical_Error->GetBinContent(iTer-1,jTer-1))/(2*(nBinPT2D-iTer)-1));
            Histogram_2D_Statistical_Error_Projection->Fill(iTer,fabs(Histogram_2D_Statistical_Error->GetBinContent(jTer-1,iTer-1))/(2*(nBinPT2D-iTer)-1));
        }
    }
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )
    {
        Histogram_2D_Statistical_Error_Projection_Normalised->SetBinError(iTer,Histogram_2D_Statistical_Error_Projection->GetBinContent(iTer));
        Histogram_2D_Statistical_Error_Projection_Normalised->SetBinContent(iTer,0.);
    }
    
    
    // Output File for Fit Check
    TFile*  outFileFi2  =   new TFile("./result/Syst2.root","recreate");
    
    TH2F   *Histogram_2D_Statistical_Systematic_Overlap =   new TH2F    ("Syst_Stat",       "Syst_Stat",    nBinPT2D-2, 2,  nBinPT2D,   nBinPT2D-2, 2,  nBinPT2D);
    
    Histogram_2D_Statistical_Systematic_Overlap         ->Add(Histogram_2D_Statistical_Error,Histogram_2D_Systematics_Error,1.,-1.);
    
    TCanvas    *Canvas_Satistical_Systematic_Overlap    =   new TCanvas();
    TLegend    *Legend_Satistical_Systematic_Overlap    =   new TLegend(0.25,0.75,0.55,0.9);
    
    gStyle->SetOptStat(0);
    
    Histogram_1D_Statistical_Error_Normalised   ->SetFillColorAlpha(kGray,0.5);
    Histogram_1D_Statistical_Error_Normalised   ->Draw("E3 same");
    Histogram_1D_Statistical_Error_Normalised   ->GetXaxis()->SetTitle("p_{T} bin");
    Histogram_1D_Statistical_Error_Normalised   ->GetYaxis()->SetTitle("Error (%)");
    Histogram_1D_Statistical_Error_Normalised   ->SetTitle("");
    Histogram_1D_Systematics_Error_Normalised   ->SetFillColorAlpha(kRed,0.5);
    Histogram_1D_Systematics_Error_Normalised   ->Draw("E3 same");
    Legend_Satistical_Systematic_Overlap        ->AddEntry(Histogram_1D_Systematics_Error_Normalised,   "Systematics",  "F");
    Legend_Satistical_Systematic_Overlap        ->AddEntry(Histogram_1D_Statistical_Error_Normalised,   "Statistical",  "F");
    Legend_Satistical_Systematic_Overlap        ->Draw("same");
    
    Canvas_Satistical_Systematic_Overlap        ->Write();
    Canvas_Satistical_Systematic_Overlap        ->SaveAs("./graphs/Canvas_Satistical_Systematic_Overlap.pdf");
    Canvas_Satistical_Systematic_Overlap        ->SaveAs("./graphs/Canvas_Satistical_Systematic_Overlap.png");
    
    Legend_Satistical_Systematic_Overlap        ->Clear();
    delete Canvas_Satistical_Systematic_Overlap;
    
    TCanvas    *Canvas_Satistical_Systematic_Overla2    =   new TCanvas();
    Canvas_Satistical_Systematic_Overla2->Divide(2,5);
    TH1D  **hProjection =   new TH1D   *[nBinPT2D];
    TH1D  **hProjectio2 =   new TH1D   *[nBinPT2D];
    for ( Int_t iTer = 0; iTer < nBinPT2D-2; iTer++ )
    {
        Canvas_Satistical_Systematic_Overla2->cd(iTer+1);
        hProjection[iTer]   =   new TH1D(*Histogram_2D_Statistical_Error_Normalised->ProjectionX(Form("%i_",iTer),iTer+1,iTer+1));
        hProjectio2[iTer]   =   new TH1D(*Histogram_2D_Systematics_Error_Normalised->ProjectionX(Form("%i_",iTer),iTer+1,iTer+1));
        hProjection[iTer]   ->SetFillColorAlpha(kGray,0.5);
        if ( iTer == 0 )    hProjection[iTer]->GetXaxis()->SetRangeUser(2,11);
        if ( iTer == 9 )    hProjection[iTer]->GetXaxis()->SetRangeUser(3,12);
        hProjection[iTer]   ->Draw("E3 same");
        hProjection[iTer]   ->GetXaxis()->SetTitle("p_{T} bin");
        hProjection[iTer]   ->GetYaxis()->SetTitle("Error (%)");
        hProjection[iTer]   ->SetTitle("");
        hProjectio2[iTer]   ->SetFillColorAlpha(kRed,0.5);
        hProjectio2[iTer]   ->Draw("E3 same");
        if (iTer == 0)
        {
            Legend_Satistical_Systematic_Overlap        ->AddEntry(hProjectio2[iTer],   "Systematics",  "F");
            Legend_Satistical_Systematic_Overlap        ->AddEntry(hProjection[iTer],   "Statistical",  "F");
        }
        Legend_Satistical_Systematic_Overlap        ->Draw("same");
        TCanvas * c22 = new TCanvas ();
        hProjection[iTer]   ->Draw("E3 same");
        hProjectio2[iTer]   ->Draw("E3 same");
        Legend_Satistical_Systematic_Overlap        ->Draw("same");
        c22 ->Write();
        c22 ->SaveAs(Form("./graphs/Canvas_Satistical_Systematic_Overlap_2D_%i.pdf",iTer+3));
        c22 ->SaveAs(Form("./graphs/Canvas_Satistical_Systematic_Overlap_2D_%i.png",iTer+3));
        delete c22;
    }
    
    Canvas_Satistical_Systematic_Overla2        ->Write();
    Canvas_Satistical_Systematic_Overla2        ->SaveAs("./graphs/Canvas_Satistical_Systematic_Overlap_2D.pdf");
    Canvas_Satistical_Systematic_Overla2        ->SaveAs("./graphs/Canvas_Satistical_Systematic_Overlap_2D.png");
    
    delete Canvas_Satistical_Systematic_Overla2;
    
    
    Histogram_1D_Statistical_Error                          ->Write();
    Histogram_1D_Statistical_Error_Normalised               ->Write();
    Histogram_2D_Statistical_Error                          ->Write();
    Histogram_2D_Statistical_Error_Normalised               ->Write();
    Histogram_2D_Statistical_Error_Projection               ->Write();
    Histogram_2D_Statistical_Error_Projection_Normalised    ->Write();
    Histogram_1D_Systematics_Error                          ->Write();
    Histogram_1D_Systematics_Error_Normalised               ->Write();
    Histogram_2D_Systematics_Error                          ->Write();
    Histogram_2D_Systematics_Error_Normalised               ->Write();
    Histogram_2D_Systematics_Error_Projection               ->Write();
    Histogram_2D_Systematics_Error_Projection_Normalised    ->Write();
    Histogram_1D_Systematics_Error_Mean                     ->Write();
    Histogram_1D_Systematics_Error_Widt                     ->Write();
    Histogram_2D_Statistical_Systematic_Overlap             ->Write();
    
    fh1D_FnlTT_  -> Write();
    fh2D_FnlTT_  -> Write();
    fh1D_FnlTTM  -> Write();
    fh2D_FnlTTM  -> Write();
    fh1D_FnlTTW  -> Write();
    fh2D_FnlTTW  -> Write();
    fh2D_FnlTT1  -> Write();
    fh2D_FnlTM1  -> Write();
    fh2D_FnlTW1  -> Write();
    
    // Closing Files
    for ( Int_t iTer = 0; iTer < nOption2+nOptions*2+2; iTer++ )
    {
        insFileHst[iTer]->Close();
    }
    /*
    TH1F *  h1D_RXT     = new TH1F (*h1D_Rec);
    TH1F *  h1D_R0T     = new TH1F (*h1D_Raw);
    TH1F *  h1D_R1T     = new TH1F (*h1D_Ra1);
    TH1F *  h1D_R2T     = new TH1F (*h1D_Ra2);
    TH1F *  h1D_R3T     = new TH1F (*h1D_Ra3);
    TH1F *  h1D_R4T     = new TH1F (*h1D_Ra4);
    TH1F *  h1D_R5T     = new TH1F (*h1D_Ra5);
    TH1F *  h1D_R6T     = new TH1F (*h1D_Ra6);
    
    h1D_RXT->Divide(h1D_RXT,h1D_Raw);
    h1D_R0T->Divide(h1D_R0T,h1D_Raw);
    h1D_R1T->Divide(h1D_R1T,h1D_Raw);
    h1D_R2T->Divide(h1D_R2T,h1D_Raw);
    h1D_R3T->Divide(h1D_R3T,h1D_Raw);
    h1D_R4T->Divide(h1D_R4T,h1D_Raw);
    h1D_R5T->Divide(h1D_R5T,h1D_Raw);
    h1D_R6T->Divide(h1D_R6T,h1D_Raw);
    
    h1D_RXT->SetNameTitle("True",               "True");
    h1D_R0T->SetNameTitle("Fit",                "Fit");
    h1D_R1T->SetNameTitle("FMass",              "FMass");
    h1D_R2T->SetNameTitle("FWidth",             "FWidth");
    h1D_R3T->SetNameTitle("CH3",                "CH3");
    h1D_R4T->SetNameTitle("FMassWidth",         "FMassWidth");
    h1D_R5T->SetNameTitle("FMassWidth_CH3",     "FMassWidth_CH3");
    h1D_R6T->SetNameTitle("CH2",                "CH2");
    
    for (Int_t iPoint = 0; iPoint < nBinPT1D; iPoint++ )
    {
        if ( fArrPT1D[iPoint+1] <= 0.41 ) continue;
        h1D_R0T->SetBinError(iPoint+1,h1D_Raw->GetBinError(iPoint+1)/h1D_Raw->GetBinContent(iPoint+1));
        h1D_R1T->SetBinError(iPoint+1,0);
        h1D_R2T->SetBinError(iPoint+1,0);
        h1D_R3T->SetBinError(iPoint+1,0);
        h1D_R4T->SetBinError(iPoint+1,0);
        h1D_R5T->SetBinError(iPoint+1,0);
        h1D_R6T->SetBinError(iPoint+1,0);
    }
    
    TLegend * fAllLegend    =   new TLegend();
    //fAllLegend->AddEntry(h1D_RXT,"True MC","P");
    fAllLegend->AddEntry(h1D_R0T,"Baseline Fit","P");
    fAllLegend->AddEntry(h1D_R1T,"Chebychev 3","P");
    fAllLegend->AddEntry(h1D_R2T,"Chebychev 5","P");
    fAllLegend->AddEntry(h1D_R3T,"Fixed Mass","P");
    fAllLegend->AddEntry(h1D_R4T,"Fixed Width","P");
    fAllLegend->AddEntry(h1D_R5T,"Range1 Fit","P");
    fAllLegend->AddEntry(h1D_R6T,"Range 2 Fit","P");
    
    
    TCanvas * fAllSyst = new TCanvas();
    gStyle->SetOptStat(0);
    h1D_R0T->GetXaxis()->SetRangeUser(0.4,10.);
    h1D_R0T->GetYaxis()->SetRangeUser(0.9,1.1);
    h1D_R0T->SetFillColor(kGray);
    h1D_R0T->SetMarkerColor(kBlack);
    h1D_R0T->SetMarkerStyle(7);
    h1D_R0T->Draw("E3");
    h1D_R0T->Draw("same hist P ");
    h1D_R1T->SetMarkerColor(kBlue);
    h1D_R1T->SetMarkerStyle(24);
    h1D_R1T->Draw("same hist P ");
    h1D_R2T->SetMarkerColor(kGreen);
    h1D_R2T->SetMarkerStyle(25);
    h1D_R2T->Draw("same hist P ");
    h1D_R3T->SetMarkerColor(kRed);
    h1D_R3T->SetMarkerStyle(26);
    h1D_R3T->Draw("same hist P ");
    h1D_R4T->SetMarkerColor(46);
    h1D_R4T->SetMarkerStyle(30);
    h1D_R4T->Draw("same hist P ");
    h1D_R5T->SetMarkerColor(38);
    h1D_R5T->SetMarkerStyle(32);
    h1D_R5T->Draw("same hist P ");
    h1D_R6T->SetMarkerColor(4);
    h1D_R6T->SetMarkerStyle(40);
    h1D_R6T->Draw("same hist P ");
    h1D_RXT->SetMarkerColor(kBlack);
    h1D_RXT->SetMarkerStyle(4);
    //h1D_RXT->Draw("same hist P ");
    fAllLegend->Draw("same");
    fAllLegend->Write();
    fAllSyst->Write();
    
    h1D_Raw->Write();
    h1D_Rec->Write();
    h1D_Ra1->Write();
    h1D_Ra2->Write();
    h1D_Ra3->Write();
    h1D_Ra4->Write();
    h1D_Ra5->Write();
    h1D_Ra6->Write();
    h1D_RXT->Write();
    h1D_R0T->Write();
    h1D_R1T->Write();
    h1D_R2T->Write();
    h1D_R3T->Write();
    h1D_R4T->Write();
    h1D_R5T->Write();
    h1D_R6T->Write();
    h2D_Raw->Write();
     */
    
    outFileFit->Close();
    outFileFi2->Close();
}
