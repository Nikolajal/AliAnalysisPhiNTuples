// File for 1-Dimensional Analysis:
// !TODO: All Set!


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"
#include "RooMsgService.h"

void CorrectionsPTFit ( bool fSilent = false )
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
    
    // Opening Data and Efficiencies File
    TFile*  insFile_DT              =   new TFile   (fFitResults);
    TFile*  insFile_EF              =   new TFile   (fEfficiHist);
    TFile*  insCheck_               =   new TFile   ("./result/HEPData-ins1182213-v1-root.root");
    TFile*  insFile_MC_P6           =   new TFile   ("./result/MCTruth_Pythia6.root");
    TFile*  insFile_MC_P8           =   new TFile   ("./result/MCTruth_Pythia8.root");
    
    // Fit Variables for Roofit
    RooRealVar  xTransverseMom1D("xTransverseMom1D","xTransverseMom1D", fMinPT1D,fMaxPT1D);
    RooRealVar  xTransverseMom2D("xTransverseMom2D","xTransverseMom2D", fMinPT2D,fMaxPT2D);
    RooRealVar  yTransverseMom2D("yTransverseMom2D","yTransverseMom2D", fMinPT2D,fMaxPT2D);
    
    //-// Recovering Data histograms
    
    //------ 1D Histograms Recovery ------//
    
    hName                           =   "h1D_Raw";                              // Name of the histogram in the preprocessed file
    TH1F *  h1D_Raw                 =   (TH1F*)(insFile_DT->Get(hName));
    hName                           =   "hNP_1D_Eff_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Eff                 =   (TH1F*)(insFile_EF->Get(hName));
    hName                           =   "Entry_MC_DT";                          // Name of the histogram in the preprocessed file
    TH1F *  h1D_Eve                 =   (TH1F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Tru_P6              =   (TH1F*)(insFile_MC_P6->Get(hName));
    hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Tru_P8              =   (TH1F*)(insFile_MC_P8->Get(hName));
    
    //------ 2D Histograms Recovery ------//
    
    hName                           =   "h2D_Raw";                              // Name of the histogram in the preprocessed file
    TH2F *  h2D_Raw                 =   (TH2F*)(insFile_DT->Get(hName));
    hName                           =   "hNP_2D_Eff_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Eff                 =   (TH2F*)(insFile_EF->Get(hName));
    
    //------ TGraphAsymmErrors Rec ------//
    hName                           =   "Table 2/Graph1D_y1";                   // Name of the histogram in the preprocessed file
    TGraphAsymmErrors *  gCheck_    =   (TGraphAsymmErrors*)(insCheck_->Get(hName));
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Creating the histograms-------------------------------------------------------------------------------
    
    //--- Generating the binning array ---//
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    //---------- 1D Histograms ----------//
   
    hName   = "h1D_Res";
    hTitle  = "Number of #phi in |y|<5";
    TH1F*   h1D_Res     =   new TH1F (hName,hTitle,nBinPT1D,fArrPT1D);
    TH1_SetDescription (h1D_Res);
    
    //---------- 2D Histograms ----------//
    
    hName   = "h2D_Res";
    hTitle  = "Number of #phi in |y|<5";
    TH2F*   h2D_Res     =   new TH2F (hName,hTitle,nBinPT2D,fArrPT2D,nBinPT2D,fArrPT2D);
    TH2_SetDescription (h2D_Res);
    
    
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    // Output File for Fit Check
    TFile*  outFile_FT  =   new TFile(fAnlResHist,"recreate");
    
    //------- Histograms Scaling  -------//
    
    //         N_raw    1     1    Eff    1     1     1     1     1
    // N_res = ----- X --- X --- X --- X --- X --- X --- X --- X ---
    //          eff    DpT   Dy    Evn   BR    Vrt   Sig   Trk   PID
    
    // Error propagation
    h1D_Res->Sumw2();
    h2D_Res->Sumw2();
    
    // Producing of N_res by applying corrections to N_raw
    h1D_Res->Divide(h1D_Raw,h1D_Eff,1.,1.,"");
    h2D_Res->Divide(h2D_Raw,h2D_Eff,1.,1.,"");
    
    // Scaling in pT
    h1D_Res->Scale(1.,"width");
    h2D_Res->Scale(1.,"width");
    
    // Scaling in Rapidity Interval
    h1D_Res->Scale(1./kRapidityInterval);
    h2D_Res->Scale(1./kRapidityInterval);
    
    // Scaling for events corrected for efficiency
    h1D_Res->Scale(kEventEfficiency_/(h1D_Eve->GetBinContent(2)));
    h2D_Res->Scale(kEventEfficiency_/(h1D_Eve->GetBinContent(2)));
    
    // Scaling for Branching Ratio
    h1D_Res->Scale(1./kBranchingRatio__);
    h2D_Res->Scale(1./(kBranchingRatio__*kBranchingRatio__));
    
    // Scaling for Vertex Selection Efficiency
    h1D_Res->Scale(1./kVertexEfficiency);
    h2D_Res->Scale(1./kVertexEfficiency);
    
    // Scaling for Signal Extraction Efficiency
    h1D_Res->Scale(1./kSignalExtraction);
    h2D_Res->Scale(1./kSignalExtraction);

    // Scaling for Tracking Efficiency
    h1D_Res->Scale(1./kSignalExtraction);
    h2D_Res->Scale(1./kSignalExtraction);
    
    // Scaling for PID Efficiency
    h1D_Res->Scale(1./kParticleIdentifi);
    h2D_Res->Scale(1./kParticleIdentifi);
    
    if ( bPythiaTest )
    {
        // Scaling for MC biases
        h1D_Res->Scale(1./(kEventEfficiency_*kPythia1DEfficien));
        h2D_Res->Scale(1./(kEventEfficiency_*kPythia2DEfficien));
    }
    
    //------- Histograms Fitting  -------//
    Float_t fInt1D, E1Dhig, E1Dlow;
    Float_t** fInt2D;   Float_t** E2Dhig;   Float_t** E2Dlow;
    
    TGraphAsymmErrors * g1D_Res = SetSystErrors(h1D_Res);
    
    SetLevyTsalPT1D();
    //fInt1D = EvlIntegral(g1D_Res,E1Dhig,E1Dlow);
    fInt1D = EvlMeanPT__(g1D_Res,E1Dhig,E1Dlow);
    cout << " PT:" << fInt1D << " +-" << E1Dhig << endl;
    fInt1D = EvlMeanPT__(gCheck_,E1Dhig,E1Dlow);
    cout << " PT:" << fInt1D << " +-" << E1Dhig << endl;
    //fInt2D = EvlIntegral(h2D_Res,E2Dhig,E2Dlow);
    
    //------- Graphics of Fitting -------//
    TGraphAsymmErrors * gMeanPT = new TGraphAsymmErrors();
    
    //---- Graphics of Presentation -----//
    TGrCompare1D(g1D_Res,gCheck_,h1D_Tru_P6,h1D_Tru_P8);
    
    //---------------------//
    // Output and wrap up  //-------------------------------------------------------------------------------
    //---------------------//
    
    // Output File for Results
    TFile*  outFile_RS  =   new TFile(fAnlResults,"recreate");
    
    // Writing Results to file
    h1D_Raw->Write();
    h2D_Raw->Write();
    h1D_Res->Write();
    h2D_Res->Write();
    
    // Closing opened files
    insFile_DT->Close();
    insFile_EF->Close();
    outFile_FT->Close();
    outFile_RS->Close();
}
