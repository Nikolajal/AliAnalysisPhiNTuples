#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

void PresentationGraphs ( const char* const DTFile = fFileMCForm, const char* const MCFile = fFileMCForm )
{
    // Retrieving Analysis 1D and 2D files
    TFile*  insFile_PP  =   new TFile   (fInvMasHist);
    TFile*  insFile_EF  =   new TFile   (fEfficiHist);
    TFile*  insFile_1D  =   new TFile   (fFitResHist);
    TFile*  insFile_FR  =   new TFile   (fFitResults);
    TFile*  insFiZZ_1D  =   new TFile   (fFileZZ1D1D);
    TFile*  insFiZZ_2D  =   new TFile   (fFileZZ2D1D);
    TFile*  insFile_PT  =   new TFile   (fFileZZ1DPT);
    //TFile*  insFile_MC  =   new TFile   (MCFile);
    //TFile*  insFile_DT;
    //if ( bPythiaTest )  insFile_DT  =   insFile_MC;
    //else                insFile_DT  =   new TFile   (DTFile);
    TFile*  insFile_DT  =   new TFile   ("./result/LHC10.root");
    TFile*  insFile_MC  =   new TFile   ("./result/LHC14j4.root");
    
    // Retrieving Analysis 1D and 2D files from PYTHIA
    TFile*  insMCle_1D  =   new TFile   ("./result_Pythia/DTAnalysis1D.root");
    TFile*  insMCle_2D  =   new TFile   ("./result_Pythia/DTAnalysis2D.root");
    //TFile*  insMCle_PT  =   new TFile   ("./MC_PYTHIA_result/DTAnalysis1D.root");
    
    //Montecarlo Reference
    Int_t nMCcompare = 2;
    TFile*  insFile_MC_P6  =   new TFile   ("./result/MCTruth_Pythia6.root");
    TFile*  insFile_MC_P8  =   new TFile   ("./result/MCTruth_Pythia8.root");
    TFile*  insFile_M2_P6  =   new TFile   ("./result/MCTruth_Profile_Pythia6.root");
    TFile*  insFile_M2_P8  =   new TFile   ("./result/MCTruth_Profile_Pythia8.root");
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    //Recovering histograms
    TH2F *** hModel2D           = new TH2F **   [nBinPT2D];
    TH2F *** hData_2D           = new TH2F **   [nBinPT2D];
    TH1F **  hModel1D           = new TH1F *    [nBinPT1D];
    TH1F **  hData_1D           = new TH1F *    [nBinPT1D];
    TH1D **  hSlcPTFitX_2DMC    = new TH1D *    [nBinPT2D];
    TH1D **  hSlcPTFitY_2DMC    = new TH1D *    [nBinPT2D];
    TH1D **  hSlcPTFitX_2DDT    = new TH1D *    [nBinPT2D];
    TH1D **  hSlcPTFitY_2DDT    = new TH1D *    [nBinPT2D];
    TF1  **  hSlcPTFncX_2DMC    = new TF1  *    [nBinPT2D];
    TF1  **  hSlcPTFncY_2DMC    = new TF1  *    [nBinPT2D];
    TF1  **  hSlcPTFncX_2DDT    = new TF1  *    [nBinPT2D];
    TF1  **  hSlcPTFncY_2DDT    = new TF1  *    [nBinPT2D];
    TF1  *   fSlcPTFncX_DT;
    TF1  *   fSlcPTFncY_DT;
    TF1  *   fSlcPTFncX_MC;
    TF1  *   fSlcPTFncY_MC;
    TH1F *   hSlcPTFncX_DT;
    TH1F *   hSlcPTFncY_DT;
    TH1F *   hSlcPTFncX_MC;
    TH1F *   hSlcPTFncY_MC;
    TH1F *   hTOFPID;
    TH1F *   hTPCPID;
    TH1F *   hTOFPID_1;
    TH1F *   hTPCPID_1;
    TH1F *   hTOFPID_2;
    TH1F *   hTPCPID_2;
    TH1F *   hEff1D;
    TH2F *   hEff2D;
    TH2F *   hEff2D1D;
    TH1F *   hTru1D;
    TH2F *   hTru2D;
    TH1F *   hRec1D;
    TH2F *   hRec2D;
    TH1F *   hRes1D;
    TH2F *   hRes2D;
    TH1F *   hRaw1D;
    TH2F *   hRaw2D;
    TH2F *   hSliceFullX;
    TH2F *   hSliceFullY;
    
    /*
    //TList * hHistoPid = (TList*)(insFile_DT->Get("MyOutputContainerList"));
    TList * hHistoPid = (TList*)(insFile_DT->Get("MyOutputContainerListUTL"));
    
    hName = Form("fHistTOFPID0");
    hTOFPID     =   (TH1F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTOFPID1");
    hTOFPID_1   =   (TH1F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTOFPID2");
    hTOFPID_2   =   (TH1F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTPCPID0");
    hTPCPID     =   (TH1F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTPCPID1");
    hTPCPID_1   =   (TH1F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTPCPID2");
    hTPCPID_2   =   (TH1F*)(hHistoPid->FindObject(hName));
*/
    hName = Form("hNP_1D_Eff_PT_S");
    hEff1D      =   (TH1F*)(insFile_EF->Get(hName));
    hName = Form("hNP_2D_Eff_PT_S");
    hEff2D      =   (TH2F*)(insFile_EF->Get(hName));
    hName = Form("hNP_2D_Eff_X2_S");
    hEff2D1D    =   (TH2F*)(insFile_EF->Get(hName));
    
    hName = Form("hNP_1D_Tru_PT_S");
    hTru1D = (TH1F*)(insFile_EF->Get(hName));
    hName = Form("hNP_2D_Tru_PT_S");
    hTru2D = (TH2F*)(insFile_EF->Get(hName));
    hName = Form("hNP_1D_Rec_PT_S");
    hRec1D = (TH1F*)(insFile_EF->Get(hName));
    hName = Form("hNP_2D_Rec_PT_S");
    hRec2D = (TH2F*)(insFile_EF->Get(hName));
    hName = Form("hNP_1D_Res_PT_S");
    hRes1D = (TH1F*)(insFile_1D->Get(hName));
    hName = Form("hNP_2D_Res_PT_S");
    hRes2D = (TH2F*)(insFile_1D->Get(hName));
    hName = Form("hNP_1D_Raw_PT_S");
    hRaw1D = (TH1F*)(insFile_1D->Get(hName));
    hName = Form("hNP_2D_Raw_PT_S");
    hRaw2D = (TH2F*)(insFile_1D->Get(hName));
    
    hName = Form("FIT_XProjection_PT_Full_DT");
    fSlcPTFncX_DT = (TF1*)(insFile_PT->Get(hName));
    hName = Form("FIT_YProjection_PT_Full_DT");
    fSlcPTFncY_DT = (TF1*)(insFile_PT->Get(hName));
    hName = Form("FIT_XProjection_PT_Full_MC");
    fSlcPTFncX_MC = (TF1*)(insFile_PT->Get(hName));
    hName = Form("FIT_YProjection_PT_Full_MC");
    fSlcPTFncY_MC = (TF1*)(insFile_PT->Get(hName));
    hName = Form("THF_XProjection_PT_Full_DT");
    hSlcPTFncX_DT = (TH1F*)(insFile_PT->Get(hName));
    hName = Form("THF_YProjection_PT_Full_DT");
    hSlcPTFncY_DT = (TH1F*)(insFile_PT->Get(hName));
    hName = Form("THF_XProjection_PT_Full_MC");
    hSlcPTFncX_MC = (TH1F*)(insFile_PT->Get(hName));
    hName = Form("THF_YProjection_PT_Full_MC");
    hSlcPTFncY_MC = (TH1F*)(insFile_PT->Get(hName));
    
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        // Importing 1D Histograms
        hName = Form("hIM_1D_Rec_PT_B_S_%i",iHisto);
        hData_1D[iHisto]        =   (TH1F*)(insFile_PP->Get(hName));
        
    }
    
    auto fName = "";
    auto iHistoDT = 0;
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        // 2D set-up
        hModel2D[iHisto]        =   new TH2F * [nBinPT2D];
        hData_2D[iHisto]        =   new TH2F * [nBinPT2D];
        
        if ( fArrPT2D[iHisto+1] <= 0.41 ) continue;
        
        fName = "DT";
        hName = Form("THF_XProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFitX_2DDT[iHistoDT]   =   (TH1D*)(insFile_PT->Get(hName));
        
        hName = Form("THF_YProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFitY_2DDT[iHistoDT]   =   (TH1D*)(insFile_PT->Get(hName));
        
        hName = Form("FIT_XProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFncX_2DDT[iHistoDT]   =   (TF1*)(insFile_PT->Get(hName));
        
        hName = Form("FIT_YProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFncY_2DDT[iHistoDT]   =   (TF1*)(insFile_PT->Get(hName));
        
        fName = "MC";
        hName = Form("THF_XProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFitX_2DMC[iHistoDT]   =   (TH1D*)(insFile_PT->Get(hName));
        
        hName = Form("THF_YProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFitY_2DMC[iHistoDT]   =   (TH1D*)(insFile_PT->Get(hName));
        
        hName = Form("FIT_XProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFncX_2DMC[iHistoDT]   =   (TF1*)(insFile_PT->Get(hName));
        
        hName = Form("FIT_YProjection_PT_%.1f_%.1f_%s",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fName);
        hSlcPTFncY_2DMC[iHistoDT]   =   (TF1*)(insFile_PT->Get(hName));
        
        iHistoDT++;
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            // Importing 2D Histograms
            hName = Form("hIM_2D_Rec_PT_PT_BB_BS_SB_SS_%i_%i",iHisto,jHisto);
            hData_2D[iHisto][jHisto]    =   (TH2F*)(insFile_PP->Get(hName));
            
            hName = Form("Model2D_%.1f_%.1f_%.1f_%.1f",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hModel2D[iHisto][jHisto]    =   (TH2F*)(insFile_FR->Get(hName));
        }
    }
    //*/
    /*--------------*/
    /*  PRODUCTION  */
    /*--------------*/
    TFile*  outFile_CK  =   new TFile(fFileZZ2DPT,"recreate");
    gStyle->SetOptStat(0);
    /*
    fSlcPTFncX_MC->Write();
    fSlcPTFncY_MC->Write();
    hSlcPTFncX_MC->Write();
    hSlcPTFncY_MC->Write();
    
    hName = "hTru2D";
    hTru2D->Scale(1./nEntries_MC);
    fMosaicCanvas (hTru2D,"",hName,3,3,true,true,true,1e-3,1e-9);
    hName = "hRec2D";
    hRec2D->Scale(1./nEntries_MC);
    fMosaicCanvas (hRec2D,"",hName,3,3,true,true,true,1e-3,1e-9);
    hName = "hRaw2D";
    hRaw2D->Scale(1./nEntries_DT);
    fMosaicCanvas (hRaw2D,"",hName,3,3,true,true,true,1e-3,1e-9);
    hName = "hRes2D";
    hRes2D->Scale(1./nEntries_DT);
    fMosaicCanvas (hRes2D,"",hName,3,3,true,true,true,1e-3,1e-9);
    
    hName = "DT_X_2D_Slices";
    fMosaicCanvas (hSlcPTFitX_2DDT,"",hName,3,3,true,true,hSlcPTFncX_2DDT,"");
    hName = "DT_Y_2D_Slices";
    fMosaicCanvas (hSlcPTFitY_2DDT,"",hName,3,3,true,true,hSlcPTFncY_2DDT,"");
    hName = "MC_X_2D_Slices";
    fMosaicCanvas (hSlcPTFitX_2DMC,"",hName,3,3,true,true,hSlcPTFncX_2DMC,"");
    hName = "MC_Y_2D_Slices";
    fMosaicCanvas (hSlcPTFitY_2DMC,"",hName,3,3,true,true,hSlcPTFncY_2DMC,"");
    
    
    hName = "DC_X_2D_Slices";
    fMosaicCanvas (hSlcPTFitX_2DDT,"",hName,3,3,true,true,hSlcPTFitX_2DMC,"same");
    hName = "DC_Y_2D_Slices";
    fMosaicCanvas (hSlcPTFitX_2DDT,"",hName,3,3,true,true,hSlcPTFitY_2DMC,"same");
    
    TCanvas * cSlcPTFncX_DT = new TCanvas("fSlcPTFncX_DT","fSlcPTFncX_DT");
    hSlcPTFncX_DT->Draw("same");
    fSlcPTFncX_DT->Draw("same");
    cSlcPTFncX_DT->Write();
    cSlcPTFncX_DT->SaveAs("fSlcPTFncX_DT.pdf");
    cSlcPTFncX_DT->SaveAs("fSlcPTFncX_DT.png");
    delete cSlcPTFncX_DT;
    
    TCanvas * cSlcPTFncY_DT = new TCanvas("fSlcPTFncY_DT","fSlcPTFncY_DT");
    hSlcPTFncY_DT->Draw("same");
    fSlcPTFncY_DT->Draw("same");
    cSlcPTFncY_DT->Write();
    cSlcPTFncY_DT->SaveAs("fSlcPTFncY_DT.pdf");
    cSlcPTFncY_DT->SaveAs("fSlcPTFncY_DT.png");
    delete cSlcPTFncX_DT;
    
    TCanvas * cSlcPTFncX_MC = new TCanvas("fSlcPTFncX_MC","fSlcPTFncX_MC");
    hSlcPTFncX_MC->Draw("same");
    fSlcPTFncX_MC->Draw("same");
    cSlcPTFncX_MC->Write();
    cSlcPTFncX_MC->SaveAs("fSlcPTFncX_MC.pdf");
    cSlcPTFncX_MC->SaveAs("fSlcPTFncX_MC.png");
    delete cSlcPTFncX_MC;
    
    TCanvas * cSlcPTFncY_MC = new TCanvas("fSlcPTFncY_MC","fSlcPTFncY_MC");
    hSlcPTFncY_MC->Draw("same");;
    fSlcPTFncY_MC->Draw("same");
    cSlcPTFncY_MC->Write();
    cSlcPTFncY_MC->SaveAs("fSlcPTFncY_MC.pdf");
    cSlcPTFncY_MC->SaveAs("fSlcPTFncY_MC.png");
    delete cSlcPTFncY_MC;
    
    TCanvas * chEff1D = new TCanvas("chEff1D","chEff1D");
    hEff1D->Draw("same");
    chEff1D->Write();
    chEff1D->SaveAs("hEff1D.pdf");
    chEff1D->SaveAs("hEff1D.png");
    delete chEff1D;
    
    TCanvas * chEff2D = new TCanvas("chEff2D","chEff2D");
    hEff2D->Draw("lego");
    chEff2D->Write();
    chEff2D->SaveAs("hEff2D.pdf");
    chEff2D->SaveAs("hEff2D.png");
    delete chEff2D;
    
    TCanvas * chEff2D1D = new TCanvas("chEff2D1D","chEff2D1D");
    hEff2D1D->Draw("lego");
    chEff2D1D->Write();
    chEff2D1D->SaveAs("hEff2D.pdf");
    chEff2D1D->SaveAs("hEff2D.png");
    delete chEff2D1D;
    
    TCanvas * cCompare_hEff2D = new TCanvas("cCompare_hEff2D","cCompare_hEff2D");
    hEff2D1D->Draw("surf2");
    hEff2D->Draw("surf same");
    cCompare_hEff2D->Write();
    cCompare_hEff2D->SaveAs("cCompare_hEff2D.pdf");
    cCompare_hEff2D->SaveAs("cCompare_hEff2D.png");
    delete cCompare_hEff2D;
    
    TCanvas * cCompare_TOFPID = new TCanvas("cCompare_TOFPID","cCompare_TOFPID");
    gPad->SetLogz();
    hTOFPID->Draw("colz");
    hTOFPID_1->SetLineColor(2);
    hTOFPID_2->SetLineColor(3);
    hTOFPID_1->Draw("box same");
    hTOFPID_2->Draw("box same");
    cCompare_TOFPID->Write();
    cCompare_TOFPID->SaveAs("cCompare_TOFPID.pdf");
    cCompare_TOFPID->SaveAs("cCompare_TOFPID.png");
    delete cCompare_TOFPID;
    
    TCanvas * cCompare_TPCPID = new TCanvas("cCompare_TPCPID","cCompare_TPCPID");
    gPad->SetLogz();
    hTPCPID->Draw("colz");
    hTPCPID_1->SetLineColor(2);
    hTPCPID_2->SetLineColor(3);
    hTPCPID_1->Draw("box same");
    hTPCPID_2->Draw("box same");
    cCompare_TPCPID->Write();
    cCompare_TPCPID->SaveAs("cCompare_TPCPID.pdf");
    cCompare_TPCPID->SaveAs("cCompare_TPCPID.png");
    delete cCompare_TPCPID;
    */
    
    /*-----------------*/
    /*  MC Comparison  */
    /*-----------------*/
    
    // Retrieving Histograms
    hName   = "hNP_1D_Tru_PT_S";
    TH1F *  MC1D_Curve_P6   =   (TH1F*)(insFile_MC_P6->Get(hName));
    TH1F *  MC1D_Curve_P8   =   (TH1F*)(insFile_MC_P8->Get(hName));
    hName   = "hNP_2D_Tru_PT_S";
    TH2F *  MC2D_Curve_P6   =   (TH2F*)(insFile_MC_P6->Get(hName));
    TH2F *  MC2D_Curve_P8   =   (TH2F*)(insFile_MC_P8->Get(hName));
    
    // Creating the histograms-------------------------------------------------------------------------------
    
    hName   = "THF_XProjection_PT_Full_DT";
    TH1D *  DTFX_Curve_DT   =   (TH1D*)(insFile_PT->Get(hName));
    hName   = "THF_XProjection_PT_Full_MC";
    TH1F *  MCFX_Curve_P6   =   new TH1F ("","",nBinPT2D,fArrPT2D);
    TH1F *  MCFX_Curve_P8   =   new TH1F ("","",nBinPT2D,fArrPT2D);
    
    hName   = "THF_YProjection_PT_Full_DT";
    TH1D *  DTFY_Curve_DT   =   (TH1D*)(insFile_PT->Get(hName));
    hName   = "THF_XProjection_PT_Full_MC";
    TH1F *  MCFY_Curve_P6   =   new TH1F ("","",nBinPT2D,fArrPT2D);
    TH1F *  MCFY_Curve_P8   =   new TH1F ("","",nBinPT2D,fArrPT2D);
    
    TH1D ***MC2D_UtilX_TRU  = new TH1D** [nMCcompare+1];
    TH1D ***MC2D_UtilY_TRU  = new TH1D** [nMCcompare+1];
    for ( Int_t iHisto = 0; iHisto < nMCcompare+1; iHisto++ )
    {
        MC2D_UtilX_TRU[iHisto] = new TH1D* [nBinPT2D];
        MC2D_UtilY_TRU[iHisto] = new TH1D* [nBinPT2D];
    }
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        hName   = Form("XProjection_PT_%.1f_%.1f_DT",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        MC2D_UtilX_TRU[0][iHisto] = (static_cast<TH1D*>(hRes2D->ProjectionX(hName,iHisto+1,iHisto+1)));
        setMarker(MC2D_UtilX_TRU[0][iHisto]);
        hName   = Form("XProjection_PT_%.1f_%.1f_MC_P6",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        MC2D_UtilX_TRU[1][iHisto] = (static_cast<TH1D*>(MC2D_Curve_P6->ProjectionX(hName,iHisto+1,iHisto+1)));
        setLine(MC2D_UtilX_TRU[1][iHisto],1);
        hName   = Form("XProjection_PT_%.1f_%.1f_MC_P8",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        MC2D_UtilX_TRU[2][iHisto] = (static_cast<TH1D*>(MC2D_Curve_P8->ProjectionX(hName,iHisto+1,iHisto+1)));
        setLine(MC2D_UtilX_TRU[2][iHisto],2);
        
        hName   = Form("YProjection_PT_%.1f_%.1f_DT",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        MC2D_UtilY_TRU[0][iHisto] = (static_cast<TH1D*>(hRes2D->ProjectionY(hName,iHisto+1,iHisto+1)));
        setMarker(MC2D_UtilY_TRU[0][iHisto]);
        hName   = Form("YProjection_PT_%.1f_%.1f_MC_P6",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        MC2D_UtilY_TRU[1][iHisto] = (static_cast<TH1D*>(MC2D_Curve_P6->ProjectionY(hName,iHisto+1,iHisto+1)));
        setLine(MC2D_UtilY_TRU[1][iHisto],1);
        hName   = Form("YProjection_PT_%.1f_%.1f_MC_P8",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        MC2D_UtilY_TRU[2][iHisto] = (static_cast<TH1D*>(MC2D_Curve_P8->ProjectionY(hName,iHisto+1,iHisto+1)));
        setLine(MC2D_UtilY_TRU[2][iHisto],2);
    }
    
    //----------------------//
    //                      //
    //     Production       //
    //                      //
    //----------------------//
    
    //||    1-Dimension     ||\\
    
    // Setting to 0 Option Statistics
    gStyle->SetOptStat(0);
    
    TLegend * fLegend       = new TLegend (0.55,0.60,0.88,0.88);
    fLegend                 ->SetFillColor(kWhite);
    fLegend                 ->SetLineColor(kWhite);
    
    // Setting up the Canvas
    TCanvas * c1DMCCompare = new TCanvas("1DMCCompare","1DMCCompare");
    gPad->SetLogy();                                            // Log Scale in y
    
    // Setting Marker and Line Style
    setMarker(hRes1D);
    setLine(MC1D_Curve_P6,1);
    setLine(MC1D_Curve_P8,2);
    
    // Setting Legend
    fLegend     ->AddEntry(hRes1D,    "Data",        "EP");
    fLegend     ->AddEntry(MC1D_Curve_P6,    "PYTHIA6 Perugia-2011 (350)",    "L");
    fLegend     ->AddEntry(MC1D_Curve_P8,    "PYTHIA8 Monash 2013",    "L");
    
    // Drawing on Canvas
    hRes1D->Draw("same");
    MC1D_Curve_P6->Draw("HIST L same");
    MC1D_Curve_P8->Draw("HIST L same");
    hRes1D->Draw("same");
    fLegend->Draw("same");
    
    c1DMCCompare->Write();
    delete c1DMCCompare;
    
    //||    2-Dimension     ||\\

    Double_t UTL_Err = 0;
    TCanvas * c2DMCCompare_XMosaicSlices = new TCanvas("2DMCCompare_XMosaicSlices","2DMCCompare_XMosaicSlices");
    c2DMCCompare_XMosaicSlices->Divide(3,3);
    
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        TCanvas * c2 = new TCanvas(Form("2DMCCompare_XSlice_%f.1_%f.1",fArrPT2D[iHisto],fArrPT2D[iHisto+1]),Form("2DMCCompare_XSlice_%f.1_%f.1",fArrPT2D[iHisto],fArrPT2D[iHisto+1]));
        gPad->SetLogy();
        
        // Draw on Canvas
        MC2D_UtilX_TRU[0][iHisto]->Draw("same");
        MC2D_UtilX_TRU[1][iHisto]->Draw("HIST L same");
        MC2D_UtilX_TRU[2][iHisto]->Draw("HIST L same");
        fLegend->Draw("same");
        MC2D_UtilX_TRU[0][iHisto]->Draw("same");
        
        MCFX_Curve_P6->SetBinContent(iHisto+1,MC2D_UtilX_TRU[1][iHisto]->IntegralAndError(1,nBinPT2D,UTL_Err,"width"));
        MCFX_Curve_P6->SetBinError(iHisto+1,UTL_Err);
        MCFX_Curve_P8->SetBinContent(iHisto+1,MC2D_UtilX_TRU[2][iHisto]->IntegralAndError(1,nBinPT2D,UTL_Err,"width"));
        MCFX_Curve_P8->SetBinError(iHisto+1,UTL_Err);
        
        if ( iHisto >= 2 )
        {
            c2DMCCompare_XMosaicSlices->cd(iHisto-1);
            gPad->SetLogy();
            MC2D_UtilX_TRU[0][iHisto]->Draw("same");
            MC2D_UtilX_TRU[1][iHisto]->Draw("HIST L same");
            MC2D_UtilX_TRU[2][iHisto]->Draw("HIST L same");
            fLegend->Draw("same");
            MC2D_UtilX_TRU[0][iHisto]->Draw("same");
        }
        
        c2->Write();
        delete c2;
    }
    
    TCanvas * c2DMCCompare_YMosaicSlices = new TCanvas("2DMCCompare_YMosaicSlices","2DMCCompare_YMosaicSlices");
    c2DMCCompare_YMosaicSlices->Divide(3,3);
    
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        TCanvas * c5 = new TCanvas(Form("2DMCCompare_YSlice_%f.1_%f.1",fArrPT2D[iHisto],fArrPT2D[iHisto+1]),Form("2DMCCompare_YSlice_%f.1_%f.1",fArrPT2D[iHisto],fArrPT2D[iHisto+1]));
        gPad->SetLogy();
        
        MC2D_UtilY_TRU[0][iHisto]->Draw("same");
        MC2D_UtilY_TRU[1][iHisto]->Draw("HIST L same");
        MC2D_UtilY_TRU[2][iHisto]->Draw("HIST L same");
        fLegend->Draw("same");
        MC2D_UtilY_TRU[0][iHisto]->Draw("same");
        
        MCFY_Curve_P6->SetBinContent(iHisto+1,MC2D_UtilY_TRU[1][iHisto]->IntegralAndError(1,nBinPT2D,UTL_Err,"width"));
        MCFY_Curve_P6->SetBinError(iHisto+1,UTL_Err);
        MCFY_Curve_P8->SetBinContent(iHisto+1,MC2D_UtilY_TRU[2][iHisto]->IntegralAndError(1,nBinPT2D,UTL_Err,"width"));
        MCFY_Curve_P8->SetBinError(iHisto+1,UTL_Err);
        
        if ( iHisto >= 2 )
        {
            c2DMCCompare_YMosaicSlices->cd(iHisto-1);
            gPad->SetLogy();
            MC2D_UtilY_TRU[0][iHisto]->Draw("same");
            MC2D_UtilY_TRU[1][iHisto]->Draw("HIST L same");
            MC2D_UtilY_TRU[2][iHisto]->Draw("HIST L same");
            fLegend->Draw("same");
            MC2D_UtilY_TRU[0][iHisto]->Draw("same");
        }
        
        c5->Write();
        delete c5;
    }
    
    c2DMCCompare_XMosaicSlices->Write();
    delete c2DMCCompare_XMosaicSlices;
    c2DMCCompare_YMosaicSlices->Write();
    delete c2DMCCompare_YMosaicSlices;
    
    // Setting up the Canvas
    TCanvas * c2DMCCompare_XFull = new TCanvas("2DMCCompare_XFull","2DMCCompare_XFull");
    gPad->SetLogy();                                            // Log Scale in y
    
    // Setting Marker and Line Style
    setMarker(DTFX_Curve_DT);
    setLine(MCFX_Curve_P6,1);
    setLine(MCFX_Curve_P8,2);
    
    // Drawing on Canvas
    DTFX_Curve_DT->Draw("same");
    MCFX_Curve_P6->Draw("HIST L same");
    MCFX_Curve_P8->Draw("HIST L same");
    fLegend->Draw("same");
    DTFX_Curve_DT->Draw("same");
    
    c2DMCCompare_XFull->Write();
    delete c2DMCCompare_XFull;
    
    // Setting up the Canvas
    TCanvas * c2DMCCompare_YFull = new TCanvas("2DMCCompare_YFull","2DMCCompare_YFull");
    gPad->SetLogy();                                            // Log Scale in y
    
    // Setting Marker and Line Style
    setMarker(DTFY_Curve_DT);
    setLine(MCFY_Curve_P6,1);
    setLine(MCFY_Curve_P8,2);
       
    // Drawing on Canvas
    DTFY_Curve_DT->Draw("same");
    MCFY_Curve_P6->Draw("HIST L same");
    MCFY_Curve_P8->Draw("HIST L same");
    fLegend->Draw("same");
    DTFY_Curve_DT->Draw("same");
    
    c2DMCCompare_YFull->Write();
    delete c2DMCCompare_YFull;
    
    // Closing all files
    insFile_1D      ->Close();
    insFiZZ_1D      ->Close();
    insFiZZ_2D      ->Close();
    insFile_PT      ->Close();
    insFile_EF      ->Close();
    outFile_CK      ->Close();
    
    return;
}
