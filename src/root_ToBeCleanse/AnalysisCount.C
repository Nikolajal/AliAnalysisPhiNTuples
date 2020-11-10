// File for 1-Dimensional Analysis:
// !TODO: Make Final Count 2D
// !TODO: Make the actual calculation


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

void AnalysisCount ()
{
    // Retrieving Analysis 1D and 2D files
    TFile*  insFile_1D  =   new TFile   (fFitResHist);
    TFile*  insFile_EF  =   new TFile   (fEfficiHist);
    
    //Target variables
    RooRealVar  xTransverseMom1D("xTransverseMom1D","xTransverseMom1D", fMinPT1D,fMaxPT1D);
    RooRealVar  xTransverseMom2D("xTransverseMom2D","xTransverseMom2D", fMinPT2D,fMaxPT2D);
    RooRealVar  yTransverseMom2D("yTransverseMom2D","yTransverseMom2D", fMinPT2D,fMaxPT2D);
    
    //Recovering histograms in roofit
    // // 1D
    hName                           =   "hNP_1D_Tru_PT_S";
    TH1F * hNP_1D_Tru_PT_S          =   (TH1F*)(insFile_1D->Get(hName));
    hName                           =   "hNP_1D_Res_PT_S";
    TH1F * hNP_1D_Res_PT_S          =   (TH1F*)(insFile_1D->Get(hName));
    hName                           =   "hNP_1D_Rec_PT_S";
    TH1F * hNP_1D_Rec_PT_S          =   (TH1F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_1D_Raw_PT_S";
    TH1F * hNP_1D_Raw_PT_S          =   (TH1F*)(insFile_1D->Get(hName));
    
    // // 2D
    hName                           =   "hNP_2D_Tru_PT_S";
    TH2F * hNP_2D_Tru_PT_S          =   (TH2F*)(insFile_1D->Get(hName));
    hName                           =   "hNP_2D_Res_PT_S";
    TH2F * hNP_2D_Res_PT_S          =   (TH2F*)(insFile_1D->Get(hName));
    hName                           =   "hNP_2D_Rec_PT_S";
    TH2F * hNP_2D_Rec_PT_S          =   (TH2F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_2D_Raw_PT_S";
    TH2F * hNP_2D_Raw_PT_S          =   (TH2F*)(insFile_1D->Get(hName));
    hName                           =   "hNP_2D_Eff_X2_S";
    TH2F * hNP_2D_Eff_X2_S          =   (TH2F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_2D_Eff_PT_S";
    TH2F * hNP_2D_Eff_PT_S          =   (TH2F*)(insFile_EF->Get(hName));
    
    // // Trg
    hName                           =   "hNP_XD_Trg_PT_S";
    TGraphErrors* hNP_XD_Trg_PT_S   =   (TGraphErrors*)(insFile_EF->Get(hName));
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    // Creating the histograms-------------------------------------------------------------------------------
    TGraphErrors * g1D_Res_S        = new TGraphErrors ();
    g1D_Res_S                       ->SetNameTitle("g1D_Res_S");
    
    TGraphErrors * g1D_Res_N        = new TGraphErrors ();
    g1D_Res_N                       ->SetNameTitle("g1D_Res_N");
    
    TGraphErrors * g2D_Res_S        = new TGraphErrors ();
    g2D_Res_S                       ->SetNameTitle("g2D_Res_S");
    
    TGraphErrors * g2D_Res_N        = new TGraphErrors ();
    g2D_Res_N                       ->SetNameTitle("g2D_Res_N");
    
    TGraphErrors * g1D_Trg_S        = new TGraphErrors ();
    g1D_Trg_S                       ->SetNameTitle("g1D_Trg_S");
    
    TGraphErrors * g1D_Trg_N        = new TGraphErrors ();
    g1D_Trg_N                       ->SetNameTitle("g1D_Trg_N");
    
    TGraphErrors * g2D_Trg_S        = new TGraphErrors ();
    g2D_Trg_S                       ->SetNameTitle("g2D_Trg_S");
    
    TGraphErrors * g2D_Trg_N        = new TGraphErrors ();
    g2D_Trg_N                       ->SetNameTitle("g2D_Trg_N");
    
    TGraphErrors * g1D_RMC_S        = new TGraphErrors ();
    g1D_RMC_S                       ->SetNameTitle("g1D_RMC_S");
    
    TGraphErrors * g1D_RMC_N        = new TGraphErrors ();
    g1D_RMC_N                       ->SetNameTitle("g1D_RMC_N");
    
    TGraphErrors * g2D_RMC_S        = new TGraphErrors ();
    g2D_RMC_S                       ->SetNameTitle("g2D_RMC_S");
    
    TGraphErrors * g2D_RMC_N        = new TGraphErrors ();
    g2D_RMC_N                       ->SetNameTitle("g2D_RMC_N");
    
    TGraphErrors * hCheckRec         = new TGraphErrors ();
    hCheckRec       ->SetNameTitle("hCheck1D");
    
    TGraphErrors * hCheckTru         = new TGraphErrors ();
    hCheckTru       ->SetNameTitle("hCheck2D");
    
    TGraphErrors * hCheckRecC        = new TGraphErrors ();
    hCheckRecC      ->SetNameTitle("hCheck1DC");
    
    TGraphErrors * hCheckTruC        = new TGraphErrors ();
    hCheckTruC      ->SetNameTitle("hCheck2D");
    
    TGraphErrors * hCheckBnc        = new TGraphErrors ();
    
    TGraphErrors * hCheckBn2        = new TGraphErrors ();
    
    /*------------*/
    /*  ANALYSIS  */
    /*------------*/
    
    // Get Entries
    //Retrieving Number of Events considered
    TH1F *  hUtlEntry   =   (TH1F*)(insFile_EF->Get("Entry_MC_DT"));
    Int_t   nEntries_DT =   1;//hUtlEntry->GetBinContent(2);
    Int_t   nEntries_MC =   1;//hUtlEntry->GetBinContent(1);
    
    // Normalising Factors
    double TRG_1_Y = 0;
    double TRG_2_Y = 0;
    double TRG_Nul = 0;
    double TRR_1_Y = 0;
    double TRR_2_Y = 0;
    double UTL_Er1 = 0;
    double UTL_Er2 = 0;
    hNP_XD_Trg_PT_S->GetPoint(0,TRG_Nul,TRG_1_Y);
    hNP_XD_Trg_PT_S->GetPoint(1,TRG_Nul,TRG_2_Y);
    TRR_1_Y = hNP_XD_Trg_PT_S->GetErrorY(0);
    TRR_2_Y = hNP_XD_Trg_PT_S->GetErrorY(1);
    hNP_XD_Trg_PT_S->SetPoint(1,2,1.0*TRG_2_Y);
    hNP_XD_Trg_PT_S->SetPointError(1,0,1.0*TRR_2_Y);
    
    // Output File for Fit Check
    TFile*  outFile_CK  =   new TFile(fFileZZ1DPT,"recreate");
    hCheckBnc->SetPoint(0,0.9,1);
    hCheckBnc->SetPoint(1,4.5,1);
    hCheckBnc->SetLineColor(kBlack);
    hCheckBnc->SetLineStyle(10);
    
    hCheckBn2->SetPoint(0,0.9,0.03);
    hCheckBn2->SetPoint(1,4.5,0.03);
    hCheckBn2->SetLineColor(kWhite);
    hCheckBn2->SetMarkerColor(kWhite);
    
    TCanvas * cCheckRaw = new TCanvas("cCheckRaw","cCheckRaw");
    gPad->SetLogz();
    hNP_2D_Raw_PT_S->SetLineColor(kRed);
    hNP_2D_Raw_PT_S->Draw("lego");
    hNP_2D_Rec_PT_S->Draw("lego same");
    cCheckRaw->Write();
    cCheckRaw->SaveAs("cCheckRaw.pdf");
    cCheckRaw->SaveAs("cCheckRaw.png");
    delete cCheckRaw;
    
    TCanvas * cCheckRes = new TCanvas("cCheckRes","cCheckRes");
    gPad->SetLogz();
    hNP_2D_Res_PT_S->SetLineColor(kRed);
    hNP_2D_Res_PT_S->Draw("lego");
    hNP_2D_Tru_PT_S->Draw("lego same");
    cCheckRes->Write();
    cCheckRes->SaveAs("cCheckRes.pdf");
    cCheckRes->SaveAs("cCheckRes.png");
    delete cCheckRes;
    
    // Fraction Check 1D
    TH1F Utility1 = *hNP_1D_Raw_PT_S;
    Utility1.Divide(hNP_1D_Rec_PT_S);
    Utility1.Fit(fFlatFit1D,"Q","",0.4,10.);
    Utility1.SetNameTitle("Check1D_Raw_Rec","");
    Utility1.Write();
    hCheckRec->SetPoint             (0,1.2,fFlatFit1D->GetParameter(0));
    hCheckRec->SetPointError        (0,0,fFlatFit1D->GetParError(0));
    
    TH1F Utility2 = *hNP_1D_Res_PT_S;
    Utility2.Divide(hNP_1D_Tru_PT_S);
    Utility2.Fit(fFlatFit1D,"Q","",0.4,10.);
    Utility2.SetNameTitle("Check1D_Res_Tru","");
    Utility2.Write();
    hCheckTru->SetPoint             (0,1,fFlatFit1D->GetParameter(0));
    hCheckTru->SetPointError        (0,0,fFlatFit1D->GetParError(0));
    
    cout << "1D" << endl;
    auto UTL_In1 = hNP_1D_Raw_PT_S->IntegralAndError(5,nBinPT1D,UTL_Er1,"width");
    auto UTL_In2 = hNP_1D_Rec_PT_S->IntegralAndError(5,nBinPT1D,UTL_Er2,"width");
    cout << "Integral RAW:" << UTL_In1 << " +- " << UTL_Er1 << endl;
    cout << "Integral REC:" << UTL_In2 << " +- " << UTL_Er2 << endl;
    cout << "Integral Ratio:" << UTL_In1/UTL_In2 << " +- " << (UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)) << endl;
    hCheckRecC->SetPoint            (0,1.3,UTL_In1/UTL_In2);
    hCheckRecC->SetPointError       (0,0,(UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)));

         UTL_In1 = hNP_1D_Res_PT_S->IntegralAndError(5,nBinPT1D,UTL_Er1,"width");
         UTL_In2 = hNP_1D_Tru_PT_S->IntegralAndError(5,nBinPT1D,UTL_Er2,"width");
         cout << "Integral RES:" << UTL_In1 << " +- " << UTL_Er1 << endl;
         cout << "Integral TRU:" << UTL_In2 << " +- " << UTL_Er2 << endl;
         cout << "Integral Ratio:" << UTL_In1/UTL_In2 << " +- " << (UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)) << endl;
    hCheckTruC->SetPoint            (0,1.1,UTL_In1/UTL_In2);
    hCheckTruC->SetPointError       (0,0,(UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)));
    
    // Fraction Check 2D
    TH2F Utility3 = *hNP_2D_Raw_PT_S;
    Utility3.Divide(hNP_2D_Rec_PT_S);
    Utility3.Write();
    Utility3.Fit(fFlatFit2D,"Q","",0.4,10);
    Utility3.SetNameTitle("Check2D_Raw_Rec","");
    Utility3.Write();
    hCheckRec->SetPoint             (1,2.2,fFlatFit2D->GetParameter(0));
    hCheckRec->SetPointError        (1,0,fFlatFit2D->GetParError(0));
    
    TH2F Utility4 = *hNP_2D_Res_PT_S;
    Utility4.Divide(hNP_2D_Tru_PT_S);
    Utility4.Write();
    Utility4.Fit(fFlatFit2D,"Q","",0.4,10);
    Utility4.SetNameTitle("Check2D_Res_Tru","");
    Utility4.Write();
    hCheckTru->SetPoint             (1,2,fFlatFit2D->GetParameter(0));
    hCheckTru->SetPointError        (1,0,fFlatFit2D->GetParError(0));
    
    cout << "2D" << endl;
         UTL_In1 = hNP_2D_Raw_PT_S->IntegralAndError(3,nBinPT2D,3,nBinPT2D,UTL_Er1,"width");
         UTL_In2 = hNP_2D_Rec_PT_S->IntegralAndError(3,nBinPT2D,3,nBinPT2D,UTL_Er2,"width");
    cout << "Integral RAW:" << UTL_In1 << " +- " << UTL_Er1 << endl;
    cout << "Integral REC:" << UTL_In2 << " +- " << UTL_Er2 << endl;
    cout << "Integral Ratio:" << UTL_In1/UTL_In2 << " +- " << (UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)) << endl;
    
    hCheckRecC->SetPoint            (1,2.3,UTL_In1/UTL_In2);
    hCheckRecC->SetPointError       (1,0,(UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)));

         UTL_In1 = hNP_2D_Res_PT_S->IntegralAndError(3,nBinPT2D,3,nBinPT2D,UTL_Er1,"width");
         UTL_In2 = hNP_2D_Tru_PT_S->IntegralAndError(3,nBinPT2D,3,nBinPT2D,UTL_Er2,"width");
         cout << "Integral Ratio:" << UTL_In1/UTL_In2 << " +- " << (UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)) << endl;
    cout << "Integral RES:" << UTL_In1 << " +- " << UTL_Er1 << endl;
    cout << "Integral TRU:" << UTL_In2 << " +- " << UTL_Er2 << endl;
    
    hCheckTruC->SetPoint            (1,2.1,UTL_In1/UTL_In2);
    hCheckTruC->SetPointError       (1,0,(UTL_In1/UTL_In2)*sqrt(UTL_Er1*UTL_Er1/(UTL_In1*UTL_In1)+UTL_Er2*UTL_Er2/(UTL_In2*UTL_In2)));
    
    // Second Check
    fSliceCheck(hNP_1D_Raw_PT_S,hNP_1D_Rec_PT_S,"Raw","Rec","1D_Raw_Rec","1D_Raw_Rec");
    hNP_2D_Raw_PT_S->Write();
    hNP_2D_Rec_PT_S->Write();
    fSliceCheckPT2D(hNP_2D_Raw_PT_S,hNP_2D_Rec_PT_S,"Raw","Rec");
    fSliceCheck(hNP_1D_Res_PT_S,hNP_1D_Tru_PT_S,"Res","Tru","1D_Res_Tru","1D_Res_Tru");
    hNP_2D_Res_PT_S->Write();
    hNP_2D_Tru_PT_S->Write();
    fSliceCheckPT2D(hNP_2D_Res_PT_S,hNP_2D_Tru_PT_S,"Res","Tru");
    
    // Efficency Check
    fSliceCheckPT2D(hNP_2D_Eff_X2_S,hNP_2D_Eff_PT_S,"1D","2D");
    
    TF1 aaa;
    // 1D MC
    FitModelPT1D(hNP_1D_Tru_PT_S,nEntries_MC,aaa,true);
    g1D_RMC_S->SetPoint       (0,1,fLevyFit1D->GetParameter(3));
    g1D_RMC_S->SetPointError  (0,0,fLevyFit1D->GetParError(3));
    g1D_RMC_N->SetPoint       (0,1,(fLevyFit1D->GetParameter(3))/(TRG_1_Y));
    g1D_RMC_N->SetPointError  (0,0,sqrt(TMath::Power((fLevyFit1D->GetParError(3)),2)+TMath::Power((hNP_XD_Trg_PT_S->GetErrorY(0)),2))/(TRG_1_Y));
    
    // 1D Data
    FitModelPT1D(hNP_1D_Res_PT_S,nEntries_DT,aaa,true);
    g1D_Res_S->SetPoint       (0,1,fLevyFit1D->GetParameter(3));
    g1D_Res_S->SetPointError  (0,0,fLevyFit1D->GetParError(3));
    g1D_Res_N->SetPoint       (0,1,(fLevyFit1D->GetParameter(3))/(TRG_1_Y));
    g1D_Res_N->SetPointError  (0,0,sqrt(TMath::Power((fLevyFit1D->GetParError(3)),2)+TMath::Power((hNP_XD_Trg_PT_S->GetErrorY(0)),2))/(TRG_1_Y));
    
    // 2D MC
    FitModelPT2D(hNP_2D_Tru_PT_S,nEntries_MC,xTransverseMom1D,true,"MC");
    g2D_RMC_S->SetPoint       (0,2,fLevyFit1D->GetParameter(3));
    g2D_RMC_S->SetPointError  (0,0,fLevyFit1D->GetParError(3));
    g2D_RMC_N->SetPoint       (0,2,(fLevyFit1D->GetParameter(3))/(TRG_2_Y));
    g2D_RMC_N->SetPointError  (0,0,sqrt(TMath::Power((fLevyFit1D->GetParError(3)),2)+TMath::Power((hNP_XD_Trg_PT_S->GetErrorY(1)),2))/(TRG_2_Y));
    
    // 2D Data
    FitModelPT2D(hNP_2D_Res_PT_S,nEntries_DT,xTransverseMom1D,true,"DT");
    g2D_Res_S->SetPoint       (0,2,fLevyFit1D->GetParameter(3));
    g2D_Res_S->SetPointError  (0,0,fLevyFit1D->GetParError(3));
    g2D_Res_N->SetPoint       (0,2,(fLevyFit1D->GetParameter(3))/(TRG_2_Y));
    g2D_Res_N->SetPointError  (0,0,sqrt(TMath::Power((fLevyFit1D->GetParError(3)),2)+TMath::Power((hNP_XD_Trg_PT_S->GetErrorY(1)),2))/(TRG_2_Y));
    
    // Normalised MC
    g1D_Trg_N->SetPoint       (0,1,1);
    g1D_Trg_N->SetPointError  (0,0,(hNP_XD_Trg_PT_S->GetErrorY(0))/TRG_1_Y);
    g2D_Trg_N->SetPoint       (0,2,1);
    g2D_Trg_N->SetPointError  (0,0,(hNP_XD_Trg_PT_S->GetErrorY(1))/TRG_2_Y);
    g1D_Trg_S->SetPoint       (0,1,(TRG_1_Y));
    g1D_Trg_S->SetPointError  (0,0,(hNP_XD_Trg_PT_S->GetErrorY(0)));
    g2D_Trg_S->SetPoint       (0,2,(TRG_2_Y));
    g2D_Trg_S->SetPointError  (0,0,(hNP_XD_Trg_PT_S->GetErrorY(1)));
    
    
    g1D_Res_S   ->SetMarkerStyle(22);
    g2D_Res_S   ->SetMarkerStyle(22);
    g1D_Res_N   ->SetMarkerStyle(22);
    g2D_Res_N   ->SetMarkerStyle(22);
    g1D_Res_S   ->SetMarkerColor(02);
    g2D_Res_S   ->SetMarkerColor(02);
    g1D_Res_N   ->SetMarkerColor(02);
    g2D_Res_N   ->SetMarkerColor(02);
    
    g1D_Trg_S   ->SetMarkerStyle(22);
    g2D_Trg_S   ->SetMarkerStyle(22);
    g1D_Trg_N   ->SetMarkerStyle(22);
    g2D_Trg_N   ->SetMarkerStyle(22);
    g1D_Trg_S   ->SetMarkerColor(04);
    g2D_Trg_S   ->SetMarkerColor(04);
    g1D_Trg_N   ->SetMarkerColor(04);
    g2D_Trg_N   ->SetMarkerColor(04);
    
    g1D_RMC_S   ->SetMarkerStyle(22);
    g2D_RMC_S   ->SetMarkerStyle(22);
    g1D_RMC_N   ->SetMarkerStyle(22);
    g2D_RMC_N   ->SetMarkerStyle(22);
    g1D_RMC_S   ->SetMarkerColor(8);
    g2D_RMC_S   ->SetMarkerColor(8);
    g1D_RMC_N   ->SetMarkerColor(8);
    g2D_RMC_N   ->SetMarkerColor(8);
    
    hCheckRec           ->SetMarkerStyle(22);
    hCheckRec           ->SetMarkerColor(02);
    hCheckTru           ->SetMarkerStyle(23);
    hCheckTru           ->SetMarkerColor(04);
    hCheckRecC          ->SetMarkerStyle(24);
    hCheckRecC          ->SetMarkerColor(02);
    hCheckTruC          ->SetMarkerStyle(25);
    hCheckTruC          ->SetMarkerColor(04);
    hCheckRecC->Write();
    hCheckTruC->Write();

    TCanvas * cFinal_        = new TCanvas   ("gFinal_","gFinal_");
    cFinal_->SetBorderSize(0);
    cFinal_->SetFillColor(kWhite);

    TLegend * fLegend       = new TLegend   (0.35,0.30,0.65,0.60);
    fLegend->SetTextFont(22);
    fLegend->SetBorderSize(2);
    fLegend->SetLineColor(kBlack);
    fLegend->SetFillColor(kWhite);
    
    TMultiGraph * m1D   = new TMultiGraph();
    m1D->Add(g1D_Res_S);
    fLegend->AddEntry(g1D_Res_S,"MC Recreate","lp");
    m1D->Add(g1D_Trg_S);
    fLegend->AddEntry(g1D_Trg_S,"MC Truth","lp");
    m1D->Add(g1D_RMC_S);
    fLegend->AddEntry(g1D_RMC_S,"MC Fit","lp");
    
    TMultiGraph * m2D   = new TMultiGraph();
    m2D->Add(g2D_Res_S);
    m2D->Add(g2D_Trg_S);
    m2D->Add(g2D_RMC_S);

    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    pad1->SetFillStyle(4000);
    pad1->SetFillColor(kWhite);
    pad1->SetFrameBorderMode(0);
    pad1->Draw();
    pad1->cd();
    
    m1D->GetXaxis()->SetLimits(0.85,2.15);
    m1D->SetMaximum(0.034);
    m1D->SetMinimum(0.027);
    m1D->Draw("ALP");
    
    cFinal_->cd();
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000); //will be transparent
    pad2->SetFillColor(kWhite);
    pad2->SetFrameBorderMode(0);
    pad2->SetFrameFillStyle(0);
    pad2->Draw();
    pad2->cd();
    
    m2D->GetXaxis()->SetLimits(0.85,2.15);
    m2D->SetMaximum(0.0019);
    m2D->SetMinimum(0.001);
    (m2D->GetHistogram())->Draw("Y+");
    m2D->Draw("LP");
    
    fLegend->Draw();
    cFinal_->cd();
    cFinal_->Write();
    //delete cFinal_;
    
    TCanvas * cFinalN        = new TCanvas   ("gFinalN","gFinalN");
    cFinalN->SetBorderSize(0);
    cFinalN->SetFillColor(kWhite);
    
    TMultiGraph * m1DN  = new TMultiGraph();
    m1DN->Add(g1D_Res_N);
    m1DN->Add(g1D_Trg_N);
    m1DN->Add(g1D_RMC_N);
    m1DN->Add(g2D_Res_N);
    m1DN->Add(g2D_Trg_N);
    m1DN->Add(g2D_RMC_N);
    m1DN->Draw("ALP");
    
    fLegend->Draw();
    cFinalN->cd();
    cFinalN->Write();
    //delete cFinalN;
    
    /*
    TCanvas * c1        = new TCanvas   ("gFinal_","gFinal_");
    TLegend * fLegend1  = new TLegend   (0.36,0.27,0.64,0.73);
    TMultiGraph * m1D   = new TMultiGraph();
    TMultiGraph * m2D   = new TMultiGraph();

    m1D->Add(g1D_Res_S);
    fLegend1->AddEntry(g1D_Res_S,"Strip coverage","lp");
    m1D->Add(g1D_Trg_S);
    fLegend1->AddEntry(g1D_Trg_S,"Strip covera_ge","lp");
    m1D->Add(g1D_RMC_S);
    fLegend1->AddEntry(g1D_RMC_S,"Strip cov_erage","lp");
    m2D->Add(g2D_Res_S);
    m2D->Add(g2D_Trg_S);
    m2D->Add(g2D_RMC_S);
    
    c1->cd();
    TPad *pad1          = new TPad      ("pad1","",0,0,1,1);
    pad1->Draw("same");
    pad1->cd();
    m1D->Draw();
    
    c1->cd();
    TPad *pad2          = new TPad      ("pad2","",0,0,1,1);
    pad2->Draw("same");
    pad2->cd();
    m2D->Draw("X+Y+");

    fLegend1->Draw("same");
    c1->cd();
    
    
    TMultiGraph * m1    = new TMultiGraph();
    hName               = "hNP_XD_XXX_PT_S";
    hTitle              = "Multidimensional #phi production statistics";
    m1                  ->SetNameTitle(hName,hTitle);
    m1                  ->GetXaxis()->SetTitle("");
    m1                  ->GetYaxis()->SetTitle("#frac{dN_{#phi}}{dy}");
    m1                  ->Add(hNP_XD_RMC_PT_S);
    m1                  ->Add(hCheckBn2);
    m1                  ->Add(hNP_XD_Res_PT_S);
    m1                  ->Add(hNP_XD_Trg_PT_S);
    TLegend * fLegend1  = new TLegend   (0.36,0.27,0.64,0.73);
    fLegend1            ->SetFillColor(kWhite);
    fLegend1            ->SetLineColor(kWhite);
    fLegend1            ->AddEntry(hNP_XD_Res_PT_S,"Data Analysis","leP");
    fLegend1            ->AddEntry(hNP_XD_RMC_PT_S,"MC Analysis","lep");
    fLegend1            ->AddEntry(hNP_XD_Trg_PT_S,"MCTruth","lep");
    //fLegend1            ->AddEntry(hNP_2D_RMC_PT_S,"2D Fit (MC)","lep");
    m1                  ->Draw("ap");
    fLegend1            ->Draw("same");
    gStyle->SetOptStat(0);
    c1                  ->Write();
    c1                  ->SaveAs("gFinal_.pdf");
    c1                  ->SaveAs("gFinal_.png");
    delete c1;
    
    TCanvas * c2        = new TCanvas("gFinalN","gFinalN");
    
    TMultiGraph * m2    = new TMultiGraph();
    hName               = "hNP_XD_XX2_PT_S";
    hTitle              = "Multidimensional #phi production statistics normalised to MC";
    m2                  ->SetNameTitle(hName,hTitle);
    m2                  ->GetXaxis()->SetTitle("N-Tuples");
    m2                  ->GetYaxis()->SetTitle("#frac{dN_{#phi}}{dy}");
    m2                  ->Add(hNP_XD_RMC_PT_N);
    m2                  ->Add(hNP_XD_Res_PT_N);
    m2                  ->Add(hNP_XD_Trg_PT_N);
    m2                  ->Add(hCheckBnc,"L");
    TLegend * fLegend3  = new TLegend   (0.6,0.4,0.9,0.9);
    fLegend3            ->SetFillColor(kWhite);
    fLegend3            ->SetLineColor(kBlack);
    fLegend3            ->AddEntry(hNP_XD_Res_PT_S,"Data Analysis","leP");
    fLegend3            ->AddEntry(hNP_XD_RMC_PT_S,"MC Analysis","lep");
    fLegend3            ->AddEntry(hNP_XD_Trg_PT_S,"MCTruth","lep");
    m2                  ->Draw("ap");
    fLegend3            ->Draw("same");
    gStyle->SetOptStat(0);
    c2                  ->Write();
    c2                  ->SaveAs("gFinalN.pdf");
    c2                  ->SaveAs("gFinalN.png");
    
    
    TCanvas * c3        = new TCanvas();
    TMultiGraph * m3    = new TMultiGraph();
    hName               = "hCheck";
    hTitle              = "Check Values";
    m3                  ->SetNameTitle(hName,hTitle);
    m3                  ->GetXaxis()->SetTitle("N-Tuples");
    m3                  ->GetYaxis()->SetTitle("");
    m3                  ->Add(hCheckRec);
    m3                  ->Add(hCheckRecC);
    m3                  ->Add(hCheckBnc,"L");
    TLegend * fLegend2  = new TLegend   (0.62,0.43,0.89,0.89);
    fLegend2            ->SetFillColor(kWhite);
    fLegend2            ->SetLineColor(kWhite);
    fLegend2            ->AddEntry(hCheckRec,"N^{Raw}/N^{Rec} (Fit)","lep");
    fLegend2            ->AddEntry(hCheckRecC,"N^{Raw}/N^{Rec} (Int)","lep");
    m3                  ->Draw("ap");
    fLegend2            ->Draw("same");
    c3                  ->Write();
    delete c3;
    
    TCanvas * c4        = new TCanvas();
    m3                  ->Add(hCheckTru);
    m3                  ->Add(hCheckTruC);
    fLegend2            ->AddEntry(hCheckTru,"N^{Res}/N^{Tru} (Fit)","lep");
    fLegend2            ->AddEntry(hCheckTruC,"N^{Res}/N^{Tru} (Int)","lep");
    m3                  ->Draw("ap");
    fLegend2            ->Draw("same");
    c4                  ->Write();
    delete c4;
    
    //Writing Results on File
    hNP_XD_Res_PT_S ->Write();
    hNP_XD_RMC_PT_S ->Write();
    hNP_XD_Trg_PT_S ->Write();
    hNP_XD_Res_PT_N ->Write();
    hNP_XD_RMC_PT_N ->Write();
    */
    
    insFile_1D      ->Close();
    outFile_CK      ->Close();
    
    return;
}
