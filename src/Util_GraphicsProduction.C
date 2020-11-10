#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

void Util_GraphicsProduction ()
{
    // Retrieving Analysis 1D and 2D files
    TFile   *in_Efficiency_Pythia8          =   new TFile   (fEfficiHist);
    TFile   *in_RawYieldHist_Pythia8        =   new TFile   (fFitResults);
    TFile   *in_ResYieldHist_Pythia8        =   new TFile   (fAnlResHist);
    
    //------------------//
    //  Importing input //-----------------------------------------------------------------------------------
    //------------------//
    
    // Generating the binning array
    
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    
    //Fetching histograms
    
    // - // All Efficiencies
    
    // - // - // 1D
    hName   =   "hNP_1D_Eff_PT_S";
    TH1F    *fh_Efficiency_1D_1Dbin =   (TH1F*)(in_Efficiency_Pythia8->Get(hName));
    
    hName   =   "hNP_1D_Eff_X2_S";
    TH1F    *fh_Efficiency_1D_2Dbin =   (TH1F*)(in_Efficiency_Pythia8->Get(hName));
    
    // - // - // 2D
    hName   =   "hNP_2D_Eff_PT_S";
    TH2F    *fh_Efficiency_2D_2Dbin =   (TH2F*)(in_Efficiency_Pythia8->Get(hName));
    
    hName   =   "hNP_2D_Eff_X2_S";
    TH2F    *fh_Efficiency_2D_1Dsqr =   (TH2F*)(in_Efficiency_Pythia8->Get(hName));

    
    // - // All Montecarlo Truths
    
    // - // - // 1D
    hName   =   "hNP_1D_Rec_PT_S";
    TH1F    *fh_Rec_1D              =   (TH1F*)(in_Efficiency_Pythia8->Get(hName));
    
    hName   =   "hNP_1D_Gen_PT_S";
    TH1F    *fh_Gen_1D              =   (TH1F*)(in_Efficiency_Pythia8->Get(hName));
    
    hName   =   "hNP_1D_Tru_PT_S";
    TH1F    *fh_Tru_1D              =   (TH1F*)(in_Efficiency_Pythia8->Get(hName));
    
    // - // - // 2D
    hName   =   "hNP_2D_Rec_PT_S";
    TH2F    *fh_Rec_2D              =   (TH2F*)(in_Efficiency_Pythia8->Get(hName));
    
    hName   =   "hNP_2D_Gen_PT_S";
    TH2F    *fh_Gen_2D              =   (TH2F*)(in_Efficiency_Pythia8->Get(hName));
    
    hName   =   "hNP_2D_Tru_PT_S";
    TH2F    *fh_Tru_2D              =   (TH2F*)(in_Efficiency_Pythia8->Get(hName));
    
    
    // - // All Data Yields
    
    // - // - // 1D
    hName   =   "h1D_Raw";
    TH1F    *fh_Raw_1D              =   (TH1F*)(in_RawYieldHist_Pythia8->Get(hName));
    
    hName   =   "Entry_DT";
    TH1F    *fh_Ent_1D              =   (TH1F*)(in_RawYieldHist_Pythia8->Get(hName));
    
    hName   =   "hSliceFX_";
    TH1F    *fh_Res_FX              =   (TH1F*)(in_ResYieldHist_Pythia8->Get(hName));
    
    hName   =   "hSliceFY_";
    TH1F    *fh_Res_FY              =   (TH1F*)(in_ResYieldHist_Pythia8->Get(hName));
    
    // - // - // 2D
    hName   =   "h2D_Raw";
    TH2F    *fh_Raw_2D              =   (TH2F*)(in_RawYieldHist_Pythia8->Get(hName));
    
    //---------------------//
    //  Setting up output  //--------------------------------------------------------------------------------
    //---------------------//
    
    // output File
    
    TFile   *ot_Graphical   =   new TFile   ("./graphs/graphs.root","RECREATE");
    
    // Scaling RAW
    
    fh_Raw_1D->Scale(1./(fh_Ent_1D->GetBinContent(1)*kPythia1DEfficien),"width");
    fh_Raw_2D->Scale(1./(fh_Ent_1D->GetBinContent(1)*kPythia2DEfficien),"width");
    
    
    
    // Utility (some to circumvent stupid root)
    
    TH1D       *hUtility1X,   *hUtility2X,    *hUtility1Y,    *hUtility2Y;
    TH1D       *hLogEffic_  =   new TH1D    ("EF1D","EF1D",nBinPT1D,fArrPT1D);
    TH1D      **hLogEfficX  =   new TH1D   *[nBinPT2D];
    TH1D      **hLogEfficY  =   new TH1D   *[nBinPT2D];
    TLegend    *fLegend1    =   new TLegend(0.2,0.2,0.8,0.8);
    TLegend    *fLegend2    =   new TLegend(0.2,0.2,0.8,0.8);
    TLegend    *fLegend3    =   new TLegend(0.68,0.68,0.88,0.88);
    TLegend    *fLegend4    =   new TLegend(0.68,0.68,0.88,0.88);
    TLegend    *fLegend5    =   new TLegend(0.68,0.68,0.88,0.88);
    TLegend    *fLegend6    =   new TLegend(0.5,0.5,0.88,0.88);
    TLegend    *fLegend7    =   new TLegend(0.5,0.5,0.88,0.88);
    
    //---------------------//
    //  Generating output  //--------------------------------------------------------------------------------
    //---------------------//
    
    // Efficiencies
    
    // - // 1D Graphics
    
    TCanvas     *fc_Efficiency  =   new TCanvas("","");
    gStyle->SetOptStat(0);
    hLogEffic_  ->  SetTitle(Form("1D Efficiency"));
    hLogEffic_  ->  SetTitle(Form(""));
    hLogEffic_  ->  GetXaxis()  ->  SetTitle("p_{T}#phi (Gev/c)");
    hLogEffic_  ->  GetYaxis()  ->  SetTitle("#varepsilon");
    hLogEffic_  ->  SetMaximum(1.);
    hLogEffic_  ->  SetMinimum(0.);
    hLogEffic_  ->  Draw("");
    fh_Efficiency_1D_1Dbin  ->  SetMarkerColor(kBlue);
    fh_Efficiency_1D_1Dbin  ->  SetLineColor(kBlue);
    fh_Efficiency_1D_1Dbin  ->  SetMarkerStyle(26);
    fh_Efficiency_1D_1Dbin  ->  SetMarkerSize(1);
    fh_Efficiency_1D_1Dbin  ->  Draw("EP SAME");
    
    // - // 2D Graphics
    
    TCanvas     *fc_Efficiency_XSliced  =   new TCanvas("","",855,1440);
    fc_Efficiency_XSliced->Divide(3,4);
    TCanvas     *fc_Efficiency_YSliced  =   new TCanvas("","",855,1440);
    fc_Efficiency_YSliced->Divide(3,4);
    
    for ( int iChk = 0; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        fc_Efficiency_XSliced->cd(iChk+1);
        gStyle->SetOptStat(0);
        hLogEfficX[iChk]    =   new TH1D    (Form("EFX_%i",iChk),Form("EFX_%i",iChk),nBinPT2D,fArrPT2D);
        hLogEfficX[iChk]    ->  SetTitle(Form("2D Efficiency for %.1f < p_{T}#phi_{1} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hLogEfficX[iChk]    ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{2} (Gev/c)");
        hLogEfficX[iChk]    ->  GetYaxis()  ->  SetTitle("#varepsilon");
        hLogEfficX[iChk]    ->  Draw("");
        hUtility1X  =   (fh_Efficiency_2D_2Dbin->ProjectionX(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1X  ->  SetMarkerColor(kBlue);
        hUtility1X  ->  SetLineColor(kBlue);
        hUtility1X  ->  SetMarkerStyle(26);
        hUtility1X  ->  SetMarkerSize(1);
        hUtility1X  ->  Draw("EP SAME");
        hUtility2X  =   (fh_Efficiency_2D_1Dsqr->ProjectionX(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2X  ->  SetMarkerColor(kRed);
        hUtility2X  ->  SetLineColor(kRed);
        hUtility2X  ->  SetMarkerStyle(32);
        hUtility2X  ->  SetMarkerSize(1);
        hUtility2X  ->  Draw("EP SAME");
        if ( iChk == 0 )
        {
            fLegend1    ->  SetTextSize(0.06);
            fLegend1    ->  AddEntry(hUtility1X,"2DEff measured","EP");
            fLegend1    ->  AddEntry(hUtility2X,"2DEff as #varepsilon_{1D}#times#varepsilon_{1D}","EP");
            fLegend1    ->  Draw("SAME");
        }
        
        // Y-Projection
        fc_Efficiency_YSliced->cd(iChk+1);
        gStyle->SetOptStat(0);
        hLogEfficY[iChk]    =   new TH1D    (Form("EFY_%i",iChk),Form("EFY_%i",iChk),nBinPT2D,fArrPT2D);
        hLogEfficY[iChk]    ->  SetTitle(Form("2D Efficiency for %.1f < p_{T}#phi_{2} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hLogEfficY[iChk]    ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{1} (Gev/c)");
        hLogEfficY[iChk]    ->  GetYaxis()  ->  SetTitle("#varepsilon");
        hLogEfficY[iChk]    ->  Draw("");
        hUtility1Y  =   (fh_Efficiency_2D_2Dbin->ProjectionY(Form("2DEfficiency_YSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1Y  ->  SetMarkerColor(kBlue);
        hUtility1Y  ->  SetLineColor(kBlue);
        hUtility1Y  ->  SetMarkerStyle(26);
        hUtility1Y  ->  SetMarkerSize(1);
        hUtility1Y  ->  Draw("EP SAME");
        hUtility2Y  =   (fh_Efficiency_2D_1Dsqr->ProjectionY(Form("2DEfficiency_YSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2Y  ->  SetMarkerColor(kRed);
        hUtility2Y  ->  SetLineColor(kRed);
        hUtility2Y  ->  SetMarkerStyle(32);
        hUtility2Y  ->  SetMarkerSize(1);
        hUtility2Y  ->  Draw("EP SAME");
        if ( iChk == 0 )
        {
            fLegend2    ->  SetTextSize(0.06);
            fLegend2    ->  AddEntry(hUtility1Y,"2DEff measured","EP");
            fLegend2    ->  AddEntry(hUtility2Y,"2DEff as #varepsilon_{1D}#times#varepsilon_{1D}","EP");
            fLegend2    ->  Draw("SAME");
        }
    }
    
    // Saving Efficiency Graphics
    fc_Efficiency->Write();
    fc_Efficiency->SaveAs("./graphs/1DEfficiency.pdf");
    fc_Efficiency->SaveAs("./graphs/1DEfficiency.png");
    
    fc_Efficiency_XSliced->Write();
    fc_Efficiency_XSliced->SaveAs("./graphs/2DEfficiencySliceX.pdf");
    fc_Efficiency_XSliced->SaveAs("./graphs/2DEfficiencySliceX.png");
    
    fc_Efficiency_YSliced->Write();
    fc_Efficiency_YSliced->SaveAs("./graphs/2DEfficiencySliceY.pdf");
    fc_Efficiency_YSliced->SaveAs("./graphs/2DEfficiencySliceY.png");
    
    // p$_T$ Spectra
    
    // - // 1D Graphics
    
    TCanvas     *fc_1DDataYieldRaw     =   new TCanvas("","");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    fh_Raw_1D   ->  SetMarkerColor(kRed);
    fh_Raw_1D   ->  SetMarkerStyle(25);
    fh_Raw_1D   ->  SetMarkerSize(1);
    fh_Raw_1D   ->  SetTitle(Form("Raw yield compared to MC Truth"));
    fh_Raw_1D   ->  GetXaxis()  ->  SetTitle("p_{T}#phi (Gev/c)");
    fh_Raw_1D   ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{2}}{dp_{T}dy} (GeV/c)^{-1}");
    fh_Raw_1D   ->  Draw("EP");
    fh_Rec_1D   ->  SetLineColor(kBlue+3);
    fh_Rec_1D   ->  Draw("HIST L SAME");
    
    fLegend3    ->  AddEntry(fh_Raw_1D,"Raw measured yield","EP");
    fLegend3    ->  AddEntry(fh_Rec_1D,"Pythia6","L");
    fLegend3    ->  Draw("SAME");
    
    TH1F   *fh_Rw2_1D   =   new TH1F    (*fh_Raw_1D);
    fh_Rw2_1D   ->  Divide(fh_Efficiency_1D_1Dbin);
    
    TCanvas     *fc_1DDataYieldRw2     =   new TCanvas("","");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    fh_Rw2_1D   ->  SetMarkerColor(kRed);
    fh_Rw2_1D   ->  SetMarkerStyle(25);
    fh_Rw2_1D   ->  SetMarkerSize(1);
    fh_Rw2_1D   ->  SetTitle(Form("Raw yield, corrected for efficiency, compared to MC Truth"));
    fh_Rw2_1D   ->  GetXaxis()  ->  SetTitle("p_{T}#phi (Gev/c)");
    fh_Rw2_1D   ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{2}}{dp_{T}dy} (GeV/c)^{-1}");
    fh_Rw2_1D   ->  Draw("EP");
    fh_Gen_1D   ->  SetLineColor(kBlue+3);
    fh_Gen_1D   ->  Draw("HIST L SAME");
    
    fLegend4    ->  AddEntry(fh_Rw2_1D,"Raw / #varepsilon","EP");
    fLegend4    ->  AddEntry(fh_Gen_1D,"MC Truth","L");
    fLegend4    ->  Draw("SAME");
    
    TH1F   *fh_Res_1D   =   new TH1F    (*fh_Rw2_1D);
    fh_Res_1D   ->  Scale(1./kBranchingRatio__);
    
    TCanvas     *fc_1DDataYieldRes     =   new TCanvas("","");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    fh_Res_1D   ->  SetMarkerColor(kRed);
    fh_Res_1D   ->  SetMarkerStyle(25);
    fh_Res_1D   ->  SetMarkerSize(1);
    fh_Res_1D   ->  SetTitle(Form("Final measured yield MC Truth"));
    fh_Res_1D   ->  GetXaxis()  ->  SetTitle("p_{T}#phi (Gev/c)");
    fh_Res_1D   ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{2}}{dp_{T}dy} (GeV/c)^{-1}");
    fh_Res_1D   ->  Draw("EP");
    fh_Tru_1D   ->  SetLineColor(kBlue+3);
    fh_Tru_1D   ->  Draw("HIST L SAME");
    
    fLegend5    ->  AddEntry(fh_Rw2_1D,"Final Result","EP");
    fLegend5    ->  AddEntry(fh_Gen_1D,"Pythia6","L");
    fLegend5    ->  Draw("SAME");
    
    // - // 2D Graphics
    
    TCanvas     *fc_1DDataYieldRaw_XSliced  =   new TCanvas("","",855,1440);
    fc_1DDataYieldRaw_XSliced->Divide(2,5);
    TCanvas     *fc_1DDataYieldRaw_YSliced  =   new TCanvas("","",855,1440);
    fc_1DDataYieldRaw_YSliced->Divide(2,5);
    
    for ( int iChk = 2; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        fc_1DDataYieldRaw_XSliced->cd(iChk-1);
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        hUtility1X  =   (fh_Raw_2D->ProjectionX(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1X  ->  SetMarkerColor(kRed);
        hUtility1X  ->  SetMarkerStyle(25);
        hUtility1X  ->  SetMarkerSize(1);
        hUtility1X  ->  SetTitle(Form("Raw measured yield for %.1f < p_{T}#phi_{2} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hUtility1X  ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{1} (Gev/c)");
        hUtility1X  ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{3}}{dp_{T}#phi_{1}dp_{T}#phi_{2}dy} (GeV/c)^{-1}");
        hUtility1X  ->  Draw("EP SAME");
        hUtility2X  =   (fh_Rec_2D->ProjectionX(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2X  ->  SetLineColor(kBlue+3);
        hUtility2X  ->  Draw("HIST L SAME");
        if ( iChk == 2 )
        {
            fLegend6    ->  SetTextSize(0.1);
            fLegend6    ->  AddEntry(hUtility1X,"Raw","EP");
            fLegend6    ->  AddEntry(hUtility2X,"MC","L");
            fLegend6    ->  Draw("SAME");
        }
        
        // Y-Projection
        fc_1DDataYieldRaw_YSliced->cd(iChk-1);
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        hUtility1Y  =   (fh_Raw_2D->ProjectionY(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1Y  ->  SetMarkerColor(kRed);
        hUtility1Y  ->  SetMarkerStyle(25);
        hUtility1Y  ->  SetMarkerSize(1);
        hUtility1Y  ->  SetTitle(Form("Raw measured yield for %.1f < p_{T}#phi_{1} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hUtility1Y  ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{2} (Gev/c)");
        hUtility1Y  ->  GetYaxis()  ->  SetTitle("1/N_{events}#timesdN^{3}/dp_{T}#phi_{1}dp_{T}#phi_{2}dy (GeV/c)^{-1}");
        hUtility1Y  ->  Draw("EP SAME");
        hUtility2Y  =   (fh_Rec_2D->ProjectionY(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2Y  ->  SetLineColor(kBlue+3);
        hUtility2Y  ->  Draw("HIST L SAME");
        if ( iChk == 2 )
        {
            fLegend7    ->  SetTextSize(0.1);
            fLegend7    ->  AddEntry(hUtility1Y,"Raw","EP");
            fLegend7    ->  AddEntry(hUtility2Y,"MC","L");
            fLegend7    ->  Draw("SAME");
        }
    }
    
    fc_1DDataYieldRaw_XSliced->Write();
    fc_1DDataYieldRaw_XSliced->SaveAs("./graphs/2DDataYieldRaw_XSliced.pdf");
    fc_1DDataYieldRaw_XSliced->SaveAs("./graphs/2DDataYieldRaw_XSliced.png");
    
    fc_1DDataYieldRaw_YSliced->Write();
    fc_1DDataYieldRaw_YSliced->SaveAs("./graphs/2DDataYieldRaw_YSliced.pdf");
    fc_1DDataYieldRaw_YSliced->SaveAs("./graphs/2DDataYieldRaw_YSliced.png");
    
    TH2F   *fh_Rw2_2D   =   new TH2F    (*fh_Raw_2D);
    fh_Rw2_2D   ->  Divide(fh_Efficiency_2D_1Dsqr);
    
    for ( int iChk = 2; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        fc_1DDataYieldRaw_XSliced->cd(iChk-1);
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        hUtility1X  =   (fh_Rw2_2D->ProjectionX(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1X  ->  SetMarkerColor(kRed);
        hUtility1X  ->  SetMarkerStyle(25);
        hUtility1X  ->  SetMarkerSize(1);
        hUtility1X  ->  SetTitle(Form("Raw measured yield for %.1f < p_{T}#phi_{2} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hUtility1X  ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{1} (Gev/c)");
        hUtility1X  ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{3}}{dp_{T}#phi_{1}dp_{T}#phi_{2}dy} (GeV/c)^{-1}");
        hUtility1X  ->  Draw("EP SAME");
        hUtility2X  =   (fh_Gen_2D->ProjectionX(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2X  ->  SetLineColor(kBlue+3);
        hUtility2X  ->  Draw("HIST L SAME");
        if ( iChk == 2 )
        {
            fLegend6    ->  Clear();
            fLegend6    ->  SetTextSize(0.1);
            fLegend6    ->  AddEntry(hUtility1X,"Raw/#varepsilon","EP");
            fLegend6    ->  AddEntry(hUtility2X,"MC","L");
            fLegend6    ->  Draw("SAME");
        }
        
        // Y-Projection
        fc_1DDataYieldRaw_YSliced->cd(iChk-1);
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        hUtility1Y  =   (fh_Rw2_2D->ProjectionY(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1Y  ->  SetMarkerColor(kRed);
        hUtility1Y  ->  SetMarkerStyle(25);
        hUtility1Y  ->  SetMarkerSize(1);
        hUtility1Y  ->  SetTitle(Form("Raw measured yield for %.1f < p_{T}#phi_{1} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hUtility1Y  ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{2} (Gev/c)");
        hUtility1Y  ->  GetYaxis()  ->  SetTitle("1/N_{events}#timesdN^{3}/dp_{T}#phi_{1}dp_{T}#phi_{2}dy (GeV/c)^{-1}");
        hUtility1Y  ->  Draw("EP SAME");
        hUtility2Y  =   (fh_Gen_2D->ProjectionY(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2Y  ->  SetLineColor(kBlue+3);
        hUtility2Y  ->  Draw("HIST L SAME");
        if ( iChk == 2 )
        {
            fLegend7    ->  Clear();
            fLegend7    ->  SetTextSize(0.1);
            fLegend7    ->  AddEntry(hUtility1Y,"Raw/#varepsilon","EP");
            fLegend7    ->  AddEntry(hUtility2Y,"MC","L");
            fLegend7    ->  Draw("SAME");
        }
    }
    
    fc_1DDataYieldRaw_XSliced->Write();
    fc_1DDataYieldRaw_XSliced->SaveAs("./graphs/2DDataYieldRw2_XSliced.pdf");
    fc_1DDataYieldRaw_XSliced->SaveAs("./graphs/2DDataYieldRw2_XSliced.png");
    
    fc_1DDataYieldRaw_YSliced->Write();
    fc_1DDataYieldRaw_YSliced->SaveAs("./graphs/2DDataYieldRw2_YSliced.pdf");
    fc_1DDataYieldRaw_YSliced->SaveAs("./graphs/2DDataYieldRw2_YSliced.png");
    
    TH2F   *fh_Res_2D   =   new TH2F    (*fh_Rw2_2D);
    fh_Res_2D   ->  Scale(1./(kBranchingRatio__*kBranchingRatio__));
    
    for ( int iChk = 2; iChk < nBinPT2D; iChk++ )
    {
        // X-Projection
        fc_1DDataYieldRaw_XSliced->cd(iChk-1);
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        hUtility1X  =   (fh_Res_2D->ProjectionX(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1X  ->  SetMarkerColor(kRed);
        hUtility1X  ->  SetMarkerStyle(25);
        hUtility1X  ->  SetMarkerSize(1);
        hUtility1X  ->  SetTitle(Form("Raw measured yield for %.1f < p_{T}#phi_{2} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hUtility1X  ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{1} (Gev/c)");
        hUtility1X  ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{3}}{dp_{T}#phi_{1}dp_{T}#phi_{2}dy} (GeV/c)^{-1}");
        hUtility1X  ->  Draw("EP SAME");
        hUtility2X  =   (fh_Tru_2D->ProjectionX(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2X  ->  SetLineColor(kBlue+3);
        hUtility2X  ->  Draw("HIST L SAME");
        if ( iChk == 2 )
        {
            fLegend6    ->  Clear();
            fLegend6    ->  SetTextSize(0.1);
            fLegend6    ->  AddEntry(hUtility1X,"Res","EP");
            fLegend6    ->  AddEntry(hUtility2X,"MC","L");
            fLegend6    ->  Draw("SAME");
        }
        
        // Y-Projection
        fc_1DDataYieldRaw_YSliced->cd(iChk-1);
        gStyle->SetOptStat(0);
        gPad->SetLogy();
        hUtility1Y  =   (fh_Res_2D->ProjectionY(Form("2DEfficiency_XSlice_%i",iChk),iChk+1,iChk+1));
        hUtility1Y  ->  SetMarkerColor(kRed);
        hUtility1Y  ->  SetMarkerStyle(25);
        hUtility1Y  ->  SetMarkerSize(1);
        hUtility1Y  ->  SetTitle(Form("Raw measured yield for %.1f < p_{T}#phi_{1} < %.1f",fArrPT2D[iChk],fArrPT2D[iChk+1]));
        hUtility1Y  ->  GetXaxis()  ->  SetTitle("p_{T}#phi_{2} (Gev/c)");
        hUtility1Y  ->  GetYaxis()  ->  SetTitle("1/N_{events}#timesdN^{3}/dp_{T}#phi_{1}dp_{T}#phi_{2}dy (GeV/c)^{-1}");
        hUtility1Y  ->  Draw("EP SAME");
        hUtility2Y  =   (fh_Tru_2D->ProjectionY(Form("2DEfficiency_XSlice_%i_1Dx",iChk),iChk+1,iChk+1));
        hUtility2Y  ->  SetLineColor(kBlue+3);
        hUtility2Y  ->  Draw("HIST L SAME");
        if ( iChk == 2 )
        {
            fLegend7    ->  Clear();
            fLegend7    ->  SetTextSize(0.1);
            fLegend7    ->  AddEntry(hUtility1Y,"Res","EP");
            fLegend7    ->  AddEntry(hUtility2Y,"MC","L");
            fLegend7    ->  Draw("SAME");
        }
    }
    
    fc_1DDataYieldRaw_XSliced->Write();
    fc_1DDataYieldRaw_XSliced->SaveAs("./graphs/2DDataYieldRes_XSliced.pdf");
    fc_1DDataYieldRaw_XSliced->SaveAs("./graphs/2DDataYieldRes_XSliced.png");
    
    fc_1DDataYieldRaw_YSliced->Write();
    fc_1DDataYieldRaw_YSliced->SaveAs("./graphs/2DDataYieldRes_YSliced.pdf");
    fc_1DDataYieldRaw_YSliced->SaveAs("./graphs/2DDataYieldRes_YSliced.png");
    
    // Saving Yields Graphics
    fc_1DDataYieldRaw->Write();
    fc_1DDataYieldRaw->SaveAs("./graphs/1DDataYieldRaw.pdf");
    fc_1DDataYieldRaw->SaveAs("./graphs/1DDataYieldRaw.png");
    
    fc_1DDataYieldRw2->Write();
    fc_1DDataYieldRw2->SaveAs("./graphs/1DDataYieldRw2.pdf");
    fc_1DDataYieldRw2->SaveAs("./graphs/1DDataYieldRw2.png");
       
    fc_1DDataYieldRes->Write();
    fc_1DDataYieldRes->SaveAs("./graphs/1DDataYieldRes.pdf");
    fc_1DDataYieldRes->SaveAs("./graphs/1DDataYieldRes.png");
    
    hUtility1X  ->  Clear();
    hUtility1Y  ->  Clear();
    hUtility2X  ->  Clear();
    hUtility2Y  ->  Clear();
    
    TCanvas     *fc_2D_Raw_Rec_Overlap  =   new TCanvas("fc_2D_Raw_Rec_Overlap","fc_2D_Raw_Rec_Overlap",800,600);
    gPad->SetLogz();
    gStyle->SetOptStat(0);
    fh_Raw_2D   ->SetLineColor(kRed);
    fh_Raw_2D   ->SetMaximum(1.e-3);
    fh_Raw_2D   ->SetMinimum(1.e-8);
    fh_Raw_2D   ->GetXaxis()->SetTitleOffset(2.);
    fh_Raw_2D   ->GetYaxis()->SetTitleOffset(2.);
    fh_Raw_2D   ->Draw("lego");
    fh_Rec_2D   ->SetLineColor(kBlue);
    fh_Rec_2D   ->SetMaximum(1.e-3);
    fh_Rec_2D   ->SetMinimum(1.e-8);
    fh_Rec_2D   ->Draw("lego same");
    TLegend    *fLegendRawRec   =   new TLegend();
    fLegendRawRec->AddEntry(fh_Raw_2D,"Raw Yield Result",   "EP");
    fLegendRawRec->AddEntry(fh_Rec_2D,"Montecarlo Truth",   "EP");
    fLegendRawRec->Draw("same");
    
    fc_2D_Raw_Rec_Overlap->Write();
    fc_2D_Raw_Rec_Overlap->SaveAs("./graphs/fc_2D_Raw_Rec_Overlap.pdf");
    fc_2D_Raw_Rec_Overlap->SaveAs("./graphs/fc_2D_Raw_Rec_Overlap.png");
    delete fc_2D_Raw_Rec_Overlap;
    
    TCanvas     *fc_2D_Raw_Re2_Overlap  =   new TCanvas("fc_2D_Raw_Re2_Overlap","fc_2D_Raw_Re2_Overlap",800,600);
    gPad->SetLogz();
    gStyle->SetOptStat(0);
    fh_Rw2_2D   ->SetLineColor(kRed);
    fh_Rw2_2D   ->SetMaximum(1.e-3);
    fh_Rw2_2D   ->SetMinimum(1.e-8);
    fh_Rw2_2D   ->GetXaxis()->SetTitleOffset(2.);
    fh_Rw2_2D   ->GetYaxis()->SetTitleOffset(2.);
    fh_Rw2_2D   ->Draw("lego");
    fh_Gen_2D   ->SetLineColor(kBlue);
    fh_Gen_2D   ->SetMaximum(1.e-3);
    fh_Gen_2D   ->SetMinimum(1.e-8);
    fh_Gen_2D   ->Draw("lego same");
    
    fc_2D_Raw_Re2_Overlap->Write();
    fc_2D_Raw_Re2_Overlap->SaveAs("./graphs/fc_2D_Raw_Re2_Overlap.pdf");
    fc_2D_Raw_Re2_Overlap->SaveAs("./graphs/fc_2D_Raw_Re2_Overlap.png");
    delete fc_2D_Raw_Re2_Overlap;
    
    TCanvas     *fc_2D_Res_Tru_Overlap  =   new TCanvas("fc_2D_Res_Tru_Overlap","fc_2D_Res_Tru_Overlap",800,600);
    gPad->SetLogz();
    gStyle->SetOptStat(0);
    fh_Res_2D   ->SetLineColor(kRed);
    fh_Res_2D   ->SetMaximum(1.e-3);
    fh_Res_2D   ->SetMinimum(1.e-8);
    fh_Res_2D   ->GetXaxis()->SetTitleOffset(2.);
    fh_Res_2D   ->GetYaxis()->SetTitleOffset(2.);
    fh_Res_2D   ->Draw("lego");
    fh_Tru_2D   ->SetLineColor(kBlue);
    fh_Tru_2D   ->SetMaximum(1.e-3);
    fh_Tru_2D   ->SetMinimum(1.e-8);
    fh_Tru_2D   ->Draw("lego same");
    TLegend    *fLegendResTru   =   new TLegend();
    fLegendResTru->AddEntry(fh_Res_2D,"Final Yield Result", "EP");
    fLegendResTru->AddEntry(fh_Tru_2D,"Montecarlo Truth",   "EP");
    fLegendResTru->Draw("same");
    
    fc_2D_Res_Tru_Overlap->Write();
    fc_2D_Res_Tru_Overlap->SaveAs("./graphs/fc_2D_Res_Tru_Overlap.pdf");
    
    
    /*
    // 2D PT Eval
    Double_t***     fResults2D =   new Double_t**  [nBinPT2D+1];
    for ( Int_t iTer = 0; iTer < nBinPT2D+1; iTer++ )
    {
        fResults2D[iTer]        =   new Double_t*   [2];
        fResults2D[iTer][0]     =   new Double_t    [6];
        fResults2D[iTer][1]     =   new Double_t    [6];
    }
    
    fResults2D = ExtrapolateVl( fh_Res_2D );
    for ( int iChk = 0; iChk < nBinPT2D; iChk++ )
    {
        hUtility1X->SetBinContent(iChk,fh_Tru_2D->ProfileX(Form("x_%i",iChk),iChk+1,iChk+1)->Integral("width"));
        hUtility1Y->SetBinContent(iChk,fh_Tru_2D->ProfileY(Form("y_%i",iChk),iChk+1,iChk+1)->Integral("width"));
        
    }
    for ( int iChk = 2; iChk < nBinPT2D; iChk++ )
    {
        hUtility2X->SetBinContent   (iChk,fResults2D[iChk][0][0]);
        hUtility2Y->SetBinContent   (iChk,fResults2D[iChk][1][0]);
        hUtility2X->SetBinError     (iChk,fResults2D[iChk][0][1]);
        hUtility2Y->SetBinError     (iChk,fResults2D[iChk][1][1]);
    }
    
    TCanvas     *fc_1DDataYieldReX     =   new TCanvas("","");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    hUtility2X  ->  SetMarkerColor(kRed);
    hUtility2X  ->  SetMarkerStyle(25);
    hUtility2X  ->  SetMarkerSize(1);
    hUtility2X  ->  SetTitle(Form("Raw yield compared to MC Truth"));
    hUtility2X  ->  GetXaxis()  ->  SetTitle("p_{T}#phi (Gev/c)");
    hUtility2X  ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{2}}{dp_{T}dy} (GeV/c)^{-1}");
    hUtility2X  ->  Draw("EP");
    hUtility1X  ->  SetLineColor(kBlue+3);
    hUtility1X  ->  Draw("HIST L SAME");
    
    fLegend3    ->  Clear();
    fLegend3    ->  AddEntry(fh_Res_FX,"Raw measured yield","EP");
    fLegend3    ->  AddEntry(hUtility1X,"MC Truth","L");
    fLegend3    ->  Draw("SAME");
    
    TCanvas     *fc_1DDataYieldReY      =   new TCanvas("","");
    gStyle->SetOptStat(0);
    gPad->SetLogy();
    fh_Res_FY  ->  SetMarkerColor(kRed);
    fh_Res_FY  ->  SetMarkerStyle(25);
    fh_Res_FY  ->  SetMarkerSize(1);
    fh_Res_FY  ->  SetTitle(Form("Raw yield compared to MC Truth"));
    fh_Res_FY  ->  GetXaxis()  ->  SetTitle("p_{T}#phi (Gev/c)");
    fh_Res_FY  ->  GetYaxis()  ->  SetTitle("#frac{1}{N_{events}}#times#frac{dN^{2}}{dp_{T}dy} (GeV/c)^{-1}");
    fh_Res_FY  ->  Draw("EP");
    hUtility1Y  ->  SetLineColor(kBlue+3);
    hUtility1Y  ->  Draw("HIST L SAME");
    
    fLegend4    ->  Clear();
    fLegend4    ->  AddEntry(fh_Res_FY,"Raw measured yield","EP");
    fLegend4    ->  AddEntry(hUtility1Y,"MC Truth","L");
    fLegend4    ->  Draw("SAME");
    
    fc_1DDataYieldReX->Write();
    fc_1DDataYieldReX->SaveAs("./graphs/1DDataYieldReX.pdf");
    fc_1DDataYieldReX->SaveAs("./graphs/1DDataYieldReX.png");
    
    fc_1DDataYieldReY->Write();
    fc_1DDataYieldReY->SaveAs("./graphs/1DDataYieldReY.pdf");
    fc_1DDataYieldReY->SaveAs("./graphs/1DDataYieldReY.png");
    */
    
    
    TFile*  insFile_DT  =   new TFile   ("./result/LHC10.root");
    TFile*  insFile_MC  =   new TFile   ("./result/LHC14j4.root");
    
    TList * hHistoPid = (TList*)(insFile_DT->Get("MyOutputContainerListUTL"));
    hName = Form("fHistTOFPID0");
    TH2F*hTOFPID     =   (TH2F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTOFPID1");
    TH2F*hTOFPID_1   =   (TH2F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTOFPID2");
    TH2F*hTOFPID_2   =   (TH2F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTPCPID0");
    TH2F*hTPCPID     =   (TH2F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTPCPID1");
    TH2F*hTPCPID_1   =   (TH2F*)(hHistoPid->FindObject(hName));
    hName = Form("fHistTPCPID2");
    TH2F*hTPCPID_2   =   (TH2F*)(hHistoPid->FindObject(hName));
    
    TList * hHistoPid_MC = (TList*)(insFile_MC->Get("MyOutputContainerListUTL"));
    hName = Form("fHistTOFPID0");
    TH2F*hTOFPID_MC     =   (TH2F*)(hHistoPid_MC->FindObject(hName));
    hName = Form("fHistTOFPID1");
    TH2F*hTOFPID_1_MC   =   (TH2F*)(hHistoPid_MC->FindObject(hName));
    hName = Form("fHistTOFPID2");
    TH2F*hTOFPID_2_MC   =   (TH2F*)(hHistoPid_MC->FindObject(hName));
    hName = Form("fHistTPCPID0");
    TH2F*hTPCPID_MC     =   (TH2F*)(hHistoPid_MC->FindObject(hName));
    hName = Form("fHistTPCPID1");
    TH2F*hTPCPID_1_MC   =   (TH2F*)(hHistoPid_MC->FindObject(hName));
    hName = Form("fHistTPCPID2");
    TH2F*hTPCPID_2_MC   =   (TH2F*)(hHistoPid_MC->FindObject(hName));
    
    TFile * hout = new TFile ("quick.root","recreate");
    
    hTOFPID  ->Write();
    hTOFPID_MC  ->Write();
    hTPCPID  ->Write();
    hTPCPID_MC  ->Write();
    
    TH1F   *DataTOFcount  =   new TH1F ("DataTOFcount","DataTOFcount",100, 0, 4);
    TH1F   *DataTPCcount  =   new TH1F ("DataTPCcount","DataTPCcount",100, 0, 4);
    TH1F   *MC__TOFcount  =   new TH1F ("MC__TOFcount","MC__TOFcount",100, 0, 4);
    TH1F   *MC__TPCcount  =   new TH1F ("MC__TPCcount","MC__TPCcount",100, 0, 4);
    
    for ( Int_t iBin = 1; iBin <= 100; iBin++ )
    {
        double error;
        if ( iBin <= 8 ) continue;
        auto DataTOFslice = hTOFPID->ProjectionX(Form("DT_TOF_%i",iBin),iBin,iBin);
        DataTOFcount->SetBinContent (iBin, DataTOFslice->IntegralAndError(-1,200,error));
        DataTOFcount->SetBinError   (iBin, error);
        DataTOFslice->Write();
        auto DataTPCslice = hTPCPID->ProjectionX(Form("DT_TPC_%i",iBin),iBin,iBin);
        DataTPCcount->SetBinContent (iBin, DataTPCslice->IntegralAndError(-1,200,error));
        DataTPCcount->SetBinError   (iBin, error);
        DataTPCslice->Write();
        auto DataTOFslice_MC = hTOFPID_MC->ProjectionX(Form("MC_TOF_%i",iBin),iBin,iBin);
        MC__TOFcount->SetBinContent (iBin, DataTOFslice_MC->IntegralAndError(-1,200,error));
        MC__TOFcount->SetBinError   (iBin, error);
        DataTOFslice_MC->Write();
        auto DataTPCslice_MC = hTPCPID_MC->ProjectionX(Form("MC_TPC_%i",iBin),iBin,iBin);
        MC__TPCcount->SetBinContent (iBin, DataTPCslice_MC->IntegralAndError(-1,200,error));
        MC__TPCcount->SetBinError   (iBin, error);
        DataTPCslice_MC->Write();
    }
    
    TH1F   *DataEfficiency  =   new TH1F ("e_Data","e_Data",100, 0, 4);
    TH1F   *MC__Efficiency  =   new TH1F ("e_MC","e_MC",100, 0, 4);
    DataEfficiency->Divide(DataTOFcount,DataTPCcount,1.,1.,"B");
    MC__Efficiency->Divide(MC__TOFcount,MC__TPCcount,1.,1.,"B");
    
    TH1F   *hRatioEfficiency= new TH1F ("hRatioEfficiency","hRatioEfficiency",100, 0, 4);
    hRatioEfficiency->Divide(DataEfficiency,MC__Efficiency);
    
    DataTOFcount->Write();
    DataTPCcount->Write();
    MC__TOFcount->Write();
    MC__TPCcount->Write();
    DataEfficiency->Write();
    MC__Efficiency->Write();
    hRatioEfficiency->Write();
    
    hout->Close();
    insFile_DT->Close();
    insFile_MC->Close();
    
    //---------------------//
    //  Wrapping up to end //--------------------------------------------------------------------------------
    //---------------------//
    
    delete  fc_Efficiency;
    delete  fc_Efficiency_XSliced;
    delete  fc_Efficiency_YSliced;
    delete  fc_1DDataYieldRaw;
    delete  fc_1DDataYieldRw2;
    delete  fc_1DDataYieldRes;
    delete  fc_1DDataYieldRaw_XSliced;
    delete  fc_1DDataYieldRaw_YSliced;
    //delete  fc_1DDataYieldReX;
    //delete  fc_1DDataYieldReY;
    ot_Graphical            ->Close();
    in_Efficiency_Pythia8   ->Close();
    
    return;
}
