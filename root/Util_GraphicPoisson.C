#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"

void Util_GraphicPoisson ()
{
    TFile   *in_Efficiency_Pythia6          =   new TFile   ("./result/LHC14j4.root");
    TFile   *in_Efficiency_Pythia8          =   new TFile   ("/Volumes/[HD][Nikolajal]_WD/Rubini/MCZGenerated.root");
    
    TTree *PTreePTr6    =   (TTree*)in_Efficiency_Pythia6->Get(fTreeTruName);
    TTree *PTreePTr8    =   (TTree*)in_Efficiency_Pythia8->Get(fTreeTruName);
    
    EVPHI               evPhi_Tr6;
    PTreePTr6  ->SetBranchAddress    ("nPhi",           &evPhi_Tr6.nPhi);
    PTreePTr6  ->SetBranchAddress    ("bEta",           &evPhi_Tr6.bEta);
    
    EVPHI               evPhi_Tr8;
    PTreePTr8  ->SetBranchAddress    ("evPhi.nPhi",           &evPhi_Tr8.nPhi);
    PTreePTr8  ->SetBranchAddress    ("evPhi.bEta",           &evPhi_Tr8.bEta);
    
    TH1D * hPythia6 = new TH1D ("hPythia6","hPythia6",11,-0.5,10.5);
    hPythia6->SetLineWidth(3);
    hPythia6->SetLineColor(kBlue);
    hPythia6->SetLineStyle(9);
    hPythia6->GetXaxis()->SetTitle("Number of #phi-meson produced in event");
    hPythia6->GetYaxis()->SetTitle("Frequency");
    
    TH1D * hPythia6_Y = new TH1D ("hPythia6_Y","hPythia6_Y",11,-0.5,10.5);
    hPythia6_Y->SetLineWidth(3);
    hPythia6_Y->SetLineColor(kBlue);
    hPythia6_Y->SetLineStyle(9);
    hPythia6_Y->GetXaxis()->SetTitle("N-tuples of #phi-meson produced in event");
    hPythia6_Y->GetYaxis()->SetTitle("Frequency");
    
    TH1D * hPythia8 = new TH1D ("hPythia8","hPythia8",11,-0.5,10.5);
    hPythia8->SetLineWidth(3);
    hPythia8->SetLineColor(kBlue+3);
    hPythia8->SetLineStyle(10);
    
    TH1D * hPythia8_Y = new TH1D ("hPythia8_Y","hPythia8_Y",11,-0.5,10.5);
    hPythia8_Y->SetLineWidth(3);
    hPythia8_Y->SetLineColor(kBlue+3);
    hPythia8_Y->SetLineStyle(10);
     
    TH1D * hPoisson = new TH1D ("hPoisson","hPoisson",11,-0.5,10.5);
    hPoisson->SetLineWidth(3);
    hPoisson->SetLineColor(kBlack);
    hPoisson->SetLineStyle(1);
    
    TH1D * hPoisson_Y = new TH1D ("hPoisson_Y","hPoisson_Y",11,-0.5,10.5);
    hPoisson_Y->SetLineWidth(3);
    hPoisson_Y->SetLineColor(kBlack);
    hPoisson_Y->SetLineStyle(1);
    
    TH1D * hData = new TH1D ("hData","hData",11,-0.5,10.5);
    SetGraphicStyle(hData,"DT");
    
    Int_t nEntrie6 = PTreePTr6->GetEntries();
    for ( Int_t iEvent = 0; iEvent < nEntrie6; iEvent++ )
    {
        PTreePTr6->GetEntry(iEvent);
        Int_t nPhiEta = 0;
        for (Int_t iPhi = 0; iPhi < evPhi_Tr6.nPhi; iPhi++ )
        {
            if ( !evPhi_Tr6.bEta[iPhi] ) continue;
            nPhiEta++;
        }
        hPythia6->Fill(nPhiEta);
    }
    
    hPythia6->Scale(1./nEntrie6);
    
    for ( Int_t iBin = 1; iBin <= 10; iBin++ )
    {
        for ( Int_t jBin = iBin; jBin < 10; jBin++ )
        {
            hPythia6_Y->Fill(iBin,TMath::Binomial(jBin,iBin)*hPythia6->GetBinContent(jBin+1));
        }
    }
    
    Int_t nEntrie8 = PTreePTr8->GetEntries();
    for ( Int_t iEvent = 0; iEvent < nEntrie8; iEvent++ )
    {
        PTreePTr8->GetEntry(iEvent);
        Int_t nPhiEta = 0;
        for (Int_t iPhi = 0; iPhi < evPhi_Tr8.nPhi; iPhi++ )
        {
            if ( !evPhi_Tr8.bEta[iPhi] ) continue;
            nPhiEta++;
        }
        hPythia8->Fill(nPhiEta);
    }
    
    hPythia8->Scale(1./nEntrie8);
    
    for ( Int_t iBin = 1; iBin <= 10; iBin++ )
    {
        for ( Int_t jBin = iBin; jBin < 10; jBin++ )
        {
            hPythia8_Y->Fill(iBin,TMath::Binomial(jBin,iBin)*hPythia8->GetBinContent(jBin+1));
        }
    }
    
    Int_t nEntrieP  = 1.e9;
    for ( Int_t iEvent = 0; iEvent < nEntrieP; iEvent++ )
    {
        hPoisson->Fill(fRandomGen->Poisson(0.0333));
    }

    hPoisson->Scale(1./nEntrieP);
    
    for ( Int_t iBin = 1; iBin <= 10; iBin++ )
    {
        for ( Int_t jBin = iBin; jBin < 10; jBin++ )
        {
            hPoisson_Y->Fill(iBin,TMath::Binomial(jBin,iBin)*hPoisson->GetBinContent(jBin+1));
        }
    }
    
    hData->SetBinContent(2,0.033);
    hData->SetBinContent(3,0.00144);
    hData->SetBinError(2,0.0037);
    hData->SetBinError(3,0.0003);
    
    TFile * fOutput = new TFile("./result/Pythia.root","recreate");
    hPythia6->Write();
    hPythia8->Write();
    hPoisson->Write();
    hPythia6_Y->Write();
    hPythia8_Y->Write();
    hPoisson_Y->Write();
    
    TCanvas *c1 = new TCanvas();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    
    TLegend * fLegend = new TLegend(0.6,0.75,0.9,.9);
    fLegend->AddEntry(hPythia6,"Pythia6","L");
    fLegend->AddEntry(hPythia8,"Pythia8","L");
    fLegend->AddEntry(hPoisson,"Poiss.Distr.","L");
    
    TLatex * fText = new TLatex();
    fText   ->DrawLatexNDC(0.6, 0.85, "|y| < 0.5");
    
    hPythia6->SetTitle("");
    hPythia6->GetXaxis()->SetRangeUser(0.,5.);
    hPythia6->Draw("HISTL");
    hPythia8->Draw("SAMEHISTL");
    hPoisson->Draw("SAMEHISTL");
    fLegend->Draw("SAME");
    
    c1->Write();
    c1->SaveAs("./graphs/Poiss.pdf");
    c1->SaveAs("./graphs/Poiss.png");
    
    TCanvas *c2 = new TCanvas();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    
    TLegend * fLegend_Y = new TLegend(0.6,0.75,0.9,.9);
    fLegend_Y->AddEntry(hPythia6_Y,"Pythia6","L");
    fLegend_Y->AddEntry(hPythia8_Y,"Pythia8","L");
    fLegend_Y->AddEntry(hPoisson_Y,"Poiss.Distr.","L");
    
    TLatex * fText_Y = new TLatex();
    fText_Y   ->DrawLatexNDC(0.6, 0.85, "|y| < 0.5");
    
    hPythia6_Y->SetTitle("");
    hPythia6_Y->GetXaxis()->SetRangeUser(1.,5.);
    hPythia6_Y->Draw("HISTL");
    hPythia8_Y->Draw("SAMEHISTL");
    hPoisson_Y->Draw("SAMEHISTL");
    hData->Draw("SAMEPE");
    fLegend_Y->Draw("SAME");
    
    c2->Write();
    c2->SaveAs("./graphs/Poiss2.pdf");
    c2->SaveAs("./graphs/Poiss2.png");
    
    in_Efficiency_Pythia6->Close();
    in_Efficiency_Pythia8->Close();
    fOutput->Close();
}
