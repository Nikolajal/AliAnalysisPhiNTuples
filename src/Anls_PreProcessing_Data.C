#include "../inc/SetValues.h"
// !TODO: Rebooting the selection process

void Anls_DataPreProcessing ( const char * fFileName )
{
    if ( !fFileName )
    {
        cout << "Must Specify an input root file" << endl;
        return;
    }
    
    //Retrieving Event data
    TFile *insFileDT     =   new TFile   (fFileName);
    
    //Retrieving Event data TTree
    TTree *PTreeKSig    =   (TTree*)insFileDT->Get(fTreeSigName);
    
    if ( !PTreeKSig )
    {
        cout << "Input Data Tree not found!" << endl;
        return;
    }
        
    // Define some simple data structures to Set Branch Addresses
    // Kaon +- Couples B + S
    EVKAONCOUPLE        evKaonSig;
    if ( bPythiaTest == false )
    {
        PTreeKSig  ->SetBranchAddress    ("nKaonCouple",   &evKaonSig.nKaonCouple);
        PTreeKSig  ->SetBranchAddress    ("iKaon",         &evKaonSig.iKaon);
        PTreeKSig  ->SetBranchAddress    ("jKaon",         &evKaonSig.jKaon);
        PTreeKSig  ->SetBranchAddress    ("bPhi",          &evKaonSig.bPhi);
        PTreeKSig  ->SetBranchAddress    ("bRec",          &evKaonSig.bRec);
        PTreeKSig  ->SetBranchAddress    ("bEta",          &evKaonSig.bEta);
        PTreeKSig  ->SetBranchAddress    ("InvMass",       &evKaonSig.InvMass);
        PTreeKSig  ->SetBranchAddress    ("pT",            &evKaonSig.pT);
    }
    else
    {
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.nKaonCouple",   &evKaonSig.nKaonCouple);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.iKaon",         &evKaonSig.iKaon);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.jKaon",         &evKaonSig.jKaon);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.bPhi",          &evKaonSig.bPhi);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.bRec",          &evKaonSig.bRec);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.bEta",          &evKaonSig.bEta);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.InvMass",       &evKaonSig.InvMass);
        PTreeKSig  ->SetBranchAddress    ("evKaonCouple.pT",            &evKaonSig.pT);
    }
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    Int_t       S_ArrpT[1024];
    
    // Creating the histograms-------------------------------------------------------------------------------
    // 1D
    TH1F **     hIM_1D_Rec_PT_B_S               = new TH1F *    [nBinPT1D];
    
    // 2D
    TH1F **     hIM_2D_Rec_PT_B_S               = new TH1F *    [nBinPT2D];
    TH2F ***    hIM_2D_Rec_PT_PT_BB_BS_SB_SS    = new TH2F **   [nBinPT2D];
    
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        hName = Form("hIM_1D_Rec_PT_B_S_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT1D[iHisto],fArrPT1D[iHisto+1]);
        hIM_1D_Rec_PT_B_S[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        hIM_1D_Rec_PT_B_S[iHisto]   ->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
        hIM_1D_Rec_PT_B_S[iHisto]   ->GetYaxis()->SetTitle("Entries");
        hIM_1D_Rec_PT_B_S[iHisto]   ->GetXaxis()->SetTitleOffset(1.15);
        hIM_1D_Rec_PT_B_S[iHisto]   ->GetYaxis()->SetTitleOffset(1.15);
    }
    
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        hName = Form("hIM_2D_Rec_PT_B_S_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        hIM_2D_Rec_PT_B_S[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        hIM_2D_Rec_PT_B_S[iHisto]   ->GetXaxis()->SetTitle("m_{K^{+}K^{-}} (GeV/c^{2})");
        hIM_2D_Rec_PT_B_S[iHisto]   ->GetYaxis()->SetTitle("Entries");
        hIM_2D_Rec_PT_B_S[iHisto]   ->GetXaxis()->SetTitleOffset(1.15);
        hIM_2D_Rec_PT_B_S[iHisto]   ->GetYaxis()->SetTitleOffset(1.15);
        
        hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto]    = new TH2F *    [nBinPT2D];
        
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hIM_2D_Rec_PT_PT_BB_BS_SB_SS_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto][jHisto]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto][jHisto]    ->GetXaxis()->SetTitle(Form("m_{K^{+}K^{-}} (GeV/c^{2}) in p_{T} [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1]));
            hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto][jHisto]    ->GetYaxis()->SetTitle(Form("m_{K^{+}K^{-}} (GeV/c^{2}) in p_{T} [%.1f-%.1f] GeV/c",fArrPT2D[jHisto],fArrPT2D[jHisto+1]));
            hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto][jHisto]    ->GetXaxis()->SetTitleOffset(1.5);
            hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto][jHisto]    ->GetYaxis()->SetTitleOffset(1.5);
        }
    }
    
    // Creating the Number of events Result Histogram------------------------------------------------------------------
    hName                       = "Entry_DT";
    hTitle                      = "Events in DT";
    TH1F *          hUtlEntry   = new TH1F (hName,hTitle,1,0.5,1.5);
    hUtlEntry                   ->GetXaxis()->SetTitle("");
    hUtlEntry                   ->GetYaxis()->SetTitle("Events");
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    // Evaluating entries
    Int_t nEntries = PTreeKSig->GetEntries();
    for ( Int_t iEvent = 0; iEvent < nEntries; iEvent++ )
    {
        // Recovering events
        PTreeKSig->GetEntry(iEvent);
        
        for ( Int_t iPhi = 0; iPhi < evKaonSig.nKaonCouple; iPhi++ )
        {
            S_ArrpT[iPhi] = -1;
            
            // Only Physically relevant couples of Kaons
            if ( !evKaonSig.bEta[iPhi] ) continue;
        
            // Only Phi couples of Kaons
            if ( !evKaonSig.bRec[iPhi] && bPythiaTest ) continue;
            
            // Determining pT
            S_ArrpT[iPhi] = fGetBinPT1D ( evKaonSig.pT[iPhi] );
            
            // Refusing non acceptable Couple
            if ( S_ArrpT[iPhi] == -1 ) continue;
            
            // Filling 1D Histogram
            hIM_1D_Rec_PT_B_S[S_ArrpT[iPhi]]   ->Fill(evKaonSig.InvMass[iPhi]);
            
            // Determining 2D Sig pT bin
            S_ArrpT[iPhi] = fGetBinPT2D ( evKaonSig.pT[iPhi] );
            
            // Skipping invalid pT
            if ( S_ArrpT[iPhi] != -1) hIM_2D_Rec_PT_B_S[S_ArrpT[iPhi]]   ->Fill(evKaonSig.InvMass[iPhi]);
            
        }
        for (Int_t iPhi = 0; iPhi < evKaonSig.nKaonCouple; iPhi++ )
        {
            // Only acceptable couples of Kaons
            if ( S_ArrpT[iPhi] == -1 ) continue;
            
            // Combinatorial mixing
            for (Int_t jPhi = 0; jPhi < evKaonSig.nKaonCouple; jPhi++ )
            {
                // Auto-correlation protection
                if ( iPhi == jPhi ) continue;

                // Only acceptable couples of Kaons
                if ( S_ArrpT[jPhi] == -1 || S_ArrpT[iPhi] == -1 ) continue;
                
                // Only non overlapping couples of Kaons
                if ( evKaonSig.iKaon[iPhi] == evKaonSig.iKaon[jPhi] ) continue;
                if ( evKaonSig.iKaon[iPhi] == evKaonSig.jKaon[jPhi] ) continue;
                if ( evKaonSig.jKaon[iPhi] == evKaonSig.iKaon[jPhi] ) continue;
                if ( evKaonSig.jKaon[iPhi] == evKaonSig.jKaon[jPhi] ) continue;
                
                hIM_2D_Rec_PT_PT_BB_BS_SB_SS[S_ArrpT[iPhi]][S_ArrpT[jPhi]]   ->Fill(evKaonSig.InvMass[iPhi],evKaonSig.InvMass[jPhi],0.5);
            }
        }
    }
    
    // Evaluating the Entries value
    nEntries      = PTreeKSig->GetEntries();
    hUtlEntry     ->SetBinContent(1,nEntries);
    
    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    
    TFile *outFile  =   new TFile   (fInvMasHist,"recreate");
    
    hUtlEntry->Write();
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hIM_1D_Rec_PT_B_S[iHisto]   ->Write();
    }
    
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hIM_2D_Rec_PT_B_S[iHisto]   ->Write();
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hIM_2D_Rec_PT_PT_BB_BS_SB_SS[iHisto][jHisto]->Write();
        }
    }
    
    outFile->Close();
    insFileDT->Close();
}
