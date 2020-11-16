#include "../inc/AliAnalysisPhiPair.h"
// !TODO: [INFO] About trees in input

void Anls_PreProcessing_Data ( string fFileName = "" )
{
    if ( fFileName == "" )
    {
        cout << "[WARNING] Must Specify an input root file" << endl;
        cout << "[INFO] Usage Anls_PreProcessing_Data(\"Root_file_name.root\")" << endl;
        return;
    }
    
    //Retrieving Event data
    TFile *insFileDT        =   new TFile   (fFileName.c_str());
    
    //Retrieving Event data TTree
    TTree   *TPhiCandidate  =   (TTree*)insFileDT->Get(fPhiCandidate_Tree);
    TTree   *TKaonCandidate =   (TTree*)insFileDT->Get(fKaonCandidate_Tree);
    
    if ( !TPhiCandidate && !TKaonCandidate )
    {
        cout << "Input Data Tree not found!" << endl;
        return;
    }
    /*
    if ( !TPhiCandidate )
    {
        cout << "[INFO] " << endl;
    }
    */
    
    // Define tree data structures
    Struct_PhiCandidate     evPhiCandidate;
    Struct_PhiEfficiency    evPhiEfficiency;
    Struct_KaonCandidate    evKaonCandidate;
    Struct_KaonEfficiency   evKaonEfficiency;
    
    TPhiCandidate-> SetBranchAddress    ("Multiplicity",    &evPhiCandidate.Multiplicity);
    TPhiCandidate-> SetBranchAddress    ("nPhi",            &evPhiCandidate.nPhi);
    TPhiCandidate-> SetBranchAddress    ("Px",              &evPhiCandidate.Px);
    TPhiCandidate-> SetBranchAddress    ("Py",              &evPhiCandidate.Py);
    TPhiCandidate-> SetBranchAddress    ("Pz",              &evPhiCandidate.Pz);
    TPhiCandidate-> SetBranchAddress    ("InvMass",         &evPhiCandidate.InvMass);
    TPhiCandidate-> SetBranchAddress    ("iKaon",           &evPhiCandidate.iKaon);
    TPhiCandidate-> SetBranchAddress    ("jKaon",           &evPhiCandidate.jKaon);
    
    TKaonCandidate-> SetBranchAddress   ("Multiplicity",    &evKaonCandidate.Multiplicity);
    TKaonCandidate-> SetBranchAddress   ("nKaon",           &evKaonCandidate.nKaon);
    TKaonCandidate-> SetBranchAddress   ("Px",              &evKaonCandidate.Px);
    TKaonCandidate-> SetBranchAddress   ("Py",              &evKaonCandidate.Py);
    TKaonCandidate-> SetBranchAddress   ("Pz",              &evKaonCandidate.Pz);
    TKaonCandidate-> SetBranchAddress   ("Charge",          &evKaonCandidate.Charge);
    TKaonCandidate-> SetBranchAddress   ("TOFSigma",        &evKaonCandidate.SigmaTOF);
    TKaonCandidate-> SetBranchAddress   ("TPCSigma",        &evKaonCandidate.SigmaTPC);
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    Int_t       U_AccCand[1024];
    Int_t       U_nAccept;
    
    // Creating the histograms-------------------------------------------------------------------------------
    // 1D
    TH1F      **hREC_1D_in_PT               = new TH1F     *[nBinPT1D];
    TH1F       *hREC_1D_in_Rap;
    TH1F       *hREF_1D_in_Rap;
    
    // 2D
    TH1F      **hREC_1D_in_PT_2D_bin        = new TH1F     *[nBinPT2D];
    TH2F     ***hREC_2D_in_PT               = new TH2F    **[nBinPT2D];
    
    hName = "hREC_1D_in_Rap";
    hTitle= "Rapidity difference for #phi meson candidates";
    hREC_1D_in_Rap  =   new TH1F (hName,hTitle,100,-1.,1.);
    
    hName = "hREF_1D_in_Rap";
    hTitle= "Rapidity distribution for #phi meson candidates";
    hREF_1D_in_Rap  =   new TH1F (hName,hTitle,100,-.5,.5);
    
    for ( Int_t iHisto = 0; iHisto < nBinPT1D; iHisto++ )
    {
        hName = Form("hREC_1D_in_PT_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT1D[iHisto],fArrPT1D[iHisto+1]);
        hREC_1D_in_PT[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT[iHisto],"IM 1D");
    }
    
    for ( Int_t iHisto = 0; iHisto < nBinPT2D; iHisto++ )
    {
        hName = Form("hREC_1D_in_PT_2D_bin_%i",iHisto);
        hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1]);
        hREC_1D_in_PT_2D_bin[iHisto]   = new TH1F (hName,hTitle,nBinIM1D,fArrIM1D);
        SetAxis(hREC_1D_in_PT_2D_bin[iHisto],"IM 1D");
        
        hREC_2D_in_PT[iHisto]    = new TH2F *    [nBinPT2D];
        
        for ( Int_t jHisto = 0; jHisto < nBinPT2D; jHisto++ )
        {
            hName = Form("hREC_2D_in_PT_%i_%i",iHisto,jHisto);
            hTitle= Form("m_{K^{+}K^{-}} in p_{T} range [%.1f-%.1f] GeV/c and [%.1f-%.1f] GeV/c",fArrPT2D[iHisto],fArrPT2D[iHisto+1],fArrPT2D[jHisto],fArrPT2D[jHisto+1]);
            hREC_2D_in_PT[iHisto][jHisto]    = new TH2F (hName,hTitle,nBinIM2D,fArrIM2D,nBinIM2D,fArrIM2D);
            SetAxis(hREC_2D_in_PT[iHisto][jHisto],"IM 2D");
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
    
    // Evaluating entries and saving them for later
    Int_t nEntries = TPhiCandidate->GetEntries();
    hUtlEntry     ->SetBinContent(1,nEntries);
    
    // Starting cycle
    for ( Int_t iEvent = 0; iEvent < nEntries; iEvent++ )
    {
        // Recovering events
        TPhiCandidate->GetEntry(iEvent);
        TLorentzVector  LPhi_candidate1,    LPhi_candidate2;
        U_nAccept = 0;
        
        for ( Int_t iPhi = 0; iPhi < evPhiCandidate.nPhi; iPhi++ )
        {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[iPhi],evPhiCandidate.Py[iPhi],evPhiCandidate.Pz[iPhi],evPhiCandidate.InvMass[iPhi]);
            if ( !fRapidityCut ( LPhi_candidate1.Rapidity() ) ) continue;
            if ( !fTransverseMomCut ( LPhi_candidate1.Pt() ) ) continue;
            U_AccCand[U_nAccept] = iPhi;
            U_nAccept++;
        }
        for ( Int_t iPhi = 0; iPhi < U_nAccept; iPhi++ )
        {
            LPhi_candidate1.SetXYZM(evPhiCandidate.Px[U_AccCand[iPhi]],evPhiCandidate.Py[U_AccCand[iPhi]],evPhiCandidate.Pz[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[iPhi]]);
            
            hREC_1D_in_PT[fGetBinPT1D(LPhi_candidate1.Pt())]            ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]]);
            hREC_1D_in_PT_2D_bin[fGetBinPT2D(LPhi_candidate1.Pt())]     ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]]);
            hREF_1D_in_Rap                                              ->  Fill(LPhi_candidate1.Rapidity());
            for ( Int_t jPhi = 0; jPhi < U_nAccept; jPhi++ )
            {
                LPhi_candidate2.SetXYZM(evPhiCandidate.Px[U_AccCand[jPhi]],evPhiCandidate.Py[U_AccCand[jPhi]],evPhiCandidate.Pz[U_AccCand[jPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]]);
                
                // Only non overlapping couples of Kaons
                if ( evPhiCandidate.iKaon[U_AccCand[iPhi]] == evPhiCandidate.iKaon[U_AccCand[jPhi]] ) continue;
                if ( evPhiCandidate.jKaon[U_AccCand[iPhi]] == evPhiCandidate.jKaon[U_AccCand[jPhi]] ) continue;
                if ( evPhiCandidate.iKaon[U_AccCand[iPhi]] == evPhiCandidate.jKaon[U_AccCand[jPhi]] ) continue;
                if ( evPhiCandidate.jKaon[U_AccCand[iPhi]] == evPhiCandidate.iKaon[U_AccCand[jPhi]] ) continue;
                
                hREC_2D_in_PT[fGetBinPT2D(LPhi_candidate1.Pt())][fGetBinPT2D(LPhi_candidate2.Pt())] ->  Fill(evPhiCandidate.InvMass[U_AccCand[iPhi]],evPhiCandidate.InvMass[U_AccCand[jPhi]],0.5);
                
                hREC_1D_in_Rap->    Fill(LPhi_candidate1.Rapidity()-LPhi_candidate2.Rapidity(),0.5);
            }
        }
    }

    //--------------------------//
    //  Printing output objects //
    //--------------------------//
    
    TFile *outFile  =   new TFile   (fInvMasHist,"recreate");
    
    hUtlEntry->Write();
    hREF_1D_in_Rap->Write();
    hREC_1D_in_Rap->Write();
    for (int iHisto = 0; iHisto < nBinPT1D; iHisto++)
    {
        hREC_1D_in_PT[iHisto]   ->Write();
    }
    
    for (int iHisto = 0; iHisto < nBinPT2D; iHisto++)
    {
        hREC_1D_in_PT_2D_bin[iHisto]   ->Write();
        
        for (int jHisto = 0; jHisto < nBinPT2D; jHisto++)
        {
            hREC_2D_in_PT[iHisto][jHisto]->Write();
        }
    }
    
    outFile->Close();
    insFileDT->Close();
}
