#include "../inc/SetValues.h"
// !TODO: Rebooting the selection process

int FractionCount ()
{
    //Retrieving Event data
    TFile *insFile       =   new TFile   ("./result/LHC14j4.root");
    
    //Retrieving Event data TTree
    TTree *PTreeKSig    =   (TTree*)insFile->Get(fTreeSigName);
    TTree *PTreePTru    =   (TTree*)insFile->Get(fTreeTruName);
    
    if ( !PTreeKSig || !PTreePTru )
    {
        cout << "Input Data Tree not found!" << endl;
        return 0;
    }
        
    // Define some simple data structures to Set Branch Addresses
    // Kaon +- Couples B + S
    EVKAONCOUPLE        evKaonSig;
    PTreeKSig  ->SetBranchAddress    ("nKaonCouple",   &evKaonSig.nKaonCouple);
    PTreeKSig  ->SetBranchAddress    ("iKaon",         &evKaonSig.iKaon);
    PTreeKSig  ->SetBranchAddress    ("jKaon",         &evKaonSig.jKaon);
    PTreeKSig  ->SetBranchAddress    ("bPhi",          &evKaonSig.bPhi);
    PTreeKSig  ->SetBranchAddress    ("bRec",          &evKaonSig.bRec);
    PTreeKSig  ->SetBranchAddress    ("bEta",          &evKaonSig.bEta);
    PTreeKSig  ->SetBranchAddress    ("InvMass",       &evKaonSig.InvMass);
    PTreeKSig  ->SetBranchAddress    ("pT",            &evKaonSig.pT);
    
    // Phi decay Kaon +- Couples S
    EVPHI               evPhi_Tru;
    PTreePTru  ->SetBranchAddress    ("nPhi",                 &evPhi_Tru.nPhi);
    PTreePTru  ->SetBranchAddress    ("bRec",                 &evPhi_Tru.bRec);
    PTreePTru  ->SetBranchAddress    ("bEta",                 &evPhi_Tru.bEta);
    PTreePTru  ->SetBranchAddress    ("bKdc",                 &evPhi_Tru.bKdc);
    PTreePTru  ->SetBranchAddress    ("pT",                   &evPhi_Tru.pT);
    
    //---------------------//
    //  Setting up output  //
    //---------------------//
    
    // Generating the binning array--------------------------------------------------------------------------
    fSetBinPT1D();
    fSetBinIM1D();
    fSetBinPT2D();
    fSetBinIM2D();
    Int_t       S_ArrpT[1024];
    
    //-------------------------//
    //  Filling output objects //
    //-------------------------//
    
    int EVL_1D = 0;
    int TOT_1D = 0;
    int EVL_2D = 0;
    int TOT_2D = 0;
    int EVL_1D_PT = 0;
    int TOT_1D_PT = 0;
    int EVL_2D_PT = 0;
    int TOT_2D_PT = 0;
    
    // Evaluating entries
    Int_t nEntries = PTreeKSig->GetEntries();
    for ( Int_t iEvent = 0; iEvent < nEntries; iEvent++ )
    {
        // Recovering events
        PTreeKSig->GetEntry(iEvent);
        
        for (Int_t iPhi = 0; iPhi < evKaonSig.nKaonCouple; iPhi++ )
        {
            // Only Physically relevant couples of Kaons
            if ( !evKaonSig.bEta[iPhi] || !evKaonSig.bPhi[iPhi] ) continue;
            
            TOT_1D++;
            TOT_1D_PT++;
            if ( evKaonSig.InvMass[iPhi] <= fMaxIM1D && evKaonSig.InvMass[iPhi] >= fMinIM1D )
            {
                EVL_1D++;
            }
            if ( evKaonSig.pT[iPhi] <= 0.4 && evKaonSig.pT[iPhi] >= fMinPT1D )
            {
                EVL_1D_PT++;
            }
            
            // Combinatorial mixing
            for (Int_t jPhi = 0; jPhi < evKaonSig.nKaonCouple; jPhi++ )
            {
                // Auto-correlation protection
                if ( iPhi == jPhi ) continue;
                
                // Only Physically relevant couples of Kaons
                if ( !evKaonSig.bEta[jPhi] || !evKaonSig.bPhi[jPhi] ) continue;
                
                TOT_2D++;
                TOT_2D_PT++;
                if ( evKaonSig.InvMass[iPhi] <= fMaxIM1D && evKaonSig.InvMass[iPhi] >= fMinIM1D )
                {
                    if ( evKaonSig.InvMass[jPhi] <= fMaxIM1D && evKaonSig.InvMass[jPhi] >= fMinIM1D )
                    {
                        EVL_2D++;
                    }
                }
                if ( evKaonSig.pT[iPhi] <= 0.4 && evKaonSig.pT[iPhi] >= fMinPT1D )
                {
                     if ( evKaonSig.pT[jPhi] <= 0.4 && evKaonSig.pT[jPhi] >= fMinPT1D )
                    {
                        EVL_2D_PT++;
                    }
                }
            }
        }
    }
    
    cout << "InvMass 1D: " << static_cast<float>(EVL_1D)/TOT_1D << endl;
    cout << "InvMass 2D: " << static_cast<float>(EVL_2D)/TOT_2D << endl;
    cout << "TranMom 1D: " << static_cast<float>(EVL_1D_PT)/TOT_1D_PT << endl;
    cout << "TranMom 2D: " << static_cast<float>(EVL_2D_PT)/TOT_2D_PT << endl;
    
    return 0;
}
