#include "../inc/SetValues.h"
// !TODO: All set!

// Pythia
#include "Pythia8/Pythia.h"

int main (int argc, char *argv[])
{
    if (argc < 3)
    {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./GeneratorMC [filename] [nevents]" << endl;
        return -1;
    }
    
    Int_t   PEvents = atoi(argv[2]);

    // Define tree data structures
    EVKAONCOUPLE    evKaonBkg;
    EVKAONCOUPLE    evKaonSig;
    EVKAON          evKaon;
    EVPHI           evPhi;

    // Output File
    TFile * outFile     = new   TFile   (argv[1],       "recreate","",101);
    
    //Initialisation of TTree
    TTree * PTreeKSig   = new   TTree   (fTreeSigName,  "A ROOT tree for pythia MC - Kaon+- Couples");
    TTree * PTreeKBkg   = new   TTree   (fTreeBkgName,  "A ROOT tree for pythia MC - Kaon++/-- Couples");
    TTree * PTreePTru   = new   TTree   (fTreeTruName,  "A ROOT tree for pythia MC - Phi");
    
    // Setting-Up Kaon Couple TTree
    PTreeKSig->Branch    ("nKaonCouple" ,&evKaonSig.nKaonCouple,    "nKaonCouple/I");
    PTreeKSig->Branch    ("iKaon"       ,&evKaonSig.iKaon,          "iKaon[nKaonCouple]/I");
    PTreeKSig->Branch    ("jKaon"       ,&evKaonSig.jKaon,          "jKaon[nKaonCouple]/I");
    PTreeKSig->Branch    ("bPhi"        ,&evKaonSig.bPhi,           "bPhi[nKaonCouple]/O");
    PTreeKSig->Branch    ("bRec"        ,&evKaonSig.bRec,           "bRec[nKaonCouple]/O");
    PTreeKSig->Branch    ("bEta"        ,&evKaonSig.bEta,           "bEta[nKaonCouple]/O");
    PTreeKSig->Branch    ("InvMass"     ,&evKaonSig.InvMass,        "InvMass[nKaonCouple]/F");
    PTreeKSig->Branch    ("pT"          ,&evKaonSig.pT,             "pT[nKaonCouple]/F");
    
    PTreeKBkg->Branch    ("nKaonCouple" ,&evKaonBkg.nKaonCouple,    "nKaonCouple/I");
    PTreeKBkg->Branch    ("iKaon"       ,&evKaonBkg.iKaon,          "iKaon[nKaonCouple]/I");
    PTreeKBkg->Branch    ("jKaon"       ,&evKaonBkg.jKaon,          "jKaon[nKaonCouple]/I");
    PTreeKBkg->Branch    ("bPhi"        ,&evKaonBkg.bPhi,           "bPhi[nKaonCouple]/O");
    PTreeKBkg->Branch    ("bRec"        ,&evKaonBkg.bRec,           "bRec[nKaonCouple]/O");
    PTreeKBkg->Branch    ("bEta"        ,&evKaonBkg.bEta,           "bEta[nKaonCouple]/O");
    PTreeKBkg->Branch    ("InvMass"     ,&evKaonBkg.InvMass,        "InvMass[nKaonCouple]/F");
    PTreeKBkg->Branch    ("pT"          ,&evKaonBkg.pT,             "pT[nKaonCouple]/F");
    
    // Setting-Up Phi TTree
    PTreePTru->Branch    ("nPhi"        ,&evPhi.nPhi,               "nPhi/I");
    PTreePTru->Branch    ("bEta"        ,&evPhi.bEta,               "bEta[nPhi]/O");
    PTreePTru->Branch    ("bRec"        ,&evPhi.bRec,               "bRec[nPhi]/O");
    PTreePTru->Branch    ("bKdc"        ,&evPhi.bKdc,               "bKdc[nPhi]/O");
    PTreePTru->Branch    ("pT"          ,&evPhi.pT,                 "pT[nPhi]/F");
    
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    
    //Settings
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString(Form("333:mMin = %f",fMinIMMC));
    pythia.readString(Form("333:mMax = %f",fMaxIMMC));
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia.init();
    
    // Save the ID of kaons here
    int nKaon, nPhi, nRecMistake, kaonID[1024], phiID[1024], phiRec[1024], RecMistake[1024];
    bool kaonRec[1024];
    
    // Cycling through events
    for ( int iEvent = 0; iEvent < PEvents; iEvent++ )
    {
        // Next event
        pythia.next();
        
        // Resetting counters
        evPhi.nPhi              = 0;
        evKaon.nKaon            = 0;
        evKaonBkg.nKaonCouple   = 0;
        evKaonSig.nKaonCouple   = 0;
        
        // Starting cycling through event particles
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            const auto particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( particle.id() == 333 )
            {
                evPhi.ID        [evPhi.nPhi]    =   iParticle;
                evPhi.bEta      [evPhi.nPhi]    =   (fabs(particle.p().rap()) <= 0.5);
                evPhi.pT        [evPhi.nPhi]    =   particle.pT();
                evPhi.Dau1      [evPhi.nPhi]    =   particle.daughter1();
                evPhi.Dau2      [evPhi.nPhi]    =   particle.daughter2();
                auto const Dau1                 =   ( pythia.event[particle.daughter1()] );
                auto const Dau2                 =   ( pythia.event[particle.daughter2()] );
                
                evPhi.bKdc      [evPhi.nPhi]    =   ( ( particle.daughterList().size() == 2 ) &&
                                                ( Dau1.id() == -Dau2.id() ) &&
                                                ( abs(Dau1.id()) == 321 ) );
                
                evPhi.bRec      [evPhi.nPhi]    =   ( evPhi.bKdc[evPhi.nPhi] &&
                                                ( fabs(Dau1.eta()) < 0.8 ) &&
                                                ( fabs(Dau2.eta()) < 0.8 ) &&
                                                ( Dau1.pT() > 0.15 ) &&
                                                ( Dau2.pT() > 0.15 ) );
                evPhi.nPhi++;
            }
            
            //Skipping non-final particles
            if ( !particle.isFinal() )       continue;
            
            // Storing Kaons
            if ( fabs(particle.id()) == 321 )
            {
                evKaon.ID       [evKaon.nKaon]  =   iParticle;
                evKaon.bRec     [evKaon.nKaon]  =   (fabs(particle.eta()) < 0.8) && (particle.pT() > 0.15);
                evKaon.Mom1     [evKaon.nKaon]  =   particle.mother1();
                evKaon.Mom2     [evKaon.nKaon]  =   particle.mother2();
                evKaon.nKaon++;
            }
        }
        
        // Cycling through Kaons found
        for ( int iKaon = 0; iKaon < evKaon.nKaon; iKaon++ )
        {
            // Recovering the Kaon
            const auto Kaon1 = pythia.event[evKaon.ID[iKaon]];
            
            for ( int jKaon = (iKaon+1); jKaon < evKaon.nKaon; jKaon++ )
            {
                // Recovering the Kaon
                const auto Kaon2 = pythia.event[evKaon.ID[jKaon]];
                
                // Building the candidate Phi
                const auto pPhi = Kaon1.p() + Kaon2.p();
                
                //Cut on Invariant Mass not in resonance region
                if (pPhi.mCalc() < fMinIMMC) continue;
                if (pPhi.mCalc() > fMaxIMMC) continue;
                
                // Unlike Sign ( Sig + Bkg )
                if ( Kaon1.id() == -Kaon2.id() )
                {
                    evKaonSig.InvMass[evKaonSig.nKaonCouple]    =   pPhi.mCalc();
                    evKaonSig.pT[evKaonSig.nKaonCouple]         =   pPhi.pT();
                    evKaonSig.bRec[evKaonSig.nKaonCouple]       =   ( evKaon.bRec[iKaon] && evKaon.bRec[jKaon] );
                    evKaonSig.bEta[evKaonSig.nKaonCouple]       =   ( fabs(pPhi.rap()) <= 0.5 );
                    evKaonSig.iKaon[evKaonSig.nKaonCouple]      =   iKaon;
                    evKaonSig.jKaon[evKaonSig.nKaonCouple]      =   jKaon;
                    
                    evKaonSig.bPhi[evKaonSig.nKaonCouple]       =   ( evKaon.Mom1[iKaon] == evKaon.Mom1[jKaon] &&
                                                                     evKaon.Mom2[iKaon] == evKaon.Mom2[jKaon] &&
                                                                     ( pythia.event[evKaon.Mom1[iKaon]]).id() == 333 &&
                                                                     evKaon.Mom2[iKaon] == 0);
                    evKaonSig.nKaonCouple++;
                }
                
                // Like Sign ( Bkg )
                if ( Kaon1.id() == Kaon2.id() )
                {
                    if (evKaon.bRec[iKaon] && evKaon.bRec[jKaon]) continue;
                    if (fabs(pPhi.rap()) <= 0.5) continue;
                    evKaonBkg.InvMass[evKaonBkg.nKaonCouple]   =    pPhi.mCalc();
                    evKaonBkg.pT[evKaonBkg.nKaonCouple]        =    pPhi.pT();
                    evKaonBkg.nKaonCouple++;
                }
            }
        }
        PTreeKSig     ->Fill();
        PTreeKBkg     ->Fill();
        PTreePTru     ->Fill();
    }
    PTreeKSig   ->Write();
    PTreeKBkg   ->Write();
    PTreePTru   ->Write();
    outFile     ->Close();
    return 0;
}
