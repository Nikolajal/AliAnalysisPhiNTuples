#include "../inc/SetValues.h"
// !TODO: All set!

// Pythia
#include "Pythia8/Pythia.h"

int main (int argc, char *argv[])
{
    // Check everything is good
    if (argc < 3)
    {
        cout << "ERROR: Insufficient parameters given!" << endl;
        cout << "Please use as: ./GeneratorMC [filename] [nevents]" << endl;
        return -1;
    }
    
    // Definition of number of events
    Int_t   nEvents = atoi(argv[2]);
    
    // Output File
    TFile * outFile     = new   TFile   (Form("%s.root",argv[1]),   "recreate", "", 101);
    
    //Initialisation of TTree
    TTree * fPhiCandidate,  * fPhiEfficiency,   * fKaonCandidate;
    
    // Define tree data structures
    Struct_PhiCandidate     evPhiCandidate;
    Struct_PhiEfficiency    evPhiEfficiency;
    Struct_KaonCandidate    evKaonCandidate;
    
    // PhiCandidate Tree Set-Up
    fPhiCandidate = new TTree   ("PhiCandidate",    "Data Tree for Phi Candidates");
    fPhiCandidate->Branch       ("Multiplicity",    &evPhiCandidate.Multiplicity,   "Multiplicity/F");
    fPhiCandidate->Branch       ("nPhi",            &evPhiCandidate.nPhi,           "nPhi/b");
    fPhiCandidate->Branch       ("Px",              &evPhiCandidate.Px,             "Px[nPhi]/F");
    fPhiCandidate->Branch       ("Py",              &evPhiCandidate.Py,             "Py[nPhi]/F");
    fPhiCandidate->Branch       ("Pz",              &evPhiCandidate.Pz,             "Pz[nPhi]/F");
    fPhiCandidate->Branch       ("InvMass",         &evPhiCandidate.InvMass,        "InvMass[nPhi]/F");
    fPhiCandidate->Branch       ("iKaon",           &evPhiCandidate.iKaon,          "iKaon[nPhi]/b");
    fPhiCandidate->Branch       ("jKaon",           &evPhiCandidate.jKaon,          "jKaon[nPhi]/b");
    
    // KaonCandidate Tree Set-Up
    fKaonCandidate = new TTree  ("KaonCandidate",   "Data Tree for Kaon Candidates");
    fKaonCandidate->Branch      ("Multiplicity",    &evKaonCandidate.Multiplicity,  "Multiplicity/F");
    fKaonCandidate->Branch      ("nKaon",           &evKaonCandidate.nKaon,         "nKaon/b");
    fKaonCandidate->Branch      ("Px",              &evKaonCandidate.Px,            "Px[nKaon]/F");
    fKaonCandidate->Branch      ("Py",              &evKaonCandidate.Py,            "Py[nKaon]/F");
    fKaonCandidate->Branch      ("Pz",              &evKaonCandidate.Pz,            "Pz[nKaon]/F");
    fKaonCandidate->Branch      ("Charge",          &evKaonCandidate.Charge,        "Charge[nKaon]/B");
    fKaonCandidate->Branch      ("TOFSigma",        &evKaonCandidate.SigmaTOF,      "TOFSigma[nKaon]/B");
    fKaonCandidate->Branch      ("TPCSigma",        &evKaonCandidate.SigmaTPC,      "TPCSigma[nKaon]/B");
    
    // PhiEfficiency Tree Set-Up
    fPhiEfficiency = new TTree  ("PhiEfficiency",   "MC Tree for Phi Efficiency");
    fPhiEfficiency->Branch      ("nPhi",            &evPhiEfficiency.nPhi,          "nPhi/b");
    fPhiEfficiency->Branch      ("Px",              &evPhiEfficiency.Px,            "Px[nPhi]/F");
    fPhiEfficiency->Branch      ("Py",              &evPhiEfficiency.Py,            "Py[nPhi]/F");
    fPhiEfficiency->Branch      ("Pz",              &evPhiEfficiency.Pz,            "Pz[nPhi]/F");
    fPhiEfficiency->Branch      ("Selection",       &evPhiEfficiency.Selection,     "Selection[nPhi]/b");
    
    // PYTHIA INITIALISATION
    Pythia8::Pythia pythia;
    
    // Settings
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString(Form("333:mMin = %f",fMinIMMC));
    pythia.readString(Form("333:mMax = %f",fMaxIMMC));
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");
    pythia.init();
    
    // Save the ID of kaons here
    //int nKaon, nPhi, nRecMistake, kaonID[1024], phiID[1024], phiRec[1024], RecMistake[1024];
    //bool kaonRec[1024];
    
    // Cycling through events
    for ( int iEvent = 0; iEvent < nEvents; iEvent++ )
    {
        // Next event
        pythia.next();
        
        // Set counters to zero
        evPhiCandidate.nPhi     =   0;
        evPhiEfficiency.nPhi    =   0;
        evKaonCandidate.nKaon   =   0;
        
        // Starting cycling through event particles
        for ( int iParticle = 0; iParticle < pythia.event.size() ; iParticle++ )
        {
            // Saving particles
            const auto Current_Particle = pythia.event[iParticle];
            
            // Storing True Phis
            if ( Current_Particle.id() == 333 )
            {
                evPhiEfficiency.Px[evPhiEfficiency.nPhi]        =   Current_Particle.px();
                evPhiEfficiency.Py[evPhiEfficiency.nPhi]        =   Current_Particle.py();
                evPhiEfficiency.Pz[evPhiEfficiency.nPhi]        =   Current_Particle.pz();
                evPhiEfficiency.Selection[evPhiEfficiency.nPhi] =   0;
                
                
                auto const Dau1                                 =   ( pythia.event[Current_Particle.daughter1()] );
                auto const Dau2                                 =   ( pythia.event[Current_Particle.daughter2()] );
                
                if  ( ( Current_Particle.daughterList().size() == 2 ) &&
                     ( Dau1.id() == -Dau2.id() ) &&
                     ( abs(Dau1.id()) == 321 ) )                evPhiEfficiency.Selection[evPhiEfficiency.nPhi]++;
                
                if  ( evPhiEfficiency.Selection[evPhiEfficiency.nPhi] == 1 &&
                     ( fabs(Dau1.eta()) < 0.8 ) &&
                     ( fabs(Dau2.eta()) < 0.8 ) &&
                     ( Dau1.pT() > 0.15 ) &&
                     ( Dau2.pT() > 0.15 ) )                     evPhiEfficiency.Selection[evPhiEfficiency.nPhi]++;
                /*
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
                */
                
                evPhiEfficiency.nPhi++;
            }
            
            //Skipping non-final particles
            if ( !Current_Particle.isFinal() )       continue;
            
            // Storing Kaons
            if ( fabs(Current_Particle.id()) == 321 )
            {
                /*
                evKaon.ID       [evKaon.nKaon]  =   iParticle;
                evKaon.bRec     [evKaon.nKaon]  =   (fabs(particle.eta()) < 0.8) && (particle.pT() > 0.15);
                evKaon.Mom1     [evKaon.nKaon]  =   particle.mother1();
                evKaon.Mom2     [evKaon.nKaon]  =   particle.mother2();
                evKaon.nKaon++;
                 */
            }
            /*
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
             */
        }
        fKaonCandidate  ->Fill();
        fPhiCandidate   ->Fill();
        fPhiEfficiency  ->Fill();
    }
    fKaonCandidate  ->Write();
    fPhiCandidate   ->Write();
    fPhiEfficiency  ->Write();
    
    outFile     ->Close();
    return 0;
}
