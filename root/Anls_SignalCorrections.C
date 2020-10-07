// File for 1-Dimensional Analysis:
// !TODO: All Set!


#include "../inc/SetValues.h"
#include "../inc/SetFunctions.h"
#include "RooMsgService.h"

void Anls_SignalCorrections ( bool fSilent = false )
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
    hName                           =   "Entry_DT";                             // Name of the histogram in the preprocessed file
    TH1F *  h1D_EDT                 =   (TH1F*)(insFile_DT->Get(hName));
    hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Tru                 =   (TH1F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_1D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Rec                 =   (TH1F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Tru_P6              =   (TH1F*)(insFile_MC_P6->Get(hName));
    hName                           =   "hNP_1D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Tru_P8              =   (TH1F*)(insFile_MC_P8->Get(hName));
    hName                           =   "hNP_1D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Rec_P6              =   (TH1F*)(insFile_MC_P6->Get(hName));
    hName                           =   "hNP_1D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
    TH1F *  h1D_Rec_P8              =   (TH1F*)(insFile_MC_P8->Get(hName));
    
    //------ 2D Histograms Recovery ------//
    
    hName                           =   "h2D_Raw";                              // Name of the histogram in the preprocessed file
    TH2F *  h2D_Raw                 =   (TH2F*)(insFile_DT->Get(hName));
    hName                           =   "hNP_2D_Eff_X2_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Eff                 =   (TH2F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_2D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Tru                 =   (TH2F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_2D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Rec                 =   (TH2F*)(insFile_EF->Get(hName));
    hName                           =   "hNP_2D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Tru_P6              =   (TH2F*)(insFile_MC_P6->Get(hName));
    hName                           =   "hNP_2D_Tru_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Tru_P8              =   (TH2F*)(insFile_MC_P8->Get(hName));
    hName                           =   "hNP_2D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Rec_P6              =   (TH2F*)(insFile_MC_P6->Get(hName));
    hName                           =   "hNP_2D_Rec_PT_S";                      // Name of the histogram in the preprocessed file
    TH2F *  h2D_Rec_P8              =   (TH2F*)(insFile_MC_P8->Get(hName));
    
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
    h1D_Raw->Scale(1.,"width");
    h2D_Raw->Scale(1.,"width");
    h1D_Res->Scale(1.,"width");
    h2D_Res->Scale(1.,"width");
    
    // Scaling in Rapidity Interval
    h1D_Res->Scale(1./kRapidityInterval);
    h2D_Res->Scale(1./kRapidityInterval);
    
    // Scaling for events corrected for efficiency
    h1D_Raw->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
    h2D_Raw->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
    h1D_Res->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
    h2D_Res->Scale(kEventEfficiency_/(h1D_EDT->GetBinContent(1)));
    
    // Scaling for Branching Ratio
    h1D_Res->Scale(1./kBranchingRatio__);
    h2D_Res->Scale(1./(kBranchingRatio__*kBranchingRatio__));
    
    if ( bPythiaTest )
    {
        // Scaling for MC biases
        h1D_Raw->Scale(1./(kEventEfficiency_*kPythia1DEfficien));
        h2D_Raw->Scale(1./(kEventEfficiency_*kPythia2DEfficien));
        h1D_Res->Scale(1./(kEventEfficiency_*kPythia1DEfficien));
        h2D_Res->Scale(1./(kEventEfficiency_*kPythia2DEfficien));
    }
    else
    {
        // Scaling for Vertex Selection Efficiency
        h1D_Res->Scale(1./kVertexEfficiency);
        h2D_Res->Scale(1./kVertexEfficiency);
    }
    
    //------- Histograms Fitting  -------//
    
    /*
    SetLevyTsalis();
    gCheck_->Fit(fLevyFit1D,"IMRE0S","",0.4,6.0);
    cout << "M:" << fLevyFit1D->Mean(0.,10.) << endl;
    TCanvas * c133 = new TCanvas();
    //gPad->SetLogy();
    //gCheck_->Draw();
    //fLevyFit1D->Draw("same");
    h1D_Res->Divide(fLevyFit1D);
    h1D_Res->Draw();
    c133->SaveAs("./graphs/eeee.pdf");
    */
    
    // 1D PT Eval
    Double_t *  fResults1D = new Double_t [6];
    gSystem->RedirectOutput("/dev/null");
    fResults1D = ExtrapolateVl( h1D_Res );
    gSystem->RedirectOutput(0,0);
    
    // 2D PT Eval
    Double_t***     fResults2D =   new Double_t**  [nBinPT2D+1];
    for ( Int_t iTer = 0; iTer < nBinPT2D+1; iTer++ )
    {
        fResults2D[iTer]        =   new Double_t*   [2];
        fResults2D[iTer][0]     =   new Double_t    [6];
        fResults2D[iTer][1]     =   new Double_t    [6];
    }
    gSystem->RedirectOutput("/dev/null");
    fResults2D = ExtrapolateVl( h2D_Res );
    gSystem->RedirectOutput(0,0);
    
    Double_t errREC, errRAW;
    Double_t intREC = h1D_Rec->IntegralAndError(5,100,errREC,"width");
    Double_t intRAW = h1D_Raw->IntegralAndError(5,100,errRAW,"width");
    cout << "1D:" << endl;
    cout << "INT_Rec: " << intREC               << "+-" << errREC                                       << endl;
    cout << "INT_Raw: " << intRAW               << "+-" << errRAW                                       << endl;
    cout << "Raw/Rec: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
    cout << endl;
    intREC = h1D_Tru->IntegralAndError(5,100,errREC,"width");
    intRAW = h1D_Res->IntegralAndError(5,100,errRAW,"width");
    cout << "1D:" << endl;
    cout << "INT_Tru: " << intREC               << "+-" << errREC                                       << endl;
    cout << "INT_Res: " << intRAW               << "+-" << errRAW                                       << endl;
    cout << "Res/Tru: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
    cout << endl;
    intREC = h1D_Tru->IntegralAndError(-1,100,errREC,"width");
    cout << "1D:" << endl;
    cout << "INT_Tru: " << intREC               << "+-" << errREC                                                   << endl;
    cout << "INT_Res: " << fResults1D[0]        << "+-" << fResults1D[1] << " +-" << fResults1D[2]                  << endl;
    cout << "Res/Tru: " << fResults1D[0]/intREC << "+-" << (1/intREC)*(fResults1D[1]+errREC*(fResults1D[0]/intREC)) << endl;
    cout << endl;
    cout << endl;
    intREC = h2D_Rec->IntegralAndError(5,100,5,100,errREC,"width");
    intRAW = h2D_Raw->IntegralAndError(5,100,5,100,errRAW,"width");
    cout << "2D:" << endl;
    cout << "INT_Rec: " << intREC               << "+-" << errREC                                       << endl;
    cout << "INT_Raw: " << intRAW               << "+-" << errRAW                                       << endl;
    cout << "Raw/Rec: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
    cout << endl;
    intREC = h2D_Tru->IntegralAndError(5,100,5,100,errREC,"width");
    intRAW = h2D_Res->IntegralAndError(5,100,5,100,errRAW,"width");
    cout << "2D:" << endl;
    cout << "INT_Tru: " << intREC               << "+-" << errREC                                       << endl;
    cout << "INT_Res: " << intRAW               << "+-" << errRAW                                       << endl;
    cout << "Res/Tru: " << intRAW/intREC        << "+-" << (1/intREC)*(errRAW+errREC*(intRAW/intREC))   << endl;
    cout << endl;
    intREC = h2D_Tru->IntegralAndError(-1,100,-1,100,errREC,"width");
    cout << "2DX:" << endl;
    cout << "INT_Tru: " << intREC                       << "+-" << errREC                                                               << endl;
    cout << "INT_Res: " << fResults2D[0][0][0]          << "+-" << fResults2D[0][0][1]  << " +-" << fResults2D[0][0][2]                 << endl;
    cout << "Res/Tru: " << fResults2D[0][0][0]/intREC   << "+-" << (1/intREC)*(fResults2D[0][0][1]+errREC*(fResults2D[0][0][0]/intREC))<< endl;
    cout << "2DY:" << endl;
    cout << "INT_Tru: " << intREC                       << "+-" << errREC                                                               << endl;
    cout << "INT_Res: " << fResults2D[0][1][0]          << "+-" << fResults2D[0][1][1] << " +-" << fResults2D[0][1][2]                  << endl;
    cout << "Res/Tru: " << fResults2D[0][1][0]/intREC   << "+-" << (1/intREC)*(fResults2D[0][1][1]+errREC*(fResults2D[0][1][0]/intREC))<< endl;
    cout << endl;
    auto err1p  = sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERP*fResults1D[0]*kEventEfficienERP);
    auto err1m  = sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERM*fResults1D[0]*kEventEfficienERM);
    cout << "1D YIELD:" << endl;
    cout << "YIELD: " << fResults1D[0]    << "+-" <<  fResults1D[1]   << "+" <<   err1p  << "-" <<  err1m  << endl;
    cout << "MeanP: " << fResults1D[3]    << "+-" <<  fResults1D[4]   << endl;
    cout << endl;
    auto errxp  =   sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERP*fResults2D[0][0][0]*kEventEfficienERP);
    auto errxm  =   sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERM*fResults2D[0][0][0]*kEventEfficienERM);
    auto erryp  =   sqrt(fResults2D[0][1][2]*fResults2D[0][1][2]+fResults2D[0][1][0]*kEventEfficienERP*fResults2D[0][1][0]*kEventEfficienERP);
    auto errym  =   sqrt(fResults2D[0][1][2]*fResults2D[0][1][2]+fResults2D[0][1][0]*kEventEfficienERM*fResults2D[0][1][0]*kEventEfficienERM);
    cout << "2DX YIELD:" << endl;
    cout << "YIELD: " << fResults2D[0][0][0]    << "+-" <<  fResults2D[0][0][1]   << "+" <<  errxp    << "-" <<   errxm  << endl;
    cout << "MeanP: " << fResults2D[0][0][3]    << "+-" <<  fResults2D[0][0][4]   << endl;
    cout << endl;
    cout << "2DY YIELD:" << endl;
    cout << "YIELD: " << fResults2D[0][1][0]    << "+-" <<  fResults2D[0][1][1]   << "+" <<   erryp   << "-" <<   errym << endl;
    cout << "MeanP: " << fResults2D[0][1][3]    << "+-" <<  fResults2D[0][1][4]   << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "P61D: " << h1D_Tru_P6->Integral("width") << endl;
    cout << "P81D: " << h1D_Tru_P8->Integral("width") << endl;
    cout << "P62D: " << h2D_Tru_P6->Integral("width") << endl;
    cout << "P82D: " << h2D_Tru_P8->Integral("width") << endl;
    
    auto valy   =   fResults2D[0][0][0]/(fResults1D[0]*fResults1D[0]);
    auto valx   =   fResults2D[0][1][0]/(fResults1D[0]*fResults1D[0]);
    auto errx   =   sqrt(fResults2D[0][0][1]*fResults2D[0][1][1]/(fResults2D[0][0][0]*fResults2D[0][0][0])+4*(fResults1D[1]*fResults1D[1])/(fResults1D[0]*fResults1D[0]));
    auto erry   =   sqrt(fResults2D[0][1][1]*fResults2D[0][1][1]/(fResults2D[0][1][0]*fResults2D[0][1][0])+4*(fResults1D[1]*fResults1D[1])/(fResults1D[0]*fResults1D[0]));
    auto erxP   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERP,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][0][0]/fResults1D[0],2)),2)));
    auto eryP   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERP,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][1][0]/fResults1D[0],2)),2)));
    auto erxM   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERM,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][0][0]/fResults1D[0],2)),2)));
    auto eryM   =   sqrt((pow(kSyst_SigExtr2D,2)+pow(kEventEfficienERM,2)+(pow((2*fResults1D[2]),2)-pow(2*kSyst_SigExtr1D*fResults1D[0],2))/pow((fResults1D[0]),2)*pow(1+(pow(fResults2D[0][1][0]/fResults1D[0],2)),2)));
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "2DX/1D^2:" << valx << "+-" << errx << "+" << valx*erxP << "-" << valx*erxM << endl;
    cout << endl;
    cout << "2DY/1D^2:" << valy << "+-" << erry << "+" << valx*erxP << "-" << valx*erxM << endl;
    cout << endl;
    cout << "2DY/1D^2 P6:" << (h2D_Tru_P6->Integral("width"))/(pow(h1D_Tru_P6->Integral("width"),2)) << endl;
    cout << "2DY/1D^2 P8:" << (h2D_Tru_P8->Integral("width"))/(pow(h1D_Tru_P8->Integral("width"),2)) << endl;
    cout << endl;
    auto errsigstat =   sqrt(2*errx*errx + pow(1-2*fResults1D[0],2)*pow(fResults1D[1],2));
    auto errsigstatP=   sqrt(2*valx*erxP*valx*erxP + pow(1-2*fResults1D[0],2)*(pow(fResults1D[2],2)+kEventEfficienERP*kEventEfficienERP));
    auto errsigstatM=   sqrt(2*valx*erxM*valx*erxM + pow(1-2*fResults1D[0],2)*(pow(fResults1D[2],2)+kEventEfficienERM*kEventEfficienERM));
    cout << "Sigma:" << 2*fResults2D[0][0][0] + fResults1D[0] -fResults1D[0]*fResults1D[0] << "+-" << errsigstat  <<"+" << errsigstatP << "-" << errsigstatM << endl;
    cout << "Sigma6:" << 2*(h2D_Tru_P6->Integral("width")) + h1D_Tru_P6->Integral("width") - (pow(h1D_Tru_P6->Integral("width"),2)) << endl;
    cout << "Sigma8:" << 2*(h2D_Tru_P8->Integral("width")) + h1D_Tru_P8->Integral("width") - (pow(h1D_Tru_P8->Integral("width"),2)) << endl;
    cout << endl;
    cout << endl;
    cout << endl;
    auto poisson = 2*fResults2D[0][0][0]/fResults1D[0]-fResults1D[0];
    auto poissonstat =   sqrt(4*pow(poisson*fResults1D[1]/fResults1D[0],2) + 4*pow(poisson*fResults2D[0][0][1]/fResults2D[0][0][0],2) + pow(fResults1D[1],2));
    auto poissonsyst=   sqrt(poisson*poisson*(kBranchingRatio__*kBranchingRatio__+kSyst_TrackEff1D*kSyst_TrackEff1D+kSyst_PID1D*kSyst_PID1D+kSyst_SigExtr1D*kSyst_SigExtr1D+kSyst_SigExtr2D*kSyst_SigExtr2D)+pow(fResults1D[2],2));
    auto errsigst2tM=   sqrt(poisson*poisson*(kBranchingRatio__*kBranchingRatio__+kSyst_TrackEff1D*kSyst_TrackEff1D+kSyst_PID1D*kSyst_PID1D+kSyst_SigExtr1D*kSyst_SigExtr1D+kSyst_SigExtr2D*kSyst_SigExtr2D)+pow(fResults1D[2],2));
    cout << "Sigma:" <<  poisson << "+-" << poissonstat  <<"+" << poissonsyst << "-" << poissonsyst << endl;
    cout << "Sigma6:" << 2*(h2D_Tru_P6->Integral("width")/h1D_Tru_P6->Integral("width")) - h1D_Tru_P6->Integral("width") << endl;
    cout << "Sigma8:" << 2*(h2D_Tru_P8->Integral("width")/h1D_Tru_P8->Integral("width")) - h1D_Tru_P8->Integral("width") << endl;
    
    //---- Graphics of Presentation -----//
    
    TGraphAsymmErrors   *g1D_Res_Stat    =   new TGraphAsymmErrors(h1D_Res);
    TGraphAsymmErrors   *g1D_Res_Syst    =   new TGraphAsymmErrors(SetSystErrorsh(h1D_Res));
    for ( Int_t iTer = 0; iTer < 4; iTer++ )
    {
        g1D_Res_Stat->RemovePoint(0);
        g1D_Res_Syst->RemovePoint(0);
    }
    
    TGrCompare1D(g1D_Res_Stat,g1D_Res_Syst,h1D_Tru_P6,h1D_Tru_P8,gCheck_);
    TGrCompare2D(h2D_Res,h2D_Tru_P6,h2D_Tru_P8);
    TGraphAEGeneratorPT(fResults2D,h2D_Tru_P6,h2D_Tru_P8);
    
    TCanvas *c____1 = new TCanvas("c1", "Intensity of LED 1",0,  0, 800, 600);
    c____1->SetMargin(0.2,0.1,0.9,0.9);
    
    TLegend        *lLegend1D           =   new TLegend(0.6,0.45,0.9,0.6);
    TMultiGraph    *gFinal1D            =   new TMultiGraph();
    TGraphErrors   *gFinal1DDataStat    =   new TGraphErrors();
    gFinal1DDataStat->SetPoint(0,1,fResults1D[0]);
    gFinal1DDataStat->SetPointError(0,0,fResults1D[1]);
    TGraphAsymmErrors   *gFinal1DDataSyst    =   new TGraphAsymmErrors();
    gFinal1DDataSyst->SetPoint(0,1,fResults1D[0]);
    gFinal1DDataSyst->SetPointError(0,0.025,0.025,sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERM*fResults1D[0]*kEventEfficienERM),sqrt(fResults1D[2]*fResults1D[2]+fResults1D[0]*kEventEfficienERP*fResults1D[0]*kEventEfficienERP));
    TGraphErrors   *gFinal1DPythia6_    =   new TGraphErrors();
    gFinal1DPythia6_->SetPoint(0,1.05,h1D_Tru_P6->Integral("width"));
    TGraphErrors   *gFinal1DPythia8_    =   new TGraphErrors();
    gFinal1DPythia8_->SetPoint(0,1.1,h1D_Tru_P8->Integral("width"));
    
    gFinal1DDataStat->SetMarkerStyle(20);
    gFinal1DDataStat->SetMarkerColor(46);
    gFinal1DDataStat->SetLineColor(kRed);
    gFinal1DDataSyst->SetMarkerStyle(20);
    gFinal1DDataSyst->SetMarkerColor(46);
    gFinal1DDataSyst->SetLineColor(kRed);
    gFinal1DPythia6_->SetMarkerStyle(26);
    gFinal1DPythia6_->SetMarkerColor(kBlue);
    gFinal1DPythia8_->SetMarkerStyle(32);
    gFinal1DPythia8_->SetMarkerColor(kBlue+3);
    
    gFinal1D->Add(gFinal1DDataSyst,"EP5");
    gFinal1D->Add(gFinal1DDataStat,"EP");
    gFinal1D->Add(gFinal1DPythia6_,"P");
    gFinal1D->Add(gFinal1DPythia8_,"P");
    gFinal1D->GetYaxis()->SetTitle("#frac{dN_{#phi}}{dy}");
    gFinal1D->GetXaxis()->SetLimits(0.85,1.5);
    gFinal1D->SetMaximum(0.041);
    gFinal1D->SetMinimum(0.025);
    TAxis *xAxis1  =   new TAxis(*gFinal1D->GetYaxis());
    
    lLegend1D->AddEntry(gFinal1DDataStat,"Stat. Err.","EP");
    lLegend1D->AddEntry(gFinal1DDataSyst,"Syst. Err.","F");
    lLegend1D->AddEntry(gFinal1DPythia6_,"Pythia 6","P");
    lLegend1D->AddEntry(gFinal1DPythia8_,"Pythia 8","P");
    
    gFinal1D->Draw("PA");
    lLegend1D->Draw("same");
    xAxis1->Draw();
    c____1->Write();
    c____1->SaveAs("./graphs/c____1.pdf");
    c____1->SaveAs("./graphs/c____1.png");
    
    TCanvas *c____2 = new TCanvas("c22", "Intensity of LED 1",0,  0, 800, 600);
    c____2->SetMargin(0.2,0.1,0.9,0.9);
    
    TMultiGraph    *gFinal2D            =   new TMultiGraph();
    TGraphErrors   *gFinal2DDataStat    =   new TGraphErrors();
    gFinal2DDataStat->SetPoint(0,1,fResults2D[0][0][0]);
    gFinal2DDataStat->SetPointError(0,0,fResults2D[0][0][1]);
    TGraphAsymmErrors   *gFinal2DDataSyst    =   new TGraphAsymmErrors();
    gFinal2DDataSyst->SetPoint(0,1,fResults2D[0][0][0]);
    gFinal2DDataSyst->SetPointError(0,0.025,0.025,sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERM*fResults2D[0][0][0]*kEventEfficienERM),sqrt(fResults2D[0][0][2]*fResults2D[0][0][2]+fResults2D[0][0][0]*kEventEfficienERP*fResults2D[0][0][0]*kEventEfficienERP));
    TGraphErrors   *gFinal2DPythia6_    =   new TGraphErrors();
    gFinal2DPythia6_->SetPoint(0,1.05,h2D_Tru_P6->Integral("width"));
    TGraphErrors   *gFinal2DPythia8_    =   new TGraphErrors();
    gFinal2DPythia8_->SetPoint(0,1.1,h2D_Tru_P8->Integral("width"));
    
    gFinal2DDataStat->SetMarkerStyle(20);
    gFinal2DDataStat->SetMarkerColor(46);
    gFinal2DDataStat->SetLineColor(kRed);
    gFinal2DDataSyst->SetMarkerStyle(20);
    gFinal2DDataSyst->SetMarkerColor(46);
    gFinal2DDataSyst->SetLineColor(kRed);
    gFinal2DPythia6_->SetMarkerStyle(26);
    gFinal2DPythia6_->SetMarkerColor(kBlue);
    gFinal2DPythia8_->SetMarkerStyle(32);
    gFinal2DPythia8_->SetMarkerColor(kBlue+3);
    
    gFinal2D->Add(gFinal2DDataSyst,"EP5");
    gFinal2D->Add(gFinal2DDataStat,"EP");
    gFinal2D->Add(gFinal2DPythia6_,"P");
    gFinal2D->Add(gFinal2DPythia8_,"P");
    gFinal2D->GetYaxis()->SetTitle("#frac{dN_{#phi#phi}}{dy}");
    gFinal2D->GetXaxis()->SetLimits(0.85,1.5);
    gFinal2D->SetMaximum(0.0021);
    gFinal2D->SetMinimum(0.0008);
    TAxis *xAxis2  =   new TAxis(*gFinal2D->GetYaxis());
    
    gFinal2D->Draw("PA");
    lLegend1D->Draw("same");
    xAxis2->Draw();
    c____2->Write();
    c____2->SaveAs("./graphs/c____2.pdf");
    c____2->SaveAs("./graphs/c____2.png");
    
    TCanvas *c____3 = new TCanvas("c32", "Intensity of LED 1",0,  0, 800, 600);
    c____3->SetMargin(0.2,0.1,0.9,0.9);
    
    TMultiGraph    *gFinal3D            =   new TMultiGraph();
    TGraphErrors   *gFinal3DDataStat    =   new TGraphErrors();
    gFinal3DDataStat->SetPoint(0,1,poisson);
    gFinal3DDataStat->SetPointError(0,0,poissonstat);
    TGraphAsymmErrors   *gFinal3DDataSyst    =   new TGraphAsymmErrors();
    gFinal3DDataSyst->SetPoint(0,1,poisson);
    gFinal3DDataSyst->SetPointError(0,0.025,0.025,poissonsyst,poissonsyst);
    TGraphErrors   *gFinal3DPythia6_    =   new TGraphErrors();
    gFinal3DPythia6_->SetPoint(0,1.05,2*(h2D_Tru_P6->Integral("width")/h1D_Tru_P6->Integral("width")) - h1D_Tru_P6->Integral("width"));
    TGraphErrors   *gFinal3DPythia8_    =   new TGraphErrors();
    gFinal3DPythia8_->SetPoint(0,1.1,2*(h2D_Tru_P8->Integral("width")/h1D_Tru_P8->Integral("width")) - h1D_Tru_P8->Integral("width"));
    TGraph * hPoissonLine = new TGraph();
    hPoissonLine->SetPoint(0,0,0);
    hPoissonLine->SetPoint(1,3,0);
    
    lLegend1D->AddEntry(hPoissonLine,"Pois. Distr.","L");
    
    gFinal3DDataStat->SetMarkerStyle(20);
    gFinal3DDataStat->SetMarkerColor(46);
    gFinal3DDataStat->SetLineColor(kRed);
    gFinal3DDataSyst->SetMarkerStyle(20);
    gFinal3DDataSyst->SetMarkerColor(46);
    gFinal3DDataSyst->SetLineColor(kRed);
    gFinal3DPythia6_->SetMarkerStyle(26);
    gFinal3DPythia6_->SetMarkerColor(kBlue);
    gFinal3DPythia8_->SetMarkerStyle(32);
    gFinal3DPythia8_->SetMarkerColor(kBlue+3);
    
    gFinal3D->Add(gFinal3DDataSyst,"EP5");
    gFinal3D->Add(gFinal3DDataStat,"EP");
    gFinal3D->Add(gFinal3DPythia6_,"P");
    gFinal3D->Add(gFinal3DPythia8_,"P");
    gFinal3D->Add(hPoissonLine,"L");
    gFinal3D->GetYaxis()->SetTitle("#sigma^{2}_{#phi}/#mu_{#phi} -1");
    gFinal3D->GetXaxis()->SetLimits(0.85,1.5);
    gFinal3D->SetMaximum(0.1);
    gFinal3D->SetMinimum(-0.05);
    TAxis *xAxis3  =   new TAxis(*gFinal3D->GetYaxis());
    
    gFinal3D->Draw("PA");
    lLegend1D->Draw("same");
    xAxis3->Draw();
    c____3->Write();
    c____3->SaveAs("./graphs/c____3.pdf");
    c____3->SaveAs("./graphs/c____3.png");
    
    TCanvas *c____4 = new TCanvas("c3244", "Intensity of LED 1",0,  0, 800, 600);
    c____4->SetMargin(0.2,0.1,0.9,0.9);
    
    TMultiGraph    *gFinal4D            =   new TMultiGraph();
    TGraphErrors   *gFinal4DDataStat    =   new TGraphErrors();
    gFinal4DDataStat->SetPoint(0,1,valx);
    gFinal4DDataStat->SetPointError(0,0,errx);
    TGraphAsymmErrors   *gFinal4DDataSyst    =   new TGraphAsymmErrors();
    gFinal4DDataSyst->SetPoint(0,1,valx);
    gFinal4DDataSyst->SetPointError(0,0.025,0.025,erxM,erxP);
    TGraphErrors   *gFinal4DPythia6_    =   new TGraphErrors();
    gFinal4DPythia6_->SetPoint(0,1.05,(h2D_Tru_P6->Integral("width")/pow(h1D_Tru_P6->Integral("width"),2)));
    TGraphErrors   *gFinal4DPythia8_    =   new TGraphErrors();
    gFinal4DPythia8_->SetPoint(0,1.1,(h2D_Tru_P8->Integral("width")/pow(h1D_Tru_P8->Integral("width"),2)));
    TGraph * hPoissonLin2 = new TGraph();
    hPoissonLin2->SetPoint(0,0,0.5);
    hPoissonLin2->SetPoint(1,3,0.5);
    
    gFinal4DDataStat->SetMarkerStyle(20);
    gFinal4DDataStat->SetMarkerColor(46);
    gFinal4DDataStat->SetLineColor(kRed);
    gFinal4DDataSyst->SetMarkerStyle(20);
    gFinal4DDataSyst->SetMarkerColor(46);
    gFinal4DDataSyst->SetLineColor(kRed);
    gFinal4DPythia6_->SetMarkerStyle(26);
    gFinal4DPythia6_->SetMarkerColor(kBlue);
    gFinal4DPythia8_->SetMarkerStyle(32);
    gFinal4DPythia8_->SetMarkerColor(kBlue+3);
    
    gFinal4D->Add(gFinal4DDataSyst,"EP5");
    gFinal4D->Add(gFinal4DDataStat,"EP");
    gFinal4D->Add(gFinal4DPythia6_,"P");
    gFinal4D->Add(gFinal4DPythia8_,"P");
    gFinal4D->Add(hPoissonLin2,"L");
    gFinal4D->GetYaxis()->SetTitle("#frac{dN_{#phi#phi}}{dy}/(#frac{dN_{#phi}}{dy})^{2}");
    gFinal4D->GetXaxis()->SetLimits(0.85,1.5);
    gFinal4D->SetMaximum(1.6);
    gFinal4D->SetMinimum(0.4);
    TAxis *xAxis4  =   new TAxis(*gFinal4D->GetYaxis());
    
    gFinal4D->Draw("PA");
    lLegend1D->Draw("same");
    xAxis4->Draw();
    c____4->Write();
    c____4->SaveAs("./graphs/c____4.pdf");
    c____4->SaveAs("./graphs/c____4.png");
    
    ratioplot(h1D_Res,h1D_Tru_P6,h1D_Tru_P8,"1D");
    TH1D * hSlice6 = new TH1D ("hSlice6","hSlice6",nBinPT2D,fArrPT2D);
    TH1D * hSlice8 = new TH1D ("hSlice8","hSlice8",nBinPT2D,fArrPT2D);
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )
    {
        auto    hProjXPythia6   =   h2D_Tru_P6->ProjectionX(Form("P6X_%i",iTer),iTer+1,iTer+1);
        auto    hProjXPythia8   =   h2D_Tru_P8->ProjectionX(Form("P8X_%i",iTer),iTer+1,iTer+1);
        auto    hProjXDataset   =   h2D_Res->ProjectionX(Form("DTX_%i",iTer),iTer+1,iTer+1);
        auto    hProjYPythia6   =   h2D_Tru_P6->ProjectionY(Form("P6Y_%i",iTer),iTer+1,iTer+1);
        auto    hProjYPythia8   =   h2D_Tru_P8->ProjectionY(Form("P8Y_%i",iTer),iTer+1,iTer+1);
        auto    hProjYDataset   =   h2D_Res->ProjectionY(Form("DTY_%i",iTer),iTer+1,iTer+1);
        hSlice6->SetBinContent(iTer+1,hProjXPythia6->Integral("width"));
        hSlice8->SetBinContent(iTer+1,hProjXPythia8->Integral("width"));
        if ( iTer > 1 )
        {
            ratioplot(hProjXDataset,hProjXPythia6,hProjXPythia8,Form("2DX_%i",iTer+1));
            ratioplot(hProjYDataset,hProjYPythia6,hProjYPythia8,Form("2DY_%i",iTer+1));
        }
    }
    
    ratioplot(h1D_Raw,h1D_Rec_P6,h1D_Rec_P8,"1D_Raw");
    for ( Int_t iTer = 0; iTer < nBinPT2D; iTer++ )
    {
        auto    hProjXPythia6_   =   h2D_Rec_P6->ProjectionX(Form("P6X_%i_Raw",iTer),iTer+1,iTer+1);
        auto    hProjXPythia8_   =   h2D_Rec_P8->ProjectionX(Form("P8X_%i_Raw",iTer),iTer+1,iTer+1);
        auto    hProjXDataset_   =   h2D_Raw->ProjectionX(Form("DTX_%i_Raw",iTer),iTer+1,iTer+1);
        auto    hProjYPythia6_   =   h2D_Rec_P6->ProjectionY(Form("P6Y_%i_Raw",iTer),iTer+1,iTer+1);
        auto    hProjYPythia8_   =   h2D_Rec_P8->ProjectionY(Form("P8Y_%i_Raw",iTer),iTer+1,iTer+1);
        auto    hProjYDataset_   =   h2D_Raw->ProjectionY(Form("DTY_%i_Raw",iTer),iTer+1,iTer+1);
        if ( iTer > 1 )
        {
            ratioplot(hProjXDataset_,hProjXPythia6_,hProjXPythia8_,Form("2DX_%i_Raw",iTer+1));
            ratioplot(hProjYDataset_,hProjYPythia6_,hProjYPythia8_,Form("2DY_%i_Raw",iTer+1));
        }
    }
    
    
    TCanvas *cCompareFinal  =   new TCanvas();
    TH1D * hData__ = TH1DGeneratorYield(fResults2D,false);
    hData__->GetXaxis()->SetTitle("P_{T}#phi_{2} (GeV/c)");
    hData__->GetYaxis()->SetTitle("#frac{d^{2}N#phi#phi}{dydp_{T}#phi_{2}}(GeV/c)^{-1}");
    TLegend * lLegend1          =   new TLegend(0.65,0.65,0.85,0.85);
    lLegend1                    ->SetFillColor(kWhite);
    lLegend1                    ->SetLineColor(kWhite);
    lLegend1                    ->AddEntry(hSlice8,"Res","L");
    lLegend1                    ->AddEntry(hData__,"MC","EP");
    gPad->SetLogy();
    hSlice8->SetLineColor(kBlue);
    hSlice8->Draw("HIST L");
    hData__->SetMarkerStyle(25);
    hData__->SetMarkerColor(kRed);
    hData__->Draw("SAME EP");
    lLegend1->Draw("SAME");
    cCompareFinal->SaveAs("./graphs/cCompareFinal.pdf");
    ratioplot(fResults2D,hSlice6,hSlice8);
    
    //---------------------//
    // Output and wrap up  //-------------------------------------------------------------------------------
    //---------------------//
    
    // Output File for Results
    TFile*  outFile_RS  =   new TFile(fAnlResults,"recreate");
    
    h1D_Tru->Write();
    h1D_Res->Write();
    h2D_Tru->Write();
    h2D_Res->Write();
    
    // Closing opened files
    outFile_FT->Close();
    outFile_RS->Close();
    insFile_DT->Close();
    insFile_EF->Close();
}
