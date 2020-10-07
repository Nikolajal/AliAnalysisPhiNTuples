#include "TF1.h"
#include <stdio.h>

using namespace std;

Double_t            fLevyFunc1D     ( Double_t * fVar, Double_t * fParams )
{
    Double_t    fPT     = fVar[0];
    Double_t    fMass   = fParams[0];
    Double_t    fEnne   = fParams[1];
    Double_t    fSlop   = fParams[2];
    Double_t    fdNdY   = fParams[3];
    
    Double_t    fNum1   = (fEnne-1)*(fEnne-2);
    Double_t    fDen1   = fEnne*fSlop*(fEnne*fSlop+fMass*(fEnne-2));
    Double_t    fFac1   = fNum1/fDen1;
    
    Double_t    fMasT   = sqrt(fMass*fMass+fPT*fPT);
    Double_t    fNum2   = fMasT - fMass;
    Double_t    fDen2   = fEnne*fSlop;
    Double_t    fFac2   = TMath::Power((1 + fNum2/fDen2),(-fEnne));
    
    return      fPT*fdNdY*fFac1*fFac2;
}

Double_t        fLevyFunc2D ( Double_t * fVar, Double_t * fParams )
{
    Double_t * fParx = new Double_t [4];
    Double_t * fPary = new Double_t [4];
    Double_t * fVarx = new Double_t [1];
    Double_t * fVary = new Double_t [1];
    
    fVarx = &fVar[0];
    fVary = &fVar[1];
    for ( int i = 0; i < 4; i++ )
    {
        fParx[i] = fParams[i];
        fPary[i] = fParams[i+4];
    }
    return  fPary[0]*(LevyTsallis_Func(fVarx,fParx)*LevyTsallis_Func(fVary,fParx));
}

Double_t        fFlatFuncXD ( Double_t * fVar, Double_t * fParams )
{
    return fParams[0];
}

TF1 *               fLevyFit1D      = new TF1 ("fLevyFunc1D",fLevyFunc1D,fMinPT1D,fMaxPT1D,4);

TF2 * fLevyFit2D = new TF2 ("fLevyFunc2D",fLevyFunc2D,fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,8);

TF1 * fFlatFit1D = new TF1 ("fLevyFunc1D",fFlatFuncXD,fMinPT1D,fMaxPT1D,1);

TF2 * fFlatFit2D = new TF2 ("fLevyFunc2D",fFlatFuncXD,fMinPT2D,fMaxPT2D,fMinPT2D,fMaxPT2D,1);
