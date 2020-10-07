/*********************************************
 author: Roberto Preghenella (R+development)
         preghenella@bo.infn.it
*********************************************/

/*
  several function used for PbPb combined spectra
  Blast Wave is also implemented here
  further documentation will come
*/

#ifndef SpectraUtils_C
#define SpectraUtils_C

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include <Riostream.h>
#endif

Double_t integralPrecision = 1.e-12;

/*****************************************************************/
/* PT-EXPONENTIAL */
/*****************************************************************/

Double_t
PtExponential_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * TMath::Exp(-pt / T);
}

TF1 *
PtExponential(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
  
  TF1 *f = new TF1(name, PtExponential_Func, 0., 10., 3);
  f->SetParameters(mass, T, norm);
  f->SetParNames("mass", "T", "norm");
  f->FixParameter(0, mass);
  return f;
}

/*****************************************************************/
/* MT-EXPONENTIAL */
/*****************************************************************/

Double_t
MtExponential_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * TMath::Exp(-mt / T);
}

TF1 *
MtExponential(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
  
  TF1 *f = new TF1(name, MtExponential_Func, 0., 10., 3);
  f->SetParameters(mass, T, norm);
  f->SetParNames("mass", "T", "norm");
  f->FixParameter(0, mass);
  return f;
}

/*****************************************************************/
/* FERMI-DIRAC */
/*****************************************************************/

Double_t
FermiDirac_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm / (TMath::Exp(mt / T) + 1.);
}

TF1 *
FermiDirac(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
  
  TF1 *f = new TF1(name, FermiDirac_Func, 0., 10., 3);
  f->SetParameters(mass, T, norm);
  f->SetParNames("mass", "T", "norm");
  f->FixParameter(0, mass);
  return f;
}

/*****************************************************************/
/* BOLTZMANN */
/*****************************************************************/

Double_t
Boltzmann_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t T = p[1];
  Double_t norm = p[2];

  return pt * norm * mt * TMath::Exp(-mt / T);
}

TF1 *
Boltzmann(const Char_t *name, Double_t mass, Double_t T = 0.1, Double_t norm = 1.)
{
  
  TF1 *fBoltzmann = new TF1(name, Boltzmann_Func, 0., 10., 3);
  fBoltzmann->SetParameters(mass, T, norm);
  fBoltzmann->SetParNames("mass", "T", "norm");
  fBoltzmann->FixParameter(0, mass);
  return fBoltzmann;
}

/*****************************************************************/
/* LEVY-TSALLIS */
/*****************************************************************/

Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}

/*****************************************************************/
/* BOLTZMANN-GIBBS BLAST-WAVE */
/*****************************************************************/

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;
Double_t
BGBlastWave_Integrand(const Double_t *x, const Double_t *p)
{
  
  /* 
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
  */
  
  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];

  Double_t beta = beta_max * TMath::Power(r, n);
  if (beta > 0.9999999999999999) beta = 0.9999999999999999;
  Double_t rho = TMath::ATanH(beta);
  Double_t argI0 = pt * TMath::SinH(rho) * temp_1;
  if (argI0 > 700.) argI0 = 700.;
  Double_t argK1 = mt * TMath::CosH(rho) * temp_1;
  //  if (argI0 > 100 || argI0 < -100)
  //    printf("r=%f, pt=%f, beta_max=%f, temp=%f, n=%f, mt=%f, beta=%f, rho=%f, argI0=%f, argK1=%f\n", r, pt, beta_max, 1. / temp_1, n, mt, beta, rho, argI0, argK1);
  return r * mt * TMath::BesselI0(argI0) * TMath::BesselK1(argK1);
  
}

Double_t BGBlastWave_FuncN(const Double_t *x, const Double_t *p)
{
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  
  if (!fBGBlastWave_Integrand)
    fBGBlastWave_Integrand = new TF1("fBGBlastWave_Integrand", BGBlastWave_Integrand, 0., 1., 5);
  fBGBlastWave_Integrand->SetParameters(mt, pt, beta_max, temp, n);
  Double_t integral = fBGBlastWave_Integrand->Integral(0., 1., integralPrecision);
  return pt * integral;
}
  
TF1 *fBGBlastWave_FuncN = NULL;
Double_t saved_mass = -1., saved_beta_max = -1., saved_temp = -1., saved_n = -1;
Double_t saved_integral = -1.;
Double_t
BGBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t norm = p[4];
  
  if (!fBGBlastWave_FuncN)
    fBGBlastWave_FuncN = new TF1("fBGBlastWave_stuff", BGBlastWave_FuncN, 0., 10., 4);
  //  if (mass != saved_mass || beta_max != saved_beta_max || temp != saved_temp || n != saved_n) {
    fBGBlastWave_FuncN->SetParameters(mass, beta_max, temp, n);
    saved_mass = mass;
    saved_beta_max = beta_max;
    saved_temp = temp;
    saved_n = n;
    saved_integral = fBGBlastWave_FuncN->Integral(0., 10., integralPrecision);
    //  }
  Double_t value = fBGBlastWave_FuncN->Eval(pt);
  return norm * value / saved_integral;
}

TF1 *
BGBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t norm = 1.e6)
{
  
  TF1 *fBGBlastWave = new TF1(name, BGBlastWave_Func, 0., 10., 5);
  fBGBlastWave->SetParameters(mass, beta_max, temp, n, norm);
  fBGBlastWave->SetParNames("mass", "beta_max", "T", "n", "norm");
  fBGBlastWave->FixParameter(0, mass);
  fBGBlastWave->SetParLimits(1, 0., 1.);
  fBGBlastWave->SetParLimits(2, 1.e-6, 1.e6);
  fBGBlastWave->SetParLimits(3, 0.2, 5.);
  return fBGBlastWave;
}

/*****************************************************************/
/* TSALLIS BLAST-WAVE */
/*****************************************************************/

static TF1 *fTsallisBlastWave_Integrand_r = NULL;
Double_t
TsallisBlastWave_Integrand_r(const Double_t *x, const Double_t *p)
{
  /* 
     x[0] -> r (radius)
     p[0] -> mT (transverse mass)
     p[1] -> pT (transverse momentum)
     p[2] -> beta_max (surface velocity)
     p[3] -> T (freezout temperature)
     p[4] -> n (velocity profile)
     p[5] -> q
     p[6] -> y (rapidity)
     p[7] -> phi (azimuthal angle)
  */
  
  Double_t r = x[0];
  Double_t mt = p[0];
  Double_t pt = p[1];
  Double_t beta_max = p[2];
  Double_t temp_1 = 1. / p[3];
  Double_t n = p[4];
  Double_t q = p[5];
  Double_t y = p[6];
  Double_t phi = p[7];

  if (q <= 1.) return r;

  Double_t beta = beta_max * TMath::Power(r, n);
  Double_t rho = TMath::ATanH(beta);
  
  Double_t part1 = mt * TMath::CosH(y) * TMath::CosH(rho);
  Double_t part2 = pt * TMath::SinH(rho) * TMath::Cos(phi);
  Double_t part3 = part1 - part2;
  Double_t part4 = 1 + (q - 1.) * temp_1 * part3;
  Double_t expo = -1. / (q - 1.);
  //  printf("part1=%f, part2=%f, part3=%f, part4=%f, expo=%f\n", part1, part2, part3, part4, expo);
  Double_t part5 = TMath::Power(part4, expo);

  return r * part5;
}

static TF1 *fTsallisBlastWave_Integrand_phi = NULL;
Double_t
TsallisBlastWave_Integrand_phi(const Double_t *x, const Double_t *p)
{
  /* 
     x[0] -> phi (azimuthal angle)
  */
  
  Double_t phi = x[0];
  fTsallisBlastWave_Integrand_r->SetParameter(7, phi);
  Double_t integral = fTsallisBlastWave_Integrand_r->Integral(0., 1., integralPrecision);
  return integral;
}

static TF1 *fTsallisBlastWave_Integrand_y = NULL;
Double_t
TsallisBlastWave_Integrand_y(const Double_t *x, const Double_t *p)
{
  /* 
     x[0] -> y (rapidity)
  */

  Double_t y = x[0];
  fTsallisBlastWave_Integrand_r->SetParameter(6, y);
  Double_t integral = fTsallisBlastWave_Integrand_phi->Integral(-TMath::Pi(), TMath::Pi(), integralPrecision);
  return TMath::CosH(y) * integral;
}

Double_t
TsallisBlastWave_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t beta_max = p[1];
  Double_t temp = p[2];
  Double_t n = p[3];
  Double_t q = p[4];
  Double_t norm = p[5];

  if (!fTsallisBlastWave_Integrand_r)
    fTsallisBlastWave_Integrand_r = new TF1("fTsallisBlastWave_Integrand_r", TsallisBlastWave_Integrand_r, 0., 1., 8);
  if (!fTsallisBlastWave_Integrand_phi)
    fTsallisBlastWave_Integrand_phi = new TF1("fTsallisBlastWave_Integrand_phi", TsallisBlastWave_Integrand_phi, -TMath::Pi(), TMath::Pi(), 0);
  if (!fTsallisBlastWave_Integrand_y)
    fTsallisBlastWave_Integrand_y = new TF1("fTsallisBlastWave_Integrand_y", TsallisBlastWave_Integrand_y, -0.5, 0.5, 0);

  fTsallisBlastWave_Integrand_r->SetParameters(mt, pt, beta_max, temp, n, q, 0., 0.);
  Double_t integral = fTsallisBlastWave_Integrand_y->Integral(-0.5, 0.5, integralPrecision);
  return norm * pt * integral;
}

TF1 *
TsallisBlastWave(const Char_t *name, Double_t mass, Double_t beta_max = 0.9, Double_t temp = 0.1, Double_t n = 1., Double_t q = 1.1, Double_t norm = 1.e6)
{
  
  TF1 *fTsallisBlastWave = new TF1(name, TsallisBlastWave_Func, 0., 10., 6);
  fTsallisBlastWave->SetParameters(mass, beta_max, temp, n, q, norm);
  fTsallisBlastWave->SetParNames("mass", "beta_max", "T", "n", "q", "norm");
  fTsallisBlastWave->FixParameter(0, mass);
  fTsallisBlastWave->SetParLimits(1, 0.01, 0.99);
  fTsallisBlastWave->SetParLimits(2, 0.01, 1.);
  fTsallisBlastWave->SetParLimits(3, 0.1, 10.);
  fTsallisBlastWave->SetParLimits(4, 1.001, 1.2);
  return fTsallisBlastWave;
}

#endif /* SpectraUtils_C */
