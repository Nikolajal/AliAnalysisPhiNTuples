// Global Values and constants file
// !TODO: All set!

#ifndef ALIANALYSISUTILITY_H
#define ALIANALYSISUTILITY_H

//------------------------------//
//      GLOBAL VARIABLES        //
//------------------------------//

// Random Generator
TRandom *   fRandomGen      =   new TRandom();

// Title and Name for histograms
auto        hName           =   "Name";
auto        hTitle          =   "Title";

//---------------------------------------//
//- Physics is by defualt in GeV, we    -//
//- define some constants to evaluate   -//
//- results in other units              -//
//---------------------------------------//

auto const  KeV             =   1e6;
auto const  MeV             =   1e3;
auto const  GeV             =   1;
auto const  TeV             =   1e-3;

#endif
