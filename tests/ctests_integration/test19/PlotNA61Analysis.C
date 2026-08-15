#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"

int    NPoints    = 0;

float* SecMom     = 0;

float* SecSigma   = 0; 
float* SecESigma  = 0;

float* K2PiRatio  = 0;
float* K2PiERatio = 0;

const int NModels = 3;
std::string ModelName[3] = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
// int         ColorModel[2]    = { 14, 7 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel[3] = { kGreen, 7, kRed }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel[4]     = { 8, 21, 23, 25 };

const int NVersions = 4;
std::string Versions[4] = { "geant4-09-06-p03", "geant4-10-00-p01", "geant4-10-00-p02", "geant4-10-01-b01" };
// std::string Versions[4] = { "geant4-09-06-ref07", "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
// std::string Versions[2] = { "geant4-10-00-b01", "geant4-09-06-ref03" };
std::string CurrentVersion = "geant4-10-00-b01";
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };

static int isNA61UtilLoaded = 0;


void plotKPlus2PiPlusRatio(std::string beam, std::string target )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas* myc = new TCanvas("myc","",800,400);
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatio( beam, target, "20", "140" );

   myc->cd(2);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatio( beam, target, "140", "240" );

   return;

}

void plotKPlus2PiPlusRatioRegre(std::string beam, std::string target, std::string model )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas* myc = new TCanvas("myc","",800,400);
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatioRegre( beam, target, "20", "140", model );

   myc->cd(2);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatioRegre( beam, target, "140", "240", model );

   return;

}

void plotSecondarySum8( std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "0", "20" );

   myc1->cd(2);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "20", "40" );

   myc1->cd(3);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "40", "60" );

   myc1->cd(4);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "60", "100" );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "100", "140" );

   myc2->cd(2);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "140", "180" );

   myc2->cd(3);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "180", "240" );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "240", "300" );

   return;

}

void plotSecSumMC2Data( std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "0", "20" );

   myc1->cd(2);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "20", "40" );

   myc1->cd(3);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "40", "60" );

   myc1->cd(4);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "60", "100" );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "100", "140" );

   myc2->cd(2);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "140", "180" );

   myc2->cd(3);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "180", "240" );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectMC2Data( "proton", "C", secondary, "240", "300" );

   return;

}

void plotSecondarySumCombined( std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas* myc   = new TCanvas( "myc", "", 900, 900 );
   
   TPad* pad0 = new TPad( "pad0","",0.01,0.51, 0.49, 0.99 );
   TPad* pad1 = new TPad( "pad1","",0.51,0.51, 0.99, 0.99 );
   TPad* pad2 = new TPad( "pad2","",0.01,0.01, 0.49, 0.49 );
   TPad* pad3 = new TPad( "pad3","",0.51,0.01, 0.99, 0.49 );
   
   pad0->Draw();
   pad0->Divide(1.,2.,0.,0.);
   pad0->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "0", "20" );
   pad0->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "0", "20" );
   
   myc->cd();

   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "20", "40" );
   pad1->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "20", "40" );

   myc->cd();

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "40", "60" );
   pad2->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "40", "60" );

   myc->cd();
   
   pad3->Draw();
   pad3->Divide(1.,2.,0.,0.);
   pad3->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "60", "100" );
   pad3->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "60", "100" );

   TCanvas* myc1   = new TCanvas( "myc1", "", 900, 900 );

   TPad* pad4 = new TPad( "pad4","",0.01,0.51, 0.49, 0.99 );
   TPad* pad5 = new TPad( "pad5","",0.51,0.51, 0.99, 0.99 );
   TPad* pad6 = new TPad( "pad6","",0.01,0.01, 0.49, 0.49 );
   TPad* pad7 = new TPad( "pad7","",0.51,0.01, 0.99, 0.49 );
   
   pad4->Draw();
   pad4->Divide(1.,2.,0.,0.);
   pad4->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "100", "140" );
   pad4->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "100", "140" );
   
   myc1->cd();

   pad5->Draw();
   pad5->Divide(1.,2.,0.,0.);
   pad5->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "140", "180" );
   pad5->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "140", "180" );

   myc1->cd();

   pad6->Draw();
   pad6->Divide(1.,2.,0.,0.);
   pad6->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "180", "240" );
   pad6->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "180", "240" );
   
   if ( secondary == "proton" ) return;

   myc1->cd();
   
   pad7->Draw();
   pad7->Divide(1.,2.,0.,0.);
   pad7->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "240", "300" );
   pad7->cd(2); gPad->SetRightMargin(0.025);
   drawMomSpectMC2Data( "proton", "C", secondary, "240", "300" );

   return;

}


void plotSecondarySum8Regre( std::string secondary, std::string model )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA61Spectra.C");
      isNA61UtilLoaded = 1;
   }

   TCanvas *myc1 = new TCanvas("myc1","",900,900);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   gPad->SetLogy(); //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "0", "20", model );

   myc1->cd(2);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "20", "40", model );

   myc1->cd(3);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "40", "60", model );

   myc1->cd(4);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "60", "100", model );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);
   myc2->Divide(2,2);
   
   myc2->cd(1);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "100", "140", model );

   myc2->cd(2);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "140", "180", model );

   myc2->cd(3);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "180", "240", model );
   
   if ( secondary == "proton" ) return;

   myc2->cd(4);
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "240", "300", model );

   return;

}
