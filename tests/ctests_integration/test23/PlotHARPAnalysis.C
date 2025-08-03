#include <iostream>
#include <fstream>
#include <sstream>
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

#include "../test23/shared-root-macros/HARP.h"

//const int   NModels = 4;
const int   NModels = 2;
std::string ModelName[2]  = { "bertini", "ftfp" };
// std::string ModelName[4]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "qgsp_ftfp_bert" };
int         SymbModel[4]     = { 8, 21, 23, 25 };
int         ColorModel[4] = { kMagenta, 7, kRed, kBlack }; // 14 = grey, 7 = light "sky"-blue

const int NVersions = 2;
std::string Versions[2] = { "geant4-09-06-p03", "geant4-10-00-p02" }; 
int ColorVersion[6] = { kGreen, kRed, kMagenta, 7, kBlack, 14 };
int SymbVersion[5]     = { 20, 34, 21, 29, 23 };

static int isHARPUtilLoaded = 0;

void PlotHARPAnalysis( std::string beam, std::string target, std::string energy,
                       std::string secondary,
		       std::string region )
{

   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawHARPSpectra.C");
      isHARPUtilLoaded = 1;
   }

   ReadHARPData( beam, target, energy, secondary, region );
   
   TCanvas* myc = 0;
   
   if ( region == "FW" )
   {   
      myc = new TCanvas( "myc", "", 800, 800 );
      myc->Divide( 2, 2 );
   }
   else if ( region == "LA" )
   {
      myc = new TCanvas( "myc", "", 1200, 1200 );
      myc->Divide( 3, 3 );
   }

   for ( int i=0; i<NSetsHARP; i++ )
   {
      myc->cd(i+1);
      PlotHARPHisto( beam, target, energy, secondary, region, i );
   }
   
   myc->cd();
   
   return;

}

void PlotHARPAnalysisRegre( std::string beam, std::string target, std::string energy,
                       std::string secondary,
		       std::string region, 
		       std::string model )
{

   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawHARPSpectra.C");
      isHARPUtilLoaded = 1;
   }

   ReadHARPData( beam, target, energy, secondary, region );
   
   TCanvas* myc = 0;
   
   if ( region == "FW" )
   {   
      myc = new TCanvas( "myc", "", 800, 800 );
      myc->Divide( 2, 2 );
   }
   else if ( region == "LA" )
   {
      myc = new TCanvas( "myc", "", 1200, 1200 );
      myc->Divide( 3, 3 );
   }

   for ( int i=0; i<NSetsHARP; i++ )
   {
      myc->cd(i+1);
      PlotHARPHistoRegre( beam, target, energy, secondary, region, i, model );
   }
   
   myc->cd();
   
   return;

}

void PlotHARPAnaMC2Data(std::string beam, std::string target, std::string energy,
                       std::string secondary,
		       std::string region )
{

   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawHARPSpectra.C");
      isHARPUtilLoaded = 1;
   }

   ReadHARPData( beam, target, energy, secondary, region );
   
   TCanvas* myc = 0;
   
   if ( region == "FW" )
   {   
      myc = new TCanvas( "myc", "", 800, 800 );
      myc->Divide( 2, 2 );
   }
   else if ( region == "LA" )
   {
      myc = new TCanvas( "myc", "", 1200, 1200 );
      myc->Divide( 3, 3 );
   }

   for ( int i=0; i<NSetsHARP; i++ )
   {
      myc->cd(i+1);
      PlotHARPMC2Data( beam, target, energy, secondary, region, i );
   }
   
   myc->cd();
   
   return;

}

void PlotHARPAnalysisCombined( std::string beam, std::string target, std::string energy,
                               std::string secondary,
		               std::string region )
{

   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawHARPSpectra.C");
      isHARPUtilLoaded = 1;
   }

   ReadHARPData( beam, target, energy, secondary, region );
   
   TCanvas* myc = 0;
   
   float delta = 0.;
   int NX = 1;
   int NY = 1;
   
   if ( region == "FW" )
   {   
      myc = new TCanvas( "myc", "", 800, 800 );
      delta = 0.48;
      NX = 2;
      NY = 2;      
      // myc->Divide( 2, 2 );
   }
   else if ( region == "LA" )
   {
      myc = new TCanvas( "myc", "", 1200, 1200 );
      // myc->Divide( 3, 3 );
      delta = 0.31;
      NX = 3;
      NY = 3;
   }

   int iset = 0;
   for ( int i=0; i<NX; ++i )
   {
      float y2 = 1. - 0.01*(i+1) - (delta+0.01)*i;
      float y1 = y2 - delta;
      for ( int j=0; j<NY; ++j )
      {
	 std::ostringstream osCount;
	 osCount << iset;
         float x1 = 0.01 + delta*j + 0.02*j;
         float x2 = x1 + delta;
//         std::cout << " x1= " << x1 << " x2= " << x2 << " y1= " << y1 << " y2= " << y2 << std::endl;	 
         std::ostringstream osCount;
         osCount << iset;
         std::string pname = "pad" + osCount.str();
	 TPad* pad = new TPad("pname","",x1,y1,x2,y2);
         myc->cd();
	 pad->Draw();
         pad->Divide(1.,2.,0.,0.);
         pad->cd(1); gPad->SetRightMargin(0.025);
         PlotHARPHisto( beam, target, energy, secondary, region, iset, true, true );
         pad->cd(2); gPad->SetRightMargin(0.025);
         PlotHARPMC2Data( beam, target, energy, secondary, region, iset, false, false );
	 iset++;
      }   
   }
   
   
   
/*
   for ( int i=0; i<NSets; i++ )
   {
//      myc->cd(i+1);
      std::ostringstream osCount;
      osCount << i;
      std::string pname = "pad" + osCount.str();
      float x1 = 0.01*(i+1) + delta*i;
      float x2 = x1 + delta + 0.01; 
      
      std::cout << " x1= " << x1 << " x2= " << x2 << " y1= " << y1 << " y2= " << y2 << std::endl;
      
      TPad* pad = new TPad("pname","",x1,y1,x2,y2);
      pad->Draw();
      pad->Divide(1.,2.,0.,0.);
      pad->cd(1); gPad->SetRightMargin(0.025);
      PlotHARPMC2Data( beam, target, energy, secondary, region, i );
      pad->cd(2); gPad->SetRightMargin(0.025);
      PlotHARPMC2Data( beam, target, energy, secondary, region, i );

   }
*/   
   myc->cd();
   
   return;

}

