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

#include "../test23/shared-root-macros/NA49.h"
#include "../test23/shared-root-macros/ReadNA49Data.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA49.C"

// const int NModels = 2;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
// std::string ModelName[2]  = { "NuBeam-with-decays", "ftfp_bert-with-decays", "qgsp_bert-with-decays" };  
// std::string ModelName[5]  = { "NuBeam", "ftfp_bert", "qgsp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
const int NModels = 3;
std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue
int         SymbModel[4]     = { 8, 21, 23, 25 };

const int NVersions = 2;
// std::string Versions[4] = { "geant4-09-06-p02", "geant4-09-06-p03", "geant4-10-00", "geant4-10-00-p01" };
// std::string Versions[4] = { "geant4-10-00-p01", "geant4-10-00-ref03", "geant4-10-00-ref04", "geant4-10-00-ref05" };
// std::string Versions[4] = { "geant4-09-06-p03", "geant4-10-00-p01", "geant4-10-00-ref04", "geant4-10-01-beta-cand00" };
std::string Versions[2] = { "geant4-09-06-p03", "geant4-10-00-p02" }; 
int ColorVersion[6] = { kGreen, kRed, kMagenta, 7, kBlack, 14 };
int SymbVersion[5]     = { 20, 34, 21, 29, 23 };
std::string CurrentVersion = "";

static int isNA49UtilLoaded = 0;

void plot_pions_for_numi_talk()
{
     
   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }
      
   TCanvas* myc1 = new TCanvas("myc1", "", 1000, 600 );
   myc1->Divide(2,1);
   
   myc1->cd(1); 
   drawIntegratedSpectrum( "proton", "C", "piplus", "dNdxF" );
   myc1->cd(2); 
   drawIntegratedSpectrum( "proton", "C", "piminus", "dNdxF" );

   TCanvas* myc2 = new TCanvas("myc2", "", 1000, 600 );
   myc2->Divide(2,1);
   
   myc2->cd(1); 
   drawIntegratedSpectrum( "proton", "C", "piplus", "pT" );
   myc2->cd(2); 
   drawIntegratedSpectrum( "proton", "C", "piminus", "pT" );

   return;

}

void plot_protons_for_numi_talk()
{

   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas* myc = new TCanvas("myc", "", 1000, 600 );
   myc->Divide(2,1);
   
   myc->cd(1); 
   drawIntegratedSpectrum( "proton", "C", "proton", "dNdxF" );
   myc->cd(2); 
   drawIntegratedSpectrum( "proton", "C", "proton", "pT" );
   
   return;

}


void plot_dNdxF( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(3,2);

   myc->cd(1);
   drawIntegratedSpectrum( beam, target, "proton", "dNdxF" );
   
   myc->cd(2);
   drawIntegratedSpectrum( beam, target, "antiproton", "dNdxF" );

   myc->cd(3);
   drawIntegratedSpectrum( beam, target, "neutron", "dNdxF" );

   myc->cd(4);
   drawIntegratedSpectrum( beam, target, "piplus", "dNdxF" );

   myc->cd(5);
   drawIntegratedSpectrum( beam, target, "piminus", "dNdxF" );


   return;

}

void plot_pT( std::string beam, std::string target )
{

   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   TCanvas *myc = new TCanvas("myc","",800,800);
   myc->Divide(2,2);

   myc->cd(1);
   drawIntegratedSpectrum( beam, target, "proton", "pT" );
   
   myc->cd(2);
   drawIntegratedSpectrum( beam, target, "antiproton", "pT" );

   myc->cd(3);
   drawIntegratedSpectrum( beam, target, "piplus", "pT" );

   myc->cd(4);
   drawIntegratedSpectrum( beam, target, "piminus", "pT" );


   return;

}

void plot_dNdxF_pT( std::string beam, std::string target, std::string secondary )
{

   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   /*
    Here should be special provision for neutrons because only dNdxF data are available
    */

   TCanvas* myc   = new TCanvas("myc","",1000,600);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.12, 0.49, 0.99 );
   
   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, secondary, "dNdxF" ); 
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   if ( secondary == "proton" || secondary == "antiproton" )
   {
     drawIntSpectrumMC2Data( beam, target, secondary, "dNdxF", false ); // no legent for mc/data for (anti)protons
   }
   else
   {
      drawIntSpectrumMC2Data( beam, target, secondary, "dNdxF" );
   }

   myc->cd();
      
   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.12, 0.99, 0.99 );

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, secondary, "pT", false ); // no legend because it overlaps with the contents
   pad2->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawIntSpectrumMC2Data( beam, target, secondary, "pT" );
   
   myc->cd();

   TPad* pad3 = new TPad("pad3","",0.01, 0.01, 0.49, 0.11);
   pad3->Draw();
   // pad3->cd();
   myc->cd();
   TPad* pad4 = new TPad("pad4","",0.51, 0.01, 0.99, 0.11);
   pad4->Draw();
   // pad4->cd();
   for ( int m=0; m<NModels; ++m )
   {
      int NDF1=0;
      int NDF2 = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "dNdxF", ModelName[m], NDF1 );  
      std::ostringstream os1;
      os1 << (chi2/NDF1);
      std::string txt1 = " #chi^{2}/NDF = " + os1.str();
      txt1 += ( " for " + ModelName[m] + " vs NA49 Data" );
      TLatex* ltxt1 = new TLatex( 0.10, 0.75-m*0.3, txt1.c_str() );
      ltxt1->SetTextSize(0.25);
      pad3->cd();
      ltxt1->Draw();  
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "pT", ModelName[m], NDF2 );
      std::ostringstream os2;
      os2 << (chi2/NDF2);
      std::string txt2 = " #chi^{2}/NDF = " + os2.str();
      txt2 += ( " for " + ModelName[m] + " vs NA49 Data" );
      TLatex* ltxt2 = new TLatex( 0.10, 0.70-m*0.3, txt2.c_str() );
      ltxt2->SetTextSize(0.25);
      pad4->cd();
      ltxt2->Draw();
   }

   myc->cd();

   return;

} 

void plot_dNdxF_pT_Regre( std::string beam, std::string target, std::string secondary, std::string model )
{

   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("../test23/shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   /*
    Here should be special provision for neutrons because only dNdxF data are available
    */

   TCanvas* myc   = new TCanvas("myc","",1000,600);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.12, 0.49, 0.99 );
   
   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrumRegre( beam, target, secondary, "dNdxF", model ); 
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   if ( secondary == "proton" || secondary == "antiproton" )
   {
     drawIntSpectrumMC2DataRegre( beam, target, secondary, "dNdxF", model ); // , false ); // no legent for mc/data for (anti)protons
   }
   else
   {
      drawIntSpectrumMC2DataRegre( beam, target, secondary, "dNdxF", model );
   }

   myc->cd();
      
   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.12, 0.99, 0.99 );

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrumRegre( beam, target, secondary, "pT", model ); //, false ); // no legend because it overlaps with the contents
   pad2->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawIntSpectrumMC2DataRegre( beam, target, secondary, "pT", model );
   
   myc->cd();

   TPad* pad3 = new TPad("pad3","",0.01, 0.01, 0.49, 0.11);
   pad3->Draw();
   // pad3->cd();
   myc->cd();
   TPad* pad4 = new TPad("pad4","",0.51, 0.01, 0.99, 0.11);
   pad4->Draw();
   // pad4->cd();
   for ( int m=0; m<NVersions; ++m )
   {
      int NDF1=0;
      int NDF2 = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "dNdxF", model, NDF1, Versions[m] );  
      std::ostringstream os1;
      os1 << (chi2/NDF1);
      std::string txt1 = " #chi^{2}/NDF = " + os1.str();
      txt1 += ( " for " + Versions[m] + " vs NA49 Data" );
      TLatex* ltxt1 = new TLatex( 0.10, 0.75-m*0.3, txt1.c_str() );
      ltxt1->SetTextSize(0.25);
      pad3->cd();
      ltxt1->Draw();  
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "pT", model, NDF2, Versions[m] );
      std::ostringstream os2;
      os2 << (chi2/NDF2);
      std::string txt2 = " #chi^{2}/NDF = " + os2.str();
      txt2 += ( " for " + Versions[m] + " vs NA49 Data" );
      TLatex* ltxt2 = new TLatex( 0.10, 0.70-m*0.3, txt2.c_str() );
      ltxt2->SetTextSize(0.25);
      pad4->cd();
      ltxt2->Draw();
   }

   myc->cd();

   return;

} 

void plotDDiffXSec( std::string beam, std::string target, std::string secondary, int start=0, int end=22 )
{

   if ( isNA49UtilLoaded <= 0 )
   {
     //      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA49UtilLoaded = 1;
   }

   readDDiffSpectra( beam, target, secondary );

   int N1 = std::min(0,start);
   int N2 = std::max(end,NSubSets); 
   int NN = N2 - N1 ;
   int NN1 = 0;

   
   TCanvas** cnv = 0;
   if ( NN%2 == 0 )
   {
      NN1 = NN/2;
   }
   else
   {
      NN1 = NN/2 + 1;
   }

   cnv = new TCanvas*[NN1];

   for ( int i=0; i<NN1; ++i )
   {
      std::ostringstream cnt;
      cnt << i;
      std::string name = "cnv" + cnt.str();
      cnv[i] = new TCanvas( name.c_str(), "", 800, 500 );
      cnv[i]->Divide( 2, 1 );
   }
   
   for ( int i=0; i<NN; ++i )
   {
      cnv[i/2]->cd((i%2)+1);
      draw1DDiffXSec( beam, target, secondary, i );
   }

   return;

}

