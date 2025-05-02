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
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRefArray.h"
#include "TStyle.h"
#include "TGraph.h"

#include "TFileMerger.h"

#include "../test23/shared-root-macros/Chi2Calc.C"

double Chi2KESpectrumITEP( std::string, std::string, std::string, 
                           std::string, std::string,
		           std::string,
		           int&,
			   std::string,
		           bool, bool );

void plotMC2Data( std::string, std::string, std::string, 
                  std::string, std::string,
		  std::string model="bertini",
	          int NVer=4 );

void plotModelsMC2Data( std::string, std::string, std::string, 
                        std::string, std::string );

// options for std::string sec_angle are:
// principle: "59.1", "119.0"  
// also available: ...

// beam energy should be given as "1.40", or "5.00", or "7.50"... 

// read-in exp.data business
//
int NPtKE = 0;
float BeamEn, SecAng;
// there's no particular reason for having size=30 
// usually, the number of exp.points is less than that, but... just to be on the safe side...
float KE[30], ExpValue[30], YY[30], Err1[30], Err2[30];
float ErrStat[30], ErrSys[30]; 

// regression test business
//
/*
const int NVersions = 2;
std::string Versions[2] = { "geant4-09-06-p03", "geant4-10-01-b01" };
// std::string Versions[4] = { "geant4-09-04-ref10", "geant4-09-05", "geant4-09-05-p01", "geant4-09-05-ref02+MKtag"  };
// std::string Versions[4] = { "geant4-09-04-ref10", "geant4-09-05-ref01", "geant4-09-05-ref02+MKtag", "geant4-09-05-p01"  };
// int ColorVersion[4] = { kBlack, 7, kRed, kGreen }; // 7 = very light sky-blue
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };
*/

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

// model comparison business
//
const int NModels = 4;
// std::string Models[3] = { "bertini", "binary", "ftfp" };
// std::string Models[4] = { "bertini", "inclxx", "binary", "ftfp" };
std::string Models[4] = { "bertini", "ftfp", "ftfp_tune3", "fluka" };
// int ColorModel[6] = { 6, 3, 14 };
// int  ColorModel[5] = { kMagenta, kRed, kBlack, 7, 14 }; // 14 = grey, 7 = light "sky"-blue
int         ColorModel[6] = { kMagenta, kRed, kGreen, kBlue, 14, 7 }; // 14 = grey, 7 = light "sky"-blue

// --> General purspose exp.data read-in
//


#ifndef G4VAL_PLOTITEP_C
#define G4VAL_PLOTITEP_C


bool readKEData( std::string beam, std::string target, std::string energy, 
                 std::string secondary_type, std::string sec_angle,
		 bool useStatEr=true, bool useSysEr=true )
{
      
   // clear everything out
   //

   for ( int i=0; i<30; ++i )
   {
      KE[i] = 0.;
      ExpValue[i]=0.;
      YY[i] = 0.;
      Err1[i] = 0.;
      Err2[i] = 0.;
      ErrStat[i] = 0.;
      ErrSys[i] = 0.;
   }

   // read data
   //
   // try to find it first
   //
   std::string dirname = "../test47/itep";
   
   if ( gSystem->AccessPathName( dirname.c_str() ) ) // will return 1 (true) is the path does NOT exist
   {
      if ( gSystem->Getenv("TEST47DATA") )
      { 
         dirname = gSystem->Getenv("TEST47DATA");
      }
      else if ( gSystem->Getenv("G4INSTALL") )
      {
         dirname = gSystem->Getenv("G4INSTALL");
	 dirname += "/tests/test47/itep";
      }
   }
  
   std::cout << " dirname = " << dirname << std::endl;
   
   // std::string dirname += "/itep/" + beam + "/" + secondary_type + "/";
   dirname += "/" + beam + "/" + secondary_type + "/";
   std::string filename = target + energy + "GeV" + sec_angle + "deg.dat";
   std::string file = dirname + filename;
   std::cout << "reading from file: " << file << std::endl;
   
   ifstream infile;
   infile.open(file.c_str());
   
   if ( !infile.is_open() ) 
   {
      std::cout << " data file " << file << " does not exist " << std::endl;
      return false;
   }
   
   infile >> BeamEn >> SecAng >> NPtKE;
   
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      infile >> KE[ip] >> ExpValue[ip] >> ErrStat[ip] >> ErrSys[ip];
      ErrSys[ip] *= ExpValue[ip];
      YY[ip] = 1.;
      float totalerr2 = 0.;
      if ( useStatEr ) totalerr2 += ErrStat[ip]*ErrStat[ip] ;
      if ( useSysEr )  totalerr2 += ErrSys[ip]*ErrSys[ip];
      Err1[ip]  = sqrt(totalerr2);
      Err2[ip] = Err1[ip] / ExpValue[ip] ;
   }
   
   infile.close();

   return true;

}

TGraphErrors* getITEPAsGraph()
{

   double ymin = 10000.;
   double ymax = -1.;
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      if ( (ExpValue[ip]+Err1[ip]) > ymax ) ymax = ExpValue[ip]+Err1[ip];
      if ( (ExpValue[ip]-Err1[ip]) < ymin ) ymin = ExpValue[ip]-Err1[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
  
   TGraphErrors* gr = new TGraphErrors(NPtKE,KE,ExpValue,0,Err1);   
   
   gr->SetMinimum(ymin);
   gr->SetMaximum(ymax);

   return gr;

}


// chi2 calculation (integral)
//
void calcChi2KESpectrumITEP( std::string beam, std::string target, std::string energy, 
                             std::string secondary_type, 
		             std::string model,
			     std::string version = ".",
		             bool useStatEr=true, bool useSysEr=true )

{
   
   if ( (energy == "1.40" && model == "ftfp") || ( energy == "7.50" && model == "binary" ) )
   {
      std::cout << " Model " << model << " is not good/tested at " << energy << "GeV; giving up." << std::endl;
      return; 
   }

   int NDF=0;
   int tmpNDF = 0;
   double chi2 = 0.;
   
   chi2 += Chi2KESpectrumITEP( beam, target, energy, secondary_type, "59.1", model, tmpNDF, version, useStatEr, useSysEr );
   NDF += tmpNDF;
   tmpNDF=0;
   chi2 += Chi2KESpectrumITEP( beam, target, energy, secondary_type, "89.0", model, tmpNDF, version, useStatEr, useSysEr );
   NDF += tmpNDF;
   tmpNDF=0;
   chi2 += Chi2KESpectrumITEP( beam, target, energy, secondary_type, "119.0", model, tmpNDF, version, useStatEr, useSysEr );
   NDF += tmpNDF;
   tmpNDF=0;
   chi2 += Chi2KESpectrumITEP( beam, target, energy, secondary_type, "159.6", model, tmpNDF, version, useStatEr, useSysEr );
   NDF += tmpNDF;
   tmpNDF=0;
   
//   std::cout << " Chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for " << model << std::endl;
   
   return;

}

double Chi2KESpectrumITEP( std::string beam, std::string target, std::string energy, 
                           std::string secondary_type, std::string sec_angle,
		           std::string model,
		           int& NDF,
			   std::string version = ".",
		           bool useStatEr=true, bool useSysEr=true )
{

   double chi2 = 0.;

   bool status = readKEData( beam, target, energy, secondary_type, sec_angle, useStatEr, useSysEr );
   
   if ( !status )
   {
      std::cout << " Insufficient exp.data for " << beam << "+" << target << " --> " << secondary_type <<"+X" << std::endl;
      return 0.;
   } 
   
//   std::string histofile = version + "/";
//   histofile += ( beam + target + model + energy + "GeV.root" );

   std::string location = "";
   if ( version == CurrentVersion || version == "." )
   {
         location = "";
   }
   else
   {
         location = regre_test_dir + "/test47/" + version + "/";
   }
   std::string histofile = location + beam + target + model + energy + "GeV.root"; 

   std::string histoname = "KE" + secondary_type + "0" + beam + target + model + energy + "GeV";
   int counter = 5 - sec_angle.length();
   for ( int il=0; il<counter; il++ )
   {
      histoname += " "; // pad up for the fact that initially sec.angle was supposed to be char[6] - no more, no less...
   }
   histoname += sec_angle;
   
   TFile* hfile = new TFile( histofile.c_str() );
   
   std::cout << " Histo name : " << histoname << std::endl;
   TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
   
   if ( hi == 0 || NPtKE <= 0 ) 
   {
      std::cout << " Invalid case: no exp.data or simulation, or both " << std::endl;
      return 0.;
   }
   
   TGraphErrors* gdata = getITEPAsGraph();
   chi2 = Chi2( gdata, hi, NDF );
      
   // should I subtract 1DF ???
   // NDF-- ;
   
//   std::cout << " Chi2/NDF = " << Chi2 << "/" << NDF << " = " << Chi2/NDF << std::endl;
   
   return chi2;

}


void plotModelRegreSummary( std::string beam, std::string target, std::string model )
{

   double chi2 = 0.;
   int NDF = 0;
   
   std::string txt = model + " vs (ITEP) Data";
   TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
   ltxt->SetTextSize(0.2);
   
   if ( model != "ftfp" )
   {
   
      TCanvas *myc1 = new TCanvas("myc1","",800,800);

      TPad* pad1 = new TPad( "pad1", "", 0.01, 0.57, 0.49, 0.99 );
      TPad* pad2 = new TPad( "pad2", "", 0.01, 0.14, 0.49, 0.56 );
      TPad* pad3 = new TPad( "pad3", "", 0.51, 0.57, 0.99, 0.99 );
      TPad* pad4 = new TPad( "pad4", "", 0.51, 0.14, 0.99, 0.56 );

      myc1->cd();
      pad1->Draw();
      pad1->cd();
      gPad->SetLogy();
      plotMC2Data( beam, target, "1.40", "proton", "59.1", model ); //, 4 );

      myc1->cd();
      pad2->Draw();
      pad2->cd();
      gPad->SetLogy();
      plotMC2Data( beam, target, "1.40", "proton", "119.0", model ); //, 4 );
   
      myc1->cd();

      TPad* pad11 = new TPad("pad11","", 0.01, 0.01, 0.49, 0.13 );
      pad11->Draw();
      pad11->cd();
      
      //std::string txt = model + " vs (ITEP) Data";
      //TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
      //ltxt->SetTextSize(0.2);
      // --> ltxt->Draw();

      int icounter1 = 0;
      for ( int m=0; m<NVersions; ++m )
      {
         chi2 = 0.;
         NDF = 0;
         chi2 += Chi2KESpectrumITEP( beam, target, "1.40", "proton", "59.1",  model, NDF, Versions[m] );
         chi2 += Chi2KESpectrumITEP( beam, target, "1.40", "proton", "119.0", model, NDF, Versions[m] );
         std::ostringstream os;
         if ( (chi2/NDF) < 100 )
         {
            os.precision(3);
         }
         else
         {
            os.precision(4);
         }
         os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
         std::string txt1 = "#chi^{2}/NDF = ";
         txt1 += os.str();
         // ---> txt1 += ( " for " + Versions[m] );
         // ---> TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter1*0.2, txt1.c_str() );
         txt1 += ( "  " + Versions[m] );
         TLatex* ltxt1 = new TLatex(0.10, 0.8-icounter1*0.2, txt1.c_str() );
         ltxt1->SetTextSize(0.2);
	 ltxt1->SetTextColor( ColorVersion[m] );
         ltxt1->Draw();
         icounter1++;
      }

      myc1->cd();
      pad3->Draw();
      pad3->cd();
      gPad->SetLogy();
      plotMC2Data( beam, target, "1.40", "neutron", "59.1", model ); //,  4 );

      myc1->cd();
      pad4->Draw();
      pad4->cd();
      gPad->SetLogy();
      plotMC2Data( beam, target, "1.40", "neutron", "119.0", model ); //, 4 );
      
      myc1->cd();

      TPad* pad12 = new TPad("pad12","", 0.51, 0.01, 0.99, 0.13 );
      pad12->Draw();
      pad12->cd();
      
      // std::string txt = model + " vs (ITEP) Data";
      // TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
      // ltxt->SetTextSize(0.2);
      // ---> ltxt->Draw();

      icounter1 = 0;
      for ( int m=0; m<NVersions; ++m )
      {
         chi2 = 0.;
         NDF = 0;
         chi2 += Chi2KESpectrumITEP( beam, target, "1.40", "neutron", "59.1",  model, NDF, Versions[m] );
         chi2 += Chi2KESpectrumITEP( beam, target, "1.40", "neutron", "119.0", model, NDF, Versions[m] );
         if ( NDF < 2 )
         {
            TText* txt2 = new TText( 0.1, 0.4, "Insufficient dataset" );
	    txt2->SetTextSize(0.2);
	    txt2->Draw();
	    break;
         }
         std::ostringstream os;
         if ( (chi2/NDF) < 100 )
         {
            os.precision(3);
         }
         else
         {
            os.precision(4);
         }
         os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
         std::string txt1 = "#chi^{2}/NDF = ";
         txt1 += os.str();
         // ---> txt1 += ( " for " + Versions[m] );
         // ---> TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter1*0.2, txt1.c_str() );
         txt1 += ( "  " + Versions[m] );
         TLatex* ltxt1 = new TLatex(0.10, 0.8-icounter1*0.2, txt1.c_str() );
         ltxt1->SetTextSize(0.2);
	 ltxt1->SetTextColor( ColorVersion[m] );
         ltxt1->Draw();
         icounter1++;
      }

      myc1->cd();
   
      std::string output1 = beam + "-" + target + "-" + model + "-regre";
      if ( model == "bertini" || model=="inclxx" || model=="incl++" ) output1 += "-p1";
      output1 += ".gif";
   
      myc1->Print( output1.c_str() );
  
   }
      
   if ( model == "binary" ) return;

   std::string en;
   if ( beam == "proton" )
   {
      en = "7.50";
   }
   else if ( beam == "piplus" || beam == "piminus" )
   {
      en = "5.00";
   }
   else
   {
      std::cout << " Nothing available for " << beam << " beam " << std::endl;
      return;
   }

   TCanvas *myc2 = new TCanvas("myc2","",800,800);

   TPad* pad5 = new TPad( "pad5", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad6 = new TPad( "pad6", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad7 = new TPad( "pad7", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad8 = new TPad( "pad8", "", 0.51, 0.14, 0.99, 0.56 );


   myc2->cd();
   pad5->Draw();
   pad5->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "proton", "59.1", model ); //, 4 );

   myc2->cd();
   pad6->Draw();
   pad6->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "proton", "119.0", model ); //, 4 );
   
      myc2->cd();

      TPad* pad13 = new TPad("pad13","", 0.01, 0.01, 0.49, 0.13 );
      pad13->Draw();
      pad13->cd();
      
      //std::string txt11 = model + " vs (ITEP) Data";
      //TLatex* ltxt11 = new TLatex( 0.1, 0.8, txt.c_str() );
      //ltxt11->SetTextSize(0.2);
      // ---> ltxt->Draw();

      int icounter2 = 0;
      for ( int m=0; m<NVersions; ++m )
      {
         chi2 = 0.;
         NDF = 0;
         chi2 += Chi2KESpectrumITEP( beam, target, en, "proton", "59.1",  model, NDF, Versions[m] );
         chi2 += Chi2KESpectrumITEP( beam, target, en, "proton", "119.0", model, NDF, Versions[m] );
         std::ostringstream os;
         if ( (chi2/NDF) < 100 )
         {
            os.precision(3);
         }
         else
         {
            os.precision(4);
         }
         os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
         std::string txt1 = "#chi^{2}/NDF = ";
         txt1 += os.str();
         // ---> txt1 += ( " for " + Versions[m] );
         // ----> TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter2*0.2, txt1.c_str() );
         txt1 += ( "  " + Versions[m] );
         TLatex* ltxt1 = new TLatex(0.10, 0.8-icounter2*0.2, txt1.c_str() );
         ltxt1->SetTextSize(0.2);
	 ltxt1->SetTextColor( ColorVersion[m] );
         ltxt1->Draw();
         icounter2++;
      }

   myc2->cd();
   pad7->Draw();
   pad7->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "neutron", "59.1", model ); //,  4 );

   myc2->cd();
   pad8->Draw();
   pad8->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "neutron", "119.0", model ); //, 4 );

      myc2->cd();
      
      TPad* pad14 = new TPad("pad14","", 0.51, 0.01, 0.99, 0.13 );
      pad14->Draw();
      pad14->cd();
      // ---> ltxt->Draw();

      icounter2 = 0;
      for ( int m=0; m<NVersions; ++m )
      {
         chi2 = 0.;
         NDF = 0;
         chi2 += Chi2KESpectrumITEP( beam, target, en, "neutron", "59.1",  model, NDF, Versions[m] );
         chi2 += Chi2KESpectrumITEP( beam, target, en, "neutron", "119.0", model, NDF, Versions[m] );
         if ( NDF < 2 )
         {
            TText* txt2 = new TText( 0.1, 0.4, "Insufficient dataset" );
	    txt2->SetTextSize(0.2);
	    txt2->Draw();
	    break;
         }
         std::ostringstream os;
         if ( (chi2/NDF) < 100 )
         {
            os.precision(3);
         }
         else
         {
            os.precision(4);
         }
         os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
         std::string txt1 = "#chi^{2}/NDF = ";
         txt1 += os.str();
         // ---> txt1 += ( " for " + Versions[m] );
         // ---> TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter2*0.2, txt1.c_str() );
         txt1 += ( "  " + Versions[m] );
         TLatex* ltxt1 = new TLatex(0.10, 0.8-icounter2*0.2, txt1.c_str() );
         ltxt1->SetTextSize(0.2);
	 ltxt1->SetTextColor( ColorVersion[m] );
         ltxt1->Draw();
         icounter2++;
      }

      myc2->cd();

   std::string output2 = beam + "-" + target + "-" + model + "-regre";
   if ( model == "bertini" || model=="inclxx" || model=="incl++" ) output2 += "-p2";
   output2 += ".gif";
   
   myc2->Print( output2.c_str() );


   return;

}


// --> model comparison

void plotMC2Data( std::string beam, std::string target, std::string energy, 
                  std::string secondary_type, std::string sec_angle,
		  std::string model = "bertini",
	          int NVer=4 )
{
            
   bool status = readKEData( beam, target, energy, secondary_type, sec_angle );
   
   if ( !status )
   {
      TText* txt = new TText(0.25,0.6,"Insufficient exp.data");
      txt->SetTextSize(0.075);
      txt->Draw();
      return;
   }
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      if ( (YY[ip]+Err2[ip]) > ymax ) ymax = YY[ip]+Err2[ip];
      if ( (YY[ip]-Err2[ip]) < ymin ) ymin = YY[ip]-Err2[ip];
   }

   TGraph*  gr1 = new TGraphErrors(NPtKE,KE,YY,0,Err2);
   gr1->SetMarkerColor(kBlack /* 4 */ );  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   std::string xtitle = "Kinetic Energy of secondary " + secondary_type + " (GeV)";
   gr1->GetXaxis()->SetTitle( xtitle.c_str() );
   gr1->GetXaxis()->SetTitleFont(62);
   gr1->GetXaxis()->SetLabelFont(62);
   gr1->GetYaxis()->SetTitle("MC/Data");
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
   
   TGraph* gr[NVersions];
   std::string dir;
   if ( NVer > NVersions ) NVer = NVersions;
      
   for ( int iv=0; iv<NVer; iv++ )
   {

   // open G4 output root file (histo)
//   if ( iv == 0 ) 
//   { 
//      dir = "";
//   }
//   else
//   {
//      dir = Versions[iv] + "/";
//   }
//   std::string histofile = dir + beam + target + model + energy + "GeV.root";

   std::string location = "";
   if ( Versions[iv] == CurrentVersion )
   {
         location = "";
   }
   else
   {
         location = regre_test_dir + "/test47/" + Versions[iv] + "/";
   }
   std::string histofile = location + beam + target + model + energy + "GeV.root"; 

   std::string histoname = "KE" + secondary_type + "0" + beam + target + model + energy + "GeV";
   int counter = 5 - sec_angle.length();
   for ( int il=0; il<counter; il++ )
   {
      histoname += " "; // pad up for the fact that initially sec.angle was supposed to be char[6] - no more, no less...
   }
   histoname += sec_angle;
   
   TFile* hfile = new TFile( histofile.c_str() );
   
   TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
   
   //if ( iv == 1 || iv == 2 ) hi->Scale( (1./32.) );
   
   if ( hi == 0 || NPtKE <= 0 ) 
   {
      std::cout << " Invalid case for " << model << ": no exp.data or simulatiion, or both " << std::endl;
      return ;
   }
      
   
   float MC2DataX[30];
   float MC2DataY[30];
   float DX[30], DY[30];
   int np =0;
   
   int nx = hi->GetNbinsX();
   

   for ( int k=1; k <= nx; k++ )
   {        
        double xx1 = hi->GetBinLowEdge(k);
	double xx2 = hi->GetBinWidth(k);
	for (int kk=0; kk<NPtKE; kk++ )
	{
	   if ( xx1 < KE[kk] && xx1+xx2 > KE[kk] )
	   {
	      double yy = hi->GetBinContent(k);
	      MC2DataX[np] = KE[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / ExpValue[kk];
	      // also error calc here !...
	      DY[np]=Err1[kk]*MC2DataY[np]/ExpValue[kk];
	      if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	      if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	      np++;
	      break;
           }
	}
   }
   
   gr[iv] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
   if ( iv == 0 )
   {
      std::string newtitle = beam + "+" + target + " -> " + secondary_type;
      newtitle += ( " at " + energy + "GeV, theta=" + sec_angle );
      gr[iv]->SetTitle(newtitle.c_str());
      gr1->SetTitle(newtitle.c_str());
//      gr[iv]->SetTitle(hi->GetTitle());
//      gr1->SetTitle(hi->GetTitle());
   }
   else
   {
      gr[iv]->SetTitle("");
   }
   gr[iv]->GetXaxis()->SetTitle( xtitle.c_str() );
   gr[iv]->GetYaxis()->SetTitle("MC/Data");
   gr[iv]->SetMarkerColor(ColorVersion[iv]); 
   gr[iv]->SetLineColor(ColorVersion[iv]); 
   gr[iv]->SetMarkerStyle(21);
   gr[iv]->SetMarkerSize(1.6);
   
   } // end loop on versions
      
   
   if ( ymin > 0. ) ymin = 0;
   if ( ymax < 3. ) ymax = 3.;
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   // gr1->GetYaxis()->SetRangeUser(ymin, ymax); 
   gr1->GetYaxis()->SetRangeUser( 0.2, 6. ); // 0.1, 30. ); 
   
   gr1->Draw("apl");

   // ---> TLegend* leg = new TLegend(0.65, 0.65, 0.99, 0.9);   
   TLegend* leg = new TLegend(0.6, 0.85, 0.90, 0.90);   
   
   for ( int iv=0; iv<NVer; iv++ )
   {
      gr[iv]->Draw("lpsame");       
      // ---> leg->AddEntry( gr[iv], Versions[iv].c_str(), "p" );  
   }
    
   leg->AddEntry( gr1, "ITEP771 data", "p");     

   leg->Draw();
   leg->SetFillColor(kWhite);
         
   return;

}


// --> models comparison


void plotForTalk()
{
   
   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(4,2);

   myc->cd(1);
   plotModelsMC2Data( "proton", "C", "7.50", "proton", "59.1" );
   myc->cd(2);
   plotModelsMC2Data( "proton", "U", "7.50", "proton", "59.1" );
   myc->cd(1);
   plotModelsMC2Data( "piminus", "C", "5.00", "proton", "59.1" );
   myc->cd(2);
   plotModelsMC2Data( "piminus", "U", "5.00", "proton", "59.1" );

   myc->cd(5);
   plotModelsMC2Data( "proton", "C", "7.50", "proton", "119.0" );
   myc->cd(6);
   plotModelsMC2Data( "proton", "U", "7.50", "proton", "119.0" );
   myc->cd(7);
   plotModelsMC2Data( "piminus", "C", "5.00", "proton", "119.0" );
   myc->cd(8);
   plotModelsMC2Data( "piminus", "U", "5.00", "proton", "119.0" );
  
   return;
}

// --> summary plots
//
void plotModelsMC2DataSummary( std::string beam, std::string target, std::string energy )
{

   double chi2 = 0.;
   int    NDF = 0;
      
   TCanvas *myc1 = new TCanvas("myc1","",800,800);
//   myc1->Divide(2,2);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.57, 0.49, 0.99 );
   //myc1->cd(1);
   myc1->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   plotModelsMC2Data( beam, target, energy, "proton", "59.1" );
   
   TPad* pad2 = new TPad("pad2","", 0.01, 0.14, 0.49, 0.56 );
   // myc1->cd(2);   
   myc1->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   plotModelsMC2Data( beam, target, energy, "proton", "119.0" );

   myc1->cd();

   TPad* pad5 = new TPad("pad5","", 0.01, 0.01, 0.49, 0.13 );
   pad5->Draw();
   pad5->cd();
  
//   TText* txt = new TText( 0.1, 0.8, "MonteCarlo vs (ITEP) Data" );
//   txt->SetTextSize(0.2);
//   txt->Draw();

   int icounter = 0;
   for ( int m=0; m<NModels; ++m )
   {
      if ( ( energy == "1.40" && Models[m] == "ftfp" )   || 
           ( energy == "7.50" && Models[m] == "binary" ) || 
	   ( energy == "5.00" && Models[m] == "binary" ) )
      {
         continue; 
      }
      chi2 = 0.;
      NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, target, energy, "proton", "59.1",  Models[m], NDF );
      chi2 += Chi2KESpectrumITEP( beam, target, energy, "proton", "119.0", Models[m], NDF );
      std::ostringstream os;
      if ( (chi2/NDF) < 100 )
      {
         os.precision(3);
      }
      else
      {
         os.precision(4);
      }
      os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      if ( Models[m].find("flu") != std::string::npos )
      {
         txt1 += ( "  fluka.cern v4.4.0" );
      }
      else
      {
         txt1 += ( "  " + Models[m] );
      }
//      txt1 += ( " for " + Models[m] );
//      if ( Models[m].find("flu") != std::string::npos )
//      {
//         txt1 += "4.4.0";
//      }
      // --> TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter*0.2, txt1.c_str() );
      TLatex* ltxt1 = new TLatex(0.10, 0.8-icounter*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.2);
      ltxt1->SetTextColor( ColorModel[m] );
      ltxt1->Draw();
      icounter++;
   }
   
   
   TPad* pad3 = new TPad("pad3","",0.51,0.57,0.99,0.99);   
//   myc1->cd(3);
   myc1->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetLogy();
   plotModelsMC2Data( beam, target, energy, "neutron", "59.1" );
      
   TPad* pad4 = new TPad("pad4","", 0.51, 0.14, 0.99, 0.56);   
   // myc1->cd(4);
   myc1->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetLogy();
   plotModelsMC2Data( beam, target, energy, "neutron", "119.0" );

   myc1->cd();

   TPad* pad6 = new TPad("pad6","", 0.51, 0.01, 0.99, 0.13 );
   pad6->Draw();
   pad6->cd();
  
   // TText* txt = new TText( 0.1, 0.8, "MonteCarlo vs (ITEP) Data" );
   // txt->SetTextSize(0.2);
//   txt->Draw();

   icounter=0;
   for ( int m=0; m<NModels; ++m )
   {
      if ( ( energy == "1.40" && Models[m] == "ftfp" )   || 
           ( energy == "7.50" && Models[m] == "binary" ) ||
	   ( energy == "5.00" && Models[m] == "binary" ) )
      {
         continue; 
      }
      chi2 = 0.;
      NDF = 0;
      chi2 += Chi2KESpectrumITEP( beam, target, energy, "neutron", "59.1",  Models[m], NDF );
      chi2 += Chi2KESpectrumITEP( beam, target, energy, "neutron", "119.0", Models[m], NDF );
      
      if ( NDF < 2 )
      {
         TText* txt2 = new TText( 0.1, 0.4, "Insufficient dataset" );
	 txt2->SetTextSize(0.2);
	 txt2->Draw();
	 break;
      }
      
      std::ostringstream os;
      if ( (chi2/NDF) < 100 )
      {
         os.precision(3);
      }
      else
      {
         os.precision(4);
      }
      os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      if ( Models[m].find("flu") != std::string::npos )
      {
         txt1 += ( "  fluka.cern v4.4.0" );
      }
      else
      {
         txt1 += ( "  " + Models[m] );
      }
      if ( Models[m].find("flu") != std::string::npos )
      {
         txt1 += "4.4.0";
      }
      // --> TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter*0.2, txt1.c_str() );
      TLatex* ltxt1 = new TLatex(0.10, 0.8-icounter*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.2);
      ltxt1->SetTextColor(ColorModel[m]);
      ltxt1->Draw();
      icounter++;
   }
   
   myc1->cd();
   
   std::string output1 = beam + "-" + target + "-" + energy + "GeV-models.gif";
   myc1->Print( output1.c_str() );

   return;

}

// -- > basic macro
//
void plotModelsMC2Data( std::string beam, std::string target, std::string energy, 
                        std::string secondary_type, std::string sec_angle )
{
            
   bool status = readKEData( beam, target, energy, secondary_type, sec_angle );

   if ( !status )
   {
      TText* txt = new TText(0.05,0.8,"No exp. data available");
      txt->SetTextSize(0.06);
      txt->Draw();
      return;
   }
      
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      if ( (YY[ip]+Err2[ip]) > ymax ) ymax = YY[ip]+Err2[ip];
      if ( (YY[ip]-Err2[ip]) < ymin ) ymin = YY[ip]-Err2[ip];
   }

   TGraph*  gr1 = new TGraphErrors(NPtKE,KE,YY,0,Err2);
   gr1->SetMarkerColor(kBlack /* 4 */ );  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   std::string xtitle = "Kinetic Energy of secondary " + secondary_type + " (GeV)";
   gr1->GetXaxis()->SetTitle( xtitle.c_str() );
   gr1->GetXaxis()->SetTitleFont(62);
   gr1->GetXaxis()->SetLabelFont(62);
   gr1->GetYaxis()->SetTitle("MC/Data");
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
   
   TGraph* gr[NModels];

   std::string skip ="";
   if ( energy == "1.40" ) 
   {
      skip = "ftfp";
   }
   else 
   {
      skip = "binary";
   }  
   
   for ( int iv=0; iv<NModels; iv++ )
   {

      if ( Models[iv] == skip ) continue;
      
   // open G4 output root file (histo)
   std::string histofile = beam + target + Models[iv] + energy + "GeV.root";
   std::string histoname = "KE" + secondary_type + "0" + beam + target + Models[iv] + energy + "GeV";
   int counter = 5 - sec_angle.length();
   for ( int il=0; il<counter; il++ )
   {
      histoname += " "; // pad up for the fact that initially sec.angle was supposed to be char[6] - no more, no less...
   }
   histoname += sec_angle;
   
   TFile* hfile = new TFile( histofile.c_str() );
   
//   std::cout << " Histo name : " << histoname << std::endl;
   TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
   
   if ( hi == 0 || NPtKE <= 0 ) 
   {
      std::cout << " Invalid case: no exp.data or simulatiion, or both " << std::endl;
      return ;
   }
   
   //hi->Scale( (1./32.) );
      
   float MC2DataX[30];
   float MC2DataY[30];
   float DX[30], DY[30];
   int np =0;
   
   int nx = hi->GetNbinsX();
   

   for ( int k=1; k <= nx; k++ )
   {        
        double xx1 = hi->GetBinLowEdge(k);
	double xx2 = hi->GetBinWidth(k);
	for (int kk=0; kk<NPtKE; kk++ )
	{
	   if ( xx1 < KE[kk] && xx1+xx2 > KE[kk] )
	   {
	      double yy = hi->GetBinContent(k);
	      MC2DataX[np] = KE[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / ExpValue[kk];
	      // also error calc here !...
	      DY[np]=Err1[kk]*MC2DataY[np]/ExpValue[kk];
	      if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	      if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	      np++;
	      break;
           }
	}
   }
   
   // do NOT plot if less than 2 data points
   //
   
   gr[iv] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
   if ( iv == 0 )
   {
      std::string newtitle = beam + "+" + target + " -> " + secondary_type;
      newtitle += ( " at " + energy + "GeV, theta=" + sec_angle );
      gr[iv]->SetTitle(newtitle.c_str());
      gr1->SetTitle(newtitle.c_str());
//      gr[iv]->SetTitle(hi->GetTitle());
//      gr1->SetTitle(hi->GetTitle());
   }
   else
   {
      gr[iv]->SetTitle("");
   }
   gr[iv]->GetXaxis()->SetTitle( xtitle.c_str() );
   gr[iv]->GetXaxis()->SetTitleFont(62);
   gr[iv]->GetXaxis()->SetLabelFont(62);
   gr[iv]->GetYaxis()->SetTitle("MC/Data");
   gr[iv]->GetYaxis()->SetTitleFont(62);
   gr[iv]->GetYaxis()->SetLabelFont(62);
   gr[iv]->SetMarkerColor(ColorModel[iv]); 
   gr[iv]->SetLineColor(ColorModel[iv]); 
   gr[iv]->SetMarkerStyle(21);
   gr[iv]->SetMarkerSize(1.0); // mrk size used to be 1.6
   
   } // end loop on models
      
   
   if ( ymin > 0. ) ymin = 0;
   if ( ymax < 3. ) ymax = 3.;
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   // gr1->GetYaxis()->SetRangeUser(ymin, ymax); 
   gr1->GetYaxis()->SetRangeUser( 0.2, 5.); // 0.1, 30. ); 
   
   gr1->Draw("apl");
   
   // --> TLegend* leg = new TLegend(0.65, 0.65, 0.99, 0.90);   
   TLegend* leg = new TLegend(0.6, 0.85, 0.90, 0.90);   
   
   for ( int iv=0; iv<NModels; iv++ )
   {
      if ( Models[iv] == skip ) continue;
      gr[iv]->Draw("lpsame");   
//       leg->AddEntry( gr[iv], Models[iv].c_str(), "p" );
   }     

   leg->AddEntry( gr1, "ITEP771 data", "p"); // "exp.data", "p");     

   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}

/*
void crudeMerge( std::string beam, std::string target, std::string energy, std::string model )
{

// Note: if one merges weighted/normilized histograms, in this case 
//       the resulting histogram(s) will be 32 times the statistics;
//       at the analysis stage they'd need to be scaled by 1/32 

   TFileMerger* fm = new TFileMerger();
   
   std::string output = beam + target + model + energy + "GeV.root";
   
   fm->OutputFile( output.c_str() );
   
   for ( int id=1; id<=32; id++ )
   {
      std::string filename = beam + target + model + energy + "GeV-";
      char buf[5];
      sprintf( buf, "%i", id );
      filename.append( buf ); 
      filename += ".root";    
      // std::cout << " file name = " << file << std::endl;           
      fm->AddFile( filename.c_str() );
   }
   
   fm->Merge();
   
   return;

}

void fancyMerge( std::string beam, std::string target, std::string energy, std::string model, bool doScale=false )
{
      
   std::string output = beam + target + model + energy + "GeV.root" ;
   
   targetFile = TFile::Open( output.c_str(), "RECREATE" );
   
   double scale = 1./32.;
   
   std::string input = beam + target + model + energy + "GeV-1.root";
   
   TFile* iFile1 = TFile::Open( input.c_str() );
   TIter  next( iFile1->GetListOfKeys() );
   TKey*  key = (TKey*)next();
   TH1* h = 0;
   while ( key )
   {   
         if ( !(TClass::GetClass(key->GetClassName())->InheritsFrom(TH1::Class())) ) continue;
         const char* kName = key->GetName();
         h = (TH1*)key->ReadObj();
         const char* hName = h->GetName();
         std::cout << " histoname = " << hName << std::endl;
	 TH1F* h1 = h->Clone();
	 for ( int id=2; id<=32; id++ )
	 {
	    std::string input_t = beam + target + model + energy + "GeV-" ;
            char buf[5];
            sprintf( buf, "%i", id );
            input_t.append( buf ); 
            input_t += ".root"; 
	    TFile* iFile_t = TFile::Open( input_t.c_str() );
	    TH1F* h_t = (TH1F*)iFile_t->Get( h->GetName() );
	    h1->Add( h_t );  
	    iFile_t->Close();
	 }
	 if ( doScale )
	 {
	    h1->Scale( scale );
	 }
	 targetFile->cd();
	 h1->Write();
         key = (TKey*)next();
   }
   
   targetFile->Close();
     
   return;

}
*/
#endif



