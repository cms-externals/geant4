
#include <iostream>
#include <sstream>
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
#include "TGraphErrors.h"

#include "NA61.h"
#include "ReadNA61Data.C"

#ifndef G4VAL_DRAWNA61SPECTRA_C
#define G4VAL_DRAWNA61SPECTRA_C


void drawMomSpectrum( std::string beam, std::string target, 
                      std::string secondary, 
		      std::string sec_angle_min, std::string sec_angle_max,
		      bool pub2015=false,
		      bool drawLegend=true,
                      bool largeLegend=false,
                      bool useStatEr=true, bool useSysEr=true )
{
   
   if ( pub2015 )
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max, "../test23/na61-exp-data/Pub-2015/" );
   }
   else
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   }
   
   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;
   
   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     float tmp = 0.;
     if ( useStatEr ) tmp += SecESigma[ii]*SecESigma[ii];
     if ( useSysEr )  tmp += SecESys[ii]*SecESys[ii];
      SumErr[ii] = std::sqrt( tmp );
   }

   for ( int ip=0; ip<NPointsNA61; ip++ )
   {
      if ( (SecSigma[ip]+SumErr[ip]) > ymax ) ymax = SecSigma[ip]+SumErr[ip];
      if ( (SecSigma[ip]-SumErr[ip]) < ymin ) ymin = SecSigma[ip]-SumErr[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NModels_HE];
  
   // TString YTitle = "d#sigma/dp, [mb/(GeV/c)]";
   TString YTitle = "production rate";
   
   for ( int m=0; m<NModels_HE; m++ )
   {
   
      std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = secondary;
      if ( pub2015 )
      {
         histoname += "_";
      }
      else
      {
         histoname += "Mult_"; 
      }
      histoname += ( sec_angle_min + "_" + sec_angle_max );
            
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      if ( !hi[m] )
      {
         std::cout << " can NOT find histogram = " << histoname << std::endl;
	 continue;
      }
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel_HE[m]);
      hi[m]->SetLineWidth(2);
      // FIXME !!!
      // hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	float yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 std::string htit = "31 GeV/c " + std::string( hi[m]->GetTitle() );
	 hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("p, GeV/c");
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetXaxis()->SetTitleSize(0.07);
	 if ( pub2015 )
	 {
	    //hi[m]->GetXaxis()->SetTitleSize(0.04);
	    hi[m]->GetYaxis()->SetTitle( "d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]" );
	    //hi[m]->GetYaxis()->SetTitleSize(0.04);
	    //hi[m]->GetYaxis()->SetTitleOffset(1.5);
	 }
	 else
	 {
	    hi[m]->GetYaxis()->SetTitle("production rate");
	 }
	 hi[m]->GetYaxis()->SetTitleSize(0.07);
	 hi[m]->GetYaxis()->SetTitleOffset(0.6);
//	 hi[m]->GetYaxis()->SetTitle("");
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
         hi[m]->Draw("histE1");
      }
      else hi[m]->Draw("histE1same");     
   }
   
   TLegend* leg = 0;
   if ( largeLegend )
   {   
      leg = new TLegend(0.45, 0.50, 0.99, 0.90);
      leg->SetTextSize(0.1);
   }
   else
   {
      leg = new TLegend( 0.45, 0.7, 0.99, 0.90);
      leg->SetTextSize(0.05);
//      leg = new TLegend( 0.4, 0.6, 0.99, 0.90);
//      leg->SetTextSize(0.055);
   }  
 
   for ( int m=0; m<NModels_HE; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName_HE[m].c_str(), "L" );
   }
         
//   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SecESigma );
   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SumErr );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   gr->GetYaxis()->SetTitle(YTitle);
//   gr->GetYaxis()->SetTitle("");
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   if ( drawLegend )
   {
      leg->Draw();
      leg->SetFillColor(kWhite);
   }

   delete [] SumErr;
   
   return ;

}

void drawMomSpectMC2Data( std::string beam, std::string target, 
                          std::string secondary, 
		          std::string sec_angle_min, std::string sec_angle_max,
			  bool pub2015=false,
			  bool drawLegend=true )
{

   if ( pub2015 )
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max, "../test23/na61-exp-data/Pub-2015/" );
   }
   else
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   }
   
   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;
   
   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
      float tmp = SecESigma[ii]*SecESigma[ii] + SecESys[ii]*SecESys[ii];
      SumErr[ii] = std::sqrt( tmp );
   }


   TH1F* hi[NModels_HE];
   
   float* MC2DataX = new float[NPointsNA61];
   float* MC2DataY = new float[NPointsNA61];
   float* DX = new float[NPointsNA61];
   float* DY = new float[NPointsNA61];

   int np = 0;
   
   TGraph* gr[NModels_HE];

   for ( int m=0; m<NModels_HE; m++ )
   {
   
      std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = secondary; //  + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      if ( pub2015 )
      {
         histoname += "_";
      }
      else
      {
         histoname += "Mult_"; 
      }
      histoname += ( sec_angle_min + "_" + sec_angle_max );
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );

      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {         
         float xx1 = hi[m]->GetBinLowEdge(k);
	 float xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsNA61; kk++ )
	 {
	    if ( xx1 < SecMom[kk] && xx1+xx2 > SecMom[kk] )
	    {
	       float yy  = hi[m]->GetBinContent(k);
	       float eyy = hi[m]->GetBinError(k);
	       MC2DataX[np] = SecMom[kk];
	       DX[np] = 0.;	       
	       MC2DataY[np] = yy / SecSigma[kk];
	       DY[np] = SumErr[kk]*MC2DataY[np] / SecSigma[kk];
	       if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	       if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	       np++;
	       break;
	    }
	 }

      }
      
      gr[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      int ibb = hi[m]->GetXaxis()->GetFirst();
      float xx1 = hi[m]->GetBinLowEdge(ibb);
      ibb = hi[m]->GetXaxis()->GetLast();
      float xx2 = hi[m]->GetBinLowEdge(ibb) + hi[m]->GetBinWidth(ibb);
      gr[m]->GetXaxis()->SetRangeUser(xx1,xx2);
      // gr[m]->SetTitle( hi[m]->GetTitle() );
      gr[m]->SetTitle( "" );
      gr[m]->SetMarkerColor(ColorModel_HE[m]);  
      gr[m]->SetMarkerStyle(SymbModel_HE[m]);
      gr[m]->SetMarkerSize(1.2);

   }
   
//   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SumErr );

   float* Value = new float[NPointsNA61];
   float* Error = new float[NPointsNA61];
   for ( int i=0; i<NPointsNA61; i++ )
   {
      Value[i] = 1.;
      Error[i] = SumErr[i] / SecSigma[i];
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }

   TGraph* gr1 = new TGraphErrors( NPointsNA61, SecMom, Value, 0, Error );

   TAxis* xaxis = gr1->GetXaxis();   
   xaxis->SetTitle("p, GeV/c");
   xaxis->SetTitleSize(0.07);
   xaxis->SetTitleFont(62);
   xaxis->SetLabelFont(62);
//   xaxis->SetTitleOffset(1.1);
//   xaxis->SetTitleOffset(1.);
   xaxis->SetTitleOffset(0.6);
   xaxis->CenterTitle();

   //gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   gr1->GetYaxis()->SetRangeUser( 0.1, 10. );
   gr1->SetMarkerColor(kBlack /* kBlue */ );
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
//   gr1->SetTitle(gr[0]->GetTitle());
   gr1->SetTitle("");
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetTitle("MC/Data (production rate)");
//   gr1->GetYaxis()->SetTitle("");
//   gr1->GetYaxis()->SetTitleOffset(1.1);
//   gr1->GetYaxis()->SetTitleOffset(0.9);
   gr1->GetYaxis()->SetTitleOffset(0.6);
   gr1->GetYaxis()->SetTitleSize(0.07);
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
   gr1->GetYaxis()->CenterTitle();

   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   for ( int m=0; m<NModels_HE; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }

   if ( drawLegend )
   {
      TLegend* leg = new TLegend(0.45, 0.7, 0.95, 0.95);
      leg->SetTextSize(0.045);
      for ( int m=0; m<NModels_HE; m++ )
      {
         leg->AddEntry( gr[m], ModelName_HE[m].c_str(), "p" );
      }
      leg->AddEntry( gr1, "exp.data", "p");
      leg->Draw();
      leg->SetFillColor(kWhite);
   }

   return;

}

void drawKPlus2PiPlusRatio( std::string beam, std::string target, 
		            std::string sec_angle_min, std::string sec_angle_max)
{

   readKPlus2PiPlusRatio( beam, target, sec_angle_min, sec_angle_max );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int ip=0; ip<NPointsNA61; ip++ )
   {
      if ( (K2PiRatio[ip]+K2PiERatio[ip]) > ymax ) ymax = K2PiRatio[ip]+K2PiERatio[ip];
      if ( (K2PiRatio[ip]-K2PiERatio[ip]) < ymin ) ymin = K2PiRatio[ip]-K2PiERatio[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NModels_HE];
  
   TString YTitle = "K^{+}/#pi^{+} ratio";
   
   for ( int m=0; m<NModels_HE; m++ )
   {
   
      std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = "kplus2piplus_" + sec_angle_min + "_" + sec_angle_max;
      std::cout << "histoname = " << histoname << std::endl;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel_HE[m]);
      hi[m]->SetLineWidth(2);
      // FIXME !!!
      // hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 hi[m]->GetXaxis()->SetTitle("p, GeV/c");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
         hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle(YTitle);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetTitleOffset(1.1);
         hi[m]->Draw();
      }
      else hi[m]->Draw("same");     
   }
   
   TLegend* leg = new TLegend(0.4, 0.70, 0.9, 0.9);
   leg->SetTextSize(0.04);
   
   for ( int m=0; m<NModels_HE; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName_HE[m].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, K2PiRatio, 0, K2PiERatio );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}

void drawKPlus2PiPlusRatioRegre( std::string beam, std::string target, 
		                 std::string sec_angle_min, std::string sec_angle_max,
			         std::string model )
{

/*
   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      // gROOT->LoadMacro("./shared-root-macros/DrawNA49Spectra.C");
      isNA61UtilLoaded = 1;
   }
*/   

   readKPlus2PiPlusRatio( beam, target, sec_angle_min, sec_angle_max );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int ip=0; ip<NPointsNA61; ip++ )
   {
      if ( (K2PiRatio[ip]+K2PiERatio[ip]) > ymax ) ymax = K2PiRatio[ip]+K2PiERatio[ip];
      if ( (K2PiRatio[ip]-K2PiERatio[ip]) < ymin ) ymin = K2PiRatio[ip]-K2PiERatio[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NVersions];
  
   TString YTitle = "K^{+}/#pi^{+} ratio";
   
   for ( int m=0; m<NVersions; m++ )
   {
   
//      std::string histofile = "";
//      if ( Versions[m] == CurrentVersion )
//      {
//         histofile = "./na61-histo/" + beam + target + "31.0GeV" + model;
//      }
//      else
//      {
//         histofile = Versions[m] + "/na61-histo/" + beam + target + "31.0GeV" + model;
//      }
//      histofile += ".root";

      std::string location = "";
      if ( Versions[m] == CurrentVersion )
      {
         location = ".";
      }
      else
      {
//         location = regre_test_dir + "/test19/" + Versions[m];
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[m];
      }
      std::string histofile = location + "/na61-histo/" + beam + target + "31.0GeV" + model + ".root"; 

      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = "kplus2piplus_" + sec_angle_min + "_" + sec_angle_max;
      std::cout << "histoname = " << histoname << std::endl;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorVersion[m]);
      hi[m]->SetLineWidth(2);
      // FIXME !!!
      // hi[m]->GetXaxis()->SetTitle("Kinetic energy of secondary neutron (MeV)");
      // hi[m]->GetYaxis()->SetTitle("Number of neutrons per MeV");
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 hi[m]->GetXaxis()->SetTitle("p, GeV/c");
	 hi[m]->GetXaxis()->SetTitleSize(0.07);
	 hi[m]->GetXaxis()->SetTitleOffset(0.6);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle(YTitle);
	 hi[m]->GetYaxis()->SetTitleSize(0.07);
	 hi[m]->GetYaxis()->SetTitleOffset(0.6);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetTitleOffset(1.1);
         hi[m]->Draw();
      }
      else hi[m]->Draw("same");     
   }
   
   TLegend* leg = new TLegend(0.4, 0.70, 0.9, 0.9);
   leg->SetTextSize(0.04);
   
   for ( int m=0; m<NVersions; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], Versions[m].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, K2PiRatio, 0, K2PiERatio );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);

   return ;

}

void drawMomSpectrumRegre( std::string beam, std::string target, 
                           std::string secondary, 
		           std::string sec_angle_min, std::string sec_angle_max,
			   std::string model,
			   bool pub2015=false,
			   bool drawLegend=true,
                           bool largeLegend=false,
                           bool useStatEr=true, bool useSysEr=true )
{

   // readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   
   if ( pub2015 )
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max, "../test23/na61-exp-data/Pub-2015/" );
   }
   else
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   }

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     float tmp = 0.;
     if ( useStatEr ) tmp += SecESigma[ii]*SecESigma[ii];
     if ( useSysEr )  tmp += SecESys[ii]*SecESys[ii];
      SumErr[ii] = std::sqrt( tmp );
   }

   for ( int ip=0; ip<NPointsNA61; ip++ )
   {
      if ( (SecSigma[ip]+SumErr[ip]) > ymax ) ymax = SecSigma[ip]+SumErr[ip];
      if ( (SecSigma[ip]-SumErr[ip]) < ymin ) ymin = SecSigma[ip]-SumErr[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
   
   TH1F* hi[NVersions];
  
   TString YTitle = "production rate";
   if ( pub2015 )
   {
      YTitle = "d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]";
   }

   for ( int iv=0; iv<NVersions; iv++ )
   {
      
/*
      std::string histofile = "";
      
      if ( Versions[iv] == CurrentVersion )
      {
         histofile = "./na61-histo/" + beam + target + "31.0GeV" + model;
      }
      else
      {
         histofile = Versions[iv] + "/na61-histo/" + beam + target + "31.0GeV" + model;
      }
      histofile += ".root";
      
      std::cout << " Regression, histofile: " << histofile << std::endl;
*/

      std::string location = "";
      if ( Versions[iv] == CurrentVersion )
      {
         location = ".";
      }
      else
      {
//         location = regre_test_dir + "/test19/" + Versions[iv];
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[iv];
      }
      std::string histofile = location + "/na61-histo/" + beam + target + "31.0GeV" + model + ".root"; 
      
      TFile* f = new TFile( histofile.c_str() );

      // std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      std::string histoname = secondary;
      if ( pub2015 )
      {
         histoname += "_";
      }
      else
      {
         histoname += "Mult_"; 
      }
      histoname += ( sec_angle_min + "_" + sec_angle_max );
      
      hi[iv] = (TH1F*)f->Get( histoname.c_str() );
      hi[iv]->SetStats(0);
      hi[iv]->SetLineColor(ColorVersion[iv]);
      hi[iv]->SetLineWidth(2);
      int nx = hi[iv]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[iv]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( iv == 0 ) 
      {
         std::string htit = "31 GeV/c " + std::string( hi[iv]->GetTitle() );
	 hi[iv]->SetTitle( htit.c_str() );
	 hi[iv]->Draw("histE1");
	 hi[iv]->GetXaxis()->SetTitle("p,GeV/c");
	 hi[iv]->GetXaxis()->SetTitleSize(0.07);
	 hi[iv]->GetXaxis()->SetTitleOffset(0.6);
	 hi[iv]->GetXaxis()->SetTitleFont(62);
	 hi[iv]->GetXaxis()->SetLabelFont(62);
	 hi[iv]->GetXaxis()->CenterTitle();
	 hi[iv]->GetYaxis()->SetTitle( YTitle );
	 hi[iv]->GetYaxis()->SetTitleSize(0.07);
	 hi[iv]->GetYaxis()->SetTitleOffset(0.6);
	 hi[iv]->GetYaxis()->SetTitleFont(62);
	 hi[iv]->GetYaxis()->SetLabelFont(62);
	 hi[iv]->GetYaxis()->CenterTitle();
	 hi[iv]->GetYaxis()->SetTitleOffset(0.6);
      }
      else hi[iv]->Draw("histE1same");     

   }
   
//   TLegend* leg = new TLegend(0.45, 0.70, 0.95, 0.9);
//   leg->SetTextSize(0.04);

   TLegend* leg = 0;
   if ( largeLegend )
   {   
      leg = new TLegend(0.45, 0.50, 0.99, 0.90);
      leg->SetTextSize(0.1);
   }
   else
   {
      leg = new TLegend( 0.45, 0.7, 0.99, 0.90);
      leg->SetTextSize(0.05);
   }  
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      hi[iv]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // hi[m]->SetTitle("");
      leg->AddEntry( hi[iv], Versions[iv].c_str(), "L" );
   }
      
   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SumErr );
   // gr->GetYaxis()->SetRangeUser( 0., 0.02 );
   // gr->SetRangeUser( ymin, ymax*1.5 );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");
   
   leg->AddEntry( gr, "exp.data", "p");

   if ( drawLegend )
   {
      leg->Draw();
      leg->SetFillColor(kWhite);
   }

   return ;

}

void drawMomSpectMC2DataRegre( std::string beam, std::string target, 
                               std::string secondary, 
		               std::string sec_angle_min, std::string sec_angle_max,
			       std::string model,
			       bool pub2015=false,
			       bool drawLegend=true,
                               bool largeLegend=false,
                               bool useStatEr=true, bool useSysEr=true )
{

   // readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );

   if ( pub2015 )
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max, "../test23/na61-exp-data/Pub-2015/" );
   }
   else
   {
      readMomSpectrum( beam, target, secondary, sec_angle_min, sec_angle_max );
   }
   
   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;
   
   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     float tmp = 0.;
     if ( useStatEr ) tmp += SecESigma[ii]*SecESigma[ii];
     if ( useSysEr )  tmp += SecESys[ii]*SecESys[ii];
      SumErr[ii] = std::sqrt( tmp );
   }

   TH1F* hi[NVersions];
   
   float* MC2DataX = new float[NPointsNA61];
   float* MC2DataY = new float[NPointsNA61];
   float* DX = new float[NPointsNA61];
   float* DY = new float[NPointsNA61];

   int np = 0;
   
   TGraph* gr[NVersions];

   for ( int m=0; m<NVersions; m++ )
   {
   
      std::string location = "";
      if ( Versions[m] == CurrentVersion )
      {
         location = ".";
      }
      else
      {
//         location = regre_test_dir + "/test19/" + Versions[m];
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[m];
      }
      std::string histofile = location + "/na61-histo/" + beam + target + "31.0GeV" + model + ".root"; 
      TFile* f = new TFile( histofile.c_str() );
      
      // std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      std::string histoname = secondary;
      if ( pub2015 )
      {
         histoname += "_";
      }
      else
      {
         histoname += "Mult_"; 
      }
      histoname += ( sec_angle_min + "_" + sec_angle_max );
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );

      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {         
         float xx1 = hi[m]->GetBinLowEdge(k);
	 float xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsNA61; kk++ )
	 {
	    if ( xx1 < SecMom[kk] && xx1+xx2 > SecMom[kk] )
	    {
	       float yy  = hi[m]->GetBinContent(k);
	       float eyy = hi[m]->GetBinError(k);
	       MC2DataX[np] = SecMom[kk];
	       DX[np] = 0.;	       
	       MC2DataY[np] = yy / SecSigma[kk];
	       DY[np] = SumErr[kk]*MC2DataY[np] / SecSigma[kk];
	       if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	       if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	       np++;
	       break;
	    }
	 }
	 
      }
      
      gr[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      int ibb = hi[m]->GetXaxis()->GetFirst();
      float xx1 = hi[m]->GetBinLowEdge(ibb);
      ibb = hi[m]->GetXaxis()->GetLast();
      float xx2 = hi[m]->GetBinLowEdge(ibb) + hi[m]->GetBinWidth(ibb);
      gr[m]->GetXaxis()->SetRangeUser(xx1,xx2);
      gr[m]->SetTitle( hi[m]->GetTitle() );
      gr[m]->SetMarkerColor(ColorVersion[m]);  
      gr[m]->SetMarkerStyle(SymbVersion[m]);
      gr[m]->SetMarkerSize(1.2);

   }
   
//   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SumErr );

   float* Value = new float[NPointsNA61];
   float* Error = new float[NPointsNA61];
   for ( int i=0; i<NPointsNA61; i++ )
   {
      Value[i] = 1.;
      Error[i] = SumErr[i] / SecSigma[i];
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }

   TGraph* gr1 = new TGraphErrors( NPointsNA61, SecMom, Value, 0, Error );

   TAxis* xaxis = gr1->GetXaxis();   
   xaxis->SetTitle("p, GeV/c");
   xaxis->SetTitleSize(0.07);
   xaxis->SetTitleFont(62);
   xaxis->SetLabelFont(62);
   xaxis->SetTitleOffset(0.6);
   xaxis->CenterTitle();

   //gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   gr1->GetYaxis()->SetRangeUser( 0.1, 10. );
   gr1->SetMarkerColor(kBlack /* kBlue */ );
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
   gr1->SetTitle(gr[0]->GetTitle());
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetTitle("MC/Data");
   gr1->GetYaxis()->SetTitleOffset(0.6);
   gr1->GetYaxis()->SetTitleSize(0.07);
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
   gr1->GetYaxis()->CenterTitle();

   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   for ( int m=0; m<NVersions; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }


   if ( drawLegend )
   {
      TLegend* leg = new TLegend(0.45, 0.70, 0.95, 0.9);
      leg->SetTextSize(0.04);
      for ( int m=0; m<NModels_HE; m++ )
      {
         leg->AddEntry( gr[m], ModelName_HE[m].c_str(), "p" );
      }
      leg->AddEntry( gr1, "exp.data", "p");
      leg->Draw();
      leg->SetFillColor(kWhite);
   }

   return;


}
#endif
