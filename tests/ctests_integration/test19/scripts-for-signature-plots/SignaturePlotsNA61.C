// =======================
//
// Julia Yarba, FNAL
// Apr.5, 2015 
//
// As of now, we chosse the following to be one of our signature Geant4 plots:
// FTF(P) simulated results for 158GeV/c p+C -> pi+ compared vs NA49 data on <pT> vs xF   
//
// =======================

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

std::string TEST_NAME="test19";

#include "../test23/shared-root-macros/NA61.h"
#include "../test23/shared-root-macros/ReadNA61Data.C"
//
#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_HighEnergy.h"
//
// #include "Chi2Calc.C"
// #include "Chi2CalcNA61.C"
//
#include "../test23/shared-root-macros/DrawNA61Spectra.C"

void DrawMomSpect( std::string, std::string, std::string );
void DrawMC2Data( std::string, std::string, std::string, bool leg=false );


// NOTE: One can (re-)scale exp.data and MC results with 229.3mb
//       and preesnt results in [mb]

void SignaturePlotsNA61()
{

   TCanvas* myc1 = new TCanvas( "myc1", "", 900., 700. );
   
   readMomSpectrum( "proton", "C", "piplus", "0", "20" );
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     SecSigma[ii] *= 229.3;
     SecESigma[ii] *= 229.3;
     SecESys[ii] *= 229.3;
   }

   TPad* pad1 = new TPad( "pad1", "", 0.025, 0.44, 0.5, 0.9 );
   myc1->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetRightMargin(0.025);
   // drawMomSpectrumRegre( "proton", "C", "piplus", "0", "20", "ftfp", false);
   DrawMomSpect( "piplus", "0", "20" );
   
// -->    TLatex* theta1 = new TLatex( 1., 0.0053, "0 < #theta_{#pi^{+}} < 20 [mrad]" );
   TLatex* theta1 = new TLatex( 1., 1.2, "0 < #theta_{#pi^{+}} < 20 [mrad]" );
   theta1->SetTextSize(0.05);
   pad1->cd();
   theta1->Draw(); 
   
   TPad* pad2 = new TPad( "pad2", "", 0.025, 0.02, 0.5, 0.44 );
   myc1->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   DrawMC2Data( "piplus", "0", "20" );

   readMomSpectrum( "proton", "C", "piplus", "100", "140" );
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     SecSigma[ii] *= 229.3;
     SecESigma[ii] *= 229.3;
     SecESys[ii] *= 229.3;
   }

   TPad* pad3 = new TPad( "pad3", "", 0.525, 0.44, 1.0, 0.9 );
   myc1->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetRightMargin(0.025);
   DrawMomSpect( "piplus", "100", "140" );

//   TLatex* theta3 = new TLatex( 3.5, 0.060, "100 < #theta_{#pi^{+}} < 140 [mrad]" );
   TLatex* theta3 = new TLatex( 3.5, 13.8, "100 < #theta_{#pi^{+}} < 140 [mrad]" );
   theta1->SetTextSize(0.05);
   pad3->cd();
   theta3->Draw(); 
   
   TPad* pad4 = new TPad( "pad4", "", 0.525, 0.02, 1.0, 0.44 );
   myc1->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   DrawMC2Data( "piplus", "100", "140", true );   
   
   std::string data = "Data: N.Abgrall et al., Phys.Rev. C84, 034604 (2011)";
   TText* txt = new TText( 0.28, 0.91, data.c_str() );
   txt->SetTextSize(0.025);

   std::string htitle = "31 GeV/c proton + C #rightarrow #pi^{+} + X";
   TLatex* ltxt = new TLatex( 0.35, 0.95, htitle.c_str() );
   ltxt->SetTextSize(0.025);
   myc1->cd(0);
   ltxt->Draw();
   txt->Draw();
   
/*
   std::string cap = "Mean p_{T} as a function of x_{F} for #pi^{+} produced in p+C interactions at 158GeV/c.";   
   TLatex* lcap = new TLatex( 0.02, 0.065, cap.c_str() );
   lcap->SetTextSize(0.028);
   TText* cap1 = new TText( 0.05, 0.03, 
                 "Geant4 (FTF) results are compared with experimental data from NA49.");
   cap1->SetTextSize(0.028);
   
   TLatex* xtitle = new TLatex( 0.525, 0.085, "x_{F}" );
   xtitle->SetTextSize(0.034);
   
   myc1->cd();
   // lcap->Draw();
   // cap1->Draw();
   xtitle->Draw();
*/

   myc1->Print("proton-C-31GeV-piplus-ftfp-regre.gif");

}

void DrawMomSpect( std::string secondary, 
                   std::string sec_angle_min,
		   std::string sec_angle_max )
{

// -->   readMomSpectrum( "proton", "C", secondary, sec_angle_min, sec_angle_max );
   
   TGaxis::SetMaxDigits(2);
   
   std::string YTitle = "d#sigma/dp, [mb/(GeV/c)]";
//   std::string YTitle = "production rate";

   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     float tmp = 0.;
     tmp += SecESigma[ii]*SecESigma[ii];
     tmp += SecESys[ii]*SecESys[ii];
     SumErr[ii] = std::sqrt( tmp );
   }

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int ip=0; ip<NPointsNA61; ip++ )
   {
      if ( (SecSigma[ip]+SumErr[ip]) > ymax ) ymax = SecSigma[ip]+2.*SumErr[ip];
      if ( (SecSigma[ip]-SumErr[ip]) < ymin ) ymin = SecSigma[ip]-2.*SumErr[ip];
      if ( ymin < 0. ) ymin = 0.;
   }

   TH1F* hi[NVersions];
   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      std::string location = "";
      if ( Versions[iv] == CurrentVersion )
      {
         location = ".";
      }
      else
      {
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[iv];
      }
      
      std::string histofile = location + "/na61-histo/protonC31.0GeVftfp.root";       
      TFile* f = new TFile( histofile.c_str() );
 
      std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      hi[iv] = (TH1F*)f->Get( histoname.c_str() );

      hi[iv]->SetStats(0);
      hi[iv]->Scale( 229.3 );
      // hi[iv]->SetTitle("158 GeV/c proton on Carbon #rightarrow #pi^{+} + X" );
      hi[iv]->SetTitle("");
      hi[iv]->SetLineColor(ColorVersion[iv]);
      hi[iv]->SetLineWidth(6-iv);
      //hi[iv]->GetXaxis()->SetTitle("x_{F}");
      hi[iv]->GetXaxis()->SetTitle("");
      hi[iv]->GetXaxis()->SetLabelSize(0.04);
      hi[iv]->GetXaxis()->SetLabelFont(62);
      hi[iv]->GetYaxis()->SetTitle( YTitle.c_str() );
      hi[iv]->GetYaxis()->SetTitleSize(0.05);
      hi[iv]->GetYaxis()->SetTitleOffset(1.);
      hi[iv]->GetYaxis()->SetTitleFont(62);
      hi[iv]->GetYaxis()->SetLabelSize(0.04);
      hi[iv]->GetYaxis()->SetLabelFont(62);
      hi[iv]->GetYaxis()->CenterTitle();
      hi[iv]->GetYaxis()->SetRangeUser( ymin, ymax );

      if ( iv == 0 )
      {
         hi[iv]->Draw(); // we don't need to worry about "histE1" option because it's a profile histogram
      }
      else hi[iv]->Draw("same");
      
   }

   TGraph* gr = new TGraphErrors( NPointsNA61, SecMom, SecSigma, 0, SumErr );

   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
   
   gr->Draw("p");

   return;

}

void DrawMC2Data( std::string secondary, 
                  std::string sec_angle_min,
		  std::string sec_angle_max,
		  bool leg=false )
{

//   readMomSpectrum( "proton", "C", secondary, sec_angle_min, sec_angle_max );
   
   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;
   
   float* SumErr =  new float[NPointsNA61];
   for ( int ii=0; ii<NPointsNA61; ++ii )
   {
     float tmp = 0.;
     tmp += SecESigma[ii]*SecESigma[ii];
     tmp += SecESys[ii]*SecESys[ii];
     SumErr[ii] = std::sqrt( tmp );
   }

   TH1F* hi[NVersions];
   
   float* MC2DataX = new float[NPointsNA61];
   float* MC2DataY = new float[NPointsNA61];
   float* DX = new float[NPointsNA61];
   float* DY = new float[NPointsNA61];

   int np = 0;
   
   TGraph* gr[NVersions];
   
   float xx1 = -1.;
   float xx2 = -1.;

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
      std::string histofile = location + "/na61-histo/protonC31.0GeVftfp.root";       
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      hi[m]->Scale( 229.3 );

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
      xx1 = hi[m]->GetBinLowEdge(ibb);
      ibb = hi[m]->GetXaxis()->GetLast();
      xx2 = hi[m]->GetBinLowEdge(ibb) + hi[m]->GetBinWidth(ibb);
      gr[m]->GetXaxis()->SetRangeUser(xx1,xx2);
      gr[m]->SetTitle( hi[m]->GetTitle() );
      gr[m]->SetMarkerColor(ColorVersion[m]);  
      gr[m]->SetMarkerStyle(SymbVersion[m]);
      gr[m]->SetMarkerSize(1.2);

   }

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

   TGraphErrors* gr1 = new TGraphErrors( NPointsNA61, SecMom, Value, 0, Error );
   gr1->GetXaxis()->SetRangeUser( xx1, xx2 );

   TAxis* xaxis = gr1->GetXaxis();   
   xaxis->SetTitle("momentum of secondary #pi^{+} [GeV/c]");
   xaxis->SetTitleSize(0.05);
   xaxis->SetTitleFont(62);
   xaxis->SetLabelFont(62);
   xaxis->SetLabelSize(0.04);
   // xaxis->SetTitleOffset(0.6);
   xaxis->CenterTitle();

   gr1->GetYaxis()->SetRangeUser( 0.4, 2.5 );
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
   // gr1->SetTitle(gr[0]->GetTitle());
   gr1->SetTitle("");
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetTitle("MC/Data");
   // gr1->GetYaxis()->SetTitleOffset(0.6);
   gr1->GetYaxis()->SetTitleSize(0.05);
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
   gr1->GetYaxis()->CenterTitle();

   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");

   for ( int m=0; m<NVersions; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }

   if ( leg )
   {
      TLegend* legend = new TLegend( 0.45, 0.65, 0.99, 0.99 );
      legend->SetTextSize(0.05);
      legend->SetTextFont(62);
      for ( int m=0; m<NVersions; ++m )
      {
         legend->AddEntry( gr[m], Versions[m].c_str(), "p" );
      }
      legend->AddEntry( gr1, "exp.data", "p" );
      legend->Draw();
      legend->SetFillColor(kWhite);
   }

   return;

}
