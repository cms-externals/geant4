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
#include "TGraph.h"

#include "NA49.h"
#include "ReadNA49Data.C"

// #include "DbReader.h" // <--- will appear in NA49Models.C
#include "TGraphAsymmErrors.h"

//const int NModels = 2;
// std::string ModelName[4]  = { "qgsp_ftfp_bert", "ftfp_bert", "ftfp", "qgsp" };  
//std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
// std::string ModelName[5]  = { "NuBeam", "qgsp_bert", "ftfp_bert", "NuBeam-with-res-decays", "qgsp-g4lund-str-fragm"};  
// std::string ModelName[4]  = { "ftfp", "qgsp", "ftfp_bert", "qgsp_ftfp_bert" };  
// const int NModels = 3;
// std::string ModelName[3]  = { "ftfp", "qgsp", "qgsp-g4lund-str-fragm" };
//int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

// static int isReadNA49Loaded = 0;

#include "REGRESSION_TEST.h"


#ifndef G4VAL_DRAWNA49SPECTRA_C
#define G4VAL_DRAWNA49SPECTRA_C

void drawIntegratedSpectrum( std::string beam, std::string target,  
                             std::string secondary, std::string histo,
                             bool drawLegend=true,
			     bool NuERange=false )
{
       
   int recid = 14649;
   if ( secondary == "piminus" ) recid = 14651;
   if ( secondary == "proton" ) recid = 14653;
   if ( secondary == "antiproton" ) recid = 14655; // FIXME !!! There's something goofy with pbar exp.data !!!
   if ( secondary == "neutron" ) recid = 14657;
   
   if ( histo == "pT" ) recid += 1;   
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
/* remove DBReader, revert back to using local ASCII files 
   DbReader dbr;
   TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( dbr.getDataByID( recid ) );

   int NBins = gr1->GetN();
   double* YY = gr1->GetY();
   for ( int i=0; i<NBins; ++i )
   {
      if ( (YY[i]+gr1->GetErrorYhigh(i)) > ymax ) ymax = YY[i]+gr1->GetErrorYhigh(i);
      if ( (YY[i]-gr1->GetErrorYlow(i)) < ymin ) ymin = YY[i]-gr1->GetErrorYlow(i);  
   }
*/

   readIntegratedSpectra( beam, target, secondary );

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
         Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
         Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
         Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }

   TH1F* hi[NModels_HE];
   std::string YTitle;
   if ( histo == "dNdxF" )
   {
      YTitle = "dN/dxF";
   }
   else if ( histo == "pT" )
   {
      YTitle = "Average pT, GeV/c";
   }
   else if ( histo == "pT2" )
   {
      YTitle = "dpT^2/dxF, (GeV/c)^2";
   }
   
   for ( int m=0; m<NModels_HE; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      if ( secondary == "piplus" || secondary == "piminus"  || secondary == "antiproton" )
      {
         hi[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
      }
      
/*
      if ( histo == "pT" && m<3 )
      {
         hi[m]->Scale(32.);
      }
*/
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel_HE[m]);
      hi[m]->SetLineWidth(2);
      if ( m == 0 ) hi[m]->SetLineWidth(3.5);
      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
         if ( histo == "dNdxF" )
	 {
	    hi[m]->Draw("histE1");
	 }
	 else
	 {
	    hi[m]->Draw("e");
	 }
	 // hi[m]->Draw();
	 hi[m]->GetXaxis()->SetTitle("xF");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->SetTitleOffset(0.9);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle( YTitle.c_str() );
//	 hi[m]->GetYaxis()->SetTitleOffset(1.5);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
      }
      else
      {
         if ( histo == "dNdxF" )
	 {
	    hi[m]->Draw("histE1same");
	 }
	 else
	 {
	    hi[m]->Draw("esame");
	 } 
      }   

      if ( NuERange )
      {
         for ( int ir=0; ir<4; ++ir )
	 {
	    std::ostringstream osCount;
            // osCount.clear();
            osCount << ir;
            std::string hname_tmp = histoname + "_" + osCount.str();
	    TH1F* htmp = (TH1F*)f->Get( hname_tmp.c_str() );
	    htmp->SetLineColor(kBlue /* kBlack */ );
	    htmp->SetLineWidth(2);
	    htmp->SetLineStyle(ir);
	    htmp->Draw("same");
         }
      }
   
   }

   TLegend* leg = 0;
   if ( secondary == "neutron" )
   {
// for neutrons
      leg = new TLegend(0.6, 0.62, 0.97, 0.92);
//      leg->SetTextSize(0.035);
      leg->SetTextSize(0.06);
//      leg->SetTextSize(0.065);
   }
   else
   {
// for other secondaries   
      leg = new TLegend(0.65, 0.67, 0.99, 0.92);
      leg->SetTextSize(0.053);
   }
   
   for ( int m=0; m<NModels_HE; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName_HE[m].c_str(), "L" );
   }
     
   TGraph* gr1 = new TGraphErrors( NPointsNA49, xF, Value, 0, Error );

   // gr->GetYaxis()->SetRangeUser( 0., 2.5 );
   // gr->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   // gr->SetRangeUser( 0.5, 2. );
   gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   int nbh = hi[0]->GetNbinsX();
   gr1->GetXaxis()->SetRangeUser( hi[0]->GetBinLowEdge(1), hi[0]->GetBinLowEdge(nbh)+hi[0]->GetBinWidth(nbh) );
   gr1->SetMarkerColor(kBlack /* kBlue */ );
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
    
   gr1->Draw("p");
      
   leg->AddEntry( gr1, "exp.data", "p");

   if ( drawLegend )
   {
      leg->Draw();
      leg->SetFillColor(kWhite);   
   }

   return;

}

void drawIntSpectrumMC2Data( std::string beam, std::string target, std::string secondary, std::string histo, bool drawLegend=true )
{

   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
	 Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
	 Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
	 Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
   }

   TH1F* hi[NModels_HE];
   
   float* MC2DataX = new float[NPointsNA49];
   float* MC2DataY = new float[NPointsNA49];
   float* DX = new float[NPointsNA49];
   float* DY = new float[NPointsNA49];
   int np = 0;
   
   TGraph* gr[NModels_HE];

   for ( int m=0; m<NModels_HE; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {
         double xx1 = hi[m]->GetBinLowEdge(k);
	 double xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsNA49; kk++ )
	 {
	    if ( xx1 < xF[kk] && xx1+xx2 > xF[kk] )
	    {
	       double yy  = hi[m]->GetBinContent(k);
	       double eyy = hi[m]->GetBinError(k);
	       MC2DataX[np] = xF[kk];
	       DX[np] = 0.;
	       
	       MC2DataY[np] = yy / Value[kk];
	       // also need error calc here !...
	       // DY[np]=0.;
	       DY[np] = Error[kk]*MC2DataY[np] / Value[kk];
	       //
	       // the right way to do it is this:
	       // enew2 = (e0**2 * c1**2 + e1**2 * c0**2) / c1**4
	       // enew = sqrt( enew2 )
	       //
	       // but in our case of 1M beam events (~330K inelastic interactions)
	       // the sim error will be about 5% vs 0.5-1.5% in the exp.data, 
	       // thus the sim error will largely dominate... 
	       // need to get at least 10 times more the statistics !...
	       //
	       // double etmp2 = ( eyy*eyy + Error[kk]*Error[kk]*MC2DataY[np]*MC2DataY[np] );
	       // DY[np] = sqrt(etmp2) / Value[kk];
	       //
	       if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	       if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	       np++;
	       break;
	    }
	 }
      }
      
      gr[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      // gr[m]->SetTitle(hi[m]->GetTitle());
      gr[m]->SetTitle("");
//      gr[m]->GetXaxis()->SetTitle("xF");
//      gr[m]->GetYaxis()->SetTitle("MC/Data");
      gr[m]->SetMarkerColor(ColorModel_HE[m]);  
      gr[m]->SetMarkerStyle(SymbModel_HE[m]);
      gr[m]->SetMarkerSize(1.2);
      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
      {
	 gr[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
      }
      // if ( m==0 ) gr1->SetTitle(hi[m]->GetTitle());   
   }
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      Error[i] /= Value[i];
      Value[i] = 1.;
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
   TGraph* gr1 = new TGraphErrors( NPointsNA49, xF, Value, 0, Error );
   if ( secondary == "proton" )
   {
      gr1->GetYaxis()->SetRangeUser( 0.2, 5.0 );
   }
   else
   {
      gr1->GetYaxis()->SetRangeUser( 0.5, 2.0 );
   }
   // gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   // gr1->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   gr1->SetMarkerColor(kBlack /* kBlue */ );
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
   gr1->SetTitle("");
   
   TAxis* xaxis = gr1->GetXaxis();
   
   xaxis->SetTitle("xF");
   xaxis->SetTitleSize(0.05);
   xaxis->SetTitleFont(62);
   xaxis->SetLabelFont(62);
   xaxis->SetTitleOffset(0.9);
   xaxis->CenterTitle();
   
   if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
   {
         xaxis->SetLimits( -0.55, 0.55 );
   }
   
   if ( histo == "dNdxF" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (dN/dxF)");
   }
   else if ( histo == "pT" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (Average pT, GeV/c)" );
   }
   gr1->GetYaxis()->SetTitleSize(0.05);
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
//   gr1->GetYaxis()->SetTitleOffset(1.5);
   gr1->GetYaxis()->CenterTitle();

   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   for ( int m=0; m<NModels_HE; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = 0;
// for neutrons
   if ( secondary == "neutron" )
   {
      leg = new TLegend(0.15, 0.1, 0.55, 0.35); 
//      leg->SetTextSize(0.035);
        leg->SetTextSize(0.055);  
//      leg->SetTextSize(0.065);
   }  
   else if ( secondary == "proton" )
   {
     leg = new TLegend( 0.1, 0.65, 0.47, 0.95);
      leg->SetTextSize(0.05);
   }
// for other secondaries
   else
   {
      leg = new TLegend(0.1, 0.15, 0.47, 0.4); 
      leg->SetTextSize(0.05);
   }  
   
   for ( int m=0; m<NModels_HE; m++ )
   {
      leg->AddEntry( gr[m], ModelName_HE[m].c_str(), "p" );
   }

   leg->AddEntry( gr1, "exp.data", "p");

   if ( drawLegend )
   {
      leg->Draw();
      leg->SetFillColor(kWhite);
   }
   
   return;

}

void drawIntegratedSpectrumRegre( std::string beam, std::string target, std::string secondary, 
                                   std::string histo, std::string model,
				   bool drawLegend=true )
{

   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
         Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
         Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
         Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
   
   TH1F* hi[NVersions];
   std::string YTitle;
   if ( histo == "dNdxF" )
   {
      YTitle = "dN/dxF";
   }
   else if ( histo == "pT" )
   {
      YTitle = "d<pT>/dxF, GeV/c";
   }
   else if ( histo == "pT2" )
   {
      YTitle = "d<pT>^2/dxF, (GeV/c)^2";
   }
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      
/*
      std::string histofile = "";      
      if ( Versions[iv] == CurrentVersion )
      {
         histofile = "./na49-histo/" + beam + target + "158.0GeV" + model;
      }
      else
      {
         histofile = Versions[iv] + "/na49-histo/" + beam + target + "158.0GeV" + model;
      }
      histofile += ".root";
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
      std::string histofile = location + "/na49-histo/" + beam + target + "158.0GeV" + model + ".root"; 
      
 //     std::cout << " Regression, histofile: " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );
      std::string histoname = secondary + "_" + histo ;
      
      hi[iv] = (TH1F*)f->Get( histoname.c_str() );
      
      hi[iv]->SetStats(0);
      hi[iv]->SetLineColor(ColorVersion[iv]);
      hi[iv]->SetLineWidth(6-iv);
      int nx = hi[iv]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[iv]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( iv == 0 ) 
      {
         hi[iv]->Draw();
	 hi[iv]->GetXaxis()->SetTitle("xF");
	 hi[iv]->GetXaxis()->SetTitleSize(0.05);
	 hi[iv]->GetXaxis()->SetTitleFont(62);
	 hi[iv]->GetXaxis()->SetTitleOffset(0.9);
	 hi[iv]->GetXaxis()->CenterTitle();
	 hi[iv]->GetXaxis()->SetLabelFont(62);
	 hi[iv]->GetYaxis()->SetTitle( YTitle.c_str() );
	 // hi[iv]->GetYaxis()->SetTitleOffset(1.5);
         hi[iv]->GetYaxis()->SetTitleSize(0.05);
	 hi[iv]->GetYaxis()->SetTitleFont(62);
         hi[iv]->GetYaxis()->CenterTitle();
	 hi[iv]->GetYaxis()->SetLabelFont(62);
      }
      else hi[iv]->Draw("same");     

   }

   TLegend* leg = new TLegend(0.6, 0.6, 0.98, 0.9);
//   TLegend* leg = new TLegend(0.65, 0.70, 0.98, 0.95);
   if ( secondary == "neutron" )
   {
// for neutrons
//      leg->SetTextSize(0.028);
      leg->SetTextSize(0.06);
   }
   else
   {
// for other particles
      leg->SetTextSize(0.045);
   }
   
   for ( int iv=0; iv<NVersions; iv++ )
   {
      hi[iv]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      // std::string htitle = "pi- on " + target + ", " + model;
      // hi[iv]->SetTitle( htitle.c_str() );
      // leg->AddEntry( hi[iv], entry.c_str(), "L" );
      leg->AddEntry( hi[iv], Versions[iv].c_str(), "L" );
   }

   TGraph* gr = new TGraphErrors( NPointsNA49, xF, Value, 0, Error );

   // gr->GetYaxis()->SetRangeUser( 0., 2.5 );
   // gr->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   // gr->SetRangeUser( 0., 2.5 );
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

   return;

}


void drawIntSpectrumMC2DataRegre( std::string beam, std::string target, std::string secondary, 
                                  std::string histo, std::string model,
				  bool drawLegend=true )
{

   readIntegratedSpectra( beam, target, secondary );

   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
	 Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
	 Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
	 Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
   }

   TH1F* hi[NVersions];
   
   float* MC2DataX = new float[NPointsNA49];
   float* MC2DataY = new float[NPointsNA49];
   float* DX = new float[NPointsNA49];
   float* DY = new float[NPointsNA49];
   int np = 0;
   
   TGraph* gr[NVersions];

   for ( int m=0; m<NVersions; m++ )
   {

//      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + model + ".root"; 
/*
      std::string histofile = "";      
      if ( Versions[m] == CurrentVersion )
      {
         histofile = "./na49-histo/" + beam + target + "158.0GeV" + model;
      }
      else
      {
         histofile = Versions[m] + "/na49-histo/" + beam + target + "158.0GeV" + model;
      }
      histofile += ".root";
*/

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
      std::string histofile = location + "/na49-histo/" + beam + target + "158.0GeV" + model + ".root"; 

      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {
         double xx1 = hi[m]->GetBinLowEdge(k);
	 double xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsNA49; kk++ )
	 {
	    if ( xx1 < xF[kk] && xx1+xx2 > xF[kk] )
	    {
	       double yy  = hi[m]->GetBinContent(k);
	       double eyy = hi[m]->GetBinError(k);
	       MC2DataX[np] = xF[kk];
	       DX[np] = 0.;
	       
	       MC2DataY[np] = yy / Value[kk];
	       // also need error calc here !...
	       // DY[np]=0.;
	       DY[np] = Error[kk]*MC2DataY[np] / Value[kk];
	       //
	       // the right way to do it is this:
	       // enew2 = (e0**2 * c1**2 + e1**2 * c0**2) / c1**4
	       // enew = sqrt( enew2 )
	       //
	       // but in our case of 1M beam events (~330K inelastic interactions)
	       // the sim error will be about 5% vs 0.5-1.5% in the exp.data, 
	       // thus the sim error will largely dominate... 
	       // need to get at least 10 times more the statistics !...
	       //
	       // double etmp2 = ( eyy*eyy + Error[kk]*Error[kk]*MC2DataY[np]*MC2DataY[np] );
	       // DY[np] = sqrt(etmp2) / Value[kk];
	       //
	       if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	       if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	       np++;
	       break;
	    }
	 }
      }
      
      gr[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      // gr[m]->SetTitle(hi[m]->GetTitle());
      gr[m]->SetTitle("");
//      gr[m]->GetXaxis()->SetTitle("xF");
//      gr[m]->GetYaxis()->SetTitle("MC/Data");
      gr[m]->SetMarkerColor(ColorVersion[m]);  
      gr[m]->SetMarkerStyle(SymbVersion[m]);
//      gr[m]->SetMarkerStyle(22);
      gr[m]->SetMarkerSize(1.2);
      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
      {
	 gr[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
      }
      
      // if ( m==0 ) gr1->SetTitle(hi[m]->GetTitle());
   
   }
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      Error[i] /= Value[i];
      Value[i] = 1.;
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
   TGraph* gr1 = new TGraphErrors( NPointsNA49, xF, Value, 0, Error );
   if ( secondary == "proton" )
   {
      gr1->GetYaxis()->SetRangeUser( 0.2, 5. );
   }
   else
   {
      gr1->GetYaxis()->SetRangeUser( 0.5, 2. );
   }
   // gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   // gr1->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   gr1->SetMarkerColor(kBlack /* kBlue */ );
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
   gr1->SetTitle("");
   
   TAxis* xaxis = gr1->GetXaxis();
   
   xaxis->SetTitle("xF");
   xaxis->SetTitleSize(0.05);
   xaxis->SetTitleOffset(0.9);
   xaxis->SetTitleFont(62);
   xaxis->SetLabelFont(62);
   xaxis->CenterTitle();
   
   if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
   {
         xaxis->SetLimits( -0.55, 0.55 );
   }
   if ( histo == "dNdxF" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (dN/dxF)");
//      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
//      {
//         xaxis->SetLimits( -0.55, 0.55 );
//      }
   }
   else if ( histo == "pT" )
   {
      gr1->GetYaxis()->SetTitle("MC/Data (Average pT, GeV/c)" );
   }
   gr1->GetYaxis()->SetTitleSize(0.05);
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
//   gr1->GetYaxis()->SetTitleOffset(1.5);
   gr1->GetYaxis()->CenterTitle();



   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   for ( int m=0; m<NVersions; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      gr[m]->Draw("lpsame");
   }
  
   TLegend* leg = 0;
   if ( secondary == "neutron" )
   {
      leg = new TLegend(0.1, 0.65, 0.55, 0.9); 
   }
   else
   { 
      leg = new TLegend(0.1, 0.65, 0.5, 0.95);  
   }
   leg->SetTextSize(0.045); 
   
   for ( int m=0; m<NVersions; m++ )
   {
      leg->AddEntry( gr[m], Versions[m].c_str(), "p" );
   }

   leg->AddEntry( gr1, "exp.data", "p");

   if ( drawLegend )
   {
      leg->Draw();
      leg->SetFillColor(kWhite);
   }
   
   return;

}

void drawStdDiffGraph( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );

   //double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   //double ymax = -1. ;
   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
	 Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
	 Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
	 Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
   }

   TH1F* hi[NModels_HE];
   
   float* StdDiffX = new float[NPointsNA49];
   float* StdDiffY = new float[NPointsNA49];
   int np = 0;
   
   TGraph* gr[NModels_HE];

   for ( int m=0; m<NModels_HE; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName_HE[m] + ".root"; 

      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {
         double xx1 = hi[m]->GetBinLowEdge(k);
	 double xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsNA49; kk++ )
	 {
	    if ( xx1 < xF[kk] && xx1+xx2 > xF[kk] )
	    {	       
	       double yy = hi[m]->GetBinContent(k);
	       double yerror = hi[m]->GetBinError(k);
	       StdDiffX[np] = xF[kk];
	       StdDiffY[np] = ( yy - Value[kk]) / sqrt( Error[kk]*Error[kk] + yerror*yerror );
	       if ( StdDiffY[np] > ymax ) ymax = StdDiffY[np];
	       if ( StdDiffY[np] < ymin ) ymin = StdDiffY[np];
	       np++;
	       break;
	    }
	 }
      }
            
      gr[m] = new TGraph( np, StdDiffX, StdDiffY );
      gr[m]->SetTitle(hi[m]->GetTitle());
      // gr[m]->SetTitle("");
      gr[m]->GetXaxis()->SetTitle("xF");
      gr[m]->GetYaxis()->SetTitle("SSMD");
      gr[m]->GetYaxis()->SetTitleOffset(1.5);
      gr[m]->SetMarkerColor(ColorModel_HE[m]);  
      gr[m]->SetMarkerStyle(SymbModel_HE[m]);
      gr[m]->SetMarkerSize(1.6);
      if ( secondary == "piplus" || secondary == "piminus" || secondary == "antiproton" )
      {
	 // gr[m]->GetXaxis()->SetRangeUser(-0.55,0.55);
	 gr[m]->GetXaxis()->SetLimits(-0.55,0.55);
      }
      
      // if ( m==0 ) gr1->SetTitle(hi[m]->GetTitle());
   
   }
   
   int iymin = int(ymin) - 1;
   int iymax = int(ymax) + 1;
      
   TLine* line = new TLine();
         
   for ( int m=0; m<NModels_HE; m++ )
   {
      if ( m == 0 )
      {
         // gr[m]->GetYaxis()->SetRangeUser( -2., 1.0 );
         gr[m]->GetYaxis()->SetRangeUser( float(iymin), float(iymax) );
	 /// gr[m]->GetYaxis()->SetLimits( ymin, ymax );
	 gr[m]->Draw("ap");
	 line->DrawLine(-0.55,0.,0.55,0.);
      }
      else
      {
         // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
         gr[m]->Draw("psame");
      }
   }
  
   TLegend* leg = new TLegend(0.1, 0.70, 0.43, 0.9);   
   leg->SetTextSize(0.025);
   
   for ( int m=0; m<NModels_HE; m++ )
   {
      leg->AddEntry( gr[m], ModelName_HE[m].c_str(), "p" );
   }

   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}

void drawStdDiffHisto( std::string beam, std::string target, std::string secondary, std::string histo )
{

   readIntegratedSpectra( beam, target, secondary );

   float* Value = new float[NPointsNA49];
   float* Error = new float[NPointsNA49];
   
   for ( int i=0; i<NPointsNA49; i++ )
   {
      if ( histo == "dNdxF" )
      {
	 Value[i] = dNdxF[i];
	 Error[i] = err_dNdxF[i];
      }
      else if ( histo == "pT" )
      {
	 Value[i] = pT[i];
	 Error[i] = err_pT[i];
      }
      else if ( histo == "pT2" )
      {
	 Value[i] = pT2[i];
	 Error[i] = err_pT2[i];
      }
   }

   TH1F* hi[NModels_HE];
   
   float** StdDiffY = new float*[NModels_HE];
   for ( int m=0; m<NModels_HE; m++ )
   {
      StdDiffY[m] = new float[NPointsNA49];
   }
   
   int np = 0;
   
   TH1F* hst[NModels_HE];

   float ymax = -1.;
   float ymin = 1000000.;
   
   for ( int m=0; m<NModels_HE; m++ )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + histo ;
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );

      int nx = hi[m]->GetNbinsX();
      
      np=0;
      for ( int k=1; k<=nx; k++ ) 
      {
         float xx1 = hi[m]->GetBinLowEdge(k);
	 float xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsNA49; kk++ )
	 {
	    if ( xx1 < xF[kk] && xx1+xx2 > xF[kk] )
	    {	       
	       float yy = hi[m]->GetBinContent(k);
	       float yerror = hi[m]->GetBinError(k);
	       StdDiffY[m][np] = ( yy - Value[kk]) / sqrt( Error[kk]*Error[kk] + yerror*yerror );
	       if ( StdDiffY[m][np] > ymax ) ymax = StdDiffY[m][np];
	       if ( StdDiffY[m][np] < ymin ) ymin = StdDiffY[m][np];
	       np++;
	       break;
	    }
	 }
      } 

   }
   
   float hstxmax = std::max( std::fabs(ymin), std::fabs(ymax) );            
   int ihstxmax = int(hstxmax) + 1;
   ymax = float(ihstxmax);
   ymin = - ymax ;
   int nbins = (2.*float(ihstxmax) / 0.1);

   float hstymax = -1;

   for ( int m=0; m<NModels_HE; m++ )
   {
      
      std::ostringstream osnum;
      osnum << m;
      std::string hstname = "hst" + osnum.str(); 
                   
      hst[m] = new TH1F( hstname.c_str(), hi[m]->GetTitle(), nbins, ymin, ymax );
      std::string xtitle = "SSMD (from " + histo + " spectra)";
      hst[m]->GetXaxis()->SetTitle(xtitle.c_str());
      
      for ( int k=0; k<NPointsNA49; k++ )
      {
         hst[m]->Fill( StdDiffY[m][k] );
      }
      
      // now find out glocal ymax (to set plot's limit later)
      int nx = hst[m]->GetNbinsX();
      for ( int k=1; k<=nx; k++ )
      {
         if ( (hst[m]->GetBinContent(k)+hst[m]->GetBinError(k)) > hstymax ) hstymax = hst[m]->GetBinContent(k)+hst[m]->GetBinError(k);
      }                     
   }
   
   TLegend* leg = new TLegend(0.1, 0.70, 0.45, 0.9);  
   
   for ( int m=0; m<NModels_HE; ++m )
   {

      leg->AddEntry( hst[m], ModelName_HE[m].c_str(), "L" );
      
      hst[m]->SetStats(0);
      hst[m]->SetLineColor( ColorModel_HE[m] );
      // hst[m]->SetLineStyle(NModels_HE-m);
      hst[m]->SetLineWidth(2.);       
      if ( m == 0 )
      {
         hst[m]->GetYaxis()->SetRangeUser( 0., hstymax );
	 hst[m]->Draw();
      }
      else
      {
         hst[m]->Draw("same");
      }
   }
   
   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}

// this routine will draw a ** single ** plot for a selected xF-bin 
// (specified by icount, from 0 to 21)
//
void draw1DDiffXSec( std::string beam, std::string target, std::string secondary, int icount )
{

//   TLegend* leg = new TLegend(0.55, 0.7, 0.99, 0.9);
//   leg->SetTextSize(0.035);
   TLegend* leg = new TLegend(0.45, 0.7, 0.95, 0.9);
   leg->SetTextSize(0.042);
   
   for ( int m=0; m<NModels_HE; ++m )
   {

      std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + ModelName_HE[m] + ".root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::ostringstream cnt;
      cnt << (icount+7);
      std::string hname = "pT";
      if ( secondary == "piplus" )
      {
	    hname +="pip";
      }
      else if ( secondary == "piminus" )
      {
	    hname += "pim";
      }
      else if ( secondary == "proton" )
      {
            hname += "pro";
      }
      hname += cnt.str();
      TH1F* h = (TH1F*)f->Get( hname.c_str() );
      h->Scale(226.);
      // h->Scale( 251. ); // is it G4's xsec p+c at 158GeV ???
      h->SetStats(0);
      h->SetLineColor(ColorModel_HE[m]);
      h->SetLineWidth(2);
      h->SetTitle(0);
      h->GetYaxis()->SetRangeUser( 0.001, 700. );
      leg->AddEntry( h, ModelName_HE[m].c_str(), "L" );
      if ( m == 0 ) 
      {
	    gPad->SetLogy();
	    std::string newname = "p+C ->" + secondary + " + X at 158GeV/c, xF=";
            std::ostringstream xf;
	    xf << DDiff_xF[icount];
	    newname += xf.str();
	    h->SetTitle( newname.c_str() );
//	    h->SetTitleSize(0.055);
	    h->GetXaxis()->SetTitleOffset(1.1);
	    h->GetXaxis()->SetTitleSize(0.045);
	    h->GetXaxis()->SetTitleFont(62);
	    h->GetXaxis()->CenterTitle();
	    h->GetXaxis()->SetTitle( "p_{T}, GeV/c" );
	    h->GetYaxis()->SetTitleOffset(0.9);
	    h->GetYaxis()->SetTitle( "f(x_{F},p_{T}), mb/(GeV^{2}/c^{3})" );
//	    h->GetYaxis()->SetTitle( "" );
	    h->GetYaxis()->SetTitleSize(0.045);
	    h->GetYaxis()->SetTitleFont(62);
	    h->GetYaxis()->CenterTitle();
	    h->GetYaxis()->SetRangeUser(0.01, 1000. );
	    h->Draw();
      }
      else
      { 
	    h->Draw("same");
      }
   } 
   
   int NPt = DDiffDataHolder[icount].GetEntries();
   double*    tmpPT    = new double[NPt];
   double*    tmpXSec  = new double[NPt];
   double*    tmpEXSec = new double[NPt];
   
   // std::cout << " icount = " << icount << " xF = " << DDiff_xF[icount] << std::endl;
   
   for ( int i=0; i<NPt; ++i )
   {
      
       TVector3* tmpVec = (TVector3*)DDiffDataHolder[icount].At(i);
      tmpPT[i]    = tmpVec->x();
      tmpXSec[i]  = tmpVec->y();
      tmpEXSec[i] = tmpVec->z();
      // std::cout << " pt, xsec, exsec = " << tmpPT[i] << " " << tmpXSec[i] << " " << tmpEXSec[i] << std::endl;
   }
   
   TGraph* gr = new TGraphErrors( NPt, tmpPT, tmpXSec, 0, tmpEXSec );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->Draw("psame");
   gr->SetMarkerSize(1.0);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   leg->Draw();
   leg->SetFillColor(kWhite);  
   
   delete [] tmpPT;
   delete [] tmpXSec;
   delete [] tmpEXSec; 
               
   return;

}

#endif
