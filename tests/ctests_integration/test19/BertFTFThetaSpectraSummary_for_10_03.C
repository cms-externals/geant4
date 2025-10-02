#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <math.h>
#include <algorithm>
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
#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test23/shared-root-macros/HARP.h"
#include "../test23/shared-root-macros/ReadHARPData.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcHARP.C"

const int NNN = 4;
std::string model[5] = { "bertini", "ftfp", "ftfp", "ftfp", "ftfp" };
std::string version[5] = { "geant4-10-03-cand-00", "geant4-10-03-cand-00", "geant4-10-02-ref-10", "geant4-10-01-patch-03", "geant4-09-06-patch-04" };
int color[5] = { kMagenta, kRed, kGreen, kBlack, 7 };


// forward declarations
//
void DrawBertFTF_ThetaSpectrum_for_10_03( std::string, std::string, std::string,
                                          std::string,
			                  float, float );


void BertFTF_ThetaSpectra_for_10_03( std::string beam, std::string target, std::string energy,
                                     std::string secondary )
{

   ReadHARPData( beam, target, energy, secondary, "LA" );
   
   TCanvas* myc1 = new TCanvas( "myc1", "", 900, 800 );
   myc1->cd();
   
   std::string txt1 = energy + " GeV/c " + beam + " + " + target + " #rightarrow " + secondary;
   txt1 += " + X";
   TLatex* ltxt1 = new TLatex( 0.35, 0.95, txt1.c_str() );
   ltxt1->SetTextSize(0.025);
   myc1->cd();
   ltxt1->Draw();

   TPad* pad1 = new TPad("pad1", "", 0.01, 0.02, 0.99, 0.92 );
   myc1->cd();
   pad1->Draw();
   pad1->Divide(3.,3.,0.,0.);
   
   for ( int i=0; i<8; ++i )
   {
      float pmin = 0.1 + 0.05*i;
      float pmax = pmin + 0.05;
      pad1->cd(i+1);
      // ---> gPad->SetTopMargin(0.1);
      gPad->SetRightMargin(0.1);
      gPad->SetLeftMargin(0.1);
      DrawBertFTF_ThetaSpectrum_for_10_03( beam, target, energy, secondary, pmin, pmax );      
   }

   TLegend* leg = new TLegend(0.02,0.2,0.98,0.8);
   TLegendEntry* entry = 0;
   for ( int m=0; m<NNN; ++m )
   {
      std::string info = model[m] + " " + version[m];
      entry = leg->AddEntry( "", info.c_str(),"L" );
      entry->SetLineColor( color[m] );
      entry->SetLineWidth(3);
      entry->SetTextFont(62);
   }
   entry = leg->AddEntry("","exp.data","p");
   entry->SetMarkerStyle(22);
   entry->SetMarkerColor(kBlue);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(62);
   leg->SetFillColor(kWhite);
   
   pad1->cd(9);
   leg->Draw();
   
   std::string output = beam + "-" + target + "-" + energy + "GeV" + "-" + secondary;
   output += "-theta-spectra.gif";
   myc1->Print( output.c_str() );

   return;

}

void DrawBertFTF_ThetaSpectrum_for_10_03( std::string beam, std::string target, std::string energy,
                                          std::string secondary,
			                  float pmin, float pmax )
{
  
   int ibin = findHARPMomBin( pmin, pmax );
   if ( ibin < 0 ) 
   {
      std::cout << " INVALID ibin = " << ibin << std::endl;
      return;
   }

   TGraphErrors* gr = getHARPAsThetaGraph( pmin, pmax );

   float ymin = 1000000.;
   float ymax = -1.; 
   ymin = gr->GetYaxis()->GetXmin();
   ymax = gr->GetYaxis()->GetXmax();
   
   TH1F* hi[NNN]; // bertini for 10.2.ref07 + ftf for 10.2.ref07 & 10.2.p02
      
   std::string location[5]; 
   location[0] = "."; 
   location[1] = ".";
   location[2] = regre_test_dir + "/" + TEST_NAME + "/" + version[2];
   location[3] = regre_test_dir + "/" + TEST_NAME + "/" + version[3];
   location[4] = regre_test_dir + "/" + TEST_NAME + "/" + version[4];
     
   for ( int m=0; m<NNN; ++m )
   {

      std::string histofile = location[m] + "/harp-histo/" + beam + target + energy + "GeV" + model[m] + ".root";
// ==>      std::cout << " histofile = " << histofile << std::endl;     
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_mom_" + std::to_string(ibin+1); // MC spectra start 1 mom bin lower than the data
// ==>      std::cout << " histoname = " << histoname << std::endl;    
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      if ( !hi[m] )
      {
         std::cout << " NULL histo " << histoname << " from file " << histofile << std::endl;
	 continue;
      }

      hi[m]->SetStats(0);
      hi[m]->SetLineColor(color[m]);
      hi[m]->SetLineWidth(2);

      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }

      if ( m == 0 ) 
      {
	 // std::string htit = energy + "GeV/c  " + std::string( hi[m]->GetTitle() );
	 // hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("#Theta (rad)");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
	 hi[m]->GetYaxis()->SetTitleOffset(1.0);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
         hi[m]->Draw("histE1");
      }
      else hi[m]->Draw("histE1same"); 

   } 

   for ( int m=0; m<NNN; m++ )
   {
      if ( hi[m] ) hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
   }
   
   gr->Draw( "psame" );  

   return;

}

void BertFTFThetaSpectraSummary_for_10_03( std::string beam, std::string target, std::string energy )
{

   BertFTF_ThetaSpectra_for_10_03( beam, target,energy,"piplus");
   BertFTF_ThetaSpectra_for_10_03( beam, target,energy,"piminus");

   return;  

}

void BertFTFThetaSpectraSummary( std::string secondary )
{

   BertFTF_ThetaSpectra_for_10_03("proton","C","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("proton","C","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("proton","C","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("proton","Be","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("proton","Be","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("proton","Be","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("proton","Ta","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("proton","Ta","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("proton","Ta","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("piplus","C","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piplus","C","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piplus","C","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("piplus","Be","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piplus","Be","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piplus","Be","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("piplus","Ta","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piplus","Ta","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piplus","Ta","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("piminus","C","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piminus","C","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piminus","C","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("piminus","Be","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piminus","Be","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piminus","Be","8.0",secondary );

   BertFTF_ThetaSpectra_for_10_03("piminus","Ta","3.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piminus","Ta","5.0",secondary );
   BertFTF_ThetaSpectra_for_10_03("piminus","Ta","8.0",secondary );

   return;

}
