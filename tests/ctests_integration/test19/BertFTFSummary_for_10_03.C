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
void DrawBertFTF_for_10_03( std::string, std::string, std::string,
                            std::string,
		            std::string,
		            int ); 
			    
void BertFTF_for_10_03( std::string beam, std::string target, std::string energy,
                        std::string secondary,
		        std::string region ) 
{

   ReadHARPData( beam, target, energy, secondary, region );

   TCanvas* myc  = 0;
   TPad*    pad = 0;

   if ( region == "FW" )
   {   
      myc = new TCanvas( "myc", "", 700, 700 );
      std::string txt1 = energy + " GeV/c " + beam + " + " + target + " #rightarrow " + secondary;
      txt1 += " + X (FW)";
      TLatex* ltxt1 = new TLatex( 0.25, 0.96, txt1.c_str() );
      ltxt1->SetTextSize(0.025);
      myc->cd();
      ltxt1->Draw();
      // myc->Divide( 2, 2 );
      pad = new TPad("pad", "", 0.01, 0.13, 0.99, 0.95 );
      myc->cd();
      pad->Draw();
      pad->Divide(2.,2.,0.,0.);
   }
   else if ( region == "LA" )
   {
      myc = new TCanvas( "myc", "", 1000, 1000 );
      std::string txt2 = energy + " GeV/c " + beam + " + " + target + " #rightarrow " + secondary;
      txt2 += " + X (LA)";
      TLatex* ltxt2 = new TLatex( 0.25, 0.96, txt2.c_str() );
      ltxt2->SetTextSize(0.025);
      ltxt2->Draw();
      myc->cd();
      pad = new TPad("pad", "", 0.01, 0.13, 0.99, 0.95 );
      myc->cd();
      pad->Draw();
      pad->Divide(3.,3.,0.,0.);
   }
   
   for ( int i=0; i<NSetsHARP; i++ )
   {
      // myc->cd(i+1);
      myc->cd();
      pad->cd(i+1);
      gPad->SetRightMargin(0.1);
      gPad->SetLeftMargin(0.1);
      DrawBertFTF_for_10_03( beam, target, energy, secondary, region, i );
   }

   TLegend* leg = new TLegend(0.7,0.01,0.99,0.12);
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
   myc->cd();
   leg->Draw();

   std::string general = "#chi^{2}/NDF calculated over ";
   general += ( region + " theta bins" );

   TLatex* ltxt = new TLatex(0.1, 0.10, general.c_str() );
   ltxt->SetTextSize(0.02);
   myc->cd();
   ltxt->Draw();
   
   for ( int m=0; m<NNN; ++m )
   {

      double chi2 = 0.;
      int NDF = 0;
      std::ostringstream os;

      chi2 = Chi2MomSpectrumHARP( beam, target, energy, secondary, region, model[m], NDF, version[m] );

      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = " + os.str();
      txt1 += ( " for " + model[m] + " " + version[m] );
      TLatex* ltxt1 = new TLatex( 0.1, 0.08-0.02*m, txt1.c_str() );
      ltxt1->SetTextSize(0.02);
      myc->cd();
      ltxt1->Draw();
      
      std::cout << " Geant4 model/physlis = " << model[m] << " " << version[m] << std::endl;
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << std::endl;

   }
   
   std::string output = beam + "-" + target + "-" + energy + "GeV" + "-" + secondary;
   output += ( "-" + region + ".gif" );
   myc->Print( output.c_str() );

   return;

}

void DrawBertFTF_for_10_03( std::string beam, std::string target, std::string energy,
                            std::string secondary,
		            std::string region,
		            int ibin )
{
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   // this assumes that ReadHARPData has been called already...
   //
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      if ( (YHARP[ibin][i]+EYHARP[ibin][i]) > ymax ) ymax = YHARP[ibin][i]+EYHARP[ibin][i];
      if ( (YHARP[ibin][i]-EYHARP[ibin][i]) < ymin && (YHARP[ibin][i]-EYHARP[ibin][i]) > 0. ) ymin = (YHARP[ibin][i]-EYHARP[ibin][i]);
   }
   
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
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = secondary + "_" + region + "_" + std::to_string(ibin); 
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      if ( hi[m] == NULL )
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
/*
	 std::string htit = std::to_string(AngleBinHARP[0][ibin]) 
	                  + " < #Theta < " 
			  + std::to_string(AngleBinHARP[1][ibin])
			  + " [rad]"; 
*/
	 std::ostringstream os1, os2;
	 os1 << AngleBinHARP[0][ibin];
	 os2 << AngleBinHARP[1][ibin];
	 std::string htit = os1.str() + " < #Theta < " + os2.str() + " [rad]";
	 hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("momentum (GeV/c)");
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
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
   }

   float* X = new float[NPointsHARP[ibin]];
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      X[i] = 0.5 * (XMinHARP[ibin][i]+XMaxHARP[ibin][i]);
   }

   TGraph* gr = new TGraphErrors( NPointsHARP[ibin], X, YHARP[ibin], 0, EYHARP[ibin] );
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");
   
   delete [] X;
   
   return;

}

void BertFTFSummary_for_10_03( std::string beam, std::string target, std::string energy )
{

   BertFTF_for_10_03( beam, target,energy,"piplus","FW");
   BertFTF_for_10_03( beam, target,energy,"piplus","LA");
   BertFTF_for_10_03( beam, target,energy,"piminus","FW");
   BertFTF_for_10_03( beam, target,energy,"piminus","LA");

   return;
   

}

