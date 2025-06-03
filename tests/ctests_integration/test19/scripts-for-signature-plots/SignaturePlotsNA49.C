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

#include "../test23/shared-root-macros/NA49.h"
#include "../test23/shared-root-macros/ReadNA49Data.C"
// #include "Chi2Calc.C"
// #include "Chi2CalcNA49.C"
#include "../test23/shared-root-macros/DrawNA49Spectra.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_HighEnergy.h"

std::string TEST_NAME="test19";

void SignaturePlotsNA49()
{

   TCanvas* myc1 = new TCanvas( "myc1", "", 600., 700. );
   
   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.47, 0.99, 1.0 );
   myc1->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetRightMargin(0.025);
   // drawIntegratedSpectrumRegre( "proton", "C", "piplus", "pT", "ftfp", false ); //, false ); // no legend 
   DrawPTvsXF();
   
   std::string data = "Data: C. Alt et al., Eur.Phys.J.C49:897-917,2007";
   TText* txt = new TText( -0.45, 0.6, data.c_str() );
   txt->SetTextSize(0.045);
   
   std::string htitle = "158 GeV/c p+C #rightarrow #pi^{+} + X";
   TLatex* ltxt = new TLatex( -0.45, 0.64, htitle.c_str() );
   ltxt->SetTextSize(0.045);
   pad1->cd(0);
   ltxt->Draw();
   txt->Draw();
   
   
   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.1, 0.99, 0.47 );
   myc1->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetRightMargin(0.025);
   // gPad->SetLogy();
   drawIntSpectrumMC2DataRegre( "proton", "C", "piplus", "pT", "ftfp" );
   
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

}

void DrawPTvsXF()
{

   readIntegratedSpectra( "proton", "C", "piplus" );
   
   std::string YTitle = "< p_{T} >[GeV/c]";
   
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
      
      std::string histofile = location + "/na49-histo/protonC158.0GeVftfp.root";       
      TFile* f = new TFile( histofile.c_str() );
      
      std::string histoname = "piplus_pT";
      hi[iv] = (TH1F*)f->Get( histoname.c_str() );
      
      hi[iv]->SetStats(0);
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
      hi[iv]->GetYaxis()->SetTitleFont(62);
      hi[iv]->GetYaxis()->SetLabelSize(0.04);
      hi[iv]->GetYaxis()->SetLabelFont(62);
      hi[iv]->GetYaxis()->CenterTitle();
      hi[iv]->GetYaxis()->SetRangeUser( 0.2, 0.7 );

      if ( iv == 0 )
      {
         hi[iv]->Draw(); // we don't need to worry about "histE1" option because it's a profile histogram
      }
      else hi[iv]->Draw("same");
      
   }

   TGraph* gr = new TGraphErrors( NPointsNA49, xF, pT, 0, err_pT );

   gr->SetMarkerColor(kBlue);
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");

   return;

}
