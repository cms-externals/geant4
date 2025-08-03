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

// #include "../test23/shared-root-macros/Chi2Calc.C"

#include "../test47/PlotITEPAnalysis.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

void SignaturePlotsITEP771()
{

   TCanvas* myc1 = new TCanvas( "myc1", "", 600., 700. );
   
   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.52, 0.99, 0.98 );
   myc1->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetRightMargin(0.025);
   gPad->SetLogy();

   bool status = readKEData( "piminus", "Cu", "5.00", "neutron", "59.1" );
   
   if ( !status )
   {
      TText* txt = new TText(0.25,0.6,"Insufficient exp.data");
      txt->SetTextSize(0.075);
      txt->Draw();
      return;
   }

   // NOTE: requires calling readKEData first ! (see above)
   //
   TGraphErrors* gr1 = getITEPAsGraph();
   
   TH1F* hi[NVersions];
   
   for ( int iv=0; iv<NVersions; ++iv )
   {
      std::string location = "";
      if ( Versions[iv] == CurrentVersion )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test47/" + Versions[iv] + "/";
      }
      std::string histofile = location + "piminusCubertini5.00GeV.root"; 
      std::string histoname = "KEneutron0piminusCubertini5.00GeV 59.1";
   
      TFile* hfile = new TFile( histofile.c_str() );
   
      hi[iv] = (TH1F*)hfile->Get( histoname.c_str() );
      
      hi[iv]->SetTitle("");
      hi[iv]->SetStats(0);
      hi[iv]->SetLineColor(ColorVersion[iv]);
      hi[iv]->SetLineWidth(6-iv);
      hi[iv]->GetXaxis()->SetRangeUser(0.01,0.2);
      if ( iv == 0 ) 
      {
         // hi[iv]->Draw("histE1");
	 hi[iv]->Draw();
	 hi[iv]->GetXaxis()->SetTitle("");
	 hi[iv]->GetXaxis()->SetTitleSize(0.05);
	 hi[iv]->GetXaxis()->SetTitleFont(62);
	 hi[iv]->GetXaxis()->SetTitleOffset(0.9);
	 hi[iv]->GetXaxis()->SetLabelSize(0.05);
	 hi[iv]->GetXaxis()->SetLabelFont(62);
	 hi[iv]->GetXaxis()->CenterTitle();
	 hi[iv]->GetYaxis()->SetTitle( "E#frac{d^{3}#sigma}{dp^{3}} [mb/GeV^{2}]" );
	 // hi[iv]->GetYaxis()->SetTitleOffset(1.5);
         hi[iv]->GetYaxis()->SetTitleSize(0.055);
	 hi[iv]->GetYaxis()->SetTitleOffset(0.75);
	 hi[iv]->GetYaxis()->SetTitleFont(62);
	 hi[iv]->GetYaxis()->SetLabelSize(0.05);
	 hi[iv]->GetYaxis()->SetLabelFont(62);
         hi[iv]->GetYaxis()->CenterTitle();
      }
      else 
      {
         //hi[iv]->Draw("histE1same");
	 hi[iv]->Draw("same");
      }     
      
   }
   
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerColor(4);
   gr1->SetMarkerSize(1.6);
   gr1->Draw("p");
   
   TLatex* reac = new TLatex( 0.03, 70000., "5 GeV/c #pi^{-} + Cu #rightarrow n + X");
   reac->SetTextSize(0.055);
   TText* data = new TText( 0.03, 40000, "Data: Yu. D. Bayukov et al., Sov.J.Nucl.Phys. 42 116, 1985" );
   data->SetTextSize(0.055);
   pad1->cd();
   reac->Draw(); 
   data->Draw();
   

   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.12, 0.99, 0.52 );
   myc1->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetRightMargin(0.025);
   //gPad->SetLogy();
   plotMC2Data( "piminus", "Cu", "5.00", "neutron", "59.1", "bertini" ); //, 4 );
   
   TText* cap = new TText( 0.03, 0.11, 
                           "Invariant cross section as a function of kinetic energy for neutron" ); 
   cap->SetTextSize(0.028);
   TLatex* cap1 = new TLatex( 0.03, 0.08, "produced in #pi^{-} + Cu interactions at 5 GeV/c. Geant4 (Bertini) results are");
   cap1->SetTextSize(0.028);
   TText* cap2 = new TText( 0.03, 0.05, "compared with the data from ITEP-771.");
   cap2->SetTextSize(0.028);
   
   TLatex* xtitle = new TLatex( 0.25, 0.1, "Kinetic energy of secondary neutron [GeV]");
   xtitle->SetTextSize(0.028);
   
   myc1->cd();
   // cap->Draw();
   // cap1->Draw();
   // cap2->Draw();
   xtitle->Draw();

   return;

}
