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

std::string TEST_NAME="test47";

#include "../test47/PlotITEPAnalysis.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

const int NNN = 4;
std::string model[4] = { "bertini", "ftfp", "ftfp", "ftfp" };
std::string version[4] = { "geant4-10-03-ref-07", "geant4-10-03-ref-07", "geant4-10-03-patch-01", "geant4-10-02-patch-03" };
int color[4] = { kMagenta, kRed, kGreen, kBlack };


// forward declarations
//
void DrawBertFTF_for_10_03( std::string, std::string, std::string,
                            std::string,
		            std::string );
			    
void BertFTF_for_10_03( std::string beam, std::string target, std::string energy )
{

   TCanvas* myc  = 0;
   TPad*    pad = 0;

   myc = new TCanvas( "myc", "", 800, 800 );

   TLegend* leg = new TLegend(0.55,0.875,0.94,0.99);
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

   myc->cd();
   // myc->Divide( 2, 2 );
   pad = new TPad("pad", "", 0.01, 0.14, 0.99, 0.85 );
   myc->cd();
   pad->Draw();
   pad->Divide(2.,2.,0.,0.);

   std::string txt1 = energy + " GeV/c " + beam + " + " + target + " #rightarrow proton + X";;
   TLatex* ltxt1 = new TLatex( 0.075, 0.852, txt1.c_str() );
   ltxt1->SetTextSize(0.02);
   myc->cd();
   ltxt1->Draw();
   std::string txt2 = energy + " GeV/c " + beam + " + " + target + " #rightarrow neutron + X";;
   TLatex* ltxt2 = new TLatex( 0.55, 0.852, txt2.c_str() );
   ltxt2->SetTextSize(0.02);
   myc->cd();
   ltxt2->Draw();
   
   TLatex* gt1 = new TLatex( 0.1, 4.5, "#Theta(proton) = 59.1 [deg]" );
   gt1->SetTextSize(0.075);   
   myc->cd();
   pad->cd(1);
   gPad->SetBottomMargin(0.1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   DrawBertFTF_for_10_03( beam, target, energy, "proton", "59.1" );
   gt1->Draw();

   TLatex* gt2 = new TLatex( 0.04, 4.5, "#Theta(neutron) = 59.1 [deg]" );
   gt2->SetTextSize(0.075);
   myc->cd();
   pad->cd(2);
   gPad->SetBottomMargin(0.1);
   gPad->SetRightMargin(0.1);
   gPad->SetLeftMargin(0.1);
   DrawBertFTF_for_10_03( beam, target, energy, "neutron", "59.1" );
   gt2->Draw();

   TLatex* gt3 = new TLatex( 0.12, 4.5, "#Theta(proton) = 119.0 [deg]" );
   gt3->SetTextSize(0.065);
   myc->cd();
   pad->cd(3);
   gPad->SetTopMargin(0.1);
   gPad->SetRightMargin(0.15);
   gPad->SetLeftMargin(0.15);
   DrawBertFTF_for_10_03( beam, target, energy, "proton", "119.0" );
   gt3->Draw();
  
   TLatex* gt4 = new TLatex( 0.04, 4.5, "#Theta(neutron) = 119.0 [deg]" );
   gt4->SetTextSize(0.065);
   myc->cd();
   pad->cd(4);
   gPad->SetTopMargin(0.1);
   gPad->SetRightMargin(0.1);
   gPad->SetLeftMargin(0.1);
   DrawBertFTF_for_10_03( beam, target, energy, "neutron", "119.0" );
   gt4->Draw();
  
   std::string general = "#chi^{2}/NDF calculated vs ITEP771 data";

   TLatex* ltxt3 = new TLatex(0.3, 0.11, general.c_str() );
   ltxt3->SetTextSize(0.02);
   myc->cd();
   ltxt3->Draw();
   
   for ( int m=0; m<NNN; ++m )
   {

      double chi2 = 0.;
      int NDF = 0;
      std::ostringstream os;

      chi2 += Chi2KESpectrumITEP( beam, target, energy, "proton", "59.1",  model[m], NDF, version[m] );
      chi2 += Chi2KESpectrumITEP( beam, target, energy, "proton", "119.0", model[m], NDF, version[m] );

      os << (chi2/NDF);
      std::string txt4 = "#chi^{2}/NDF = " + os.str();
      txt4 += ( " for " + model[m] + " " + version[m] );
      TLatex* ltxt4 = new TLatex( 0.05, 0.08-0.02*m, txt4.c_str() );
      ltxt4->SetTextSize(0.02);
      myc->cd();
      ltxt4->Draw();
      
      std::cout << " Geant4 model/physlis = " << model[m] << " " << version[m] << std::endl;
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for proton" << std::endl;

   }
   for ( int m=0; m<NNN; ++m )
   {

      double chi2 = 0.;
      int NDF = 0;
      std::ostringstream os;

      chi2 += Chi2KESpectrumITEP( beam, target, energy, "neutron", "59.1",  model[m], NDF, version[m] );
      chi2 += Chi2KESpectrumITEP( beam, target, energy, "neutron", "119.0", model[m], NDF, version[m] );

      os << (chi2/NDF);
      std::string txt5 = "#chi^{2}/NDF = " + os.str();
      txt5 += ( " for " + model[m] + " " + version[m] );
      TLatex* ltxt5 = new TLatex( 0.5, 0.08-0.02*m, txt5.c_str() );
      ltxt5->SetTextSize(0.02);
      myc->cd();
      ltxt5->Draw();
      
      std::cout << " Geant4 model/physlis = " << model[m] << " " << version[m] << std::endl;
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for neutron" << std::endl;

   }
  
   std::string output = beam + "-" + target + "-" + energy + "GeV.gif";
   myc->Print( output.c_str() );

   return;

}

void DrawBertFTF_for_10_03( std::string beam, std::string target, std::string energy,
                            std::string secondary,
		            std::string sec_angle )
{

   bool status = readKEData( beam, target, energy, secondary, sec_angle );
   
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
/*
   std::string gtitle = "#Theta(" + secondary + ") = ";
   std::ostringstream os;
   os << sec_angle;
   gtitle += ( os.str() + " [deg]" );
   gr1->SetTitle( gtitle.c_str() );
*/
   gr1->SetTitle("");
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   std::string xtitle = "Kinetic Energy of secondary " + secondary + " [GeV]";
   gr1->GetXaxis()->SetTitle( xtitle.c_str() );
   gr1->GetXaxis()->SetTitleSize(0.05);
   gr1->GetXaxis()->SetTitleFont(62);
   gr1->GetXaxis()->SetLabelFont(62);
   gr1->GetXaxis()->CenterTitle();
   // gr1->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
   gr1->GetYaxis()->SetTitle( "MC/Data" );
   gr1->GetYaxis()->SetTitleOffset(1.0);
   gr1->GetYaxis()->SetTitleSize(0.05);
   gr1->GetYaxis()->SetTitleFont(62);
   gr1->GetYaxis()->SetLabelFont(62);
   gr1->GetYaxis()->CenterTitle();
         
   std::string location[4]; 
   location[0] = "."; 
   location[1] = ".";
   location[2] = regre_test_dir + "/" + TEST_NAME + "/" + version[2];
   location[3] = regre_test_dir + "/" + TEST_NAME + "/" + version[3];

   TGraph* grMC[NNN];
     
   for ( int m=0; m<NNN; ++m )
   {
      
      std::string histofile = location[m] + "/" + beam + target + model[m] + energy + "GeV.root"; 
      TFile* f = new TFile( histofile.c_str() );

      std::string histoname = "KE" + secondary + "0" + beam + target + model[m] + energy + "GeV";
      int counter = 5 - sec_angle.length();
      for ( int il=0; il<counter; il++ )
      {
         histoname += " "; // pad up for the fact that initially sec.angle was supposed to be char[6] - no more, no less...
      }
      histoname += sec_angle;
      TH1F* hi = (TH1F*)f->Get( histoname.c_str() );
      
      if ( hi == NULL || NPtKE <= 0 ) 
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

      grMC[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      
      grMC[m]->SetMarkerColor(color[m]);  
      grMC[m]->SetMarkerStyle(21);
      grMC[m]->SetMarkerSize(1.0); // mrk size used to be 1.6

   }
    
   // gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   gr1->GetYaxis()->SetRangeUser( 0.2, 5. );
   gr1->Draw("apl");
   
   for ( int m=0; m<NNN; ++m )
   {
      grMC[m]->Draw( "lpsame" );
   }
     
   return;

}

