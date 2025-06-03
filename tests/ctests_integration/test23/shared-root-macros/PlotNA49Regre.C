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

#include "NA49.h"
#include "ReadNA49Data.C"
#include "Chi2Calc.C"
#include "Chi2CalcNA49.C"
#include "DrawNA49Spectra.C"

#ifndef G4VAL_PLOTNA49REGRE
#define G4VAL_PLOTNA49REGRE


void plot_dNdxF_pT_Regre( std::string beam, std::string target, 
                          std::string secondary, 
			  std::string model )
{

   // This part has NOT changed
   //
   /*
    Here should be special provision for neutrons because only dNdxF data are available
    */

   if ( secondary == "neutron" )
   {
   /*
    Here is special provision for neutrons because only dNdxF data are available
    */
      TCanvas* myc = new TCanvas( "myc", "", 600, 600 );
      TPad* pad1 = new TPad( "pad1", "", 0.01, 0.12, 0.99, 0.99 );
      pad1->Draw();
      pad1->Divide(1.,2.,0.,0.);
      pad1->cd(1); gPad->SetRightMargin(0.025);
      drawIntegratedSpectrumRegre( beam, target, secondary, "dNdxF", model ); 
      pad1->cd(2); gPad->SetRightMargin(0.025);
      gPad->SetLogy();
      drawIntSpectrumMC2DataRegre( beam, target, secondary, "dNdxF", model, false );
      myc->cd();
      for ( int m=0; m<NVersions; ++m )
      {
         int NDF1=0;
         double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "dNdxF", model, NDF1, Versions[m] );  
         std::ostringstream os1;
         os1 << (chi2/NDF1);
         std::string txt1 = " #chi^{2}/NDF = " + os1.str();
         txt1 += ( " for " + Versions[m] + " vs NA49 Data" );
         TLatex* ltxt1 = new TLatex( 0.10, 0.10-m*0.025, txt1.c_str() );
         ltxt1->SetTextSize(0.025);
         ltxt1->Draw();  
      } 
      myc->cd();
      std::string output = "proton-C-158GeV-neutron-" + model + "-regre.gif";
      myc->Print( output.c_str() );
      return;      
   }

   TCanvas* myc   = new TCanvas("myc","",1200,600);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.01, 0.4, 0.99 );

   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrumRegre( beam, target, secondary, "dNdxF", model, false ); // no legend 
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawIntSpectrumMC2DataRegre( beam, target, secondary, "dNdxF", model, false ); // no legend

   myc->cd();
      
   TPad* pad2 = new TPad( "pad2", "", 0.41, 0.01, 0.8, 0.99 );

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025); // should the margin be smaller ?
   drawIntegratedSpectrumRegre( beam, target, secondary, "pT", model, false ); // no legend
   pad2->cd(2); gPad->SetRightMargin(0.025); // should the margin be smaller ?
   gPad->SetLogy();
   drawIntSpectrumMC2DataRegre( beam, target, secondary, "pT", model, false ); // no legend
   
   myc->cd();

   TLegend* leg = new TLegend( 0.78, 0.10, 0.99, 0.9 ); // 0.78-0.8 seems like an overlap 
                                                                  // with the pad but it works
   TLegendEntry* lentry = 0; 

   leg->SetFillColor(kWhite);
   leg->SetBorderSize(0);

   lentry = leg->AddEntry("","NA49 data","p");
   lentry->SetMarkerStyle(22);
   lentry->SetMarkerColor(kBlack /* kBlue */ );
   lentry->SetMarkerSize(1.5);
   lentry->SetTextFont(62);
   lentry->SetTextColor(kBlack /* kBlue */ );
   
   lentry = leg->AddEntry("","For dN/dxF spectra","");
   lentry->SetTextFont(62);
   lentry->SetTextSize(0.025);
   lentry->SetTextColor(kBlack /* kBlue */ );

// -->
   for ( int m=0; m<NVersions; ++m )
   {
      int NDF1=0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "dNdxF", model, NDF1, Versions[m] );  
      std::ostringstream os1;
      if ( (chi2/NDF1) < 100 )
      {
         os1.precision(3);
      }
      else
      {
         os1.precision(4);
      }
      os1 << (chi2/NDF1);
      std::string txt1 = " #chi^{2}/NDF=" + os1.str();
      std::string vv = Versions[m];
      if ( vv.find("geant4-") != std::string::npos )
      {
         vv = vv.erase( 0, std::strlen("geant4-") );
      }
      txt1 += ( "   " + vv ); 
      lentry = leg->AddEntry("", txt1.c_str(),"L");
      lentry->SetLineColor( ColorVersion[m] );
      lentry->SetLineWidth(3);
      lentry->SetTextFont(62);
      lentry->SetTextSize(0.025);
      lentry->SetTextColor( ColorVersion[m] );
   }

   lentry = leg->AddEntry("","For <pT> vs xF spectra","");
   lentry->SetTextFont(62);
   lentry->SetTextSize(0.025);
   lentry->SetTextColor(kBlack /* kBlue */ );

// -->
   for ( int m=0; m<NVersions; ++m )
   {
      int NDF1=0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "pT", model, NDF1, Versions[m] );  
      std::ostringstream os1;
      if ( (chi2/NDF1) < 100 )
      {
         os1.precision(3);
      }
      else
      {
         os1.precision(4);
      }
      os1 << (chi2/NDF1);
      std::string txt1 = " #chi^{2}/NDF=" + os1.str();
      std::string vv = Versions[m];
      if ( vv.find("geant4-") != std::string::npos )
      {
         vv = vv.erase( 0, std::strlen("geant4-") );
      }
      txt1 += ( "   " + vv ); 
      lentry = leg->AddEntry("", txt1.c_str(),"L");
      lentry->SetLineColor( ColorVersion[m] );
      lentry->SetLineWidth(3);
      lentry->SetTextFont(62);
      lentry->SetTextSize(0.025);
      lentry->SetTextColor( ColorVersion[m] );
   }

   myc->cd();
   leg->Draw();

   std::string output = beam + "-" + target + "-158GeV-" + secondary + "-" + model + "-regre.gif";
   myc->Print( output.c_str() );

   return;

}

#endif
