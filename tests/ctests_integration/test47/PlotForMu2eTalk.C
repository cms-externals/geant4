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

void plotModelRegreMu2e( std::string beam, std::string target, std::string model )
{

   std::string en = "5.00";
   if ( beam == "proton" )
   {
      en = "7.50";
   }
   
   double chi2 = 0.;
   int NDF = 0;
   
   std::string txt = model + " vs (ITEP) Data";
   TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
   ltxt->SetTextSize(0.2);

   TCanvas *myc1 = new TCanvas("myc1","",800,800);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad3 = new TPad( "pad3", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad4 = new TPad( "pad4", "", 0.51, 0.14, 0.99, 0.56 );
   myc1->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "neutron", "59.1", model ); //, 4 );

   myc1->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "neutron", "89.0", model ); //, 4 );
   
   myc1->cd();

   TPad* pad11 = new TPad("pad11","", 0.01, 0.01, 0.49, 0.13 );
   pad11->Draw();
   pad11->cd();

   //std::string txt = model + " vs (ITEP) Data";
   //TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
   //ltxt->SetTextSize(0.2);
   ltxt->Draw();

   int icounter1 = 0;
   for ( int m=0; m<NVersions; ++m )
   {
         chi2 = 0.;
         NDF = 0;
         chi2 += Chi2KESpectrumITEP( beam, target, en, "neutron", "59.1",  model, NDF, Versions[m] );
         chi2 += Chi2KESpectrumITEP( beam, target, en, "neutron", "89.0", model, NDF, Versions[m] );
         std::ostringstream os;
         os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
         std::string txt1 = "#chi^{2}/NDF = ";
         txt1 += os.str();
         txt1 += ( " for " + Versions[m] );
         TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter1*0.2, txt1.c_str() );
         ltxt1->SetTextSize(0.17);
         ltxt1->Draw();
         icounter1++;
   }

   myc1->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "neutron", "119.0", model ); //,  4 );

   myc1->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetLogy();
   plotMC2Data( beam, target, en, "neutron", "159.6", model ); //, 4 );
      
   myc1->cd();

   TPad* pad12 = new TPad("pad12","", 0.51, 0.01, 0.99, 0.13 );
   pad12->Draw();
   pad12->cd();
      
   // std::string txt = model + " vs (ITEP) Data";
   // TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
   // ltxt->SetTextSize(0.2);
   ltxt->Draw();

   icounter1 = 0;
   for ( int m=0; m<NVersions; ++m )
   {
         chi2 = 0.;
         NDF = 0;
         chi2 += Chi2KESpectrumITEP( beam, target, en, "neutron", "119.0",  model, NDF, Versions[m] );
         chi2 += Chi2KESpectrumITEP( beam, target, en, "neutron", "159.6", model, NDF, Versions[m] );
         if ( NDF < 2 )
         {
            TText* txt2 = new TText( 0.1, 0.4, "Insufficient dataset" );
	    txt2->SetTextSize(0.2);
	    txt2->Draw();
	    break;
         }
         std::ostringstream os;
         os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
         std::string txt1 = "#chi^{2}/NDF = ";
         txt1 += os.str();
         txt1 += ( " for " + Versions[m] );
         TLatex* ltxt1 = new TLatex(0.10, 0.6-icounter1*0.2, txt1.c_str() );
         ltxt1->SetTextSize(0.17);
         ltxt1->Draw();
         icounter1++;
   }

   myc1->cd();

   std::string output1 = beam + "-" + target + "-" + en + "GeV-" + model; 
   output1 += "-regre.gif";
   
   myc1->Print( output1.c_str() );

   return;

}
