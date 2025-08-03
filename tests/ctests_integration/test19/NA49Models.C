
/* remove DbReader, revert to using local ASCII files
R__LOAD_LIBRARY(libDbReader.so);
*/
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

/* remove DbReader, revert to using local ASCII files
// #include "../DbReader/DbReader.h"
#include "/g4/g4p/pbs/g4-had-validation/g4-releases/geant4-10-02-patch-02/tests/DbReader/DbReader.h"
*/

std::string TEST_NAME="test19";
#include "../test19/G4MODELS_HighEnergy.h"

#include "../test23/shared-root-macros/NA49.h"
#include "../test23/shared-root-macros/ReadNA49Data.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA49.C"
#include "../test23/shared-root-macros/DrawNA49Spectra.C"

#include "../test23/shared-root-macros/PlotNA49Models.C"

void NA49Models()
{

   plot_dNdxF_pT( "proton", "C", "piplus" );
   plot_dNdxF_pT( "proton", "C", "piminus" );
   plot_dNdxF_pT( "proton", "C", "proton" );
   plot_dNdxF_pT( "proton", "C", "antiproton" );
   plot_dNdxF_pT( "proton", "C", "neutron" );

   return;

}

/*
void plot_dNdxF_pT( std::string beam, std::string target, std::string secondary )
{

   if ( secondary == "neutron" )
   {
      TCanvas* myc = new TCanvas( "myc", "", 600, 600 );
      TPad* pad1 = new TPad( "pad1", "", 0.01, 0.12, 0.99, 0.99 );
      pad1->Draw();
      pad1->Divide(1.,2.,0.,0.);
      pad1->cd(1); gPad->SetRightMargin(0.025);
      drawIntegratedSpectrum( beam, target, secondary, "dNdxF" ); 
      pad1->cd(2); gPad->SetRightMargin(0.025);
      gPad->SetLogy();
      drawIntSpectrumMC2Data( beam, target, secondary, "dNdxF" );
      myc->cd();
      TPad* pad3 = new TPad("pad3","",0.01, 0.01, 0.99, 0.11);
      pad3->Draw();
      for ( int m=0; m<NModels_HE; ++m )
      {
         int NDF1=0;
         int NDF2 = 0;
         double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "dNdxF", ModelName_HE[m], NDF1 );  
         std::ostringstream os1;
         os1 << (chi2/NDF1);
         std::string txt1 = " #chi^{2}/NDF = " + os1.str();
         txt1 += ( " for " + ModelName_HE[m] + " vs NA49 Data" );
         TLatex* ltxt1 = new TLatex( 0.10, 0.75-m*0.3, txt1.c_str() );
         ltxt1->SetTextSize(0.25);
         pad3->cd();
         ltxt1->Draw();  
      } 
      myc->cd();
      myc->Print( "proton-C-158GeV-neutron-models.gif" );
      return;      
   }
      
   TCanvas* myc   = new TCanvas("myc","",1000,600);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.12, 0.49, 0.99 );
   
   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, secondary, "dNdxF" ); 
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   if ( secondary == "proton" || secondary == "antiproton" )
   {
     drawIntSpectrumMC2Data( beam, target, secondary, "dNdxF", false ); // no legent for mc/data for (anti)protons
   }
   else
   {
      drawIntSpectrumMC2Data( beam, target, secondary, "dNdxF" );
   }

   myc->cd();
      
   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.12, 0.99, 0.99 );

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawIntegratedSpectrum( beam, target, secondary, "pT", false ); // no legend because it overlaps with the contents
   pad2->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawIntSpectrumMC2Data( beam, target, secondary, "pT" );
   
   myc->cd();

   TPad* pad3 = new TPad("pad3","",0.01, 0.01, 0.49, 0.11);
   pad3->Draw();
   // pad3->cd();
   myc->cd();
   TPad* pad4 = new TPad("pad4","",0.51, 0.01, 0.99, 0.11);
   pad4->Draw();
   // pad4->cd();
   for ( int m=0; m<NModels_HE; ++m )
   {
      int NDF1=0;
      int NDF2 = 0;
      double chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "dNdxF", ModelName_HE[m], NDF1 );  
      std::ostringstream os1;
      os1 << (chi2/NDF1);
      std::string txt1 = " #chi^{2}/NDF = " + os1.str();
      txt1 += ( " for " + ModelName_HE[m] + " vs NA49 Data" );
      TLatex* ltxt1 = new TLatex( 0.10, 0.75-m*0.3, txt1.c_str() );
      ltxt1->SetTextSize(0.25);
      pad3->cd();
      ltxt1->Draw();  
      chi2 = Chi2IntegratedSpectrumNA49( beam, target, secondary, "pT", ModelName_HE[m], NDF2 );
      std::ostringstream os2;
      os2 << (chi2/NDF2);
      std::string txt2 = " #chi^{2}/NDF = " + os2.str();
      txt2 += ( " for " + ModelName_HE[m] + " vs NA49 Data" );
      TLatex* ltxt2 = new TLatex( 0.10, 0.70-m*0.3, txt2.c_str() );
      ltxt2->SetTextSize(0.25);
      pad4->cd();
      ltxt2->Draw();
   }

   myc->cd();
   
   std::string output = beam + "-" + target + "-158GeV-" + secondary;
   output += "-models.gif";
   
   myc->Print( output.c_str() );

   return;

} 
*/
