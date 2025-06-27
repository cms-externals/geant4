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

std::string TEST_NAME="test23";

#include "../test23/G4PHYSLIST_HighEnergy.h"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#include "../test23/shared-root-macros/NA61.h"
#include "../test23/shared-root-macros/ReadNA61Data.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA61.C"
#include "../test23/shared-root-macros/DrawNA61Spectra.C"

#include "../test23/shared-root-macros/PlotNA61Models.C"


void NA61Models()
{

   plotSecondarySumCombined2pages("piplus");
   plotSecondarySumCombined2pages("piminus");
   plotSecondarySumCombined2pages("proton");
   
   plotKPlus2PiPlusRatio( "proton", "C" );
   
   return;
   
}

// this is for J.V./s talk at NBI-2014 && K.G.'s talk at G4-Okinawa-2014

void plotForTalk()
{

   TCanvas *myc = new TCanvas("myc","",800,800);

   TPad* pad1 = new TPad( "pad1", "", 0.02, 0.15, 0.50, 0.99 );
   TPad* pad2 = new TPad( "pad2", "", 0.52, 0.15, 1.00, 0.99 );
   
   myc->cd();
   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1);
   gPad->SetRightMargin(0.01);
   drawMomSpectrum( "proton", "C", "piplus", "20", "40"); // true==drawLegend && largeLegend
   pad1->cd(2); 
   gPad->SetRightMargin(0.01);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", "piplus", "20", "40" );// , false );

   myc->cd();
   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1);
   gPad->SetRightMargin(0.01);
   drawMomSpectrum( "proton", "C", "piplus", "60", "100" ); // true==drawLegend && largeLegend
   pad2->cd(2); 
   gPad->SetRightMargin(0.01);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", "piplus", "60", "100"); //, false );
   
   double chi2=0.;
   int NDF=0;
   for ( int m=0; m<NModels_HE; ++m )
   {

      chi2 = 0.;
      NDF = 0;
      std::ostringstream os1;
      readMomSpectrum( "proton", "C", "piplus", "20", "40" );
      chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "piplus", "20", "40", ModelName_HE[m], NDF );
      if ( (chi2/NDF) < 100 )
      {
         os1.precision(3);
      }
      else
      {
         os1.precision(4);
      }      
      os1 << (chi2/NDF);

//      std::cout << "chi2/ndf = " << chi2 << "/" << NDF << " = " << (chi2/NDF) << " for " << ModelName_HE[m] << std::endl; 

      std::string txt1 = "#chi^{2}/NDF = " + os1.str();
      txt1 += ( " for " + ModelName_HE[m] );
      TLatex* ltxt1 = new TLatex( 0.05, 0.12-0.025*m, txt1.c_str() );
      ltxt1->SetTextSize(0.025);
      myc->cd();
      ltxt1->Draw();

      chi2 = 0.;
      NDF = 0;
      std::ostringstream os2;
      readMomSpectrum( "proton", "C", "piplus", "60", "100" );
      chi2 = Chi2MomSpectrumNA61ThetaBin( "proton", "C", "piplus", "60", "100", ModelName_HE[m], NDF );
      if ( (chi2/NDF) < 100 )
      {
         os2.precision(3);
      }
      else
      {
         os2.precision(4);
      }      
      os2 << (chi2/NDF);

//      std::cout << "chi2/ndf = " << chi2 << "/" << NDF << " = " << (chi2/NDF) << " for " << ModelName_HE[m] << std::endl; 

      std::string txt2 = "#chi^{2}/NDF = " + os2.str();
      txt2 += ( " for " + ModelName_HE[m] );
      TLatex* ltxt2 = new TLatex( 0.55, 0.12-0.025*m, txt2.c_str() );
      ltxt2->SetTextSize(0.025);
      myc->cd();
      ltxt2->Draw();

   }

   return;

}
