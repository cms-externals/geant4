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
#include "TGraphErrors.h"

#include "../test23/shared-root-macros/NA61.h"
#include "../test23/shared-root-macros/ReadNA61Data.C"
#include "../test23/shared-root-macros/DrawNA61Spectra.C"
#include "../test23/shared-root-macros/Chi2Calc.C"
#include "../test23/shared-root-macros/Chi2CalcNA61.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test23/G4PHYSLIST_HighEnergy.h"
#include "../test23/shared-root-macros/PlotNA61Regre.C"

std::string TEST_NAME="test23";

void NA61Regre()
{

  for ( int m=0; m<NModels_HE; ++m )
  {
     plotSecondarySumCombinedRegre( "piplus", ModelName_HE[m] ); 
     plotSecondarySumCombinedRegre( "piminus", ModelName_HE[m] ); 
     plotSecondarySumCombinedRegre( "proton", ModelName_HE[m] ); 
     plotKPlus2PiPlusRatioRegre( "proton", "C", ModelName_HE[m] );
  }

}

/*
void plotKPlus2PiPlusRatioRegre(std::string beam, std::string target, std::string model )
{

   TCanvas* myc = new TCanvas("myc","",800,400);
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatioRegre( beam, target, "20", "140", model );

   myc->cd(2);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatioRegre( beam, target, "140", "240", model );

   return;

}

void plotSecondarySumCombinedRegre( std::string secondary, std::string model )
{

   TCanvas *myc1 = new TCanvas("myc1","",900,900);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad3 = new TPad( "pad3", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad4 = new TPad( "pad4", "", 0.51, 0.14, 0.99, 0.56 );
   
   myc1->cd();
   pad1->Draw();
   pad1->Divide(1,2,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
   gPad->SetLogy(); //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "0", "20", model );
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "0", "20", model );

   myc1->cd();
   pad2->Draw();
   pad2->Divide(1,2,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   gPad->SetLogy(); //gPad->SetLogx();
   drawMomSpectrumRegre( "proton", "C", secondary, "20", "40", model );
   pad2->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "20", "40", model );

   myc1->cd();
   pad3->Draw();
   pad3->Divide(1,2,0.,0.);
   pad3->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrumRegre( "proton", "C", secondary, "40", "60", model );
   pad3->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "40", "60", model );

   myc1->cd();
   pad4->Draw();
   pad4->Divide(1.,2.,0.,0.);
   pad4->cd(1); gPad->SetRightMargin(0.025);   
   drawMomSpectrumRegre( "proton", "C", secondary, "60", "100", model );
   pad4->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "60", "100", model );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);

   TPad* pad5 = new TPad( "pad5", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad6 = new TPad( "pad6", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad7 = new TPad( "pad7", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad8 = new TPad( "pad8", "", 0.51, 0.14, 0.99, 0.56 );
   
   myc2->cd();
   pad5->Draw();
   pad5->Divide(1.,2.,0.,0.);
   pad5->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrumRegre( "proton", "C", secondary, "100", "140", model );
   pad5->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "100", "140", model );

   myc2->cd();
   pad6->Draw();
   pad6->Divide(1.,2.,0.,0.);
   pad6->cd(1); gPad->SetRightMargin(0.025);
   gPad->SetLeftMargin(0.15);
   drawMomSpectrumRegre( "proton", "C", secondary, "140", "180", model );
   pad6->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "140", "180", model );

   myc2->cd();
   pad7->Draw();
   pad7->Divide(1.,2.,0.,0.);
   pad7->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrumRegre( "proton", "C", secondary, "180", "240", model );
   pad7->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2DataRegre( "proton", "C", secondary, "180", "240", model );
   
   if ( secondary != "proton" ) return;
   {
      myc2->cd();
      pad8->Draw();
      pad8->Divide(1.,2.,0.,0.);
      pad8->cd(1); gPad->SetRightMargin(0.025);
      drawMomSpectrumRegre( "proton", "C", secondary, "240", "300", model );
      pad8->cd(2); gPad->SetRightMargin(0.025);
      gPad->SetLogy();
      drawMomSpectMC2DataRegre( "proton", "C", secondary, "240", "300", model );
   }

   return;

}
*/
