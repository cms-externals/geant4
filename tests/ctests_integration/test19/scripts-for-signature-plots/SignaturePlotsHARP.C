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

#include "../test23/shared-root-macros/HARP.h"
#include "../test23/shared-root-macros/ReadHARPData.C"
// #include "Chi2Calc.C"
// #include "Chi2CalcHARP.C"
#include "../test23/shared-root-macros/DrawHARPSpectra.C"

// #include "../test23/shared-root-macros/REGRESSION_TEST.h"
#include "../test19/G4MODELS_IntermediateEnergy.h"

std::string TEST_NAME="test19";

void SignaturePlotsHARP()
{

   ReadHARPData( "piplus", "Pb", "5.0", "piplus", "FW" );

   TCanvas* myc1 = new TCanvas("myc1","",900.,700.);
      
   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.5, 0.5, 1.0 );
   myc1->cd(1);
   pad1->Draw();
   pad1->cd();
   // gPad->SetLogy();
   gPad->SetRightMargin(0.025);
   PlotHARPHisto( "piplus", "Pb", "5.0", "piplus", "FW", 3, false, false );
   TLatex* reac = new TLatex( 0.8, 800.,
           "5 GeV/c #pi^{+} + Pb #rightarrow #pi^{+} + X; 0.2 < #theta < 0.25 rad" );
   reac->SetTextSize(0.045);
   TText* data = new TText( 0.8, 750.,
           "Data: M.Apollonio et al., Nucl.Phys. A821 118, 2009" );
   data->SetTextSize(0.045);
   pad1->cd();
   reac->Draw();
   data->Draw();
   
   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.10, 0.5, 0.5 );
   myc1->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetRightMargin(0.025);
   // gPad->SetLogy();
   PlotHARPMC2Data( "piplus", "Pb", "5.0", "piplus", "FW", 3, false, true );
  
   ReadHARPData( "piplus", "Pb", "5.0", "piplus", "LA" );

   TPad* pad3 = new TPad( "pad3", "", 0.51, 0.5, 0.99, 1.0 );
   myc1->cd();
   pad3->Draw();
   pad3->cd();
   // gPad->SetLogy();
   gPad->SetRightMargin(0.025);
   PlotHARPHisto( "piplus", "Pb", "5.0", "piplus", "LA", 7, false, false );

   TLatex* reac1 = new TLatex( 0.12, 1950.,
           "5 GeV/c #pi^{+} + Pb #rightarrow #pi^{+} + X; 1.75 < #theta < 1.95 rad" );
   reac1->SetTextSize(0.045);
   TText* data1 = new TText(0.12, 1820.,"Data: M. Apollonio et al., Phys.Rev.C80 065207, 2009" );
   data1->SetTextSize(0.045);
   pad3->cd();
   reac1->Draw();
   data1->Draw();

   TPad* pad4 = new TPad( "pad4", "", 0.51, 0.10, 0.99, 0.5 );
   myc1->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetRightMargin(0.025);
   // gPad->SetLogy();
   PlotHARPMC2Data( "piplus", "Pb", "5.0", "piplus", "LA", 7, false, true );
   
   TLatex* cap = new TLatex(0.02, 0.07, 
                 "Double-differential cross-sections for #pi^{+} production in #pi^{+} + Pb interactions at 5GeV/c as a function of momentum");
   cap->SetTextSize(0.025);
   TText* cap1 = new TText(0.02, 0.04, "shown in different angular bins. Geant4.10.2.p01 (Bertini) results are compared with the HARP data.");
   cap1->SetTextSize(0.025);
   
   TLatex* xtitle1 = new TLatex( 0.12, 0.08, "momentum of secondary #pi^{+} [GeV/c]" );
   TLatex* xtitle2 = new TLatex( 0.62, 0.08, "momentum of secondary #pi^{+} [GeV/c]" );
   xtitle1->SetTextSize(0.025);
   xtitle2->SetTextSize(0.025);
   
   myc1->cd();
   // cap->Draw();
   // cap1->Draw();
   xtitle1->Draw();
   xtitle2->Draw();

   return;

}

