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

#include "NA61.h"
#include "ReadNA61Data.C"
#include "Chi2Calc.C"
#include "Chi2CalcNA61.C"
#include "DrawNA61Spectra.C"

// #include "../test19/G4MODELS_HighEnergy.h"

#ifndef G4VAL_PLOTNA61MODELS_C
#define G4VAL_PLOTNA61MODELS_C

void plotKPlus2PiPlusRatio(std::string beam, std::string target )
{

   TCanvas* myc = new TCanvas("myc","",800,400);
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatio( beam, target, "20", "140" );

   myc->cd(2);
   gPad->SetLeftMargin(0.15);
   drawKPlus2PiPlusRatio( beam, target, "140", "240" );
   
   myc->Print( "proton-C-31GeV-K-pi-ratio-models.gif" );

   return;

}

void plotSecondarySum2pages( std::string secondary )
{

   TCanvas *myc1 = new TCanvas("myc1","",900,900);

   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad3 = new TPad( "pad3", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad4 = new TPad( "pad4", "", 0.51, 0.14, 0.99, 0.56 );

   myc1->cd();
   pad1->Draw(); 
   pad1->cd();
   // gPad->SetLogy(); //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
//   drawMomSpectrum( "proton", "C", secondary, "0", "20" );
   drawMomSpectrum( "proton", "C", secondary, "0", "20", false, false );

   myc1->cd();
   pad2->Draw();
   pad2->cd();
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
//   drawMomSpectrum( "proton", "C", secondary, "20", "40" );
   drawMomSpectrum( "proton", "C", secondary, "20", "40", false, false );

   myc1->cd();
   pad3->Draw();
   pad3->cd();
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
//   drawMomSpectrum( "proton", "C", secondary, "40", "60" );
   drawMomSpectrum( "proton", "C", secondary, "40", "60", false, false );

   myc1->cd();
   pad4->Draw();
   pad4->cd();
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
//   drawMomSpectrum( "proton", "C", secondary, "60", "100" );
   drawMomSpectrum( "proton", "C", secondary, "60", "100", false, false );

   TCanvas *myc2 = new TCanvas("myc2","",900,900);

   TPad* pad5 = new TPad( "pad5", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad6 = new TPad( "pad6", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad7 = new TPad( "pad7", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad8 = new TPad( "pad8", "", 0.51, 0.14, 0.99, 0.56 );
   
   myc2->cd();
   pad5->Draw();
   pad5->cd();
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "100", "140" );

   myc2->cd();
   pad6->Draw();
   pad6->cd();
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "140", "180" );

   myc2->cd();
   pad7->Draw();
   pad7->cd();
   //gPad->SetLogx();
   gPad->SetLeftMargin(0.15);
   drawMomSpectrum( "proton", "C", secondary, "180", "240" );
   
   if ( secondary != "proton" ) 
   {
      myc2->cd();
      pad8->Draw();
      pad8->cd();
      //gPad->SetLogx();
      gPad->SetLeftMargin(0.15);
      drawMomSpectrum( "proton", "C", secondary, "240", "300" );
   }

   TLegend* leg = new TLegend(0.7,0.01,0.99,0.11);
   TLegendEntry* entry = 0;
   for ( int m=0; m<NModels_HE; ++m )
   {
      entry = leg->AddEntry( "", ModelName_HE[m].c_str(),"L" );
      entry->SetLineColor( ColorModel_HE[m] );
      entry->SetLineWidth(3);
      entry->SetTextFont(62);
   }
   entry = leg->AddEntry("","exp.data","p");
   entry->SetMarkerStyle(22);
   entry->SetMarkerColor(kBlack /* kBlue */ );
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(62);
   leg->SetFillColor(kWhite);
   myc1->cd();
   leg->Draw();
   myc2->cd();
   leg->Draw();

   TLatex* ltxt = new TLatex(0.1, 0.11, "MC vs NA61 Data; #chi^{2}/NDF calculated over ALL theta bins");
   ltxt->SetTextSize(0.025);
   myc1->cd();
   ltxt->Draw();
   myc2->cd();
   ltxt->Draw();
   
   for ( int m=0; m<NModels_HE; ++m )
   {

      double chi2 = 0.;
      int NDF = 0;
      std::ostringstream os;

      chi2 = Chi2MomSpectrumNA61Integral( "proton", "C", secondary, ModelName_HE[m], NDF );

      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = " + os.str();
      txt1 += ( " for " + ModelName_HE[m] );
      TLatex* ltxt1 = new TLatex( 0.1, 0.085-0.025*m, txt1.c_str() );
      ltxt1->SetTextSize(0.025);
      myc1->cd();
      ltxt1->Draw();
      myc2->cd();
      ltxt1->Draw();
      
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName_HE[m] << std::endl;

   }
   
   std::string output1 = "proton-C-31GeV-" + secondary + "-models-p1.gif";
   std::string output2 = "proton-C-31GeV-" + secondary + "-models-p2.gif";
   myc1->Print( output1.c_str() );
   myc2->Print( output2.c_str() );   

   return;

}

void plotSecondarySumCombined2pages( std::string secondary )
{

   TCanvas* myc   = new TCanvas( "myc", "", 900, 900 );
   
   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad2 = new TPad( "pad2", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad3 = new TPad( "pad3", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad4 = new TPad( "pad4", "", 0.51, 0.14, 0.99, 0.56 );
   
   pad1->Draw();
   pad1->Divide(1.,2.,0.,0.);
   pad1->cd(1); gPad->SetRightMargin(0.025);
//   drawMomSpectrum( "proton", "C", secondary, "0", "20", false, true, true ); // false==noPub2015 && true==drawLegend && largeLegend
   drawMomSpectrum( "proton", "C", secondary, "0", "20", false, false );
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "0", "20", false, false );
   
   myc->cd();

   pad2->Draw();
   pad2->Divide(1.,2.,0.,0.);
   pad2->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "20", "40", false, false );
   pad2->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "20", "40", false, false );

   myc->cd();

   pad3->Draw();
   pad3->Divide(1.,2.,0.,0.);
   pad3->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "40", "60", false, false );
   pad3->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "40", "60", false, false );

   myc->cd();
   
   pad4->Draw();
   pad4->Divide(1.,2.,0.,0.);
   pad4->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "60", "100", false, false );
   pad4->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "60", "100", false, false );

   TCanvas* myc1   = new TCanvas( "myc1", "", 900, 900 );

   TPad* pad5 = new TPad( "pad5", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad6 = new TPad( "pad6", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad7 = new TPad( "pad7", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad8 = new TPad( "pad8", "", 0.51, 0.14, 0.99, 0.56 );
   
   pad5->Draw();
   pad5->Divide(1.,2.,0.,0.);
   pad5->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "100", "140", false, false);
   pad5->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "100", "140", false, false );
   
   myc1->cd();

   pad6->Draw();
   pad6->Divide(1.,2.,0.,0.);
   pad6->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "140", "180", false, false );
   pad6->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "140", "180", false, false );

   myc1->cd();

   pad7->Draw();
   pad7->Divide(1.,2.,0.,0.);
   pad7->cd(1); gPad->SetRightMargin(0.025);
   drawMomSpectrum( "proton", "C", secondary, "180", "240", false, false );
   pad7->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "180", "240", false, false );
   
   if ( secondary != "proton" ) 
   {
      myc1->cd();
      pad8->Draw();
      pad8->Divide(1.,2.,0.,0.);
      pad8->cd(1); gPad->SetRightMargin(0.025);
      drawMomSpectrum( "proton", "C", secondary, "240", "300", false, false );
      pad8->cd(2); gPad->SetRightMargin(0.025);
      gPad->SetLogy();
      drawMomSpectMC2Data( "proton", "C", secondary, "240", "300", false, false );
   }

   TLegend* leg = new TLegend(0.7,0.01,0.99,0.11);
   TLegendEntry* entry = 0;
   for ( int m=0; m<NModels_HE; ++m )
   {
      entry = leg->AddEntry( "", ModelName_HE[m].c_str(),"L" );
      entry->SetLineColor( ColorModel_HE[m] );
      entry->SetLineWidth(3);
      entry->SetTextFont(62);
   }
   entry = leg->AddEntry("","exp.data","p");
   entry->SetMarkerStyle(22);
   entry->SetMarkerColor(kBlack /* kBlue */ );
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(62);
   leg->SetFillColor(kWhite);
   myc->cd();
   leg->Draw();
   myc1->cd();
   leg->Draw();

//   TLatex* ltxt = new TLatex(0.2, 0.11, "MC vs NA61 Data; #chi^{2}/NDF calculated over ALL theta bins");
   TLatex* ltxt = new TLatex(0.1, 0.10, "MC vs NA61 Data; #chi^{2}/NDF calculated over ALL theta bins");
//   ltxt->SetTextSize(0.025);
   ltxt->SetTextSize(0.02);
   myc->cd();
   ltxt->Draw();
   myc1->cd();
   ltxt->Draw();
   
   for ( int m=0; m<NModels_HE; ++m )
   {

      double chi2 = 0.;
      int NDF = 0;
      std::ostringstream os;

      chi2 = Chi2MomSpectrumNA61Integral( "proton", "C", secondary, ModelName_HE[m], NDF );

      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = " + os.str();
      txt1 += ( " for " + ModelName_HE[m] );
//      TLatex* ltxt1 = new TLatex( 0.2, 0.085-0.025*m, txt1.c_str() );
      TLatex* ltxt1 = new TLatex( 0.1, 0.08-0.02*m, txt1.c_str() );
//      ltxt1->SetTextSize(0.025);
      ltxt1->SetTextSize(0.02);
      myc->cd();
      ltxt1->Draw();
      myc1->cd();
      ltxt1->Draw();
      
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName_HE[m] << std::endl;

   }

   std::string output1 = "proton-C-31GeV-" + secondary + "-models-p1.gif";
   std::string output2 = "proton-C-31GeV-" + secondary + "-models-p2.gif";
   myc->Print( output1.c_str() );
   myc1->Print( output2.c_str() );

   return;

}

void plotModels2Data2015( std::string secondary )
{

   // updated data from run 2009 (pub.2015/2016)
   //

   // 1st page 
   //   
   TCanvas* myc1 = new TCanvas( "myc1", "", 900, 900 );
   
   // 2nd page 
   //   
   TCanvas* myc2 = new TCanvas( "myc2", "", 900, 900 );
   
   // 3rd page 
   //   
   TCanvas* myc3 = new TCanvas( "myc3", "", 900, 900 );

   TLatex* ltxt = new TLatex(0.1, 0.11, "MC vs NA61 Data; #chi^{2}/NDF calculated over ALL theta bins");
   ltxt->SetTextSize(0.02);
   myc1->cd();
   ltxt->Draw();
   myc2->cd();
   ltxt->Draw();
   myc3->cd();
   ltxt->Draw();

   for ( int m=0; m<NModels_HE; ++m )
   {
      
      double chi2 = 0.;
      int NDF = 0;
      std::ostringstream os;

      chi2 = Chi2NA61IntegralData2015( "proton", "C", secondary, ModelName_HE[m], NDF );

      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = " + os.str();
      txt1 += ( " for " + ModelName_HE[m] );

      TLatex* ltxt1 = new TLatex( 0.1, 0.085-0.025*m, txt1.c_str() );
      ltxt1->SetTextSize(0.025);
      myc1->cd();
      ltxt1->Draw();
      myc2->cd();
      ltxt1->Draw();
      myc3->cd();
      ltxt1->Draw();

      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for " << ModelName_HE[m] << std::endl;

   }


   TPad* pad11 = new TPad( "pad11", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad12 = new TPad( "pad12", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad13 = new TPad( "pad13", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad14 = new TPad( "pad14", "", 0.51, 0.14, 0.99, 0.56 );
   
   myc1->cd();   
   pad11->Draw();
   pad11->Divide(1.,2.,0.,0.);
   pad11->cd(1); gPad->SetRightMargin(0.025);
   // ---> for pions & protons ---> gPad->SetLogy();
   pad11->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad11->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "0", "10", true, false );
      pad11->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "0", "10", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad11->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "0", "20", true, false );
      pad11->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "0", "20", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad11->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "0", "40", true, false );
      pad11->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "0", "40", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad11->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "0", "40", true, false );
      pad11->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "0", "40", true, false );
   }
/*
   pad1->cd(2); gPad->SetRightMargin(0.025);
   gPad->SetLogy();
   drawMomSpectMC2Data( "proton", "C", secondary, "0", "10", false );
*/   

   myc1->cd();   
   pad12->Draw();
   pad12->Divide(1.,2.,0.,0.);
   pad12->cd(1); gPad->SetRightMargin(0.025);
   pad12->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "10", "20", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad12->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "10", "20", true, false );
      pad12->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "10", "20", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad12->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "20", "40", true, false );
      pad12->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "20", "40", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad12->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "40", "60", true, false );
      pad12->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "40", "60", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad12->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "40", "60", true, false );
      pad12->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "40", "60", true, false );
   }

   myc1->cd();   
   pad13->Draw();
   pad13->Divide(1.,2.,0.,0.);
   pad13->cd(1); gPad->SetRightMargin(0.025);
   pad13->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "20", "40", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad13->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "20", "40", true, false );
      pad13->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "20", "40", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad13->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "40", "60", true, false );
      pad13->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "40", "60", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad13->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "60", "100", true, false );
      pad13->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "60", "100", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad13->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "60", "100", true, false );
      pad13->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "60", "100", true, false );
   }

   myc1->cd();
   pad14->Draw();
   pad14->Divide(1.,2.,0.,0.);
   pad14->cd(1); gPad->SetRightMargin(0.025);
   pad14->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "40", "60", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad14->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "40", "60", true, false );
      pad14->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "40", "60", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad14->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "60", "100", true, false );
      pad14->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "60", "100", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad14->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "100", "140", true, false );
      pad14->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "100", "140", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad14->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "100", "140", true, false );
      pad14->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "100", "140", true, false );
   }
   
   TLegend* leg = new TLegend(0.7,0.01,0.99,0.11);
   TLegendEntry* entry = 0;
   for ( int m=0; m<NModels_HE; ++m )
   {
      entry = leg->AddEntry( "", ModelName_HE[m].c_str(),"L" );
      entry->SetLineColor( ColorModel_HE[m] );
      entry->SetLineWidth(3);
      entry->SetTextFont(62);
   }
   entry = leg->AddEntry("","exp.data","p");
   entry->SetMarkerStyle(22);
   entry->SetMarkerColor(kBlack /* kBlue */ );
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(62);
   leg->SetFillColor(kWhite);

   myc1->cd();
   leg->Draw();
   
   std::string output1 = "proton-C-31GeV-" + secondary + "-Pub2015-models-p1.gif";
   myc1->Print( output1.c_str() );

   // 2nd page 
   //
   
// ---> moved up   TCanvas* myc2 = new TCanvas( "myc2", "", 900, 900 );

   TPad* pad21 = new TPad( "pad21", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad22 = new TPad( "pad22", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad23 = new TPad( "pad23", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad24 = new TPad( "pad24", "", 0.51, 0.14, 0.99, 0.56 );
   
   myc2->cd();   
   pad21->Draw();
   pad21->Divide(1.,2.,0.,0.);
   pad21->cd(1); gPad->SetRightMargin(0.025);
   pad21->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "60", "100", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad21->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "60", "100", true, false );
      pad21->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "60", "100", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad21->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "100", "140", true, false );
      pad21->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "100", "140", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad21->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "140", "180", true, false );
      pad21->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "140", "180", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad21->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "140", "180", true, false );
      pad21->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "140", "180", true, false );
   }

   myc2->cd();   
   pad22->Draw();
   pad22->Divide(1.,2.,0.,0.);
   pad22->cd(1); gPad->SetRightMargin(0.025);
   pad22->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "100", "140", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad22->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "100", "140", true, false );
      pad22->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "100", "140", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad22->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "140", "180", true, false );
      pad22->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "140", "180", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad22->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "180", "240", true, false );
      pad22->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "180", "240", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad22->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "180", "240", true, false );
      pad22->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "180", "240", true, false );
   }

   myc2->cd();   
   pad23->Draw();
   pad23->Divide(1.,2.,0.,0.);
   pad23->cd(1); gPad->SetRightMargin(0.025);
   pad23->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "140", "180", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad23->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "140", "180", true, false );
      pad23->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "140", "180", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad23->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "180", "240", true, false );
      pad23->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "180", "240", true, false );
   }
   else if ( secondary == "lambda" )
   {
      pad23->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "240", "300", true, false );
      pad23->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "240", "300", true, false );
   }
   else if ( secondary == "k0s" )
   {
      pad23->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "240", "300", true, false );
      pad23->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "240", "300", true, false );
   }

   myc2->cd();
   pad24->Draw();
   pad24->Divide(1.,2.,0.,0.);
   pad24->cd(1); gPad->SetRightMargin(0.025);
   pad24->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   // drawMomSpectrum( "proton", "C", secondary, "180", "240", true, false );
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad24->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "180", "240", true, false );
      pad24->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "180", "240", true, false );
   }
   else if ( secondary == "kplus" || secondary == "kminus" )
   {
      pad24->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "240", "300", true, false );
      pad24->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "240", "300", true, false );
   }
   
   myc2->cd();
   leg->Draw();
   
   std::string output2 = "proton-C-31GeV-" + secondary + "-Pub2015-models-p2.gif";
   myc2->Print( output2.c_str() );
   
   if ( secondary == "kplus" || secondary == "kminus" || 
        secondary == "k0s"   || secondary == "lambda" ) return; 
   
   // 3rd page 
   //
   
// ---> moved up   TCanvas* myc3 = new TCanvas( "myc3", "", 900, 900 );

   TPad* pad31 = new TPad( "pad31", "", 0.01, 0.57, 0.49, 0.99 );
   TPad* pad32 = new TPad( "pad32", "", 0.01, 0.14, 0.49, 0.56 );
   TPad* pad33 = new TPad( "pad33", "", 0.51, 0.57, 0.99, 0.99 );
   TPad* pad34 = new TPad( "pad34", "", 0.51, 0.14, 0.99, 0.56 );
   
   myc3->cd();   
   pad31->Draw();
   pad31->Divide(1.,2.,0.,0.);
   pad31->cd(1); gPad->SetRightMargin(0.025);
   pad31->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad31->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "240", "300", true, false );
      pad31->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "240", "300", true, false );
   }
   
   myc3->cd();   
   pad32->Draw();
   pad32->Divide(1.,2.,0.,0.);
   pad32->cd(1); gPad->SetRightMargin(0.025);
   pad32->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   if ( secondary == "proton" || secondary == "piplus" || secondary == "piminus" )
   {
      pad32->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "300", "360", true, false );
      pad32->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "300", "360", true, false );
   }

   myc3->cd();   
   pad33->Draw();
   pad33->Divide(1.,2.,0.,0.);
   pad33->cd(1); gPad->SetRightMargin(0.025);
   pad33->cd(2); gPad->SetRightMargin(0.025); gPad->SetLogy();
   if ( secondary == "piplus" || secondary == "piminus" )
   {
      pad33->cd(1);
      drawMomSpectrum( "proton", "C", secondary, "360", "420", true, false );
      pad33->cd(2);
      drawMomSpectMC2Data( "proton", "C", secondary, "360", "420", true, false );
   }

   myc3->cd();
   leg->Draw();
   
   std::string output3 = "proton-C-31GeV-" + secondary + "-Pub2015-models-p3.gif";
   myc3->Print( output3.c_str() );
   
   return;

}

#endif
