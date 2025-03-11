#include "TLegend.h"
#include <string>

#include "../test23/shared-root-macros/Chi2Calc.C"

#include "../test23/shared-root-macros/REGRESSION_TEST.h"

#ifndef G4VAL_PLOTGAMMAN_C
#define G4VAL_PLOTGAMMAN_C

#include "../test75/GammaData/gammaCu_300MeV.C"
#include "../test75/GammaData/gamma_400_750MeV.C"

void PlotGammaNuclear( std::string, std::string, std::string, std::string );
void PlotGammaNRegression( std::string, std::string, std::string, std::string );
double calcChi2GammaNuclear( std::string, std::string, 
                             std::string, std::string, int&, 
			     std::string version="." );

void OverlayData(const TString&, const TString&, TLegend* leg=0 );
TGraphErrors* GetOverlayCu( const TString& ); 
TGraphErrors* GetOverlayPb( const TString& );

void PlotGamma300Cu() 
{
   gStyle->SetOptStat(0);
   
   int NDF = 0;   
   double chi2 = 0.;

   TCanvas* myc = new TCanvas( "myc", "", 800, 800 );   
   // myc->Divide(2,2);
   
   // myc->cd(1);
   TPad* pad1 = new TPad("pad1","", 0.01, 0.55, 0.49, 0.99);
   myc->cd();
   pad1->Draw();   
   pad1->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "45" );
   chi2 += calcChi2GammaNuclear( "300", "Cu", "proton", "45", NDF );

   // myc->cd(2);
   TPad* pad2 = new TPad("pad2","",0.51,0.55,0.99,0.99);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "90" );
   chi2 += calcChi2GammaNuclear( "300", "Cu", "proton", "90", NDF );

   //myc->cd(3);
   TPad* pad3 = new TPad("pad3","", 0.01, 0.10, 0.49, 0.54 );
   myc->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "135" );
   chi2 += calcChi2GammaNuclear( "300", "Cu", "proton", "135", NDF );

   // myc->cd(4);
   TPad* pad4 = new TPad("pad4","", 0.51, 0.10, 0.99, 0.54);
   myc->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "150" );
   chi2 += calcChi2GammaNuclear( "300", "Cu", "proton", "150", NDF );

   TPad* pad5 = new TPad("pad5","", 0.01, 0.01, 0.99, 0.09 );
   myc->cd();
   pad5->Draw();
   pad5->cd();
   std::ostringstream os;
   os << (chi2/NDF);
   std::string txt1 = "Bertini vs Data (R. Schumacher et al., Phys. Rev. C 25, 2269 (1982))";
   std::string txt2 = "Integral over all agnular bins: #chi^{2}/NDF = ";
   txt2 += os.str();
   TLatex* ltxt1 = new TLatex(0.10, 0.7, txt1.c_str() );
   ltxt1->SetTextSize(0.25);
   ltxt1->Draw();
   TLatex* ltxt2 = new TLatex(0.1, 0.4, txt2.c_str() );
   ltxt2->SetTextSize(0.25);
   ltxt2->Draw();

   myc->Print("gamma-300-Cu-p-models.gif");


   return;
   
}

void PlotGamma300CuRegre() 
{
   gStyle->SetOptStat(0);
   
   TCanvas* myc = new TCanvas( "myc", "", 800, 800 );   
   TPad* pad1 = new TPad("pad1","", 0.01, 0.55, 0.49, 0.99);
   TPad* pad2 = new TPad("pad2","", 0.51, 0.55, 0.99, 0.99);
   TPad* pad3 = new TPad("pad3","", 0.01, 0.10, 0.49, 0.54 );
   TPad* pad4 = new TPad("pad4","", 0.51, 0.10, 0.99, 0.54);

   // myc->Divide(2,2);
   
   // myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "300", "Cu", "proton", "45" );

   // myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "300", "Cu", "proton", "90" );

   // myc->cd(3);
   myc->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "300", "Cu", "proton", "135" );

   // myc->cd(4);
   myc->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "300", "Cu", "proton", "150" );

   TPad* pad5 = new TPad("pad5","", 0.01, 0.01, 0.99, 0.09 );
   myc->cd();
   pad5->Draw();
   pad5->cd();

   std::string txt = "Bertini vs Data (R. Schumacher et al., Phys. Rev. C 25, 2269 (1982))";
   TLatex* ltxt = new TLatex( 0.1, 0.8, txt.c_str() );
   ltxt->SetTextSize(0.25);
   ltxt->Draw();
   for ( int m=0; m<NVersions; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += calcChi2GammaNuclear( "300", "Cu", "proton",  "45", NDF, Versions[m] );
      chi2 += calcChi2GammaNuclear( "300", "Cu", "proton",  "90", NDF, Versions[m] );
      chi2 += calcChi2GammaNuclear( "300", "Cu", "proton", "135", NDF, Versions[m] );
      chi2 += calcChi2GammaNuclear( "300", "Cu", "proton", "150", NDF, Versions[m] );
      std::ostringstream os;
      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Versions[m] );
      TLatex* ltxt1 = new TLatex(0.10, 0.6-m*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.25);
      ltxt1->Draw();
   }

   myc->Print("gamma-300-Cu-p-bert-regre.gif");

   return;
   
}

void PlotGamma668Cu()
{

   gStyle->SetOptStat(0);

   double chi2 = 0.;
   int NDF = 0;

   TCanvas* myc = new TCanvas( "myc", "", 800, 800 );   
   TPad* pad1 = new TPad("pad1","", 0.01, 0.55, 0.49, 0.99);
   TPad* pad2 = new TPad("pad2","", 0.01, 0.10, 0.49, 0.54 );
   TPad* pad3 = new TPad("pad3","", 0.51, 0.55, 0.99, 0.99);
   TPad* pad4 = new TPad("pad4","", 0.51, 0.10, 0.99, 0.54);

   // myc->Divide(2,2);
   
   // myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi-", "28" );
   chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-", "28", NDF );

   // myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi-", "44" );   
   chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-", "44", NDF );

   TPad* pad5 = new TPad("pad5","", 0.01, 0.01, 0.49, 0.09 );
   myc->cd();
   pad5->Draw();
   pad5->cd();
   std::ostringstream os;
   os << (chi2/NDF);
   std::string txt1 = "Bertini vs Data (K. Baba et al., Nucl. Phys. A322, 349 (1979))";
   std::string txt2 = "Integral over agnular bins: #chi^{2}/NDF = ";
   txt2 += os.str();
   TLatex* ltxt1 = new TLatex(0.01, 0.7, txt1.c_str() );
   ltxt1->SetTextSize(0.22);
   ltxt1->Draw();
   TLatex* ltxt2 = new TLatex(0.01, 0.4, txt2.c_str() );
   ltxt2->SetTextSize(0.22);
   ltxt2->Draw();

   
   chi2 = 0.;
   NDF=0;
   // myc->cd(3);
   myc->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi+", "28" );
   chi2 += calcChi2GammaNuclear( "668", "Cu", "pi+", "28", NDF );

   // myc->cd(4);
   myc->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi+", "44" );
   chi2 += calcChi2GammaNuclear( "668", "Cu", "pi+", "44", NDF );

   TPad* pad6 = new TPad("pad6","", 0.51, 0.01, 0.99, 0.09 );
   myc->cd();
   pad6->Draw();
   pad6->cd();
   std::ostringstream os1;
   os1 << (chi2/NDF);
//   std::string txt1 = "Bertini cascade vs Data (R. Schumacher et al., Phys. Rev. C 25, 2269 (1982))";
   std::string txt3 = "Integral over all agnular bins: #chi^{2}/NDF = ";
   txt3 += os.str();
//   TLatex* ltxt1 = new TLatex(0.10, 0.7, txt1.c_str() );
//   ltxt1->SetTextSize(0.25);
   ltxt1->Draw();
   TLatex* ltxt3 = new TLatex(0.01, 0.4, txt3.c_str() );
   ltxt3->SetTextSize(0.22);
   ltxt3->Draw();

   myc->Print("gamma-668-Cu-pions-models.gif" );

   return;

}

void PlotGamma668CuRegre()
{

   gStyle->SetOptStat(0);
   
   double chi2=0.;
   int NDF=0;

   TCanvas* myc = new TCanvas( "myc", "", 800, 800 ); 
     
   TPad* pad1 = new TPad("pad1","", 0.01, 0.55, 0.49, 0.99);
   TPad* pad2 = new TPad("pad2","", 0.01, 0.10, 0.49, 0.54 );
   TPad* pad3 = new TPad("pad3","", 0.51, 0.55, 0.99, 0.99);
   TPad* pad4 = new TPad("pad4","", 0.51, 0.10, 0.99, 0.54);

//   myc->Divide(2,2);
   
   // myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "668", "Cu", "pi-", "28" );

   // myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "668", "Cu", "pi-", "44" );   

   TPad* pad5 = new TPad("pad5","", 0.01, 0.01, 0.49, 0.09 );
   myc->cd();
   pad5->Draw();
   pad5->cd();

   std::string txt = "Bertini vs Data (K. Baba et al., Nucl. Phys. A322, 349 (1979))";
   TLatex* ltxt = new TLatex( 0.01, 0.8, txt.c_str() );
   ltxt->SetTextSize(0.22);
   ltxt->Draw();
   for ( int m=0; m<NVersions; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-",  "28", NDF, Versions[m] );
      chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-",  "44", NDF, Versions[m] );
      std::ostringstream os;
      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Versions[m] );
      TLatex* ltxt1 = new TLatex(0.10, 0.6-m*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.22);
      ltxt1->Draw();
   }

   // myc->cd(3);
   myc->cd();
   pad3->Draw();
   pad3->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "668", "Cu", "pi+", "28" );

   // myc->cd(4);
   myc->cd();
   pad4->Draw();
   pad4->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "668", "Cu", "pi+", "44" );

   TPad* pad6 = new TPad("pad5","", 0.51, 0.01, 0.99, 0.09 );
   myc->cd();
   pad6->Draw();
   pad6->cd();

//   std::string txt = "Bertini vs Data (K. Baba et al., Nucl. Phys. A322, 349 (1979))";
//   TLatex* ltxt = new TLatex( 0.01, 0.8, txt.c_str() );
//   ltxt->SetTextSize(0.2);
   ltxt->Draw();
   for ( int m=0; m<NVersions; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-",  "28", NDF, Versions[m] );
      chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-",  "44", NDF, Versions[m] );
      std::ostringstream os;
      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Versions[m] );
      TLatex* ltxt1 = new TLatex(0.10, 0.6-m*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.22);
      ltxt1->Draw();
   }

   myc->Print("gamma-668-Cu-pions-bert-regre.gif");

   return;


}

void PlotGamma668Pb()
{

   gStyle->SetOptStat(0);

   double chi2 = 0.;
   int NDF = 0;
   
   std::string txt = "Bertini vs Data (K. Baba et al., Nucl. Phys. A322, 349 (1979))";
   TLatex* ltxt = new TLatex( 0.01, 0.8, txt.c_str() );
   ltxt->SetTextSize(0.25);

   TCanvas* myc = new TCanvas( "myc", "", 800, 500 );   
   TPad* pad1 = new TPad("pad1","", 0.01, 0.14, 0.49, 0.99);
   TPad* pad2 = new TPad("pad2","", 0.51, 0.14, 0.99, 0.99);

   //myc->Divide(2,1);
   
   // myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Pb", "pi-", "44" );
   chi2 = calcChi2GammaNuclear( "668", "Cu", "pi-",  "44", NDF );
   
   
   TPad* pad3 = new TPad("pad3","", 0.01, 0.01, 0.49, 0.13 );
   myc->cd();
   pad3->Draw();
   pad3->cd();

   ltxt->Draw();

   std::ostringstream os1;
   os1 << (chi2/NDF);
   std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
   txt1 += os1.str();
   TLatex* ltxt1 = new TLatex(0.1, 0.5, txt1.c_str() );
   ltxt1->SetTextSize(0.25);
   ltxt1->Draw();
   
   NDF = 0;
   chi2 = 0.;
   // myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Pb", "pi+", "44" );   
   chi2 = calcChi2GammaNuclear( "668", "Cu", "pi+",  "44", NDF );

   TPad* pad4 = new TPad("pad4","", 0.51, 0.01, 0.99, 0.13 );
   myc->cd();
   pad4->Draw();
   pad4->cd();

   ltxt->Draw();
   
   std::ostringstream os2;
   os2 << (chi2/NDF);
   std::string txt2 = "Integral over all agnular bins: #chi^{2}/NDF = ";
   txt2 += os2.str();
   os2 << (chi2/NDF);
   TLatex* ltxt2 = new TLatex(0.1, 0.5, txt2.c_str() );
   ltxt2->SetTextSize(0.25);
   ltxt2->Draw();
   
   myc->Print("gamma-668-Pb-pions-models.gif" );

   return;

}

void PlotGamma668PbRegre()
{

   gStyle->SetOptStat(0);

   std::string txt = "Bertini vs Data (K. Baba et al., Nucl. Phys. A322, 349 (1979))";
   TLatex* ltxt = new TLatex( 0.01, 0.85, txt.c_str() );
   ltxt->SetTextSize(0.2);

   TCanvas* myc = new TCanvas( "myc", "", 800, 450 );   
   TPad* pad1 = new TPad("pad1","", 0.01, 0.17, 0.49, 0.99);
   TPad* pad2 = new TPad("pad2","", 0.51, 0.17, 0.99, 0.99);
// calcChi2GammaNuclear( "668", "Cu", "pi-",  "44", NDF, Versions[m] );
   // myc->Divide(2,1);
   
   // myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "668", "Pb", "pi-", "44" );

   TPad* pad3 = new TPad("pad3","", 0.01, 0.01, 0.49, 0.17 );
   myc->cd();
   pad3->Draw();
   pad3->cd();
   ltxt->Draw();
   for ( int m=0; m<NVersions; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += calcChi2GammaNuclear( "668", "Cu", "pi-",  "44", NDF, Versions[m] );
      std::ostringstream os;
      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Versions[m] );
      TLatex* ltxt1 = new TLatex(0.1, 0.65-m*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.2);
      ltxt1->Draw();
   }

   // myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   PlotGammaNRegression( "668", "Pb", "pi+", "44" );   

   TPad* pad4 = new TPad("pad4","", 0.51, 0.01, 0.99, 0.17 );
   myc->cd();
   pad4->Draw();
   pad4->cd();
   ltxt->Draw();
   for ( int m=0; m<NVersions; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      chi2 += calcChi2GammaNuclear( "668", "Cu", "pi+",  "44", NDF, Versions[m] );
      std::ostringstream os;
      os << (chi2/NDF);
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Versions[m] );
      TLatex* ltxt1 = new TLatex(0.1, 0.65-m*0.2, txt1.c_str() );
      ltxt1->SetTextSize(0.2);
      ltxt1->Draw();
   }

   myc->Print("gamma-668-Pb-pions-bert-regre.gif");

   return;

}

void PlotGammaNRegression( std::string energy, std::string target, std::string secondary, std::string angle )
{

   gStyle->SetOptStat(0);
   
   std::string fName = "gamma"+ energy + "MeV" + target + "Bertini.root"; 
   
   std::string hName = "";

//   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   TLegend* leg1 = new TLegend(0.25, 0.15, 0.55, 0.35);
   
   if ( secondary == "proton" )
   {
      hName = "p";
   } 
   else if ( secondary == "pi-" )
   {
      hName = "pim";
   }
   else if (secondary == "pi+" )
   {
      hName = "pip";
   }
   hName += angle ;
   hName +="deg";

   std::string dir = "";
   
   for ( int iv=0; iv<NVersions; iv++ )
   {

//      dir = "";
//      if ( Versions[iv] != CurrentVersion ) 
//      {
//         dir = regre_test_dir + Versions[iv] + "/";
//      }
//      std::string hFileName = dir + fName ;

   std::string location = "";
   if ( Versions[iv] == CurrentVersion || Versions[iv] == "." )
   {
         location = "";
   }
   else
   {
         location = regre_test_dir + "/test75/" + Versions[iv] + "/";
   }
   std::string hFileName = location + fName; 


      TFile* hfile = new TFile( hFileName.c_str() );
      TH1F* histo1 = (TH1F*)hfile->Get( hName.c_str() );   
      histo1->SetLineColor(ColorVersion[iv]);
/*      
      if ( iv == 0 )
      {
      histo1->SetLineWidth(6);
      }
      else
      {
         histo1->SetLineWidth(3);
      }
*/
      histo1->SetLineWidth(6-iv);
      TString hTit1 = "gamma + " + target;
      TString hTit2 = " X + " + secondary + " (" + angle + "deg)";
      histo1->SetTitle( hTit1 + " #rightarrow " + hTit2 );
      histo1->GetXaxis()->SetTitleSize(0.04);
      histo1->GetXaxis()->CenterTitle();
      if ( secondary == "proton" )
      {
         histo1->GetXaxis()->SetTitle("kinetic energy of secondary proton (MeV)");
         histo1->GetYaxis()->SetTitle("d^{2}#sigma / dE d#Omega  [#mub MeV^{-1} sr^{-1}]");
      }
      else if ( secondary == "pi-" || secondary == "pi+" )
      {
         histo1->GetYaxis()->SetTitle("d^{2}#sigma / dp d#Omega  [#mub MeV^{-1} sr^{-1}]"); 
	 histo1->GetYaxis()->SetRangeUser(0.001, 10. );  
      }
      histo1->GetYaxis()->SetTitleOffset(0.9);
      histo1->GetYaxis()->SetTitleSize(0.05);
      if ( iv == 0 )
      {
	 histo1->Draw();
      }
      else
      {
         histo1->Draw("same");
      }
      leg1->AddEntry( histo1, Versions[iv].c_str(), "L" );
   }

   OverlayData( hName, target, leg1 );
   leg1->Draw();
   leg1->SetFillColor(kWhite);

   return;

} 

double calcChi2GammaNuclear( std::string energy, std::string target, 
                             std::string secondary, std::string angle, int& NDF, 
			     std::string version="." )
{

   double chi2 = 0.;

   std::string fName = "";
//   if ( version != CurrentVersion )
//   {
//      fName = regre_test_dir + version + "/";
//   }

   std::string location = "";
   if ( version == CurrentVersion || version == "." )
   {
         location = "";
   }
   else
   {
         location = regre_test_dir + "/test75/" + version + "/";
   }
   fName = location + "gamma" + energy + "MeV" + target;
//   fName += "gamma" + energy + "MeV" + target;
   fName += "Bertini.root";
   TFile* ifile = new TFile( fName.c_str() );

   std::string hName = "";
   if ( secondary == "proton" )
   {
      hName = "p";
   } 
   else if ( secondary == "pi-" )
   {
      hName = "pim";
   }
   else if (secondary == "pi+" )
   {
      hName = "pip";
   }
   hName += angle ;
   hName +="deg";
   
   TH1F* histo1 = (TH1F*)ifile->Get( hName.c_str() );

   TGraphErrors* gdata = 0;

   if (  target == "Cu" )
   {
      gdata = GetOverlayCu( hName );
   }
   else
   {
      gdata = GetOverlayPb(hName);
   }
   if ( !gdata ) return 0.;

   chi2 += Chi2( gdata, histo1, NDF );

   return chi2;

}


void PlotGammaNuclear( std::string energy, std::string target, std::string secondary, std::string angle )
{
   
// FIXME !!!
// Will need an update In case of several models
//

   std::string fName = "gamma"+ energy + "MeV" + target;
   
   std::string fName1 = fName + "Bertini.root";
   TFile* ifile1 = new TFile( fName1.c_str() );
   
// In case of CHIPS...
// -->    std::string fName2 = fName + "CHIPS.root";
// -->   TFile* ifile2 = new TFile( fName2.c_str() );
   
//   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   TLegend* leg1 = new TLegend(0.25, 0.15, 0.55, 0.35);
   
   std::string hName = "";
   
   if ( secondary == "proton" )
   {
      hName = "p";
   } 
   else if ( secondary == "pi-" )
   {
      hName = "pim";
   }
   else if (secondary == "pi+" )
   {
      hName = "pip";
   }
   hName += angle ;
   hName +="deg";
      
   TH1F* histo1 = (TH1F*)ifile1->Get( hName.c_str() );   
   histo1->SetLineColor(kRed);
   histo1->SetLineWidth(3);
   TString hTit1 = "gamma + " + target;
   TString hTit2 = " X + " + secondary + " (" + angle + "deg)";
   histo1->SetTitle( hTit1 + " #rightarrow " + hTit2 );
   
   histo1->GetXaxis()->SetTitleSize(0.04);
   histo1->GetXaxis()->CenterTitle();
   
   if ( secondary == "proton" )
   {
      histo1->GetXaxis()->SetTitle("kinetic energy of secondary proton (MeV)");
      histo1->GetYaxis()->SetTitle("d^{2}#sigma / dE d#Omega  [#mub MeV^{-1} sr^{-1}]");
   }
   else if ( secondary == "pi-" || secondary == "pi+" )
   {
      histo1->GetYaxis()->SetTitle("d^{2}#sigma / dp d#Omega  [#mub MeV^{-1} sr^{-1}]");   
   }
   histo1->GetYaxis()->SetTitleOffset(0.9);
   histo1->GetYaxis()->SetTitleSize(0.05);
   // histo1->GetYaxis()->CenterTitle();
   histo1->Draw();
   leg1->AddEntry( histo1, "Bertini", "L" );
   OverlayData( hName, target, leg1 );
   leg1->Draw();
   leg1->SetFillColor(kWhite);
   
   return;

}

TGraphErrors* GetOverlayCu( const TString& hname ) 
{

    if (hname == "p45deg")  return gCup45();
    if (hname == "p90deg")  return gCup90();
    if (hname == "p135deg") return gCup135();
    if (hname == "p150deg") return gCup150();
    
    if ( hname == "pim28deg" ) return gGam668Cu_PiMin_28();
    if ( hname == "pip28deg" ) return gGam668Cu_PiPls_28();
    if ( hname == "pim44deg" ) return gGam668Cu_PiMin_44();
    if ( hname == "pip44deg" ) return gGam668Cu_PiPls_44();

    return 0;
    
}

TGraphErrors* GetOverlayPb( const TString& hname ) 
{

   if ( hname == "pim44deg" ) return gGam668Pb_PiMin_44();
   if ( hname == "pip44deg" ) return gGam668Pb_PiPls_44();

   return 0;

}

void OverlayData(const TString& hname, const TString& target, TLegend* leg=0 ) 
{

  TGraph* data = 0;
  if ( target == "Cu" )
  {
     data = GetOverlayCu(hname);
  }
  else if ( target == "Pb" )
  {
     data = GetOverlayPb(hname);
  }
  else if ( target == "C" ) // provision for future addition of gamma + C data
  {
     // data = GetOverlayC(hname);
  }
  
  if ( data )
  {
     data->SetMarkerStyle(22);
     data->SetMarkerColor(kBlue);
     data->SetMarkerSize(1.6);
     if ( leg ) leg->AddEntry( data, "exp.data", "p" );
     data->Draw("P same");
  }
  
  return;

}

#endif
