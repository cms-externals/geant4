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

#include "../test23/shared-root-macros/Chi2Calc.C"

void drawMCvsData( std::string, std::string, std::string, std::string );
void drawMCvsDataRegre( std::string, std::string, std::string, std::string, std::string );
void drawMC2DataRegre( std::string, std::string, std::string, std::string, std::string );
void drawMC2Data( std::string, std::string, std::string, std::string );
TGraphErrors* getGraphData2Data();
TGraphErrors* getGraphMC2Data( TH1F* );
void setAxisStyle( TAxis*, const char* );

int NDataPoints = 0;
const int NMaxPoints = 30; // it's larger that one would need... 
                           // but I don't want to deal with memory allocation
double Mom[30], XSec[30], Err[30]; // EStat[30], ESys[30]; 

const int NModels = 1;
std::string Models[2] = { "ftfp", "bertini" };
int ColorModel[3] = { 6, 3, 14 };

#include "../test23/shared-root-macros/SplitString.C"
#include "../test23/shared-root-macros/REGRESSION_TEST.h"

string TEST_NAME="test47";

bool readPbarProdData( std::string beam, std::string target, std::string energy,
                       std::string sec_angle,
		       bool useSysEr=true, bool useStatEr=true )
{

   
   std::string datadir = "../test47/pbar-prod";
   datadir += "/" + beam + "/anti_proton/" ;
   std::string filename = target + energy + "GeV" + sec_angle + "deg.dat";
   std::string file = datadir + filename;

   NDataPoints=0;
   for ( int i=0; i<NMaxPoints; ++i )
   {
      Mom[i]  = 0.;
      XSec[i] = 0.;
      Err[i]  = 0.;
   }
      
   ifstream infile;
   infile.open(file.c_str());

   if ( !infile.is_open() ) 
   {
      std::cout << " data file " << file << " does not exist " << std::endl;
      return false;
   }
   
   char line[256];
   for ( int i=0; i<256; ++i ) line[i] = ' ';
   std::vector<std::string> tokens;
   std::string del = " ";
   int counter = -1;
   
   while ( !infile.eof()  )
   {

      for ( int i=0; i<256; ++i ) line[i] = ' '; // cleanup the read-in buffer before next iteration

      infile.getline( line, 256 );

      std::string line1 = std::string( line ); // , 256 );

      if ( line1.find("#") == 0 ) // comment line
      {	 
	 continue;
      }
      if ( line1.find_first_not_of(del) == std::string::npos ) 
      {
	 continue;
      }
      
      tokens.clear();
      SplitString( line1, del, tokens );
      if ( tokens.size() != 7 ) // x, xmin, xmax, xsec, +estat, -estat, esys
      {
	 std::cout << " EMERGENCY !!! Wrong exp.data/format." << std::endl;
	 return false;
      }
      double x   = atof(tokens[0].c_str());
      double x1  = atof(tokens[1].c_str());
      double x2  = atof(tokens[2].c_str());
      double y   = atof(tokens[3].c_str());
      double e1  = atof(tokens[4].c_str());
      double e2  = atof(tokens[5].c_str());
      double sys = atof(tokens[6].c_str());
      
      Mom[NDataPoints]   = (x1+x2)/2.;
      XSec[NDataPoints]  = y;
      
      double err2 = 0.;
      if ( useStatEr ) err2 += ((e1-e2)*(e1-e2))/4.;
      if ( useSysEr )  err2 += sys*sys;
      Err[NDataPoints] = sqrt(err2);
      
      NDataPoints++;
      if ( NDataPoints >= NMaxPoints ) break;

   }
   
   infile.close();
      
   return true;

}

TGraphErrors* getPbarProdAsGraph()
{

   double ymin = 10000.;
   double ymax = -1.;
   for ( int ip=0; ip<NDataPoints; ++ip )
   {
      if ( (XSec[ip]+Err[ip]) > ymax ) ymax = XSec[ip]+Err[ip];
      if ( (XSec[ip]-Err[ip]) < ymin ) ymin = XSec[ip]-Err[ip];
      if ( ymin < 0. ) ymin = 0.;
   }
  
   TGraphErrors* gr = new TGraphErrors(NDataPoints,Mom,XSec,0,Err);   
   
   gr->SetMinimum(ymin);
   gr->SetMaximum(ymax);
   gr->SetMarkerColor(4);  
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.6);
   gr->GetXaxis()->SetTitle( "momentum of pbar [GeV/c]" );
   gr->GetYaxis()->SetTitle( "E d^{3}#sigma / dp^{3} [GeV mb/(GeV/c)^{3}]");

   return gr;

} 

double Chi2PbarProd( std::string beam, std::string target, std::string energy, std::string angle,
                     std::string model,
		     int& NDF,
		      std::string version = ".", 
		     bool useStatEr=true, bool useSysEr=true )
{

   // NOTE: initial version; will need to update and include features for regression testing

   bool status = readPbarProdData( beam, target, energy, angle ); // useStatEr & useSysEr are true by default
   if ( !status )
   {
      std::cout << " Insufficient exp.data for " << beam << "+" << target << " --> pbar+X, at " << angle << "deg" << std::endl;
      return 0.;
   } 
   
   std::string location = "";
   if ( version == CurrentVersion || version == "." )
   {
         location = ".";
   }
   else
   {
         location = regre_test_dir + "/" + TEST_NAME + "/" + version;
   }
   
   std::string histofile = location + "/" + beam + target + energy + "GeV" + model + ".root";
   std::string histoname = "XSecVsMomPbarAt" + angle + "deg";
   
   TFile* hfile = new TFile( histofile.c_str() );
      
   TH1F* h = (TH1F*)hfile->Get( histoname.c_str() );
      
   if ( !h || NDataPoints <= 0  ) 
   {
      std::cout << " Invalid case: no exp.data or simulation, or both " << std::endl;
      return 0.;
   }
   
   // double Chi2 = 0.;   
   TGraphErrors* gr1 = getPbarProdAsGraph();
   return Chi2( gr1, h, NDF );

}

void plotPbarProd( std::string beam, std::string target, std::string energy, std::string angle )
{

   TCanvas* myc = new TCanvas( "myc", "", 900, 600 );
//   myc->Divide( 2, 1 );
   
   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.13, 0.49, 0.99 );   
//   myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   drawMCvsData( beam, target, energy, angle );

   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.13, 0.99, 0.99 );   
//   myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   drawMC2Data( beam, target, energy, angle );

   TLegend* leg = new TLegend(0.7,0.01,0.99,0.12);
   TLegendEntry* entry = 0;
   for ( int m=0; m<NModels; ++m )
   {
      entry = leg->AddEntry( "", Models[m].c_str(),"L" );
      entry->SetLineColor( ColorModel[m] );
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
   
   int ic = 0;
   std::string txt = "Geant4 vs Data:";
   TLatex* ltxt = new TLatex( 0.1, 0.1, txt.c_str() );
   ltxt->SetTextSize(0.035);
   ltxt->Draw();
   for ( int m=0; m<NModels; ++m )
   {
      int NDF=0;
      double chi2=0.; 
      chi2 += Chi2PbarProd( beam, target, energy, angle, Models[m], NDF );  
      std::ostringstream os;
      os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Models[m] );
      TLatex* ltxt1 = new TLatex(0.10, 0.065-ic*0.035, txt1.c_str() );
      ltxt1->SetTextSize(0.035);
      ltxt1->Draw();
      ic++;
   }
   
   std::string output = beam + "-" + target + "-" + energy + "GeV-pbar-" + angle +"deg";
   output += "-models.gif";
   myc->Print( output.c_str() );
      
   return;

}


void plotPbarProdRegre( std::string beam, std::string target, std::string energy, std::string angle, 
                        std::string model )
{

   TCanvas* myc = new TCanvas( "myc", "", 900, 600 );
//   myc->Divide( 2, 1 );
   
   TPad* pad1 = new TPad( "pad1", "", 0.01, 0.15, 0.49, 0.99 );   
//   myc->cd(1);
   myc->cd();
   pad1->Draw();
   pad1->cd();
   gPad->SetLogy();
   drawMCvsDataRegre( beam, target, energy, angle, model );
   
   TPad* pad2 = new TPad( "pad2", "", 0.51, 0.15, 0.99, 0.99 );   
//   myc->cd(2);
   myc->cd();
   pad2->Draw();
   pad2->cd();
   gPad->SetLogy();
   drawMC2DataRegre( beam, target, energy, angle, model );

   TLegend* leg = new TLegend(0.7,0.01,0.99,0.14);
   TLegendEntry* entry = 0;
   for ( int m=0; m<NVersions; ++m )
   {
      entry = leg->AddEntry( "", Versions[m].c_str(),"L" );
      entry->SetLineColor( ColorVersion[m] );
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
   
   int ic = 0;
   std::string txt = "Geant4 (" + model + ") vs Data";
   TLatex* ltxt = new TLatex( 0.1, 0.12, txt.c_str() );
   ltxt->SetTextSize(0.028);
   ltxt->Draw();
   for ( int m=0; m<NVersions; ++m )
   {
      int NDF=0;
      double chi2=0.; 
      chi2 += Chi2PbarProd( beam, target, energy, angle, model, NDF, Versions[m] );  
      std::ostringstream os;
      os << (chi2/NDF);
//      std::string txt1 = "Integral over all agnular bins: #chi^{2}/NDF = ";
      std::string txt1 = "#chi^{2}/NDF = ";
      txt1 += os.str();
      txt1 += ( " for " + Versions[m] );
      TLatex* ltxt1 = new TLatex(0.10, 0.092-ic*0.028, txt1.c_str() );
      ltxt1->SetTextSize(0.028);
      ltxt1->Draw();
      ic++;
   }   
   std::string output = beam + "-" + target + "-" + energy + "GeV-pbar-" + angle +"deg-";
   output += model;
   output += "-regre.gif";
   myc->Print( output.c_str() );
   return;

}

void setAxisStyle ( TAxis* axis, const char* title )
{

   axis->SetTitle( title );
   axis->SetTitleOffset(0.8);
   axis->SetTitleSize(0.05);
   axis->SetTitleFont(62);
   axis->SetLabelFont(62);
   axis->CenterTitle();	 

   return;

}

void drawMCvsData( std::string beam, std::string target, std::string energy, 
                   // std::string secondary_type, // --> anti_proton
		   std::string sec_angle )
{

   bool status = readPbarProdData( beam, target, energy, sec_angle ); // useStatEr & useSysEr are true by default
   if ( !status )
   {
      TText* txt = new TText(0.25,0.6,"Insufficient exp.data");
      txt->SetTextSize(0.075);
      txt->Draw();
      return;
   }
   
   TGraphErrors* gr1 = getPbarProdAsGraph();
   
   double ymin = gr1->GetMinimum();
   double ymax = gr1->GetMaximum();
   
   TH1F* h[NModels];
   
   for ( int m=0; m<NModels; ++m )
   {
      std::string histofile = beam + target + energy + "GeV" + Models[m] + ".root";
      std::string histoname = "XSecVsMomPbarAt" + sec_angle + "deg";
   
      TFile* hfile = new TFile( histofile.c_str() );
      
      h[m] = (TH1F*)hfile->Get( histoname.c_str() );
      
      if ( h[m] == 0 || NDataPoints <= 0  ) 
      {
         std::cout << " Invalid case: no exp.data or simulation, or both " << std::endl;
	 return;
      }
      h[m]->SetStats( 0 );
      h[m]->SetLineColor( ColorModel[m] );
      h[m]->SetLineWidth( 2 );
      
      int nx = h[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	float yy = h[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }

      if ( m == 0 ) 
      {
         setAxisStyle( h[m]->GetXaxis(), gr1->GetXaxis()->GetTitle() );
         setAxisStyle( h[m]->GetYaxis(), gr1->GetYaxis()->GetTitle() );
	 h[m]->Draw( "histE1" );
      }
      else
      {
         h[m]->Draw( "histE1same" );
      }
   
   }
   
   for ( int m=0; m<NModels; ++m )
   {
      h[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // h[m]->SetTitle("");
// -->      leg->AddEntry( h[m], Models[m].c_str(), "L" );
   }
   
   gr1->Draw( "p" );
      
   return;

}

void drawMCvsDataRegre( std::string beam, std::string target, std::string energy, 
                        // std::string secondary_type, --> anti_proton
			std::string sec_angle,
			std::string model )
{

   bool status = readPbarProdData( beam, target, energy, sec_angle ); // useStatEr & useSysEr are true by default
   if ( !status )
   {
      TText* txt = new TText(0.25,0.6,"Insufficient exp.data");
      txt->SetTextSize(0.075);
      txt->Draw();
      return;
   }

   TGraphErrors* gr1 = getPbarProdAsGraph();
   
   double ymin = gr1->GetMinimum();
   double ymax = gr1->GetMaximum();
   
   TH1F* h[NVersions];
   
   for ( int m=0; m<NVersions; ++m )
   {

      std::string location = "";
      if ( Versions[m] == CurrentVersion )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test47/" + Versions[m] + "/";
      }
      std::string histofile = location + beam + target + energy + "GeV" + model + ".root";
      std::string histoname = "XSecVsMomPbarAt" + sec_angle + "deg";
   
      TFile* hfile = new TFile( histofile.c_str() );
      
      h[m] = (TH1F*)hfile->Get( histoname.c_str() );
      if ( h[m] == 0 || NDataPoints <= 0  ) 
      {
         std::cout << " Invalid case: no exp.data or simulation, or both " << std::endl;
	 return;
      }
      h[m]->SetStats( 0 );
      h[m]->SetLineColor( ColorVersion[m] );
      h[m]->SetLineWidth( 6-m );

      int nx = h[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) {
	float yy = h[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }

      if ( m == 0 ) 
      {
         setAxisStyle( h[m]->GetXaxis(), gr1->GetXaxis()->GetTitle() );
         setAxisStyle( h[m]->GetYaxis(), gr1->GetYaxis()->GetTitle() );
	 h[m]->Draw( "histE1" );
      }
      else
      {
         h[m]->Draw( "histE1same" );
      }
   
   } 

   for ( int m=0; m<NVersions; ++m )
   {
      h[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.5); // h[m]->SetTitle("");
// -->      leg->AddEntry( h[m], Models[m].c_str(), "L" );
   }
   
   gr1->Draw( "p" );

   return;

}

void drawMC2Data( std::string beam, std::string target, std::string energy, 
                   // std::string secondary_type, // --> anti_proton
		   std::string sec_angle )
{

   bool status = readPbarProdData( beam, target, energy, sec_angle ); // useStatEr & useSysEr are true by default
   if ( !status )
   {
      TText* txt = new TText(0.25,0.6,"Insufficient exp.data");
      txt->SetTextSize(0.075);
      txt->Draw();
      return;
   }
   
   TGraphErrors* gr1 = getGraphData2Data();

   TGraph* gr[NModels];
   
   for ( int m=0; m<NModels; ++m )
   {

      std::string histofile = beam + target + energy + "GeV" + Models[m] + ".root";

      std::string histoname = "XSecVsMomPbarAt" + sec_angle + "deg";
   
      TFile* hfile = new TFile( histofile.c_str() );
      
      TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
      if ( hi == 0 || NDataPoints <= 0  ) 
      {
         std::cout << " Invalid case: no exp.data or simulation, or both " << std::endl;
	 return;
      }
      
      gr[m] = getGraphMC2Data( hi );
      gr[m]->SetMarkerColor( ColorModel[m] );
      if ( m == 0 )
      {
         gr1->SetTitle( hi->GetTitle() );
      }
   
   }   // end loop over models
   
       
   gr1->Draw( "apl" );

   for ( int m=0; m<NModels; ++m )
   {
      gr[m]->Draw( "lpsame" );
   }

   return;

}

void drawMC2DataRegre( std::string beam, std::string target, std::string energy, 
                        // std::string secondary_type, --> anti_proton
			std::string sec_angle,
			std::string model )
{

   bool status = readPbarProdData( beam, target, energy, sec_angle ); // useStatEr & useSysEr are true by default
   if ( !status )
   {
      TText* txt = new TText(0.25,0.6,"Insufficient exp.data");
      txt->SetTextSize(0.075);
      txt->Draw();
      return;
   }

   TGraph* gr[NVersions];
   
   std::string grtitle = "";
   
   for ( int m=0; m<NVersions; ++m )
   {

      std::string location = "";
      if ( Versions[m] == CurrentVersion )
      {
         location = "";
      }
      else
      {
         location = regre_test_dir + "/test47/" + Versions[m] + "/";
      }
      std::string histofile = location + beam + target + energy + "GeV" + model + ".root";
      
      std::string histoname = "XSecVsMomPbarAt" + sec_angle + "deg";
   
      TFile* hfile = new TFile( histofile.c_str() );
      
      TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
      if ( hi == 0 || NDataPoints <= 0  ) 
      {
	 std::cout << " Invalid case: no exp.data or simulation, or both " << std::endl;
	 return;
      }
      
      gr[m] = getGraphMC2Data( hi );
      gr[m]->SetMarkerColor( ColorVersion[m] );
      if ( m == 0 ) 
      {
         grtitle = hi->GetTitle();
      }

   }
   
   TGraphErrors* gr1 = getGraphData2Data();
   gr1->SetTitle( grtitle.c_str() );
       
   gr1->Draw( "apl" );

   for ( int m=0; m<NVersions; ++m )
   {
      gr[m]->Draw( "lpsame" );
   }

   return;

}

TGraphErrors* getGraphMC2Data( TH1F* hi )
{

// this assumes that the exp.data are given via external (global) variables/arrays !!!

      double* MC2DataX = new double[NDataPoints]; 
      double* MC2DataY = new double[NDataPoints]; 
      double* DX = new double[NDataPoints];
      double* DY = new double[NDataPoints]; // [30]
      int np=0;
      
      int nx = hi->GetNbinsX();
      for ( int k=1; k <= nx; k++ )
      {        
         double xx1 = hi->GetBinLowEdge(k);
	 double xx2 = hi->GetBinWidth(k);
	 for (int kk=0; kk<NDataPoints; kk++ )
	 {
	    if ( xx1 < Mom[kk] && xx1+xx2 > Mom[kk] )
	    {
	      double yy = hi->GetBinContent(k);
	      MC2DataX[np] = Mom[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / XSec[kk];
	      // also error calc here !...
	      DY[np]=Err[kk]*MC2DataY[np]/XSec[kk];
	      np++;
	      break;
            }
	 }
      }
      
      TGraphErrors* gr = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      gr->SetMarkerStyle( 21 );
      gr->SetMarkerSize( 1.6 ); 
      
      delete [] MC2DataX;
      delete [] MC2DataY;
      delete [] DX;
      delete [] DY;
      
      return gr;

}

TGraphErrors* getGraphData2Data()
{

   double* Value = new double[NDataPoints];
   double* Error = new double[NDataPoints];
   for ( int i=0; i<NDataPoints; ++i )
   {
      Value[i] = 1.;
      Error[i] = Err[i] / XSec[i];
   }
   
   TGraphErrors* gr1 = new TGraphErrors( NDataPoints, Mom, Value, 0, Error );
   setAxisStyle( gr1->GetXaxis(), "Momentum of pbar [GeV/c]"  );
   setAxisStyle( gr1->GetYaxis(), "MC/Data (d^{3}#sigma / dp^{3} [mb/(GeV/c)^{3}])" );
//   gr1->GetYaxis()->SetRangeUser( 0.1, 10. );
   gr1->GetYaxis()->SetRangeUser( 0.1, 20. );
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
      
   delete [] Value;
   delete [] Error;

   return gr1;

}
