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

#include "HARP.h"

#include "Chi2Calc.C"

#ifndef G4VAL_DRAWHARPSPECTRA_C
#define G4VAL_DRAWHARPSPECTRA_C

double ChiSqThetaRegre[NVersions];
int NDFThetaRegre[NVersions];
double ChiSqThetaModel[NModels_IE];
int NDFThetaModel[NModels_IE];

void PlotHARPHisto( std::string beam, std::string target, std::string energy,
                    std::string secondary,
		    std::string region,
		    int ibin,
		    bool displayLeg=true,
		    bool largeLegText=false,
		    float legTextSize=0.05 )
{
   
   // this happens in the "highler level" macro 
   // ReadHARPData( beam, target, energy, secondary, region );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      if ( (YHARP[ibin][i]+EYHARP[ibin][i]) > ymax ) ymax = YHARP[ibin][i]+EYHARP[ibin][i];
      if ( (YHARP[ibin][i]-EYHARP[ibin][i]) < ymin && (YHARP[ibin][i]-EYHARP[ibin][i]) > 0. ) ymin = (YHARP[ibin][i]-EYHARP[ibin][i]);
   }
   
   TH1F* hi[NModels_IE];
   std::string YTitle;

   for ( int m=0; m<NModels_IE; m++ )
   {

      std::string histofile = "";
      
      histofile = "./harp-histo/";
      // histofile = "./harp-histo-no-res-decays/";
      // histofile = "../t23-bld/harp-histo/";
      
      // std::string histofile = "./harp-histo/" + beam + target + energy + "GeV" + ModelName_IE[m] + ".root"; 
      histofile += ( beam + target + energy + "GeV" + ModelName_IE[m] + ".root" ); 
            
      TFile* f = new TFile( histofile.c_str() );
      
      char buf[5];
      sprintf( buf, "%i", ibin );      
      std::string histoname = secondary + "_" + region + "_";
      histoname.append( buf );
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
            
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel_IE[m]);
      hi[m]->SetLineWidth(2);
      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 std::string htit = energy + "GeV/c  " + std::string( hi[m]->GetTitle() );
	 hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("momentum (GeV/c)");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
	 hi[m]->GetYaxis()->SetTitleOffset(1.1);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
         hi[m]->Draw("histE1");
      }
      else hi[m]->Draw("histE1same"); 
                        
   }
   
   TLegend* leg = 0;
   if ( largeLegText )
   {
      leg = new TLegend(0.75, 0.7, 0.99, 0.99);
//      leg->SetTextSize(0.1);
   }
   else
   {
      leg = new TLegend(0.5, 0.55, 0.9, 0.9);
      leg->SetTextSize(0.055);
   }
   leg->SetTextSize( legTextSize );
   
   for ( int m=0; m<NModels_IE; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], ModelName_IE[m].c_str(), "L" );
   }
      
   float* X = new float[NPointsHARP[ibin]];
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      X[i] = 0.5 * (XMinHARP[ibin][i]+XMaxHARP[ibin][i]);
   }
   
   TGraph* gr = new TGraphErrors( NPointsHARP[ibin], X, YHARP[ibin], 0, EYHARP[ibin] );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("psame");
      
   leg->AddEntry( gr, "exp.data", "p");

   if ( displayLeg ) 
   {
      leg->Draw();
      leg->SetFillColor(kWhite); 
   }  

   return;

}

void PlotHARPMC2Data( std::string beam, std::string target, std::string energy,
                    std::string secondary,
		    std::string region,
		    int ibin,
		    bool gtitle=true,
		    bool displayLeg=true )
{

   float ymin = 10000.; // something big... don't know if I can use FLT_MAX
   float ymax = -1. ;

   TH1F* hi[NModels_IE];
   
   float* MC2DataX = new float[NPointsHARP[ibin]];
   float* MC2DataY = new float[NPointsHARP[ibin]];
   float* DX = new float[NPointsHARP[ibin]];
   float* DY = new float[NPointsHARP[ibin]];

   int np = 0;
   
   TGraph* gr[NModels_IE];

   float xlimit1 = 0.;
   float xlimit2 = 0.;

   for ( int m=0; m<NModels_IE; m++ )
   {

      std::string histofile = "";
      
      histofile = "./harp-histo/";
      // histofile = "./harp-histo-no-res-decays/";
      // histofile = "../t23-bld/harp-histo/";
      
      // std::string histofile = "./harp-histo/" + beam + target + energy + "GeV" + ModelName_IE[m] + ".root"; 
      histofile += ( beam + target + energy + "GeV" + ModelName_IE[m] + ".root" ); 
      
//      std::cout << " histofile = " << histofile << std::endl;
      
      TFile* f = new TFile( histofile.c_str() );
     
      char buf[5];
      sprintf( buf, "%i", ibin );      
      std::string histoname = secondary + "_" + region + "_";
      histoname.append( buf );

      hi[m] = (TH1F*)f->Get( histoname.c_str() );

      int nx = hi[m]->GetNbinsX();
      
      if ( m == 0 )
      {
         TAxis* xaxis = hi[m]->GetXaxis();
	 xlimit1 = xaxis->GetBinLowEdge( xaxis->GetFirst() );
	 xlimit2 = xaxis->GetBinUpEdge(  xaxis->GetLast() ); 
      }
      
      np=0;
      for ( int k=0; k<=nx; k++ ) 
      {         
         float xx1 = hi[m]->GetBinLowEdge(k);
	 float xx2 = hi[m]->GetBinWidth(k);
	 for (int kk=0; kk<NPointsHARP[ibin]; kk++ )
	 {
	    float xdata = 0.5 * (XMinHARP[ibin][kk]+XMaxHARP[ibin][kk]);
	    if ( xx1 < xdata && xx1+xx2 > xdata )
	    {
	       if ( YHARP[ibin][kk] <=0. || EYHARP[ibin][kk] <=0. ) continue;
	       float yy  = hi[m]->GetBinContent(k);
	       float eyy = hi[m]->GetBinError(k);
	       MC2DataX[np] = xdata;
	       DX[np] = 0.;	       
	       MC2DataY[np] = yy / YHARP[ibin][kk];
	       DY[np] = EYHARP[ibin][kk]*MC2DataY[np] / YHARP[ibin][kk];
	       if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	       if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	       np++;
	       break;
            }
	 }            
      }
      
      gr[m] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
      
      int ibb = hi[m]->GetXaxis()->GetFirst();
      float xx1 = hi[m]->GetBinLowEdge(ibb);
      ibb = hi[m]->GetXaxis()->GetLast();
      float xx2 = hi[m]->GetBinLowEdge(ibb) + hi[m]->GetBinWidth(ibb);
      gr[m]->GetXaxis()->SetRangeUser(xx1,xx2);
      if ( gtitle ) 
      {
         gr[m]->SetTitle( hi[m]->GetTitle() );
      }
      else 
      {
         gr[m]->SetTitle("");
      }
      gr[m]->SetMarkerColor(ColorModel_IE[m]);  
      gr[m]->SetMarkerStyle(SymbModel_IE[m]);
      gr[m]->SetMarkerSize(1.2);
      gr[m]->GetYaxis()->SetRangeUser(0.2,5.);
//	 gr[m]->GetXaxis()->SetTitle("momentum (GeV/c)");
	 //hi[m]->GetYaxis()->SetTitle( YTitle.c_str() );
	 // hi[m]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dpd#Theta} [mb/(GeV/c/rad)]");
   }

   float* Value = new float[NPointsHARP[ibin]];
   float* Error = new float[NPointsHARP[ibin]];
   float* XX    = new float[NPointsHARP[ibin]];
//   float* Value = new float[np];
//   float* Error = new float[np];
//   float* XX    = new float[np];

   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      XX[i] = 0.5 * (XMinHARP[ibin][i]+XMaxHARP[ibin][i]);
      Value[i] = 1.;
      if ( YHARP[ibin][i] <=0. || EYHARP[ibin][i] <=0. ) 
      {
         Error[i] = 0.;
	 continue;
      }
      Error[i] = EYHARP[ibin][i] / YHARP[ibin][i];
      if ( Value[i]+Error[i] > ymax ) ymax = Value[i] + Error[i];
      if ( Value[i]-Error[i] < ymin ) ymin = Value[i] - Error[i];
      if ( ymin < 0. ) ymin = 0.; 
   }
   
//    std::cout << " ymin, ymax = " << ymin << " " << ymax << std::endl;

//   TGraph* gr1 = new TGraphErrors( NPointsHARP[ibin], XX, Value, 0, Error );
   TGraph* gr1 = new TGraphErrors( np, XX, Value, 0, Error );

   TAxis* xaxis = gr1->GetXaxis();   
   xaxis->SetTitle("momentum of secondary, GeV/c");
   xaxis->SetTitleSize(0.05);
   xaxis->SetTitleFont(62);
   xaxis->SetTitleOffset(0.9);
   xaxis->CenterTitle();

//   std::cout << " ibin= " << ibin << std::endl;
//   std::cout << " xmin= " << XMinHARP[ibin][0] << " " << " xmax= " << XMaxHARP[ibin][NPointsHARP[ibin]-1] << std::endl;
   
//   gr1->GetXaxis()->SetLimits( XMinHARP[ibin][0], XMaxHARP[ibin][NPointsHARP[ibin]-1] );
   gr1->GetXaxis()->SetLimits( xlimit1, xlimit2 );
   gr1->GetYaxis()->SetRangeUser( ymin, ymax );
   // gr1->GetXaxis()->SetRangeUser( -0.3, 0.4 );
   
   gr1->SetMarkerColor(kBlack /* kBlue */);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.0);
   // gr1->SetMarkerSize(1.5);
   if ( gtitle ) 
   {
      gr1->SetTitle(gr[0]->GetTitle());
   }
   else
   {
      gr1->SetTitle("");
   }
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetTitle("MC/Data ( d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)] )");
//   gr1->GetYaxis()->SetTitle("");
   gr1->GetYaxis()->SetTitleOffset(1.0);
   gr1->GetYaxis()->SetTitleSize(0.05);
   gr1->GetXaxis()->SetTitleFont(62);
   gr1->GetYaxis()->CenterTitle();
   gr1->GetYaxis()->SetRangeUser(0.2,5.);

   gr1->Draw("apl");
   // gr1->Draw("AC*");    
   
   TLegend* leg = new TLegend(0.6, 0.70, 0.9, 0.9);   

   for ( int m=0; m<NModels_IE; m++ )
   {
      // gr[m]->GetYaxis()->SetRangeUser( ymin-0.1, ymax+0.2 );
      leg->AddEntry( gr[m], ModelName_IE[m].c_str(), "p" );
      gr[m]->Draw("lpsame");
   }

   leg->AddEntry( gr1, "exp.data", "p");
   
   if ( displayLeg )
   {
      leg->Draw();
      leg->SetFillColor(kWhite);  
   } 

   return;

}

void PlotHARPHistoRegre( std::string beam, std::string target, std::string energy,
                    std::string secondary,
		    std::string region,
		    int ibin,
		    std::string model,
		    bool displayLeg=true,
		    bool largeLegText=false )
{
   
   // this happens in the "highler level" macro 
   // ReadHARPData( beam, target, energy, secondary, region );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      if ( (YHARP[ibin][i]+EYHARP[ibin][i]) > ymax ) ymax = YHARP[ibin][i]+EYHARP[ibin][i];
      if ( (YHARP[ibin][i]-EYHARP[ibin][i]) < ymin && (YHARP[ibin][i]-EYHARP[ibin][i]) > 0. ) ymin = (YHARP[ibin][i]-EYHARP[ibin][i]);
   }
   
   TH1F* hi[NVersions];
   std::string YTitle;

   for ( int m=0; m<NVersions; m++ )
   {

/*
      std::string histofile = "";      
      histofile = Versions[m] + "/harp-histo/";
      // histofile = "./harp-histo-no-res-decays/";
      // histofile = "../t23-bld/harp-histo/";
      histofile += ( beam + target + energy + "GeV" + model + ".root" ); 
*/
      std::string location = "";
      if ( Versions[m] == CurrentVersion )
      {
         location = ".";
      }
      else
      {
//         location = regre_test_dir + "/test19/" + Versions[m];
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[m];
      }
      std::string histofile = location + "/harp-histo/" + beam + target + energy + "GeV" + model + ".root"; 
            
      TFile* f = new TFile( histofile.c_str() );
      
      char buf[5];
      sprintf( buf, "%i", ibin );      
      std::string histoname = secondary + "_" + region + "_";
      histoname.append( buf );
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorVersion[m]);
      hi[m]->SetLineWidth(2);
      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 std::string htit = energy + "GeV/c  " + std::string( hi[m]->GetTitle() );
	 hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("momentum (GeV/c)");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
	 hi[m]->GetYaxis()->SetTitleOffset(1.0);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
         hi[m]->Draw("histE1");
      }
      else hi[m]->Draw("histE1same"); 
                  
   }
   
   TLegend* leg = 0;
   if ( largeLegText )
   {
      leg = new TLegend(0.75, 0.7, 0.99, 0.99);
      if ( largeLegText ) leg->SetTextSize(0.1);
   }
   else
   {
//      leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg = new TLegend(0.45, 0.55, 0.9, 0.9);
      leg->SetTextSize(0.045);
   }
   
   for ( int m=0; m<NVersions; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      leg->AddEntry( hi[m], Versions[m].c_str(), "L" );
   }
      
   float* X = new float[NPointsHARP[ibin]];
   for ( int i=0; i<NPointsHARP[ibin]; i++ )
   {
      X[i] = 0.5 * (XMinHARP[ibin][i]+XMaxHARP[ibin][i]);
   }
   
   TGraph* gr = new TGraphErrors( NPointsHARP[ibin], X, YHARP[ibin], 0, EYHARP[ibin] );
   gr->SetMarkerColor(kBlack /* kBlue */ );
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(1.5);
    
   gr->Draw("p");
      
   leg->AddEntry( gr, "exp.data", "p");

   if ( displayLeg ) 
   {
      leg->Draw();
      leg->SetFillColor(kWhite); 
   }  


   return;

}

void PlotHARPThetaSpectrum( std::string beam, std::string target, std::string energy,
                            std::string secondary, float pmin, float pmax )
{
      
   int ibin = findHARPMomBin( pmin, pmax );
   if ( ibin < 0 ) 
   {
      std::cout << " INVALID ibin = " << ibin << std::endl;
      return;
   }
   
   TGraphErrors* gr = getHARPAsThetaGraph( pmin, pmax );
   
   for ( int m=0; m<NModels_IE; ++m )
   {
      ChiSqThetaModel[m] = 0.;
      NDFThetaModel[m] = 0;
   }
   
   float ymin = 1000000.;
   float ymax = -1.; 
   ymin = gr->GetYaxis()->GetXmin();
   ymax = gr->GetYaxis()->GetXmax();
     
   TH1F* hi[NModels_IE];
   std::string YTitle;

   for ( int m=0; m<NModels_IE; m++ )
   {
      
      std::string histofile = "";
      
      histofile = "./harp-histo/";
      // histofile = "./harp-histo-no-res-decays/";
      // histofile = "../t23-bld/harp-histo/";
      
      // std::string histofile = "./harp-histo/" + beam + target + energy + "GeV" + ModelName_IE[m] + ".root"; 
      histofile += ( beam + target + energy + "GeV" + ModelName_IE[m] + ".root" ); 
            
      TFile* f = new TFile( histofile.c_str() );
      
      char buf[5];
      int ibin1 = ibin+1;
      sprintf( buf, "%i", ibin1 );      
      std::string histoname = secondary + "_mom_";
      histoname.append( buf );
            
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
            
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorModel_IE[m]);
      hi[m]->SetLineWidth(2);
      
      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 // std::string htit = energy + "GeV/c  " + std::string( hi[m]->GetTitle() );
	 // hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("#Theta (rad)");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
	 hi[m]->GetYaxis()->SetTitleOffset(1.0);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
         hi[m]->Draw("histE1");
      }
      else hi[m]->Draw("histE1same");
      
      ChiSqThetaModel[m] = Chi2( gr, hi[m], NDFThetaModel[m] );
   }
   

/*
   TLegend* leg = 0;
//   if ( largeLegText )
//   {
      leg = new TLegend(0.65, 0.6, 0.99, 0.9);
      leg->SetTextSize(0.04);
//   }
//   else
//   {
//      leg = new TLegend(0.5, 0.55, 0.9, 0.9);
//      leg->SetTextSize(0.055);
//   }
//   leg->SetTextSize( legTextSize );
*/   
   ymin = std::max( float(0.1), ymin );
   for ( int m=0; m<NModels_IE; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      // leg->AddEntry( hi[m], ModelName_IE[m].c_str(), "L" );
   }

   gr->Draw("psame"); 

/*
   leg->AddEntry( gr, "exp.data (HARP)", "p");
      leg->Draw();
      leg->SetFillColor(kWhite); 
*/
   
   return;

}

void PlotHARPThetaSpectrumRegre( std::string beam, std::string target, std::string energy, 
                                 std::string secondary, 
				 float pmin, float pmax, 
				 std::string model )
{
   
   for ( int m=0; m<NVersions; ++m )
   {
      ChiSqThetaRegre[m] = 0.;
      NDFThetaRegre[m] = 0;
   }
   
   int ibin = findHARPMomBin( pmin, pmax );
      
   if ( ibin < 0 ) 
   {
      std::cout << " INVALID ibin = " << ibin << std::endl;
      return;
   }
   
   TGraphErrors* gr = getHARPAsThetaGraph( pmin, pmax );
   
   float ymin = 1000000.;
   float ymax = -1.; 
   ymin = gr->GetYaxis()->GetXmin();
   ymax = gr->GetYaxis()->GetXmax();
     
   TH1F* hi[NVersions];
   std::string YTitle;
   
   for ( int m=0; m<NVersions; ++m )
   {
      std::string location = "";
      if ( Versions[m] == CurrentVersion )
      {
         location = ".";
      }
      else
      {
//         location = regre_test_dir + "/test19/" + Versions[m];
         location = regre_test_dir + "/" + TEST_NAME + "/" + Versions[m];
      }
      std::string histofile = location + "/harp-histo/" + beam + target + energy + "GeV" + model + ".root"; 
            
      TFile* f = new TFile( histofile.c_str() );
      
      char buf[5];
      int ibin1 = ibin+1;
      sprintf( buf, "%i", ibin1 );      
      std::string histoname = secondary + "_mom_";
      histoname.append( buf );
      
      hi[m] = (TH1F*)f->Get( histoname.c_str() );
      
      hi[m]->SetStats(0);
      hi[m]->SetLineColor(ColorVersion[m]);
      hi[m]->SetLineWidth(2);
      // hi[m]->SetLineWidth(5-m);

      int nx = hi[m]->GetNbinsX();
      for (int k=1; k <= nx; k++) 
      {
	double yy = hi[m]->GetBinContent(k);
	if ( yy > ymax ) ymax = yy;
	if ( yy < ymin && yy > 0. ) ymin = yy;
      }
      if ( m == 0 ) 
      {
	 // std::string htit = energy + "GeV/c  " + std::string( hi[m]->GetTitle() );
	 // hi[m]->SetTitle( htit.c_str() );
	 hi[m]->GetXaxis()->SetTitle("#Theta (rad)");
	 hi[m]->GetXaxis()->SetTitleSize(0.05);
	 hi[m]->GetXaxis()->SetTitleFont(62);
	 hi[m]->GetXaxis()->SetLabelFont(62);
	 hi[m]->GetXaxis()->CenterTitle();
	 hi[m]->GetYaxis()->SetTitle("d^{2}#sigma / dpd#Theta [mb/(GeV/c/rad)]");
	 hi[m]->GetYaxis()->SetTitleOffset(1.1);
	 hi[m]->GetYaxis()->SetTitleSize(0.05);
	 hi[m]->GetYaxis()->SetTitleFont(62);
	 hi[m]->GetYaxis()->SetLabelFont(62);
	 hi[m]->GetYaxis()->CenterTitle();
         hi[m]->Draw("histE1");
      }
      else hi[m]->Draw("histE1same"); 
      
      ChiSqThetaRegre[m] = Chi2( gr, hi[m], NDFThetaRegre[m] );

   }
   
/*
   TLegend* leg = 0;
   leg = new TLegend(0.65, 0.7, 0.99, 0.9);
   leg->SetTextSize(0.03);
*/   
   ymin = std::max( float(0.1), ymin );
   for ( int m=0; m<NVersions; m++ )
   {
      hi[m]->GetYaxis()->SetRangeUser(ymin,ymax*1.1); // hi[m]->SetTitle("");
      // leg->AddEntry( hi[m], Versions[m].c_str(), "L" );
   }

   gr->Draw("psame"); 

/*
   leg->AddEntry( gr, "exp.data (HARP)", "p");
   leg->Draw();
   leg->SetFillColor(kWhite); 
*/
   return;

}

#endif
