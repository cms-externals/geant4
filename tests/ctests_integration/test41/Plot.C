//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------------
//
//      CERN Geneva Switzerland
//
//      File name:     reader_test41
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 4 July 2007
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "TFile.h"
#include "TH1F.h"

#include "TROOT.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

void Plot()
{
  int verbose = 1;
  int idx = 0;
  const int nidx = 11;
  const string tt[nidx] = 
    {"Al  ","Be1 ","Be2 ","C   ","CH_{2} ","Fe  "," H_{2}1"," H_{2}2"," Li1 "," Li2 ","        Total"};
  const string tt1[nidx] = 
    {"al","be_1","be_2","c","ch2","fe","h2_1","h2_2","li_1","li_2","Tot"};
  const string tit[nidx] = 
    {"Al 1.5 mm","Be 0.98 mm","Be 3.73 mm","C 2.5 mm","CH_{2} 4.74 mm",
     "Fe 0.24 mm",
     "Liquid H_{2} 109 mm","Liquid H_{2} 159 mm","Li 12.78 mm","Li 6.43 mm",
     "Total"};

  string legen[5] = 
    {"Opt0 ","Opt3 ","WVI ","SS ","Opt4 "}; 

  double xi2[5] = {0.,0.,0.,0.,0.};

  // data
  const int nmax = 12;
  double x[nmax];
  double xx[nmax];
  double daty[nmax];
  double date[nmax];
  double er0[nmax];
  double err0[nmax];
  double er00[nmax] = {0};
  double y[nmax];
  double t0[nmax];
  double ert0[nmax];

  double y0[nmax];
  double y1[nmax];
  double y2[nmax];
  double y3[nmax];
  double y4[nmax];

  int i, j, k;
  string mat, s1, s2, q1, q2, q3, q4, q5, q6, q7;
  double width, z1, z2, z3;

  TGraphErrors* gr[16] = {0};
  int col[6] = {1,  2,  3,  4,  6, 13};
  int mar[6] = {21, 20, 22, 23,28,  3};
  TH1F* hh[2];

  gROOT->SetStyle("Plain");
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleOffset(0.9, "x");
  gStyle->SetTitleOffset(0.5, "y");
  gStyle->SetPadBorderMode(0);

  // -------------------------------------------------------------------
  // Control on input

  char* rl = getenv("REF");
  string rel = string(rl);

  // read data
  char* path = getenv("G4INSTALL");
  if (!path) {
    cout << "!!!ERROR: G4INSTALL is not defined"<< endl;
    exit(1); 
  }
  const int gmax = 5;


  // analysis of each plot
  for(int idx=0; idx<10; ++idx) {

    string name = tt1[idx];

    string fname = string(path) + "/tests/ctests_integration/test41/data/" + name + ".dat";
    string fin[5];
    fin[0] = name + "_opt0.log";
    fin[1] = name + "_optUB.log";
    fin[2] = name + "_optG.log";
    fin[3] = name + "_optS.log";
    fin[4] = name + "_opt4.log";

    ifstream* input = new ifstream();
    input->open(fname.c_str());
    if( !input->is_open()) {
      cout << "Input file <" << fname << "> does not exist! Exit" << endl;
      exit(1);
    }
    (*input) >> mat >> width;
    for(i=0; i<10; i++) { 
      (*input) >> z1 >> daty[i] >> date[i];
      if(1 < verbose) {
	cout << "ThetaUp= " << z1   
	     << "  p= " << daty[i] << " +-  " << date[i] << endl; 
      }
      if(i == 0) {
	x[i] = 0.5*z1;
	err0[i] = x[i];
      } else {
	x[i] = 0.5*(z1 + x[i-1] + err0[i-1]);
	err0[i] = z1 - x[i];
      }
    }
    input->close();

    for(j=0; j<gmax; j++) {
      input->open(fin[j].c_str());
      if( !input->is_open()) {
	cout << "Input file <" << fin[j] << "> does not exist! Exit" << endl;
	continue;
      }

      cout << "InputFile " << idx << " " << fin[j] << " 1 X 1" << endl;
      (*input) >> s1 >> s2;
      double xi = 0.0;
      for(i=0; i<10; i++) { 
	(*input) >> z1 >> z2 >> z3;
	if(1 < verbose) {
	  cout << "ThetaUp= " << z1
	       << "  p= " << z2 << " +-  " << z3 << endl; 
	}
        xx[i] = x[i] + 0.001*j;
	y[i]  = z2; 
	er0[i]= z3; 
	z1 = z2/daty[i];  
	t0[i] = 100.*(1.0 - z1); 
	ert0[i] = 100.*z3/daty[i]; 
        xi += (z2 - daty[i])*(z2 - daty[i])/(date[i]*date[i] + z3*z3);
      }
      xi2[j] += xi;
      if(0 == j) y0[idx] = xi*0.1; 
      if(1 == j) y1[idx] = xi*0.1; 
      if(2 == j) y2[idx] = xi*0.1; 
      if(3 == j) y3[idx] = xi*0.1; 
      if(4 == j) y4[idx] = xi*0.1; 

      delete gr[j];
      delete gr[j+5];
      gr[j] = new TGraphErrors(10,x,y,er00,er0);
      gr[j+5] = new TGraphErrors(10,xx,t0,er00,ert0);
      input->close();
    }

    gStyle->SetTitleSize(0.07, "x");
    gStyle->SetTitleSize(0.07, "y");
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetPadRightMargin(0.05);
    TCanvas* c1 = new TCanvas("c1"," ",1, 5, 800, 600);
  
    delete gr[10];
    gr[10] = new TGraphErrors(10,x,daty,err0,date);

    c1->Divide(1,2);

    c1->cd(1);

    gPad->SetLogy();
    gPad->SetGrid();
    string title = "172 MeV/c muon scattering off " + tit[idx] + ", Geant4 " + rel;  
    hh[0] = gPad->DrawFrame(0.0,0.001,0.12,100.,title.c_str());
    hh[0]->GetXaxis()->SetTitle("#theta (rad)");
    hh[0]->GetYaxis()->SetTitle("probability (%/rad)    ");
    hh[0]->Draw("AXIS SAME");

    TLegend* leg =  new TLegend(0.6, 0.5, 0.87, 0.85);
    leg->SetTextSize(0.07);

    gr[10]->SetMarkerColor(col[0]);
    gr[10]->SetMarkerStyle(mar[0]);
    gr[10]->Draw("P SAME");
    leg->AddEntry(gr[10],"Data", "p"); 

    for(i=0; i<gmax; i++) {
      gr[i]->SetMarkerColor(col[i+1]);
      gr[i]->SetMarkerStyle(mar[i+1]);
      gr[i]->Draw("P SAME");
      leg->AddEntry(gr[i],legen[i].c_str(), "p"); 
      cout << "### " << i << "  " << legen[i] << " " << i 
	   << "   Xi2/Nd= " << xi2[i]/10 << endl;
    }
    leg->Draw();
    c1->Update();

    c1->cd(2);

    gPad->SetGrid();
    hh[1] = gPad->DrawFrame(0.0,-60.,0.12,60.,"");
    hh[1]->GetXaxis()->SetTitle("#theta (rad)  ");
    hh[1]->GetYaxis()->SetTitle("(1 - Geant4/Data) (%)  ");
    hh[1]->Draw("AXIS SAME");

    TBox* b1[10];
    for(i=0; i<10; i++) {
      z1 = x[i] - err0[i];
      z2 = x[i] + err0[i];
      z3 = 100*date[i]/daty[i];
      if(z3 > 60.) z3 = 60.;
      if(i == 9) { z2 = 0.12; }
      b1[i] = new TBox(z1,-z3,z2,z3);
      b1[i]->SetFillStyle(3003);
      b1[i]->SetFillColor(15);
      b1[i]->Draw("SAME 9");
    }
    c1->Update();

    for(i=0; i<gmax; i++) {
      gr[i+5]->SetMarkerColor(col[i+1]);
      gr[i+5]->SetMarkerStyle(mar[i+1]);
      gr[i+5]->Draw("P SAME 9");
    }

    string fout = "afig_" + name + ".png";

    c1->Update();
    c1->Print(fout.c_str());
    delete c1;
  }

  // analysis of Xi2
  double x0[nmax];
  double x1[nmax];
  double x2[nmax];
  double x3[nmax];
  double x4[nmax];

  for(i=0; i<=10; ++i) { 
    x0[i] = i + 1.4; 
    x1[i] = i + 1.45; 
    x2[i] = i + 1.5; 
    x3[i] = i + 1.55; 
    x4[i] = i + 1.6; 
  }
  x0[10] += 1.0;
  x1[10] += 1.0;
  x2[10] += 1.0;
  x3[10] += 1.0;
  x4[10] += 1.0;

  y0[10] = xi2[0]*0.01;
  y1[10] = xi2[1]*0.01;
  y2[10] = xi2[2]*0.01;
  y3[10] = xi2[3]*0.01;
  y4[10] = xi2[4]*0.01;

  gStyle->SetPadLeftMargin(0.10);
  //gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
  TCanvas c2("c2"," ",1, 5, 800, 600);

  gr[11] = new TGraphErrors(11,x0,y0,er00,er00);
  gr[12] = new TGraphErrors(11,x1,y1,er00,er00);
  gr[13] = new TGraphErrors(11,x2,y2,er00,er00);
  gr[14] = new TGraphErrors(11,x3,y3,er00,er00);
  gr[15] = new TGraphErrors(11,x4,y4,er00,er00);

  gPad->SetLogy();
  gPad->SetGrid();
  string title = "172 MeV/c muon scattering - MuScat, Geant4 " + rel;  
  hh[0] = gPad->DrawFrame(0.0,0.2,14,80,title.c_str());
  hh[0]->GetXaxis()->SetNdivisions(0);
  //    hh[0]->GetXaxis()->SetTitle("");
  hh[0]->GetYaxis()->SetTitle("#Chi^{2}/N");
  hh[0]->Draw("AXIS SAME");

  TLegend* leg =  new TLegend(0.55, 0.7, 0.88, 0.88);
  leg->SetTextSize(0.03);
  TString ggg = "         ";
  for(i=0; i<11; i++) { ggg += tt[i] + "    "; }
  TLegend* leg4 =  new TLegend(0.1, 0.08, 0.95, 0.16, ggg);
  leg4->SetTextSize(0.03);

  for(i=0; i<5; i++) {
    gr[11+i]->SetMarkerColor(col[i+1]);
    gr[11+i]->SetMarkerStyle(mar[i+1]);
    gr[11+i]->SetMarkerSize(1.5);
    gr[11+i]->Draw("P SAME 9");
    int iii = int(xi2[i]);
    double qqq = iii*0.01;
    std::ostringstream os;
    os << qqq;
    string lgg = legen[i] + " #chi^{2}/N= " + os.str();
    leg->AddEntry(gr[11+i],lgg.c_str(), "p");
  } 
  
  leg->Draw("SAME");
  leg4->Draw("SAME");
  c2.Update();
  string fout = "afig_Xi2.png";
  c2.Print(fout.c_str());

  //fout = "afig_Xi2.eps";
  //c2.Print(fout.c_str());
}

