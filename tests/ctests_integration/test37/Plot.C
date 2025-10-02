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
//      File name:     reader_test37
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 5 July 2007
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "TFile.h"
#include "TH1F.h"

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

void Plot()
{
  int verbose = 1;

  // -------------------------------------------------------------------
  // Control on input


  char* rl = getenv("REF");
  string rel = string(rl);

  char* g4 = getenv("G4INSTALL");
  string path = string(g4);

  string fin[5];

  // setting
  int idx = 0;
  const int nidx = 12;
  int nmed[nidx] = {1,1,1,3,2,1,1,1,1,1,1,1}; 
  string te[nidx] = {"0.521","0.5" ,"1.0" ,"1.0" ,"0.521","0.521","1.0",
                     "0.015","0.02","0.03","0.04","0.05"};
  string tekev[nidx] = {"521","500" ,"1000" ,"1000" ,"521","521","1000",
                        "15","20","30","40","50"};
  string tt[nidx] = {"Al","Mo","Ta","AlAuAl","TaAl","Be","U",
                     "Si","Si","Si","Si","Si"};
  double hei[nidx]= {5.0,6.0,6.0,6.0,6.0,4.0,5.5,45.0,35.0,26.5,21.0,18.0};
  string legen[5] = {"Opt0","Opt4","WVI-SS","Opt3","Single Scat"}; 

  TH1D* h[5];
  string hhh[5] = {"h0","h1","h2","h3","h4"};
  double zz[5] = {0.,0.,0.,0.,0.};

  // data
  const int nmax = 120;
  double erx[nmax];
  double datx[nmax];
  double daty[nmax];
  double date[nmax];
  double x[nmax];
  double y0[nmax];
  double er0[nmax];

  int col[6] = {1, 2, 3, 4, 13, 6};
  int mar[6] = {21, 20, 22, 23, 28, 3};
  TH1F* hh;

  gROOT->SetStyle("Plain");
  gStyle->SetLabelSize(0.04, "x");
  gStyle->SetLabelSize(0.04, "y");
  gStyle->SetTitleOffset(0.9, "x");
  gStyle->SetTitleOffset(0.9, "y");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBorderMode(0);

  for(int idx=0; idx<nidx; idx++) {

    string name = tt[idx];
    string fname = string(path) + "/tests/ctests_integration/test37/Sandia_Data/" + name + ".dat";
    // open and read data file
    ifstream* input = new ifstream();
    input->open(fname.c_str());
    if( !input->is_open()) {
      cout << "Input file <" << fname << "> does not exist! Exit" << endl;
      exit(1);
    }
    cout << idx << ". file is opened: <" << fname << endl;
    if (name=="Si") {
      fin[0] = name + "_opt0_"+tekev[idx]+"keV.log";
      fin[1] = name + "_opt4_"+tekev[idx]+"keV.log";
      fin[2] = name + "_optG_"+tekev[idx]+"keV.log";
      fin[3] = name + "_opt3_"+tekev[idx]+"keV.log";
      fin[4] = name + "_optS_"+tekev[idx]+"keV.log";
    } else {
      fin[0] = name + "_opt0.log";
      fin[1] = name + "_opt4.log";
      fin[2] = name + "_optG.log";
      fin[3] = name + "_opt3.log";
      fin[4] = name + "_optS.log";    
    }

    int i, j, nm;
    string s1, s2, s3, s4, s5, s6;
    double z1 = 0.0;
    double z2;
    int nd = 0;
    do {
      (*input) >> s1 >> s2 >> s3 >> s4 >> s5 >> s6;

      // select data
      if(s2 == te[idx] && s5 == "0") {
        (*input) >> s3 >> s4;
        for(i=0; i<nmax; i++) { 
          (*input) >> datx[i] >> daty[i];
          date[i] = daty[i]*0.02;
          erx[i] = 0.0;
          if(datx[i] > 5.0) {
            z1 = 100.;
            break;
          }
          nd++;
          if(1 < verbose) 
            cout << nd << ". x= " << datx[i]   
                 << ";  y= " << daty[i] << endl; 
        }
      } else {
        (*input) >> s3 >> s4;
        for(i=0; i<nmax; i++) { 
          (*input) >> z1 >> z2;
          if(z1 > 5.0) break;
        }
      }
    } while (z1 < 15.0);

    input->close();

    TGraphErrors* gr = new TGraphErrors(nd,datx,daty,erx,date);

    TCanvas* c1 = new TCanvas("c1"," ",1, 5, 800, 600);
    gPad->SetGrid();

    double width  = 1.0;
    double height = hei[idx];

    string title = "e^{-} " + te[idx] + " MeV in " + tt[idx] + ", Geant4 " + rel;  
    hh = gPad->DrawFrame(0.0,0.0,width,height,title.c_str());
    hh->GetXaxis()->SetTitle("depth (R/R_{0})   ");
    hh->GetYaxis()->SetTitle("E_{dep} (MeV/g/cm^{2})  ");
    hh->Draw("AXIS SAME");

    TLegend* leg;
    if(idx < 5) leg = new TLegend(0.6, 0.5, 0.90, 0.85);
    else leg = new TLegend(0.75, 0.75, 0.95, 0.95);
    leg->SetTextSize(0.04);

    gr->SetMarkerColor(col[0]);
    gr->SetMarkerStyle(mar[0]);
    gr->Draw("P SAME");
    leg->AddEntry(gr,"Data", "p"); 
    
    c1->Update();

    // open and read simulation
    int nbin = 0;

    x[0] = 0.0;
    for(j=0; j<5; j++) {
      input->open(fin[j].c_str());
      if( !input->is_open()) {
	cout << "Input file <" << fin[j] << "> does not exist! Exit" << endl;
	continue;
      }
      if(1 < verbose) 
	cout << "Input file <" << fin[j] << ">" << "  Nmedia= " << nmed[idx] << endl;

      nbin = 0;
      (*input) >> s1 >> s2 >> s3;
      for(int k=0; k<nmed[idx]; k++) {
	(*input) >> s1 >> s2 >> s3 >> s4 >> nm;

	for(i=0; i<nm; i++) { 
	  (*input) >> z1 >> z2;
	  if(1 < verbose) 
	    cout << "ThetaUp= " << z1
		 << ";  p= " << z2 << endl; 
	  y0[nbin] = z2; 
	  er0[nbin] = 0.0; 
	  nbin++;
	  x[nbin] = z1;
	  if(i == nm-1) zz[k+1] = z1;
	}
      }
      input->close();
      h[j] = new TH1D( hhh[j].c_str(),"",nbin,x);
      for(i=0; i<nbin; i++) {
	h[j]->SetBinContent(i+1,y0[i]);
	h[j]->SetBinError(i+1,er0[i]);
      }
      h[j]->SetLineColor(col[j+1]);
      h[j]->SetLineWidth(2);
      h[j]->Draw("HIST C SAME");
      leg->AddEntry(h[j],legen[j].c_str(), "l"); 
    }
    TLine* line[4];
    for(j=1; j<nmed[idx]; j++) {
      z2 = zz[j];
      line[j] = new TLine(z2, 0., z2, height); 
      line[j]->SetLineColor(col[5]);
      line[j]->Draw("SAME");
    }
  
    leg->Draw();
    c1->Update();

    string fout;
    if (name=="Si") {
      fout = "Afig" + name + "_" + tekev[idx] + "keV.png";
    } else {
      fout = "Afig" + name + ".png";
    }

    c1->Print(fout.c_str());

    //fout = "Afig" + name + ".eps";
    //c1->Print(fout.c_str());
    //  fout = "Afig" + name + ".pdf";
    // c1->Print(fout.c_str());
    delete c1;
  }
}

