//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include <fstream>
#include <iostream>
#include <string>
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h" 
#include "TGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"

using namespace std;

void Plot()
{
  const int nidx = 3;
  string fnm[nidx] = {"p", "he4", "c12"};
  string tp[nidx] = {"p", "^{4}He", "^{12}C"};
  string te[nidx] = {"110", "144.3", "100"};
  string teu[nidx] = {"MeV", "MeV/u", "MeV/u"};
  string fname2[nidx] = {"H-110MeV-endep-EXP-norm-max.txt",
			 "4He-144.3MeV-endep-EXP-M03-norm-max.txt",
			 "12C100MeVen-dep-EXP-norm-max.txt"}; 

  double zmax[nidx] = {120., 180., 40.};
 
  string fname = getenv("PARTICLE");
  int idx = 0;
  for (; idx < nidx; idx++) {if (fname == fnm[idx]) break;}

  string refer = getenv("REF");

  TString finName[3];
  finName[0] = fname + "_opt0.root";
  finName[1] = fname + "_opt3.root";
  finName[2] = fname + "_opt4.root";

  TString legend[3] = {"QBBC opt0", "QBBC opt3", "QBBC opt4"};

  int n_exp[nidx] = {39, 25, 76};
  int nn = n_exp[idx];
  double *x_exp = new double[nn];
  double *y_exp = new double[nn];

  char buffer[256];

  gROOT->SetStyle("Plain");
  TCanvas *c1 = new TCanvas("c1", "c1",6,6,800,600);
  gStyle->SetOptStat(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameBorderMode(0);
  
  ifstream in;

  in.open(fname2[idx].c_str());

  if( !in.is_open()) { 
    cout << "Input file<" << fname2[idx] << "> does not exist! Exit" << endl;
    return;
  } else {
    cout << "### test44 analysis for " << fnm[idx] << "  file: "
         << finName[idx] << endl;
  }
   
  // Ignore first blank line
  in.getline(buffer,256);

  for (int i=0; i<nn; i++) {
    in >> x_exp[i] >> y_exp[i];
    x_exp[i] = 10.*x_exp[i];
  }

  string hist_title = tp[idx] + " " + te[idx] + " " + 
    teu[idx] + " " + "in Water, Geant4  " + refer;

  cout << "Data file <" << fname2[idx] << " was red " << nn << " lines" << endl;
  
  TH1F* h0 = gPad->DrawFrame(0.0,0.0,zmax[idx],1.2,hist_title.c_str());
  h0->GetXaxis()->SetTitle("z (mm)");
  h0->GetYaxis()->SetTitle("dose (relative unit)");
  h0->Draw("AXIS SAME 9");
  
  TGraph *gr = new TGraph(nn,x_exp,y_exp);
  gr->SetMarkerStyle(22);
  gr->SetMarkerSize(1.2);
  gr->Draw("P SAME 9");

  TLegend *leg = new TLegend(0.2,0.65,0.45,0.86);
  leg->SetTextFont(52);
  leg->SetTextSize(0.035);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillStyle(0);
  leg->SetMargin(0.4);
  leg->SetBorderSize(1);
  leg->AddEntry(gr,"Data","p");

  in.close();

  TH1F* hh[3];
  TFile* fq[3];

  cout << "Start loop" << endl;
  for (int j = 0; j < 3; j++) {
    cout << "File with MC <" << finName[j] << "> will be opened" << endl;
    fq[j] = new TFile(finName[j]);
    hh[j] = (TH1F*)fq[j]->Get("h1");
    if(!hh[j]) continue;
    double ymax = hh[j]->GetMaximum();
    cout << "Ymax= "<< ymax<<endl;
    hh[j]->Scale(1./ymax);
    hh[j]->SetLineColor(j+2);
    hh[j]->SetLineWidth(1);
    hh[j]->Draw("HISTO SAME 9");
  
    leg->AddEntry(hh[j], legend[j], "l");
  }
  leg->Draw("SAME 9");
  TString fout1 = "A_" + fnm[idx] + "_water.png";
  c1->Print(fout1);
  c1->Close();

  TCanvas *c2 = new TCanvas("c2","c2",6,6,800,600);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(0);
  c2->SetFrameBorderMode(0);

  gPad->SetLogy();
  TH1F* h1 = gPad->DrawFrame(0.0,0.0001,zmax[idx],10,hist_title.c_str());
  h1->GetXaxis()->SetTitle("z (mm)");
  h1->GetYaxis()->SetTitle("log(dose (relative unit))");
  h1->Draw("AXIS SAME 9");

  gr->Draw("P SAME 9");
  cout << "Start log loop" << endl;
  for (int j = 0; j < 3; j++) {
    if(!hh[j]) continue;
    hh[j]->Draw("HISTO SAME 9");
  }
  leg->Draw("SAME 9");
  TString fout2 = "B_" + fnm[idx] + "_water.png";
  c2->Print(fout2);
  c2->Close();

  delete [] x_exp;
  delete [] y_exp;
}
