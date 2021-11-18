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
#include <cassert>
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4EmDNAChemistry.hh"
#include "G4H2O.hh"
#include "G4MolecularConfiguration.hh"
#include "Randomize.hh"
#ifdef USE_ROOT
#include "TApplication.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#endif

G4DNAWaterDissociationDisplacer displacer;
G4EmDNAChemistry chemConst;

void Init(){
  chemConst.ConstructMolecule();
  chemConst.ConstructDissociationChannels();
}

//------------------------------------------------------------------------------
/*
class CounterContainer{
public:
 G4CT_COUNT_INIT(0)
 G4CT_COUNT(a)
 G4CT_COUNT(b)
 G4CT_COUNT(c)
};

class CounterContainer2: public CounterContainer{
public:
  G4CT_COUNT(d)
};

class CounterContainer3: public CounterContainer{
public:
  G4CT_COUNT(e)
};

class CounterContainer4: public CounterContainer2{
public:
  G4CT_COUNT(f)
};

class CounterContainer5: public CounterContainer3{
public:
  G4CT_COUNT(f)
};

void TestCounter(){
  G4cout << CounterContainer::a << G4endl;
  G4cout << CounterContainer::b << G4endl;
  G4cout << CounterContainer::c << G4endl;
  G4cout << CounterContainer2::d << G4endl;
  G4cout << CounterContainer3::e << G4endl;
  G4cout << CounterContainer4::f << G4endl;
  G4cout << CounterContainer5::f << G4endl;
}
*/

//------------------------------------------------------------------------------

void TestDifferentIndexes(){
  G4cout << G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay
  << G4endl;
  G4cout << G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay << G4endl;
  G4cout << G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay << G4endl;
  G4cout << G4DNAWaterDissociationDisplacer::AutoIonisation << G4endl;
  G4cout << G4DNAWaterDissociationDisplacer::DissociativeAttachment << G4endl;
  bool test1= G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay==
  G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay;
  
  bool test2= G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay==
  G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay;
  
  bool test3= G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay==
  G4DNAWaterDissociationDisplacer::AutoIonisation;
  
  assert((test1 || test2 || test3) == false);
}

//------------------------------------------------------------------------------

int TestDifferentPositions(){
  //auto h2O_conf =
  // G4H2O::Definition()->GetConfigurationWithLabel("B^1A_1");
  auto dissociationTable=G4H2O::Definition()->GetDecayTable();
  auto vecOfDissChannels=dissociationTable->GetDecayChannels("B^1A_1");
  auto autoIoniChannel = (*vecOfDissChannels)[2];
  
  int N=10000;
  
  for(int i=0; i<N;++i){
    auto vectOfDisplacements=
    displacer.GetProductsDisplacement(autoIoniChannel);
    
    bool test1=(vectOfDisplacements[0] == vectOfDisplacements[1]);
    bool test2=(vectOfDisplacements[1] == vectOfDisplacements[2]);
    bool test3=(vectOfDisplacements[0] == vectOfDisplacements[2]);
    bool finalTest= test1 || test2 || test3;
    if(finalTest == true) abort();
  }
  return 0;
}

//------------------------------------------------------------------------------
#ifdef USE_ROOT
int PlotRadialDistributionOfProducts(){
  new TCanvas();
  TH1D* h1 = new TH1D("h_p", "h_p", 25, 0, 3);

  const double Rrms=1.1; //*CLHEP::nanometer;
  int N=10000;
  
  for(int i=0; i<N;++i){
    G4ThreeVector rdVec =
    displacer.radialDistributionOfProducts(Rrms);
    // G4cout << rdVec.mag() << G4endl;
    h1->Fill(rdVec.mag());
  }
  h1->Scale(1./h1->Integral(), "width");
  h1->Draw();
  
  TF1* f1 = new TF1("f_p",
                    "sqrt(2/3.14)*pow(x,2)/pow([0],3)*exp(-0.5*pow(x,2)/pow([0], 2))", 0, 10);
  f1->SetParameter(0, Rrms/sqrt(3));
  f1->Draw("SAME");
  
  return 0;
}

//------------------------------------------------------------------------------

double electDistrib(double* x, double*){
  return G4DNAWaterDissociationDisplacer::ElectronProbaDistribution(x[0]);
}

#ifdef _WATER_DISPLACER_USE_KREIPL_
double electDistrib2(double* x, double*){
  return 4.*x[0]*exp(-2.*x[0]);
}

int PlotRadialDistributionOfElectron(){
  new TCanvas();
  TH1D* h1 = new TH1D("h_e", "h_e", 25, 0, 5);
  int N=10000;

  for(int i=0; i<N;++i){
    G4ThreeVector rdVec =
    displacer.radialDistributionOfElectron();
    h1->Fill(rdVec.mag()/CLHEP::nanometer);
  }
  h1->Scale(1./h1->Integral(), "width");
  h1->Draw();
  
//  TF1* f1 = new TF1("f_e", &electDistrib, 0, 10);
//  f1->Draw("SAME");
  
  TF1* f2 = new TF1("f_e2", &electDistrib2, 0, 10);
  f2->Draw("SAME"); // should match
  return 0;
}
#endif

#ifdef _WATER_DISPLACER_USE_TERRISOL_
double electDistrib2(double* x, double*){
  return 4*pow(x[0],2)/(sqrt(CLHEP::pi)*pow(27.22,3))
        *exp(-pow(x[0],2)/pow(27.22,2));
}

int PlotRadialDistributionOfElectron(){
  new TCanvas();
  TH1D* h1 = new TH1D("h_e", "h_e", 100, 0, 200);
  int N=10000;
  
  for(int i=0; i<N;++i){
    G4ThreeVector rdVec =
    displacer.radialDistributionOfElectron();
    h1->Fill(rdVec.mag()/CLHEP::nanometer);
  }
  h1->Scale(1./h1->Integral(), "width");
  h1->Draw();
  
  //  TF1* f1 = new TF1("f_e", &electDistrib, 0, 10);
  //  f1->Draw("SAME");
  
  TF1* f2 = new TF1("f_e2", &electDistrib2, 0, 100);
  f2->Draw("SAME"); // should match
  return 0;
}
#endif
#endif
//------------------------------------------------------------------------------

int main(int argc, char* argv[]){
  G4Random::setTheSeed(524905525242ul);
//  TestCounter(); return 0;
  G4cout << "start" << G4endl;
  Init();
  TestDifferentIndexes();
  TestDifferentPositions();
  
#ifdef USE_ROOT
//  TApplication* app=new TApplication("testDisplacer", &argc, argv);
  TRint* app=new TRint("testDisplacer", &argc, argv);
  PlotRadialDistributionOfProducts();
  PlotRadialDistributionOfElectron();
  app->Run();
#endif
  return 0;
}
