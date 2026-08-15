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
#include "HistoMIPSTest47.hh"

#include "G4UnitsTable.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "CLHEP/Units/PhysicalConstants.h"

HistoMIPSTest47::HistoMIPSTest47(G4String htitle) : TstHistoSet(htitle)
{
  xminCalorimeter = -49.5 * CLHEP::cm;
  xmaxCalorimeter = 49.5 * CLHEP::cm;
  yminCalorimeter = -49.0 * CLHEP::cm;
  ymaxCalorimeter = 49.0 * CLHEP::cm;
  rmaxCalorimeter = 45.0 * CLHEP::cm;
  zCalorimeter = 2596.0 * CLHEP::cm;
  sterad = 0.0013;

  G4cout << "HistoMIPSTest47:: Initialized with Calorimeter position " << zCalorimeter
         << " Ranges (rmax:xmin:xmax:ymin:ymax) " << rmaxCalorimeter << ":" << xminCalorimeter
         << ":" << xmaxCalorimeter << ":" << yminCalorimeter << ":" << ymaxCalorimeter << " sterad "
         << sterad << G4endl;

  std::string::size_type cnt1 =
    htitle.find("+");  // for safety should also check it it's !=std::string::npos - will do later !
  std::string::size_type cnt2 = htitle.find("-");

  G4String particle = htitle.substr(0, cnt1);
  G4String part_target = htitle.substr(0, cnt2);
  G4String target_gen = htitle.substr(cnt1 + 1, ((cnt2 - cnt1) - 1));
  G4String ene = htitle.substr(cnt2 + 1);

  std::string::size_type cnt3 = ene.find("G");

  G4double energy = atof((ene.substr(0, cnt3)).c_str());

  if (energy < 60.0)
  {
    eMin = 12.0;
  }
  else if (energy > 100.0)
  {
    eMin = 20.0;
  }
  else
  {
    eMin = 18.0;
  }

  G4cout << "HistoMIPSTest47::initialize invoked: E(eMin) " << energy << "(" << eMin << ")"
         << G4endl;

  // NOTE (JV): "GeV" is already embeded in "ene"
  //
  //  snprintf ( tag1Name, sizeof(tag1Name), "%s%s%sGeV", particle.c_str(), target_gen.c_str(),
  //  ene.c_str() );
  snprintf(tag1Name, sizeof(tag1Name), "%s%s%s", particle.c_str(), target_gen.c_str(), ene.c_str());
  snprintf(tag2Name, sizeof(tag2Name), "%s", part_target.c_str());
  snprintf(tag3Name, sizeof(tag3Name), "at %s", ene.c_str());

  //  G4cout << "HistoMIPSTest47::fileName:" << fileName.c_str()
  //	 << " Tag1:" << tag1Name << " Tag2: " << tag2Name << " Tag3: "
  //	 << tag3Name << " eMin: " << eMin << G4endl;

  char name[100], title[160];
  hiPL1.clear();
  hiPL2.clear();
  hiPL3.clear();
  hiXF1.clear();
  hiXF2.clear();
  hiXF3.clear();
  hiXW1.clear();
  hiXW2.clear();
  hiXW3.clear();
  std::string particles[8] = {"proton", "neutron", "piplus", "piminus",
                              "Kplus",  "Kminus",  "pbar",   "pizero"};

  for (unsigned int ii = 0; ii < 8; ii++)
  {
    snprintf(name, sizeof(name), "PL%s1%s", particles[ii].c_str(), tag1Name);
    snprintf(title, sizeof(title), "%s to %s %s", tag2Name, particles[ii].c_str(), tag3Name);
    hiPL1.push_back(new TH1F(name, title, 150, 0., 150.));
    hiPL1[ii]->Sumw2();
    snprintf(name, sizeof(name), "PL%s2%s", particles[ii].c_str(), tag1Name);
    hiPL2.push_back(new TH1F(name, title, 150, 0., 150.));
    hiPL2[ii]->Sumw2();
    snprintf(name, sizeof(name), "PL%s3%s", particles[ii].c_str(), tag1Name);
    hiPL3.push_back(new TH1F(name, title, 150, 0., 150.));
    hiPL3[ii]->Sumw2();
    snprintf(name, sizeof(name), "XF%s1%s", particles[ii].c_str(), tag1Name);
    hiXF1.push_back(new TH1F(name, title, 40, -1.0, 1.0));
    hiXF1[ii]->Sumw2();
    snprintf(name, sizeof(name), "XF%s2%s", particles[ii].c_str(), tag1Name);
    hiXF2.push_back(new TH1F(name, title, 40, -1.0, 1.0));
    hiXF2[ii]->Sumw2();
    snprintf(name, sizeof(name), "XF%s3%s", particles[ii].c_str(), tag1Name);
    hiXF3.push_back(new TH1F(name, title, 40, -1.0, 1.0));
    hiXF3[ii]->Sumw2();
    snprintf(name, sizeof(name), "XW%s1%s", particles[ii].c_str(), tag1Name);
    hiXW1.push_back(new TH1F(name, title, 40, -1.0, 1.0));
    hiXW1[ii]->Sumw2();
    snprintf(name, sizeof(name), "XW%s2%s", particles[ii].c_str(), tag1Name);
    hiXW2.push_back(new TH1F(name, title, 40, -1.0, 1.0));
    hiXW2[ii]->Sumw2();
    snprintf(name, sizeof(name), "XW%s3%s", particles[ii].c_str(), tag1Name);
    hiXW3.push_back(new TH1F(name, title, 40, -1.0, 1.0));
    hiXW3[ii]->Sumw2();
  }

  // FIXME !!!
  //  epTest.initialize(particle,target,energy,generator);
  //  epTest.setDebug(debug);
}

HistoMIPSTest47::~HistoMIPSTest47()
{
  for (size_t i = 0; i < hiPL1.size(); ++i)
  {
    delete hiPL1[i];
  }
  for (size_t i = 0; i < hiXF1.size(); ++i)
  {
    delete hiXF1[i];
  }
  for (size_t i = 0; i < hiXW1.size(); ++i)
  {
    delete hiXW1[i];
  }
  for (size_t i = 0; i < hiPL2.size(); ++i)
  {
    delete hiPL2[i];
  }
  for (size_t i = 0; i < hiXF2.size(); ++i)
  {
    delete hiXF2[i];
  }
  for (size_t i = 0; i < hiXW2.size(); ++i)
  {
    delete hiXW2[i];
  }
  for (size_t i = 0; i < hiPL3.size(); ++i)
  {
    delete hiPL3[i];
  }
  for (size_t i = 0; i < hiXF3.size(); ++i)
  {
    delete hiXF3[i];
  }
  for (size_t i = 0; i < hiXW3.size(); ++i)
  {
    delete hiXW3[i];
  }
}

void HistoMIPSTest47::FillEvt(G4VParticleChange* aChange, const G4LorentzVector& pinit,
                              const G4LorentzVector& labp)
{
  //   if (unInitialized) initialize();

  G4LorentzVector labv(pinit), fm;
  G4Track* sectrk = 0;
  const G4DynamicParticle* sec = 0;
  G4ParticleDefinition* pd = 0;

  G4int n = aChange->GetNumberOfSecondaries();

  G4ThreeVector mom;

  // FIXME !!
  // need position to calculate acceptance... ouch, terrible (original) code structure !
  ///
  G4double distz = zCalorimeter;  // - (aChange->GetPosition()).z();

  G4ThreeVector boostCM = labp.findBoostToCM();
  G4double roots = labp.mag();

  //  if (debug) G4cout << "HistoMIPSTest47::fill CM system " << labv << " | "
  //		    << labp << " Boost " << boostCM << " roots " << roots
  //		    << " N " << n  << G4endl;

  for (G4int j = 0; j < n; j++)
  {
    sectrk = aChange->GetSecondary(j);
    sec = sectrk->GetDynamicParticle();
    pd = sec->GetDefinition();
    //    const G4String& pname = sec->GetDefinition()->GetParticleName();
    const G4String& pname = pd->GetParticleName();
    mom = sec->GetMomentumDirection();
    G4double ke = (sec->GetKineticEnergy()) / CLHEP::GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m = (pd->GetPDGMass()) / CLHEP::GeV;
    G4double p = std::sqrt(ke * (ke + 2.0 * m));
    G4double pmax = std::sqrt(0.25 * roots * roots - m * m);
    mom *= p;
    fm = G4LorentzVector(mom, ke + m);
    labv -= fm;

    G4int type = particleType(pname);

    G4bool acc1 = false, acc2 = false;

    //    if (debug) G4cout << "Particle " << j << " type " << type << " p " << mom
    //		      << " p|m|pmax " << p << "|" << m << "|" << pmax <<G4endl;

    distz = zCalorimeter - sectrk->GetPosition().z();

    if (mom.z() > 0 && p > eMin)
    {
      G4double xCal = (distz * mom.x() / mom.z()) + sectrk->GetPosition().x();
      G4double yCal = (distz * mom.y() / mom.z()) + sectrk->GetPosition().y();
      G4double rCal = std::sqrt(xCal * xCal + yCal * yCal);
      if (rCal < rmaxCalorimeter) acc1 = true;
      if (xCal > xminCalorimeter && xCal < xmaxCalorimeter && yCal > yminCalorimeter
          && yCal < ymaxCalorimeter)
        acc2 = true;

      //      if (debug) G4cout << "x/y/r at Calorimeter " << xCal << ":" << yCal
      //			<< ":" << rCal << " Accept " << acc1 << ":" << acc2
      //			<< G4endl;
    }

    if (type >= 0 && type <= 7)
    {
      G4LorentzVector fmc = fm.boost(boostCM.x(), boostCM.y(), boostCM.z());
      G4double xf = fmc.z() / pmax;
      G4double wt = (fmc.e() * fmc.z()) / (fmc.rho() * fmc.rho() * fmc.rho() * pmax);

      //      if (debug) G4cout << "Type " << type << " CM Momentum " << fmc << " p|pt "
      //			<< fmc.rho() << "|" << fmc.perp() << "|" << pmax
      //			<< " xf " << xf << " wt " << wt << " Acc: "
      //			<< acc1 << ":" << acc2 << G4endl;

      hiPL1[type]->Fill(p);
      hiXF1[type]->Fill(xf);
      hiXW1[type]->Fill(xf, wt);
      if (acc1)
      {
        hiPL2[type]->Fill(p);
        hiXF2[type]->Fill(xf);
        hiXW2[type]->Fill(xf, wt);
      }
      if (acc2)
      {
        hiPL3[type]->Fill(p);
        hiXF3[type]->Fill(xf);
        hiXW3[type]->Fill(xf, wt);
      }
    }
  }

  // FIXME !!!
  //   epTest.fill(aChange,pinit,aPosition);

  return;
}

void HistoMIPSTest47::Write(G4int nevt, G4double cross_sec)
{
  char name[100], title[100];
  G4double xbin, scale;
  std::vector<TH1F*> hiPL0, hiXF0, hiXW0;
  std::vector<TH1F*> hiPLx, hiXFx, hiXWx;
  std::string psymbs[8] = {"p", "n", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "pbar", "#pi^{0}"};
  std::string particles[8] = {"proton", "neutron", "piplus", "piminus",
                              "Kplus",  "Kminus",  "pbar",   "pizero"};
  for (unsigned int ii = 0; ii < 8; ii++)
  {
    xbin = hiPL1[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Laboratory momentum of %s (GeV/c)", psymbs[ii].c_str());
    hiPL1[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.2f GeV/c", xbin);
    hiPL1[ii]->GetYaxis()->SetTitle(title);
    xbin = hiPL2[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Laboratory momentum of %s (GeV/c)", psymbs[ii].c_str());
    hiPL2[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.2f GeV/c", xbin);
    hiPL2[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin);
    snprintf(name, sizeof(name), "PL%s0%s", particles[ii].c_str(), tag1Name);
    hiPL0.push_back((TH1F*)hiPL2[ii]->Clone());
    hiPL0[ii]->SetName(name);
    hiPL0[ii]->Scale(scale);
    hiPL0[ii]->GetYaxis()->SetTitle("#frac{d#sigma}{dp} (mb/GeV)");
    xbin = hiPL3[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Laboratory momentum of %s (GeV/c)", psymbs[ii].c_str());
    hiPL3[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.2f GeV/c", xbin);
    hiPL3[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin);
    snprintf(name, sizeof(name), "PL%sx%s", particles[ii].c_str(), tag1Name);
    hiPLx.push_back((TH1F*)hiPL3[ii]->Clone());
    hiPLx[ii]->SetName(name);
    hiPLx[ii]->Scale(scale);
    hiPLx[ii]->GetYaxis()->SetTitle("#frac{d#sigma}{dp} (mb/GeV)");

    xbin = hiXF1[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "x_{F} of %s", psymbs[ii].c_str());
    hiXF1[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%5.2f", xbin);
    hiXF1[ii]->GetYaxis()->SetTitle(title);
    xbin = hiXF2[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "x_{F} of %s", psymbs[ii].c_str());
    hiXF2[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%5.2f", xbin);
    hiXF2[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin);
    snprintf(name, sizeof(name), "XF%s0%s", particles[ii].c_str(), tag1Name);
    hiXF0.push_back((TH1F*)hiXF2[ii]->Clone());
    hiXF0[ii]->SetName(name);
    hiXF0[ii]->Scale(scale);
    hiXF0[ii]->GetYaxis()->SetTitle(title);
    xbin = hiXF3[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "x_{F} of %s", psymbs[ii].c_str());
    hiXF3[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%5.2f", xbin);
    hiXF3[ii]->GetYaxis()->SetTitle(title);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin);
    snprintf(name, sizeof(name), "XF%sx%s", particles[ii].c_str(), tag1Name);
    hiXFx.push_back((TH1F*)hiXF3[ii]->Clone());
    hiXFx[ii]->SetName(name);
    hiXFx[ii]->Scale(scale);
    hiXFx[ii]->GetYaxis()->SetTitle(title);

    xbin = hiXW1[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "x_{F} of %s", psymbs[ii].c_str());
    hiXW1[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{E}{p^{2}})/%5.2f", xbin);
    hiXW1[ii]->GetYaxis()->SetTitle(title);
    xbin = hiXW2[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "x_{F} of %s", psymbs[ii].c_str());
    hiXW2[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{E}{p^{2}})/%5.2f", xbin);
    hiXW2[ii]->GetYaxis()->SetTitle(title);
    scale = (4. * CLHEP::pi * cross_sec) / (((double)(std::max(nevt, 1))) * xbin);
    snprintf(name, sizeof(name), "XW%s0%s", particles[ii].c_str(), tag1Name);
    hiXW0.push_back((TH1F*)hiXW2[ii]->Clone());
    hiXW0[ii]->SetName(name);
    hiXW0[ii]->Scale(scale);
    hiXW0[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
    xbin = hiXW3[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "x_{F} of %s", psymbs[ii].c_str());
    hiXW3[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{E}{p^{2}})/%5.2f", xbin);
    hiXW3[ii]->GetYaxis()->SetTitle(title);
    scale = (4. * CLHEP::pi * cross_sec) / (((double)(std::max(nevt, 1))) * xbin);
    snprintf(name, sizeof(name), "XW%sx%s", particles[ii].c_str(), tag1Name);
    hiXWx.push_back((TH1F*)hiXW3[ii]->Clone());
    hiXWx[ii]->SetName(name);
    hiXWx[ii]->Scale(scale);
    hiXWx[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  }

  for (unsigned int ii = 0; ii < 8; ii++)
  {
    hiPL0[ii]->Write();
    hiPL1[ii]->Write();
    hiPL2[ii]->Write();
    hiPL3[ii]->Write();
    hiPLx[ii]->Write();
    hiXF0[ii]->Write();
    hiXF1[ii]->Write();
    hiXF2[ii]->Write();
    hiXF3[ii]->Write();
    hiXFx[ii]->Write();
    hiXW0[ii]->Write();
    hiXW1[ii]->Write();
    hiXW2[ii]->Write();
    hiXW3[ii]->Write();
    hiXWx[ii]->Write();
  }

  // FIXME !!!
  //  epTest.write();

  return;
}

G4int HistoMIPSTest47::particleType(G4String pname)
{
  if (pname == "proton") return 0;
  if (pname == "neutron") return 1;
  if (pname == "pi+") return 2;
  if (pname == "pi-") return 3;
  if (pname == "kaon+") return 4;
  if (pname == "kaon-") return 5;
  if (pname == "anti_proton") return 6;
  if (pname == "pi0") return 7;

  return -1;
}
