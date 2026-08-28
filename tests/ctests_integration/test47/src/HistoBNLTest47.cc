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
#include "HistoBNLTest47.hh"

#include "G4UnitsTable.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "CLHEP/Units/PhysicalConstants.h"

HistoBNLTest47::HistoBNLTest47(G4String htitle) : TstHistoSet(htitle)
{
  rapidities.push_back(0.60);
  rapidities.push_back(0.80);
  rapidities.push_back(1.00);
  rapidities.push_back(1.20);
  rapidities.push_back(1.40);
  rapidities.push_back(1.60);
  rapidities.push_back(1.80);
  rapidities.push_back(2.00);
  rapidities.push_back(2.20);
  rapidities.push_back(2.40);
  rapidities.push_back(2.60);
  rapidities.push_back(2.80);
  rapidities.push_back(3.00);

  for (unsigned int ii = 0; ii < rapidities.size() - 1; ii++)
  {
    ymin.push_back(rapidities[ii]);
    ymax.push_back(rapidities[ii + 1]);
  }

  G4cout << "HistoBNLTest47:: Initialized with " << ymin.size() << " rapidity bins" << G4endl;
  for (unsigned int ii = 0; ii < ymin.size(); ii++)
    G4cout << "HistoBNLTest47::Bin " << ii << " rapidity = " << ymin[ii] << ":" << ymax[ii]
           << G4endl;

  std::string::size_type cnt1 =
    htitle.find("+");  // for safety should also check it it's !=std::string::npos - will do later !
  std::string::size_type cnt2 = htitle.find("-");

  G4String particle = htitle.substr(0, cnt1);
  G4String part_target = htitle.substr(0, cnt2);
  G4String target_gen = htitle.substr(cnt1 + 1, ((cnt2 - cnt1) - 1));
  G4String energy = htitle.substr(cnt2 + 1);

  snprintf(tag1Name, sizeof(tag1Name), "%s%s%s", particle.c_str(), target_gen.c_str(),
           energy.c_str());
  snprintf(tag2Name, sizeof(tag2Name), "%s", part_target.c_str());
  // how the hell do I do this, and what for ?
  snprintf(tag3Name, sizeof(tag3Name), "at %s", energy.c_str());

  G4cout << "HistoBNLTest47::fileName:"  // << fileName.c_str()
         << " Tag1:" << tag1Name << " Tag2: " << tag2Name << " Tag3: " << tag3Name << G4endl;

  char name[100], title[160];
  hiMT11.clear();
  hiMT12.clear();
  hiMT21.clear();
  hiMT22.clear();
  hiMT31.clear();
  hiMT32.clear();
  hiMT41.clear();
  hiMT42.clear();
  hiMT51.clear();
  hiMT52.clear();
  hiMT61.clear();
  hiMT62.clear();
  for (unsigned int ii = 0; ii < rapidities.size() - 1; ii++)
  {
    double yv = 0.5 * (rapidities[ii] + rapidities[ii + 1]);
    snprintf(name, sizeof(name), "MTproton1%s%4.2f", tag1Name, yv);
    snprintf(title, sizeof(title), "%s to p %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT11.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTproton2%s%4.2f", tag1Name, yv);
    hiMT12.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTneutron1%s%4.2f", tag1Name, yv);
    snprintf(title, sizeof(title), "%s to n %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT21.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTneutron2%s%4.2f", tag1Name, yv);
    hiMT22.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTpiplus1%s%4.2f", tag1Name, yv);
    snprintf(title, sizeof(title), "%s to #pi+ %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT31.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTpiplus2%s%4.2f", tag1Name, yv);
    hiMT32.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTpiminus1%s%4.2f", tag1Name, yv);
    snprintf(title, sizeof(title), "%s to #pi- %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT41.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTpiminus2%s%4.2f", tag1Name, yv);
    hiMT42.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTkplus1%s%4.2f", tag1Name, yv);
    snprintf(title, sizeof(title), "%s to K+ %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT51.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTkplus2%s%4.2f", tag1Name, yv);
    hiMT52.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTkminus1%s%4.2f", tag1Name, yv);
    snprintf(title, sizeof(title), "%s to K- %s (y = %8.2f)", tag2Name, tag3Name, yv);
    hiMT61.push_back(new TH1F(name, title, 800, 0., 8.));
    snprintf(name, sizeof(name), "MTkminus2%s%4.2f", tag1Name, yv);
    hiMT62.push_back(new TH1F(name, title, 800, 0., 8.));
  }

  // epTest.initialize(particle,target,energy,generator);
}

HistoBNLTest47::~HistoBNLTest47()
{
  for (size_t i = 0; i < hiMT11.size(); ++i)
  {
    delete hiMT11[i];
  }
  for (size_t i = 0; i < hiMT21.size(); ++i)
  {
    delete hiMT21[i];
  }
  for (size_t i = 0; i < hiMT31.size(); ++i)
  {
    delete hiMT31[i];
  }
  for (size_t i = 0; i < hiMT41.size(); ++i)
  {
    delete hiMT41[i];
  }
  for (size_t i = 0; i < hiMT51.size(); ++i)
  {
    delete hiMT51[i];
  }
  for (size_t i = 0; i < hiMT61.size(); ++i)
  {
    delete hiMT61[i];
  }
  for (size_t i = 0; i < hiMT12.size(); ++i)
  {
    delete hiMT12[i];
  }
  for (size_t i = 0; i < hiMT22.size(); ++i)
  {
    delete hiMT22[i];
  }
  for (size_t i = 0; i < hiMT32.size(); ++i)
  {
    delete hiMT32[i];
  }
  for (size_t i = 0; i < hiMT42.size(); ++i)
  {
    delete hiMT42[i];
  }
  for (size_t i = 0; i < hiMT52.size(); ++i)
  {
    delete hiMT52[i];
  }
  for (size_t i = 0; i < hiMT62.size(); ++i)
  {
    delete hiMT62[i];
  }
}

void HistoBNLTest47::FillEvt(G4VParticleChange* aChange, const G4LorentzVector& pinit,
                             const G4LorentzVector&)
{
  G4LorentzVector labv(pinit), fm;
  G4int n = aChange->GetNumberOfSecondaries();
  const G4DynamicParticle* sec = 0;
  G4ThreeVector mom;
  G4ParticleDefinition* pd;

  for (G4int j = 0; j < n; j++)
  {
    sec = aChange->GetSecondary(j)->GetDynamicParticle();
    pd = sec->GetDefinition();
    const G4String& pname = pd->GetParticleName();

    if (pname != "proton" && pname != "neutron" && pname != "pi+" && pname != "pi-"
        && pname != "kaon+" && pname != "kaon-")
      continue;

    mom = sec->GetMomentumDirection();
    G4double ke = (sec->GetKineticEnergy()) / CLHEP::GeV;
    if (ke < 0.0) ke = 0.0;
    G4double m = (pd->GetPDGMass()) / CLHEP::GeV;
    G4double p = std::sqrt(ke * (ke + 2.0 * m));
    G4double ee = ke + m;
    mom *= p;
    fm = G4LorentzVector(mom, ee);
    labv -= fm;
    G4double theta = mom.theta();
    G4double pt = p * std::sin(theta);
    G4double pl = p * std::cos(theta);
    G4double mt = std::sqrt(pt * pt + m * m);
    G4double mtp = (mt - std::abs(m));
    G4double yv = 0.5 * std::log((ee + pl) / (ee - pl));
    G4double wt = 1. / mt;
    for (unsigned int ii = 0; ii < ymin.size(); ii++)
    {
      if (yv > ymin[ii] && yv <= ymax[ii])
      {
        if (pname == "proton")
        {
          hiMT11[ii]->Fill(mtp);
          hiMT12[ii]->Fill(mtp, wt);
        }
        else if (pname == "neutron")
        {
          hiMT21[ii]->Fill(mtp);
          hiMT22[ii]->Fill(mtp, wt);
        }
        else if (pname == "pi+")
        {
          hiMT31[ii]->Fill(mtp);
          hiMT32[ii]->Fill(mtp, wt);
        }
        else if (pname == "pi-")
        {
          hiMT41[ii]->Fill(mtp);
          hiMT42[ii]->Fill(mtp, wt);
        }
        else if (pname == "kaon+")
        {
          hiMT51[ii]->Fill(mtp);
          hiMT52[ii]->Fill(mtp, wt);
        }
        else if (pname == "kaon-")
        {
          hiMT61[ii]->Fill(mtp);
          hiMT62[ii]->Fill(mtp, wt);
        }
      }
    }
  }

  // FIXME !!!
  // temporarily disabled !
  // will be updated shortly
  //
  // epTest.fill(aChange,pinit,aPosition);

  return;
}

void HistoBNLTest47::Write(G4int nevt, G4double cross_sec)
{
  char name[100], title[100];
  G4double xbin, dy, yv, scale;
  std::vector<TH1F*> hiMT10, hiMT20, hiMT30, hiMT40, hiMT50, hiMT60;
  for (unsigned int ii = 0; ii < ymin.size(); ii++)
  {
    dy = (ymax[ii] - ymin[ii]);
    yv = 0.5 * (ymin[ii] + ymax[ii]);

    xbin = hiMT11[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of p (GeV)");
    hiMT11[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.3f GeV", xbin);
    hiMT11[ii]->GetYaxis()->SetTitle(title);
    xbin = hiMT12[ii]->GetBinWidth(1);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin * 2. * CLHEP::pi * dy);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of p (GeV)");
    hiMT12[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{1}{m_{T}})/%6.3f GeV", xbin);
    hiMT12[ii]->GetYaxis()->SetTitle(title);
    snprintf(name, sizeof(name), "MTproton0%s%4.2f", tag1Name, yv);
    hiMT10.push_back((TH1F*)hiMT12[ii]->Clone());
    hiMT10[ii]->SetName(name);
    hiMT10[ii]->Scale(scale);
    hiMT10[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiMT21[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of n (GeV)");
    hiMT21[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.3f GeV", xbin);
    hiMT21[ii]->GetYaxis()->SetTitle(title);
    xbin = hiMT22[ii]->GetBinWidth(1);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin * 2. * CLHEP::pi * dy);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of n (GeV)");
    hiMT22[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{1}{m_{T}})/%6.3f GeV", xbin);
    hiMT22[ii]->GetYaxis()->SetTitle(title);
    snprintf(name, sizeof(name), "MTneutron0%s%4.2f", tag1Name, yv);
    hiMT20.push_back((TH1F*)hiMT22[ii]->Clone());
    hiMT20[ii]->SetName(name);
    hiMT20[ii]->Scale(scale);
    hiMT20[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiMT31[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of #pi+ (GeV)");
    hiMT31[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.3f GeV", xbin);
    hiMT31[ii]->GetYaxis()->SetTitle(title);
    xbin = hiMT32[ii]->GetBinWidth(1);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin * 2. * CLHEP::pi * dy);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of #pi+ (GeV)");
    hiMT32[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{1}{m_{T}})/%6.3f GeV", xbin);
    hiMT32[ii]->GetYaxis()->SetTitle(title);
    snprintf(name, sizeof(name), "MTpiplus0%s%4.2f", tag1Name, yv);
    hiMT30.push_back((TH1F*)hiMT32[ii]->Clone());
    hiMT30[ii]->SetName(name);
    hiMT30[ii]->Scale(scale);
    hiMT30[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiMT41[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of #pi- (GeV)");
    hiMT41[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.3f GeV", xbin);
    hiMT41[ii]->GetYaxis()->SetTitle(title);
    xbin = hiMT42[ii]->GetBinWidth(1);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin * 2. * CLHEP::pi * dy);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of #pi- (GeV)");
    hiMT42[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{1}{m_{T}})/%6.3f GeV", xbin);
    hiMT42[ii]->GetYaxis()->SetTitle(title);
    snprintf(name, sizeof(name), "MTpiminus0%s%4.2f", tag1Name, yv);
    hiMT40.push_back((TH1F*)hiMT42[ii]->Clone());
    hiMT40[ii]->SetName(name);
    hiMT40[ii]->Scale(scale);
    hiMT40[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiMT51[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of K+ (GeV)");
    hiMT51[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.3f GeV", xbin);
    hiMT51[ii]->GetYaxis()->SetTitle(title);
    xbin = hiMT52[ii]->GetBinWidth(1);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin * 2. * CLHEP::pi * dy);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of K+ (GeV)");
    hiMT52[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{1}{m_{T}})/%6.3f GeV", xbin);
    hiMT52[ii]->GetYaxis()->SetTitle(title);
    snprintf(name, sizeof(name), "MTkplus0%s%4.2f", tag1Name, yv);
    hiMT50.push_back((TH1F*)hiMT52[ii]->Clone());
    hiMT50[ii]->SetName(name);
    hiMT50[ii]->Scale(scale);
    hiMT50[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");

    xbin = hiMT61[ii]->GetBinWidth(1);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of K- (GeV)");
    hiMT61[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events/%6.3f GeV", xbin);
    hiMT61[ii]->GetYaxis()->SetTitle(title);
    xbin = hiMT62[ii]->GetBinWidth(1);
    scale = cross_sec / (((double)(std::max(nevt, 1))) * xbin * 2. * CLHEP::pi * dy);
    snprintf(title, sizeof(title), "Reduced Transverse Mass of K- (GeV)");
    hiMT62[ii]->GetXaxis()->SetTitle(title);
    snprintf(title, sizeof(title), "Events (scaled by #frac{1}{m_{T}})/%6.3f GeV", xbin);
    hiMT62[ii]->GetYaxis()->SetTitle(title);
    snprintf(name, sizeof(name), "MTkminus0%s%4.2f", tag1Name, yv);
    hiMT60.push_back((TH1F*)hiMT62[ii]->Clone());
    hiMT60[ii]->SetName(name);
    hiMT60[ii]->Scale(scale);
    hiMT60[ii]->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}} (mb/GeV^{2})");
  }

  //   TFile f(fileName.c_str(), "recreate");
  for (unsigned int ii = 0; ii < ymin.size(); ii++)
  {
    hiMT11[ii]->Write();
    hiMT10[ii]->Write();
    hiMT12[ii]->Write();
    hiMT21[ii]->Write();
    hiMT20[ii]->Write();
    hiMT22[ii]->Write();
    hiMT31[ii]->Write();
    hiMT30[ii]->Write();
    hiMT32[ii]->Write();
    hiMT41[ii]->Write();
    hiMT40[ii]->Write();
    hiMT42[ii]->Write();
    hiMT51[ii]->Write();
    hiMT50[ii]->Write();
    hiMT52[ii]->Write();
    hiMT61[ii]->Write();
    hiMT60[ii]->Write();
    hiMT62[ii]->Write();
  }

  // FIXME !
  // tmp desiabled, will be updated shortly
  //
  //   epTest.write();

  //   f.Close();

  return;
}
