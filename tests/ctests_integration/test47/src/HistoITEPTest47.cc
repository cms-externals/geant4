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
#include "HistoITEPTest47.hh"

#include "G4UnitsTable.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "CLHEP/Units/PhysicalConstants.h"

HistoITEPTest47::HistoITEPTest47(G4String htitle) : TstHistoSet(htitle)
{
  energies.push_back(0.01);
  energies.push_back(0.03);
  energies.push_back(0.05);
  energies.push_back(0.07);
  energies.push_back(0.09);
  energies.push_back(0.11);
  energies.push_back(0.13);
  energies.push_back(0.15);
  energies.push_back(0.17);
  energies.push_back(0.19);
  energies.push_back(0.21);
  energies.push_back(0.23);
  energies.push_back(0.25);
  angles.push_back(10.1 * CLHEP::deg);
  angles.push_back(15.0 * CLHEP::deg);
  angles.push_back(19.8 * CLHEP::deg);
  angles.push_back(24.8 * CLHEP::deg);
  angles.push_back(29.5 * CLHEP::deg);
  angles.push_back(34.6 * CLHEP::deg);
  angles.push_back(39.6 * CLHEP::deg);
  angles.push_back(44.3 * CLHEP::deg);
  angles.push_back(49.3 * CLHEP::deg);
  angles.push_back(54.2 * CLHEP::deg);
  angles.push_back(59.1 * CLHEP::deg);
  angles.push_back(64.1 * CLHEP::deg);
  angles.push_back(69.1 * CLHEP::deg);
  angles.push_back(74.1 * CLHEP::deg);
  angles.push_back(79.1 * CLHEP::deg);
  angles.push_back(84.1 * CLHEP::deg);
  angles.push_back(89.0 * CLHEP::deg);
  angles.push_back(98.9 * CLHEP::deg);
  angles.push_back(108.9 * CLHEP::deg);
  angles.push_back(119.0 * CLHEP::deg);
  angles.push_back(129.1 * CLHEP::deg);
  angles.push_back(139.1 * CLHEP::deg);
  angles.push_back(149.3 * CLHEP::deg);
  angles.push_back(159.6 * CLHEP::deg);
  angles.push_back(161.4 * CLHEP::deg);
  angles.push_back(165.5 * CLHEP::deg);
  angles.push_back(169.5 * CLHEP::deg);
  angles.push_back(173.5 * CLHEP::deg);
  angles.push_back(177.0 * CLHEP::deg);

  dtheta = 4.0 * CLHEP::deg;
  de = 0.02;

  for (unsigned int ii = 0; ii < angles.size(); ii++)
  {
    double cth1 = std::cos(std::min(angles[ii] + dtheta, 180 * CLHEP::deg));
    double cth2 = std::cos(std::max(angles[ii] - dtheta, 0.0 * CLHEP::deg));
    dcth.push_back(std::abs(cth1 - cth2));
    cthmin.push_back(std::min(cth1, cth2));
    cthmax.push_back(std::max(cth1, cth2));
  }

  for (unsigned int ii = 0; ii < energies.size(); ii++)
  {
    double en = energies[ii];
    emin.push_back(en - 0.5 * de);
    emax.push_back(en + 0.5 * de);
  }

  G4cout << "HistoITEPTest47:: Initialized with " << cthmin.size() << " theta bins and "
         << emin.size() << " energy bins" << G4endl;
  for (unsigned int ii = 0; ii < cthmin.size(); ii++)
    G4cout << "HistoITEPTest47::Bin " << ii << " theta " << angles[ii] / CLHEP::deg
           << " cos(theta) = " << cthmin[ii] << ":" << cthmax[ii] << " dcostheta " << dcth[ii]
           << G4endl;
  for (unsigned int ii = 0; ii < emin.size(); ii++)
    G4cout << "HistoITEPTest47::Bin " << ii << " Energy = " << emin[ii] << ":" << emax[ii]
           << G4endl;

  // OK, done with the leftover from the old code.
  //

  // Now book histogramms.
  //
  std::string tmpname, name, title, t1;

  // Now I need to decompose htitle and extract beam, target, energy(mom).
  // Then I'll have to recompose all that jazz back into name and title,
  // to better fit with older naming scheme.
  // However, some modifications will be introduced, for example, model's
  // name will be removed from the histo name and will remain part of
  // the histo file only.
  //
  std::string::size_type cnt1 =
    htitle.find("+");  // for safety should also check it it's !=std::string::npos - will do later !
  std::string::size_type cnt2 = htitle.find("-");

  // Ugh, re-composed the "core" part of the histo name
  //
  tmpname = htitle.substr(0, cnt1);
  tmpname += htitle.substr(cnt1 + 1, ((cnt2 - cnt1) - 1));
  tmpname += htitle.substr(cnt2 + 1);  // copy from cnt2+1 to the end

  // Now let's add theta bins (for KE spectra)
  //
  for (size_t ii = 0; ii < angles.size(); ++ii)
  {
    // unfortunately, I have to use old-style name formation
    // which calls for an exact number of char's...
    // ... otheriwse I'll have to redo the (Root) analysis scripts
    //
    // std::ostringstream tmp;
    G4double ang = angles[ii] / CLHEP::deg;
    // tmp << ang;
    char tmp[6];
    snprintf(tmp, sizeof(tmp), "%5.1f", ang);

    t1 = "KEproton";
    name = t1 + "0" + tmpname;  // + tmp.str();
    name.append(tmp);
    title = htitle.substr(0, cnt2) + " -> proton at " + htitle.substr(cnt2 + 1)
            + ", theta = ";  // + tmp.str();
    title.append(tmp);
    fKEproton.push_back(new TH1F(name.c_str(), title.c_str(), 800, 0., 8.));
    t1 = "KEneutron";
    name = t1 + "0" + tmpname;  // + tmp.str();
    name.append(tmp);
    title = htitle.substr(0, cnt2) + " -> neutron at " + htitle.substr(cnt2 + 1)
            + ", theta = ";  // + tmp.str();
    title.append(tmp);
    fKEneutron.push_back(new TH1F(name.c_str(), title.c_str(), 800, 0., 8.));
    t1 = "KEpiplus";
    name = t1 + "0" + tmpname;  // + tmp.str();
    name.append(tmp);
    title = htitle.substr(0, cnt2) + " -> piplus at " + htitle.substr(cnt2 + 1)
            + ", theta = ";  // + tmp.str();
    title.append(tmp);
    fKEpiplus.push_back(new TH1F(name.c_str(), title.c_str(), 800, 0., 8.));
    t1 = "KEpiminus";
    name = t1 + "0" + tmpname;  // + tmp.str();
    name.append(tmp);
    title = htitle.substr(0, cnt2) + " -> piminus at " + htitle.substr(cnt2 + 1)
            + ", theta = ";  // + tmp.str();
    title.append(tmp);
    fKEpiminus.push_back(new TH1F(name.c_str(), title.c_str(), 800, 0., 8.));
    t1 = "KEklpus";
    name = t1 + "0" + tmpname;  // + tmp.str();
    name.append(tmp);
    title = htitle.substr(0, cnt2) + " -> kplus at " + htitle.substr(cnt2 + 1)
            + ", theta = ";  // + tmp.str();
    title.append(tmp);
    fKEkplus.push_back(new TH1F(name.c_str(), title.c_str(), 800, 0., 8.));
    t1 = "KEkminus";
    name = t1 + "0" + tmpname;  // + tmp.str();
    name.append(tmp);
    title = htitle.substr(0, cnt2) + " -> kminus at " + htitle.substr(cnt2 + 1)
            + ", theta = ";  // + tmp.str();
    title.append(tmp);
    fKEkminus.push_back(new TH1F(name.c_str(), title.c_str(), 800, 0., 8.));
  }

  // Now add cenergy bins (for cos(theta) spectra)
  //
  for (unsigned int ii = 0; ii < energies.size(); ++ii)
  {
    std::ostringstream tmp;
    G4double ee = energies[ii] / CLHEP::GeV;
    tmp << ee;
    t1 = "CTproton";
    name = t1 + "0" + tmpname + tmp.str();
    title = htitle.substr(0, cnt2) + " -> proton at " + htitle.substr(cnt2 + 1)
            + ", KE = " + tmp.str() + "GeV";
    fCTproton.push_back(new TH1F(name.c_str(), title.c_str(), 80, -1., 1.));
    t1 = "CTneutron";
    name = t1 + "0" + tmpname + tmp.str();
    title = htitle.substr(0, cnt2) + " -> proton at " + htitle.substr(cnt2 + 1)
            + ", KE = " + tmp.str() + "GeV";
    fCTneutron.push_back(new TH1F(name.c_str(), title.c_str(), 80, -1., 1.));
  }
}

HistoITEPTest47::~HistoITEPTest47()
{
  // I think I might need to delete content(s) of the histo copntainers...
  // Although it hasn't been a part of the original code.

  for (size_t i = 0; i < fKEproton.size(); ++i)
  {
    delete fKEproton[i];
  }
  for (size_t i = 0; i < fKEneutron.size(); ++i)
  {
    delete fKEneutron[i];
  }
  for (size_t i = 0; i < fKEpiplus.size(); ++i)
  {
    delete fKEpiplus[i];
  }
  for (size_t i = 0; i < fKEpiminus.size(); ++i)
  {
    delete fKEpiminus[i];
  }
  for (size_t i = 0; i < fKEkplus.size(); ++i)
  {
    delete fKEkplus[i];
  }
  for (size_t i = 0; i < fKEkminus.size(); ++i)
  {
    delete fKEkminus[i];
  }
  for (size_t i = 0; i < fCTproton.size(); ++i)
  {
    delete fCTproton[i];
  }
  for (size_t i = 0; i < fCTneutron.size(); ++i)
  {
    delete fCTneutron[i];
  }
}

void HistoITEPTest47::FillEvt(G4VParticleChange* aChange, const G4LorentzVector&,
                              const G4LorentzVector&)
{
  G4int n = aChange->GetNumberOfSecondaries();

  const G4DynamicParticle* sec = 0;

  for (G4int j = 0; j < n; j++)  // loop over secondaries
  {
    sec = aChange->GetSecondary(j)->GetDynamicParticle();

    const G4String& pname = sec->GetDefinition()->GetParticleName();

    if (pname != "proton" && pname != "neutron" && pname != "pi+" && pname != "pi-"
        && pname != "kaon+" && pname != "kaon-")
      continue;

    G4ThreeVector mom = sec->GetMomentum() / CLHEP::GeV;
    G4double ke = (sec->GetKineticEnergy()) / CLHEP::GeV;
    if (ke < 0.0) ke = 0.0;

    G4double theta = mom.theta();
    G4double cth = std::cos(theta);
    G4double wt = 1. / (mom.mag());

    for (unsigned int ii = 0; ii < angles.size(); ii++)
    {
      if (cth > cthmin[ii] && cth <= cthmax[ii])
      {
        if (pname == "proton")
        {
          fKEproton[ii]->Fill(ke, wt);
        }
        else if (pname == "neutron")
        {
          fKEneutron[ii]->Fill(ke, wt);
        }
        else if (pname == "pi+")
        {
          fKEpiplus[ii]->Fill(ke, wt);
        }
        else if (pname == "pi-")
        {
          fKEpiminus[ii]->Fill(ke, wt);
        }
        else if (pname == "kaon+")
        {
          fKEkplus[ii]->Fill(ke, wt);
        }
        else if (pname == "kaon-")
        {
          fKEkminus[ii]->Fill(ke, wt);
        }
      }
    }

    for (unsigned int ii = 0; ii < energies.size(); ii++)
    {
      if (ke > emin[ii] && ke <= emax[ii])
      {
        if (pname == "proton")
        {
          fCTproton[ii]->Fill(cth, wt);
        }
        else if (pname == "neutron")
        {
          fCTneutron[ii]->Fill(cth, wt);
        }
      }
    }

  }  // end loop over secondaries

  // FIXME !!!
  //
  // --> restore later !!!   epTest.fill(aChange,pinit,aPosition);

  return;
}

void HistoITEPTest47::Write(G4int nevt /* stat */, G4double xsec)
{
  G4cout << " xsec = " << xsec << G4endl;

  for (size_t i = 0; i < fKEproton.size(); ++i)
  {
    double xbin = fKEproton[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * dcth[i]);
    G4cout << " xbin = " << xbin << " scale = " << scale << G4endl;
    fKEproton[i]->Scale(scale);
    fKEproton[i]->Write();
  }
  for (size_t i = 0; i < fKEneutron.size(); ++i)
  {
    double xbin = fKEneutron[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * dcth[i]);
    fKEneutron[i]->Scale(scale);
    fKEneutron[i]->Write();
  }
  for (size_t i = 0; i < fKEpiplus.size(); ++i)
  {
    double xbin = fKEpiplus[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * dcth[i]);
    fKEpiplus[i]->Scale(scale);
    fKEpiplus[i]->Write();
  }
  for (size_t i = 0; i < fKEpiminus.size(); ++i)
  {
    double xbin = fKEpiminus[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * dcth[i]);
    fKEpiminus[i]->Scale(scale);
    fKEpiminus[i]->Write();
  }
  for (size_t i = 0; i < fKEkplus.size(); ++i)
  {
    double xbin = fKEkplus[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * dcth[i]);
    fKEkplus[i]->Scale(scale);
    fKEkplus[i]->Write();
  }
  for (size_t i = 0; i < fKEkminus.size(); ++i)
  {
    double xbin = fKEkminus[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * dcth[i]);
    fKEkminus[i]->Scale(scale);
    fKEkminus[i]->Write();
  }

  for (size_t i = 0; i < fCTproton.size(); ++i)
  {
    double xbin = fCTproton[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * de);
    fCTproton[i]->Scale(scale);
    fCTproton[i]->Write();
  }
  for (size_t i = 0; i < fCTneutron.size(); ++i)
  {
    double xbin = fCTneutron[i]->GetBinWidth(1);
    double scale = xsec / ((double)nevt * xbin * 2. * CLHEP::pi * de);
    fCTneutron[i]->Scale(scale);
    fCTneutron[i]->Write();
  }

  // FIXME !!!
  // Will need to restore EPTest here...

  return;
}
