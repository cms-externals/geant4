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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"

#include "G4UnitsTable.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager(G4int nbAbsor, G4int maxAbsor) : fNbAbsor(nbAbsor), fMaxAbsor(maxAbsor)
{
  fMaxHisto = 2 * fMaxAbsor + 3;

  fHistoId.assign(fMaxHisto + 1, -1);
  fLabel.resize(fMaxHisto + 1);
  fTitle.resize(fMaxHisto + 1);
  fNbins.assign(fMaxHisto + 1, 0);
  fVmin.resize(fMaxHisto + 1);
  fVmax.resize(fMaxHisto + 1);
  fUnit.resize(fMaxHisto + 1);
  fExist.resize(fMaxHisto + 1);

  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create or get analysis manager
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  // enable inactivation of histograms

  // create absorber histogramms
  G4int nbins = 100;
  G4double vmin = 0., vmax = 10.;
  const G4String vunit_MeV = "MeV";
  G4int iCmpt = 1;
  for (G4int k = 1; k < fMaxAbsor; k++)
    SetHisto(iCmpt++, nbins, vmin, vmax, vunit_MeV);

  const G4String vunit_none = "none";
  nbins = 52;
  vmin = 0.;
  vmax = 52.;
  for (G4int k = fMaxAbsor; k < 2 * fMaxAbsor; k++)
    SetHisto(iCmpt++, nbins, vmin, vmax, vunit_none);

  nbins = 102;
  vmin = 0.;
  vmax = 102.;
  for (G4int k = 2 * fMaxAbsor; k <= fMaxHisto; k++)
    SetHisto(iCmpt++, nbins, vmin, vmax, vunit_none);

  // create selected histograms
  for (G4int k = 1; k <= fMaxHisto; k++)
  {
    bool bActive = false;
    if (fNbins[k] > 0)
    {
      bActive = true;

      if (k > fNbAbsor && k < fMaxAbsor) bActive = false;
      if (k > fMaxAbsor + fNbAbsor && k < 2 * fMaxAbsor) bActive = false;

      fHistoId[k] = analysisManager->CreateH1(fLabel[k], fTitle[k], fNbins[k], fVmin[k], fVmax[k]);
      analysisManager->SetH1Activation(fHistoId[k], bActive);
    }
    else
    {
      fNbins[k] = 10.;
      std::stringstream s;
      s << k;
      fLabel[k] = s.str();
      fTitle[k] = "Dummy";
      fVmin[k] = 0.;
      fVmax[k] = 100;
      fHistoId[k] = analysisManager->CreateH1(fLabel[k], fTitle[k], fNbins[k], fVmin[k], fVmax[k]);
      analysisManager->SetH1Activation(fHistoId[k], bActive);
    }
    fExist[k] = bActive;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih, G4int nbins, G4double valmin, G4double valmax,
                            const G4String& unit)
{
  if (ih < 1 || ih >= fMaxHisto)
  {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih << "does not exist"
           << G4endl;
    return;
  }

  // histo 1 : energy deposit in absorber 1
  // histo 2 : energy deposit in absorber 2
  // ...etc...........
  // MaxAbsor = 10 (-1)
  //
  // histo 11 : longitudinal profile of energy deposit in absorber 1 (MeV)
  // histo 12 : longitudinal profile of energy deposit in absorber 2 (MeV)
  // ...etc...........
  //
  // histo 21 : energy flow (MeV)
  // histo 22 : lateral energy leak (MeV)

  const G4String id[] = {"0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11",
                         "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"};

  G4String title;
  G4double vmin = valmin, vmax = valmax;
  G4double vunit = 1.;

  if (ih < fMaxAbsor)
  {
    title = "Edep in absorber " + id[ih] + " (" + unit + ")";
    vunit = G4UnitDefinition::GetValueOf(unit);
    vmin = valmin / vunit;
    vmax = valmax / vunit;
  }
  else if (ih > fMaxAbsor && ih < 2 * fMaxAbsor)
  {
    title = "longit. profile of Edep (MeV/event) in absorber " + id[ih - fMaxAbsor];
  }
  else if (ih == 2 * fMaxAbsor + 1)
  {
    title = "energy flow (MeV/event)";
  }
  else if (ih == 2 * fMaxAbsor + 2)
  {
    title = "lateral energy leak (MeV/event)";
  }
  else
    return;

  fLabel[ih] = id[ih];
  fTitle[ih] = title;
  fNbins[ih] = nbins;
  fVmin[ih] = vmin;
  fVmax[ih] = vmax;
  fUnit[ih] = vunit;

  /*  G4cout << "----> SetHisto " << ih << ": " << title << ";  "
         << nbins << " bins from "
         << vmin << " " << unit << " to " << vmax << " " << unit << G4endl;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= fMaxHisto)
  {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih << " xbin= " << xbin
           << " weight= " << weight << G4endl;
    return;
  }

  G4AnalysisManager::Instance()->FillH1(fHistoId[ih], xbin / fUnit[ih], weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double norm)
{
  if (ih >= fMaxHisto)
  {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << " undefined histogramm" << G4endl;
    return;
  }

  // histogramm not defined
  if (!fExist[ih]) return;
  G4AnalysisManager::Instance()->GetH1(fHistoId[ih])->scale(norm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
