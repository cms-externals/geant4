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

#include "AnalysisManager.hh"

#include "CalDataCollector.hh"

#include <cmath>
#include <fstream>
#include <iostream>

AnalysisManager* AnalysisManager::instance = 0;

AnalysisManager* AnalysisManager::Instance(G4double zLow, G4double zUp)
{
  if (instance == 0)
  {
    instance = new AnalysisManager(zLow, zUp);
  }

  return instance;
}

void AnalysisManager::Destroy()
{
  delete instance;
  instance = 0;
}

AnalysisManager::AnalysisManager(G4double zLow, G4double zUp)
  : DataManager(zLow, zUp),
    totEnergyDeposit(0),
    totEnergyDepositSumSquares(0),
    elecEnergyEnterTarget(0),
    photEnergyExitTarget(0),
    photEnergyExitTargetSumSquares(0),
    elecEnergyExitTarget(0),
    elecEnergyExitTargetSumSquares(0),
    primElecEnergyExitTarget(0),
    nmbEnergyDeposits(0),
    nmbElecEnterTarget(0),
    nmbElecExitTarget(0),
    nmbPrimElecExitTarget(0),
    nmbPhotExitTarget(0)
{}

AnalysisManager::~AnalysisManager() {}

void AnalysisManager::CreateCalorimeter(G4double pos, G4double thickn, G4double rad)
{
  G4double zLow = pos - 0.5 * thickn;
  G4double zUp = pos + 0.5 * thickn;

  AddDataCollector(new CalDataCollector(zLow, zUp, rad));

  binBoundaries.push_back(zLow);
  std::sort(binBoundaries.begin(), binBoundaries.end());
}

void AnalysisManager::ScoreParticleEnergy(G4double en, G4double x, G4double y, G4double z,
                                          const G4String& ptype)
{
  if (z >= GetLowerBound() && z < GetUpperBound())
  {
    DataManager::ScoreParticleEnergy(en, x, y, z, ptype);
  }
}

void AnalysisManager::ScoreEnergyDeposit(G4double en, G4double x, G4double y, G4double z,
                                         const G4String& ptype)
{
  totEnergyDeposit += en;
  totEnergyDepositSumSquares += en * en;
  nmbEnergyDeposits++;
  DataManager::ScoreEnergyDeposit(en, x, y, z, ptype);
  /*
  if(energyInCalorimeters == 0 && binBoundaries.size() > 1) {

     std::string histName = "EnergyVsDepth";
     std::string histTitle = "Energy deposit in calorimeter slices";

     binBoundaries.push_back(GetUpperBound());

     energyInCalorimeters =
         histogramFactory -> createHistogram1D(histName,
                                               histTitle,
                                               binBoundaries);
  }

  if(energyInCalorimeters != 0) energyInCalorimeters -> fill(z, en);
  */
}

void AnalysisManager::ScoreEnteringParticles(G4double en, G4double, G4double, G4double,
                                             const G4String& ptype)
{
  if (ptype == "electron" || ptype == "e-")
  {
    elecEnergyEnterTarget += en;
    nmbElecEnterTarget++;
  }
}

void AnalysisManager::ScoreExitingParticles(G4double en, G4double, G4double, G4double,
                                            const G4String& ptype, G4int id)
{
  if (ptype == "electron" || ptype == "e-")
  {
    elecEnergyExitTarget += en;

    if (en >= 0.0)
    {
      nmbElecExitTarget += 1;
      elecEnergyExitTargetSumSquares += en * en;
    }
    if (en < 0.0)
    {
      nmbElecExitTarget -= 1;
      elecEnergyExitTargetSumSquares -= en * en;
    }

    if (id == 1)
    {
      primElecEnergyExitTarget += en;

      if (en >= 0.0) nmbPrimElecExitTarget += 1;
      if (en < 0.0) nmbPrimElecExitTarget -= 1;
    }
  }
  if (ptype == "gamma")
  {
    photEnergyExitTarget += en;
    photEnergyExitTargetSumSquares += en * en;
    nmbPhotExitTarget += 1;
  }
}

void AnalysisManager::PrintResults()
{
  G4cout << "--------------------------------------------------------" << G4endl;
  G4cout << "TARGET: " << GetLowerBound() << " to " << GetUpperBound() << G4endl;
  G4cout << "--------------------------------------------------------" << G4endl;

  G4double totEnergyDepositErrSqu = 0.0;
  if (nmbEnergyDeposits > 0)
  {
    totEnergyDepositErrSqu = totEnergyDepositSumSquares / (totEnergyDeposit * totEnergyDeposit)
                             - 1.0 / G4double(nmbEnergyDeposits);
  }
  G4double totEnergyDepositRelErr = std::sqrt(totEnergyDepositErrSqu);

  G4double elecEnergyExitTargetErrSqu = 0.0;
  if (nmbElecExitTarget > 0)
  {
    elecEnergyExitTargetErrSqu =
      elecEnergyExitTargetSumSquares / (elecEnergyExitTarget * elecEnergyExitTarget)
      - 1.0 / G4double(nmbElecExitTarget);
  }
  G4double elecEnergyExitTargetRelErr = std::sqrt(elecEnergyExitTargetErrSqu);

  G4double photEnergyExitTargetErrSqu = 0.0;
  if (nmbPhotExitTarget > 0)
  {
    photEnergyExitTargetErrSqu =
      photEnergyExitTargetSumSquares / (photEnergyExitTarget * photEnergyExitTarget)
      - 1.0 / G4double(nmbPhotExitTarget);
  }
  G4double photEnergyExitTargetRelErr = std::sqrt(photEnergyExitTargetErrSqu);

  G4cout << "  Total Energy Deposit: " << totEnergyDeposit
         << "   Rel. Error: " << totEnergyDepositRelErr << std::endl
         << "  Electron Energy/Number Entering Target: " << elecEnergyEnterTarget << " "
         << nmbElecEnterTarget << G4endl
         << "  Electron Energy/Number Exiting Target:  " << elecEnergyExitTarget << " "
         << nmbElecExitTarget << "  Rel. Error Energy: " << elecEnergyExitTargetRelErr << G4endl
         << "  Primary Electron Energy/Number Exiting Target:  " << primElecEnergyExitTarget << " "
         << nmbPrimElecExitTarget << G4endl
         << "  Photon Energy/Number Exiting Target:    " << photEnergyExitTarget << " "
         << nmbPhotExitTarget << "  Rel. Error Energy: " << photEnergyExitTargetRelErr << G4endl;

  G4double backScEnergy = 1.0 - ((totEnergyDeposit + photEnergyExitTarget) / elecEnergyEnterTarget);

  G4cout << "  Fraction of incid. electr. energy escaping as photons: "
         << photEnergyExitTarget / elecEnergyEnterTarget << G4endl
         << "  Fraction of electron energy backscattered (indir. calcul.): " << backScEnergy
         << G4endl << "  Fraction of electron energy backscattered (direct calcul.): "
         << elecEnergyExitTarget / elecEnergyEnterTarget << G4endl
         << "  Fraction of electron energy backscattered, only primaries (direct calcul.): "
         << primElecEnergyExitTarget / elecEnergyEnterTarget << G4endl
         << "  Fraction of electrons backscattered (direct calcul.): "
         << G4double(nmbElecExitTarget) / G4double(nmbElecEnterTarget) << G4endl
         << "  Fraction of electrons backscattered, only primaries "
         << "(direct calcul.): " << G4double(nmbPrimElecExitTarget) / G4double(nmbElecEnterTarget)
         << G4endl;

  DataManager::PrintResults();
  G4cout << "--------------------------------------------------------" << G4endl;
}
