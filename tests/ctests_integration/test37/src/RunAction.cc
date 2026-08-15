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

#include "RunAction.hh"

#include "G4ProcessManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"

#include <fstream>
#include <iomanip>
#include <map>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
  : detector(det), primary(kin)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {}

static const std::map<std::string, G4double> density = {
  {"G4_Ta", 16.654 * g / cm3},
  {"G4_U", 18.95 * g / cm3},
  {"G4_Mo", 10.22 * g / cm3},
  {"G4_Al", 2.699 * g / cm3},
  {"G4_Be", 1.848 * g / cm3},
  {"G4_Au", 19.32 * g / cm3},
  {"G4_Si", 2.33 * g / cm3},
  // 1.7 for graphite
  // 2.0 AMORPHOUS C
  // https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=906
  {"G4_C", 2. * g / cm3},
  {"G4_Fe", 7.874 * g / cm3},
  {"G4_Cu", 8.96 * g / cm3}};

// metarial, [(energy, mfp), ...]
static const std::map<std::string, std::vector<std::pair<G4double, G4double>>> mfp = {
  {"G4_Ta",
   {{0.300 * MeV, 0.160 * g / cm2},
    {0.314 * MeV, 0.167 * g / cm2},
    {0.500 * MeV, 0.325 * g / cm2},
    {0.521 * MeV, 0.339 * g / cm2},
    {1.000 * MeV, 0.763 * g / cm2},
    {1.033 * MeV, 0.788 * g / cm2}}},
  {"G4_U",
   {{0.300 * MeV, 0.174 * g / cm2},
    {0.500 * MeV, 0.349 * g / cm2},
    {1.000 * MeV, 0.809 * g / cm2}}},
  {"G4_Mo",
   {{0.100 * MeV, 0.0246 * g / cm2},
    {0.300 * MeV, 0.137 * g / cm2},
    {0.500 * MeV, 0.281 * g / cm2},
    {1.000 * MeV, 0.673 * g / cm2}}},
  {"G4_Al",
   {{0.058 * MeV, 0.00663 * g / cm2},
    {0.109 * MeV, 0.0203 * g / cm2},
    {0.314 * MeV, 0.113 * g / cm2},
    {0.521 * MeV, 0.234 * g / cm2},
    {1.000 * MeV, 0.551 * g / cm2},
    {1.033 * MeV, 0.569 * g / cm2}}},
  {"G4_Be",
   {{0.058 * MeV, 0.00604 * g / cm2},
    {0.109 * MeV, 0.0189 * g / cm2},
    {0.314 * MeV, 0.107 * g / cm2},
    {0.521 * MeV, 0.224 * g / cm2},
    {1.000 * MeV, 0.536 * g / cm2},
    {1.033 * MeV, 0.554 * g / cm2}}},
  {"G4_Au", {{1.000 * MeV, 0.772 * g / cm2}}},
  {"G4_Si",
   {
     {15.00 * keV, 0.69394 * mg / cm2},
     {20.00 * keV, 1.1409 * mg / cm2},
     {30.00 * keV, 2.3056 * mg / cm2},
     {40.00 * keV, 3.7980 * mg / cm2},
     {50.00 * keV, 5.5866 * mg / cm2},
   }},
  {"G4_C", {{1.0 * MeV, 4.89e-1 * g / cm2}}},
  {"G4_Fe",
   {{0.3 * MeV, 1.21e-1 * g / cm2},
    {0.5 * MeV, 2.49e-1 * g / cm2},
    {1.0 * MeV, 6.06e-1 * g / cm2}}},
  {"G4_Cu",
   {{0.3 * MeV, 1.25e-1 * g / cm2}, {0.5 * MeV, 2.58e-1 * g / cm2}, {1.0 * MeV, 6.25e-1 * g / cm2}}}

};

G4double getDensity(const G4String& material)
{
  if (density.count(material) == 1) return density.at(material);
  G4cout << "Try get undefined material " << material << " from density... Return 0.0" << G4endl;
  return 0.;
}

G4double getMFP(const G4String& material, const G4double energy, const G4double eps = 0.001 * MeV)
{
  if (mfp.count(material) == 0)
  {
    G4cout << "Try get undefined material " << material << " from mfp... Return 1e10" << G4endl;
    return 1e+10;
  }
  for (auto& list : mfp.at(material))
    if (std::abs(energy - list.first) < eps) return list.second;
  G4cout << "energy " << energy << " not found" << G4endl;
  throw std::runtime_error("energy " + std::to_string(energy) + " not found for material "
                           + material + ".");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // get particle energy
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();
  G4cout << "### Run " << aRun->GetRunID() << " start; E(MeV)= " << energy / MeV
         << "   diff= " << energy - 1.033 * MeV << G4endl;

  // get name of material 1, 2 & 3
  matName1 = detector->GetAbsorber1Material()->GetName();
  matName2 = detector->GetAbsorber2Material()->GetName();
  matName3 = detector->GetAbsorber3Material()->GetName();

  MFP1 = 1.e+10;
  MFP2 = 1.e+10;
  MFP3 = 1.e+10;

  density1 = 0.0;
  density2 = 0.0;
  density3 = 0.0;

  // the CSDA range values (i.e. MFPx) were taken from the SANDIA REPORT
  // (G.J.Lockwood et al. Sandia report SAND79-0414.UC-34a, Febr. 1987) for Ta,
  // U, Mo, Al, Be and U while they were computed in case of Si at the
  // corresponding experimental e- primary energies (U. Werner et al. J. Phys.
  // D: Appl. Phys. 21 (1988) 116)

  //
  MFP1 = getMFP(matName1, energy);
  MFP2 = getMFP(matName2, energy);
  MFP3 = getMFP(matName3, energy);

  density1 = getDensity(matName1);
  density2 = getDensity(matName2);
  density3 = getDensity(matName3);

  // initialize EnergyDeposit per Layer
  for (G4int k = 0; k < 110; k++)
  {
    energyDeposit1[k] = 0.0;
    energyDeposit2[k] = 0.0;
    energyDeposit3[k] = 0.0;
  }

  // initialize EnergyDeposit per Medium
  energyDepositRun1 = 0.;
  energyDepositRun2 = 0.;
  energyDepositRun3 = 0.;

  // counters
  n_steps = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4double NumbrOfEvents = G4double(aRun->GetNumberOfEvent());
  G4double totalAbsEnergy =
    (energyDepositRun1 / MeV) + (energyDepositRun2 / MeV) + (energyDepositRun3 / MeV);

  G4int n12 = detector->GetNbOfLayersOfMedium1();
  G4double thick11 = detector->GetAbsorber1Thickness();
  G4double thick12 = G4double(n12);
  G4double LayerTh1 = thick11 / thick12;
  G4double deltaX1 = LayerTh1 * density1 / MFP1;
  // G4cout << "#1 n= " << n12 << " density= " << density1 << " MFP= " << MFP1 << G4endl;

  G4int n22 = detector->GetNbOfLayersOfMedium2();
  G4double thick21 = detector->GetAbsorber2Thickness();
  G4double thick22 = G4double(n22);
  G4double LayerTh2 = thick21 / thick22;
  G4double deltaX2 = LayerTh2 * density2 / MFP2;
  // G4cout << "#2 n= " << n22 << " density= " << density2 << " MFP= " << MFP2 << G4endl;

  G4int n32 = detector->GetNbOfLayersOfMedium3();
  G4double thick31 = detector->GetAbsorber3Thickness();
  G4double thick32 = G4double(n32);
  G4double LayerTh3 = thick31 / thick32;
  G4double deltaX3 = LayerTh3 * density3 / MFP3;
  // G4cout << "#3 n= " << n32 << " density= " << density3 << " MFP= " << MFP3 << G4endl;

  G4cout << "Thicknesses(mm)= " << thick21 / mm << "  " << thick22 / mm << "  " << thick32 / mm
         << "  " << G4endl;
  G4cout << "Bins(R/R0)=      " << deltaX1 << "  " << deltaX2 << "  " << deltaX3 << "  " << G4endl;

  G4cout << " ----------------------------------------------------------" << G4endl;
  G4cout << " ----------------  RUN SUMMARY ----------------------------" << G4endl;
  G4cout << " ----------------------------------------------------------" << G4endl;
  G4cout << "       Mean number of steps " << G4double(n_steps) / NumbrOfEvents << G4endl;

  asciiFileName = "Sandia.out";
  std::ofstream asciiFile(asciiFileName);
  if (asciiFile.is_open())
  {
    asciiFile << " FMR(z/r0)      ||       J(MeV/g/cm2)" << G4endl;
  }
  else
  {
    G4cout << "ERROR file <" << asciiFileName << "> is not opened" << G4endl;
    return;
  }

  G4double pos = 0.0;
  G4double norm = g / (NumbrOfEvents * deltaX1 * MFP1 * MeV * cm2);

  if (matName1 != "G4_Galactic")
  {
    asciiFile << "             Medium 1 ==>   " << matName1 << "  " << n12 << G4endl;
    G4cout << "             Medium 1 ==>   " << matName1 << G4endl;
    G4cout << "\t FMR (z/r0)    ||       J  " << G4endl;
    // Write to file
    for (G4int k = 1; k <= n12; k++)
    {
      pos += deltaX1;
      normalizedvalue1[k] = energyDeposit1[k] * norm;
      G4cout << "\t   " << pos << "\t          " << normalizedvalue1[k] << G4endl;
      asciiFile << std::setiosflags(std::ios::fixed) << std::setprecision(5)
                << std::setiosflags(std::ios::right) << std::setw(10);
      asciiFile << pos;
      asciiFile << "           ";
      asciiFile << std::setiosflags(std::ios::fixed) << std::setprecision(5)
                << std::setiosflags(std::ios::right) << std::setw(10);
      asciiFile << normalizedvalue1[k] << G4endl;
    }
    G4cout << "\n Deposit energy (MeV) = " << (energyDepositRun1 / MeV) / NumbrOfEvents << G4endl;
  }

  G4cout << " ----------------------------------------------------------" << G4endl;
  if (matName2 != "G4_Galactic")
  {
    G4cout << "               Medium 2 ==>   " << matName2 << G4endl;
    G4cout << "\t FMR (z/r0)    ||       J  " << G4endl;
    asciiFile << "       Medium 2 ==>   " << matName2 << "  " << n22 << G4endl;
    norm = g / (NumbrOfEvents * deltaX2 * MFP2 * MeV * cm2);
    for (G4int k = 1; k <= n22; k++)
    {
      pos += deltaX2;
      normalizedvalue2[k] = energyDeposit2[k] * norm;
      G4cout << "\t   " << pos << "\t          " << normalizedvalue2[k] << G4endl;
      asciiFile << std::setiosflags(std::ios::fixed) << std::setprecision(5)
                << std::setiosflags(std::ios::right) << std::setw(10);
      asciiFile << pos;
      asciiFile << "           ";
      asciiFile << std::setiosflags(std::ios::fixed) << std::setprecision(5)
                << std::setiosflags(std::ios::right) << std::setw(10);
      asciiFile << normalizedvalue2[k] << G4endl;
    }
    G4cout << "\n Deposit energy (MeV) = " << (energyDepositRun2 / MeV) / NumbrOfEvents << G4endl;
  }

  G4cout << " ----------------------------------------------------------" << G4endl;

  if (matName3 != "G4_Galactic")
  {
    G4cout << "               Medium 3 ==>   " << matName3 << G4endl;
    G4cout << "\t FMR (z/r0)    ||       J  " << G4endl;
    asciiFile << "       Medium 3 ==>   " << matName3 << "  " << n32 << G4endl;
    norm = g / (NumbrOfEvents * deltaX3 * MFP3 * MeV * cm2);
    for (G4int k = 1; k <= n32; k++)
    {
      pos += deltaX3;
      normalizedvalue3[k] = energyDeposit3[k] * norm;
      G4cout << "\t   " << pos << "\t          " << normalizedvalue3[k] << G4endl;
      asciiFile << std::setiosflags(std::ios::fixed) << std::setprecision(5)
                << std::setiosflags(std::ios::right) << std::setw(10);
      asciiFile << pos;
      asciiFile << "           ";
      asciiFile << std::setiosflags(std::ios::fixed) << std::setprecision(5)
                << std::setiosflags(std::ios::right) << std::setw(10);
      asciiFile << normalizedvalue3[k] << G4endl;
    }
    G4cout << "\n Deposit energy (MeV) = " << (energyDepositRun3 / MeV) / (NumbrOfEvents * 1.0)
           << G4endl;
  }
  asciiFile.close();

  G4cout << " ----------------------------------------------------------" << G4endl;

  G4cout << " Total Deposit energy (MeV) = " << (totalAbsEnergy / MeV) / (NumbrOfEvents * 1.0)
         << G4endl;

  G4cout << " ----------------------------------------------------------" << G4endl;
  G4cout << " ----------------  END OF RUN SUMMARY ---------------------" << G4endl;
  G4cout << " ----------------------------------------------------------" << G4endl;
}
