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
/// \file radiobiology/src/Manager.cc
/// \brief Implementation of the RadioBio::Manager class

#include "Manager.hh"

#include "G4AccumulableManager.hh"

#include <map>
#include <mutex>

namespace
{
std::mutex initializationMutex;
}

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Manager* Manager::CreateInstance()
{
  if (fInstance != nullptr) {
    G4Exception("RadioBioManager::createInstance", "RecreatingRadioBioManager", FatalException,
                "Creating another, new, instance of RadioBioManager");
    return fInstance;
  }

  fInstance = new Manager();
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Manager* Manager::GetInstance()
{
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor deletes all the quantities

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Manager::~Manager()
{
  for (const auto& q : fQuantities) {
    delete q.second;
  }
  fInstance = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::InitializeAll()
{
  std::lock_guard<std::mutex> lock(initializationMutex);
  for (const auto& q : fQuantities) {
    (q.second)->Initialize();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::ComputeAll()
{
  for (const auto& q : fQuantities) {
    (q.second)->Compute();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::ResetAll()
{
  for (const auto& q : fQuantities) {
    (q.second)->Reset();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::StoreAll()
{
  for (const auto& q : fQuantities) {
    (q.second)->Store();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VRadiobiologicalQuantity* Manager::GetQuantity(const G4String& name)
{
  const auto quantity = fQuantities.find(name);
  return quantity != fQuantities.end() ? quantity->second : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::PrintParameters()
{
  G4cout << "*******************************************" << G4endl
         << "*** right now registered quantities are ***" << G4endl;
  for (const auto& q : fQuantities) {
    G4cout << "*** " << q.first;
  }
  G4cout << "*** but their calculation might be not ****" << G4endl
         << "*** active. Ask for parameters of each ****" << G4endl
         << "*******************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::DigestAccumulables()
{
  for (const auto& q : fQuantities) {
    // Disabled quantities have no data to digest.
    if (!q.second->IsCalculationEnabled()) continue;

    // Hook in the accumulable manager the one named as the quantity (eg "Dose")
    G4VAccumulable* accumulable = G4AccumulableManager::Instance()->GetAccumulable(q.first);

    if (accumulable == nullptr) {
      G4Exception("RadioBioManager::AddFromAccumulable", "NoAccumulable", FatalException, q.first);
      continue;
    }

    // Add from the accumulable.
    q.second->AddFromAccumulable(accumulable);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Manager::Register(VRadiobiologicalQuantity* q, const G4String& name)
{
  if (q == nullptr) {
    G4Exception("RadioBioManager::Register", "RegisteringNullptr", JustWarning,
                "Asking to register a quantity with null pointer!");
    return false;
  }

  if (fQuantities.find(name) != fQuantities.end()) {
    G4Exception("RadioBioManager::Register", "RegisteringSameQuantity", FatalException,
                "Registering two radiobiological quantities with the same name!");
    return false;
  }
  fQuantities[name] = q;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
