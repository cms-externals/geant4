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
/// \file electromagnetic/TestEm11/src/ProcCounterAccumulable.cc
/// \brief Implementation of the ProcCounterAccumulable class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ProcCounterAccumulable.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ProcCounterAccumulable::ProcCounterAccumulable(const G4String& name)
  : G4VAccumulable(name), fProcCounter()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ProcCounterAccumulable::~ProcCounterAccumulable() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ProcCounterAccumulable::CountProcesses(G4String procName)
{
  std::map<G4String, G4int>::iterator it = fProcCounter.find(procName);
  if (it == fProcCounter.end())
  {
    fProcCounter[procName] = 1;
  }
  else
  {
    fProcCounter[procName]++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ProcCounterAccumulable::Merge(const G4VAccumulable& other)
{
  const ProcCounterAccumulable& otherProcCounterAccumulable =
    static_cast<const ProcCounterAccumulable&>(other);

  // map: processes count
  std::map<G4String, G4int>::const_iterator it;
  for (it = otherProcCounterAccumulable.fProcCounter.begin();
       it != otherProcCounterAccumulable.fProcCounter.end(); ++it)
  {
    G4String procName = it->first;
    G4int otherCount = it->second;
    if (fProcCounter.find(procName) == fProcCounter.end())
    {
      fProcCounter[procName] = otherCount;
    }
    else
    {
      fProcCounter[procName] += otherCount;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ProcCounterAccumulable::Reset()
{
  G4cout << "... Clearing procCounter map" << G4endl;
  fProcCounter.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ProcCounterAccumulable::Print(G4PrintOptions /*options*/) const
{
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;
  std::map<G4String, G4int>::const_iterator it;
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++)
  {
    G4String procName = it->first;
    G4int count = it->second;
    G4String space = " ";
    if (++index % 3 == 0) space = "\n";
    G4cout << " " << std::setw(20) << procName << "=" << std::setw(7) << count << space;
  }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
