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
// $Id: G4NuclearPolarizationStore.cc 97302 2016-06-01 09:30:11Z gcosmo $
//
// 23-Jan-2009 V.Ivanchenko make the class to be a singleton
// 17-Aug-2012 V.Ivanchenko added hadronic model factories

#include "G4NuclearPolarizationStore.hh"
#include "G4SystemOfUnits.hh"

G4ThreadLocal G4NuclearPolarizationStore* 
G4NuclearPolarizationStore::instance = nullptr;

G4NuclearPolarizationStore* G4NuclearPolarizationStore::GetInstance()
{
  if(nullptr == instance) {
    static G4ThreadLocalSingleton<G4NuclearPolarizationStore> inst;
    instance = inst.Instance();
  }
  return instance;
}

G4NuclearPolarizationStore::G4NuclearPolarizationStore()
{
  nuclist.reserve(10);
}

G4NuclearPolarizationStore::~G4NuclearPolarizationStore()
{
  //G4cout << "G4NuclearPolarizationStore::~G4NuclearPolarizationStore() "
  //	 << nuclist.size() << G4endl;
  for(auto nucp : nuclist) {
    G4cout << nucp << G4endl;
    delete nucp;
  }
  nuclist.clear();
}

void G4NuclearPolarizationStore::Register(G4NuclearPolarization* ptr)
{
  for (auto nucp : nuclist) {
    if(ptr == nucp) { return; }
  }
  for (auto nucp : nuclist) {
    if(nullptr == nucp) { 
      nucp = ptr;
      return; 
    }
  }
  nuclist.push_back(ptr);
  //G4cout <<"G4NuclearPolarizationStore::Register() "<< ptr <<G4endl; 
}

G4NuclearPolarization* 
G4NuclearPolarizationStore::FindOrBuild(G4int Z, G4int A, G4double Eexc)
{
  static const G4double tolerance = 10.*CLHEP::eV;
  for (auto nucp : nuclist) {
    if(nucp && Z == nucp->GetZ() && A == nucp->GetA() && 
       std::abs(Eexc - nucp->GetExcitationEnergy()) < tolerance) { 
      return nucp; 
    }
  }
  G4NuclearPolarization* ptr = new G4NuclearPolarization(Z, A, Eexc);
  Register(ptr);
  return ptr;   
}

void G4NuclearPolarizationStore::RemoveMe(G4NuclearPolarization* ptr)
{
  G4int length = nuclist.size();
  for (G4int i=0; i<length; ++i) {
    if(ptr == nuclist[i]) { 
      //G4cout <<"G4NuclearPolarizationStore::RemoveMe() "<< ptr <<G4endl; 
      delete ptr;
      if(i+1 == length) { nuclist.pop_back(); }
      else { nuclist[i] = nullptr; }
      return;
    }
  }
}
