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

#include "LimiterPhysics.hh"
#include "ElectronLimiter.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LimiterPhysics::LimiterPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LimiterPhysics::~LimiterPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LimiterPhysics::ConstructParticle()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LimiterPhysics::SetParameters(const G4String& part, const G4String& reg, 
				   G4double lim, G4double fac, G4double rms)
{
  if(reg != "EcalRegion" && reg != "HcalRegion") {
    G4cout << "LimiterPhysics::SetParameters WARNING: wrong region name: " 
	   << reg << G4endl;
    return;
  }
  G4int n = fPart.size();
  for(G4int i=0; i<n; ++i) {
    if(part == fPart[i] && reg == fRegion[i]) {
      G4cout << "LimiterPhysics::SetParameters overwrite parameters for "
	     << part << " inside " << reg << " : \n"
	     << "  Limit(MeV)= " << lim << " Factor= " << fac 
	     << " RMS= " << rms << G4endl;
      fLimit[i] = lim;
      fFactor[i] = fac;
      fRMS[i] = rms;
      return;
    }
  }
  G4cout << "LimiterPhysics::SetParameters for "
	 << part << " inside " << reg << " : \n"
	 << "  Limit(MeV)= " << lim << " Factor= " << fac 
	 << " RMS= " << rms << G4endl;
  fPart.push_back(part);
  fRegion.push_back(reg);
  fLimit.push_back(lim);
  fFactor.push_back(fac);
  fRMS.push_back(rms);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LimiterPhysics::ConstructProcess()
{
  G4int n = fPart.size();
  if(0 == n) { return; }
  G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
  G4RegionStore* rstore = G4RegionStore::GetInstance();

  std::vector<ElectronLimiter*> limiters;  

  for(G4int i=0; i<n; ++i) {
    const G4ParticleDefinition* particle = ptable->FindParticle(fPart[i]);
    if(!particle) continue;

    const G4Region* reg = rstore->GetRegion(fRegion[i]);
    if(!reg) continue;

    ElectronLimiter* ptr = nullptr;
    G4int nn = limiters.size();
   
    for(G4int j=0; j<nn; ++j) {
      if(particle == limiters[j]->GetParticle()) {
        ptr = limiters[j];
        break;
      }
    }
    if(!ptr) { 
      ptr = new ElectronLimiter(particle); 
      auto pManager = particle->GetProcessManager();
      pManager->AddDiscreteProcess(ptr);
    }
    if(fRegion[i] == "EcalRegion") {
      ptr->SetTrackingCutEcal(fLimit[i], fFactor[i], fRMS[i]);
    } else {
      ptr->SetTrackingCutHcal(fLimit[i], fFactor[i], fRMS[i]);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

