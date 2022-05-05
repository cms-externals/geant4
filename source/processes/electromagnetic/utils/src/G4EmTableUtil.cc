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
// Geant4 class G4EmTableUtil
//
// Author V.Ivanchenko 14.03.2022
//

#include "G4EmTableUtil.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4Electron.hh"
#include "G4UIcommand.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::PrepareMscProcess(G4VMultipleScattering* proc,
				      const G4ParticleDefinition& part,
				      G4EmModelManager* modelManager,
				      G4MscStepLimitType& stepLimit,
                                      G4double& facrange,
				      G4bool& latDisplacement, G4bool& master,
				      G4bool& isIon, G4bool& baseMat)
{
  auto param = G4EmParameters::Instance();
  G4int verb = (master) ? param->Verbose() : param->WorkerVerbose(); 
  proc->SetVerboseLevel(verb);

  if(part.GetPDGMass() > CLHEP::GeV ||
     part.GetParticleName() == "GenericIon") { isIon = true; }

  if(1 < verb) {
    G4cout << "### G4VMultipleScattering::PrepearPhysicsTable() for "
           << proc->GetProcessName()
           << " and particle " << part.GetParticleName()
           << " isIon: " << isIon << " isMaster: " << master
	   << G4endl;
  }

  // initialise process
  proc->InitialiseProcess(&part);

  // heavy particles 
  if(part.GetPDGMass() > CLHEP::MeV) {
    stepLimit = param->MscMuHadStepLimitType(); 
    facrange = param->MscMuHadRangeFactor(); 
    latDisplacement = param->MuHadLateralDisplacement();
  } else {
    stepLimit = param->MscStepLimitType(); 
    facrange = param->MscRangeFactor(); 
    latDisplacement = param->LateralDisplacement();
  }

  // initialisation of models
  auto numberOfModels = modelManager->NumberOfModels();
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VMscModel* msc = proc->GetModelByIndex(i);
    msc->SetIonisation(nullptr, &part);
    msc->SetMasterThread(master);
    msc->SetPolarAngleLimit(param->MscThetaLimit());
    G4double emax = std::min(msc->HighEnergyLimit(),param->MaxKinEnergy());
    msc->SetHighEnergyLimit(emax);
    msc->SetUseBaseMaterials(baseMat);
  }
  modelManager->Initialise(&part, G4Electron::Electron(), verb);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmTableUtil::BuildMscProcess(G4VMultipleScattering* proc,
                                    const G4VMultipleScattering* masterProc,
		                    const G4ParticleDefinition& part,
		                    const G4ParticleDefinition* firstPart,
				    G4int nModels, 
                                    G4bool& master, G4bool& baseMat)
{
  auto param = G4EmParameters::Instance();
  G4int verb = param->Verbose(); 

  if(!master && firstPart == &part) {
    baseMat = masterProc->UseBaseMaterial();
    // initialisation of models
    for(G4int i=0; i<nModels; ++i) {
      G4VMscModel* msc = proc->GetModelByIndex(i);
      G4VMscModel* msc0 = masterProc->GetModelByIndex(i);
      msc->SetUseBaseMaterials(baseMat);
      msc->SetCrossSectionTable(msc0->GetCrossSectionTable(), false);
      msc->InitialiseLocal(&part, msc0);
    }
  }
  if(!param->IsPrintLocked()) {
    const G4String& num = part.GetParticleName();

    // explicitly defined printout by particle name
    if(1 < verb || (0 < verb && (num == "e-" || 
		 		 num == "e+"    || num == "mu+" || 
				 num == "mu-"   || num == "proton"|| 
				 num == "pi+"   || num == "pi-" || 
				 num == "kaon+" || num == "kaon-" || 
				 num == "alpha" || num == "anti_proton" || 
				 num == "GenericIon" || num == "alpha+" || 
				 num == "alpha" ))) { 
      proc->StreamInfo(G4cout, part);
    }
  }
  if(1 < verb) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() done for "
	   << proc->GetProcessName()
	   << " and particle " << part.GetParticleName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmTableUtil::StoreMscTable(G4VMultipleScattering* proc,
                                    const G4ParticleDefinition* part,
                                    const G4String& dir,
                                    const G4int nModels, const G4int verb,
		                    const G4bool ascii)
{
  G4bool ok = true;
  for(G4int i=0; i<nModels; ++i) {
    G4VMscModel* msc = proc->GetModelByIndex(i);
    G4PhysicsTable* table = msc->GetCrossSectionTable();
    if (nullptr != table) {
      G4String ss = G4UIcommand::ConvertToString(i);
      G4String name = 
        proc->GetPhysicsTableFileName(part, dir, "LambdaMod"+ss, ascii);
      G4bool yes = table->StorePhysicsTable(name,ascii);

      if ( yes ) {
        if ( verb > 0 ) {
          G4cout << "Physics table are stored for " 
                 << part->GetParticleName()
                 << " and process " << proc->GetProcessName()
                 << " with a name <" << name << "> " << G4endl;
        }
      } else {
        G4cout << "Fail to store Physics Table for " 
               << part->GetParticleName()
               << " and process " << proc->GetProcessName()
               << " in the directory <" << dir
               << "> " << G4endl;
	ok = false;
      }
    }
  }
  return ok;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmTableUtil::StoreTable(G4VProcess* ptr,
                                 const G4ParticleDefinition* part, 
                                 G4PhysicsTable* aTable, 
                                 const G4String& dir,
                                 const G4String& tname,
                                 const G4int verb, const G4bool ascii)
{
  G4bool res = true;
  if (nullptr != aTable) {
    const G4String& name = 
      ptr->GetPhysicsTableFileName(part, dir, tname, ascii);
    if ( aTable->StorePhysicsTable(name, ascii) ) {
      if (0 < verb) G4cout << "Stored: " << name << G4endl;
    } else {
      res = false;
      G4cout << "Fail to store: " << name << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

