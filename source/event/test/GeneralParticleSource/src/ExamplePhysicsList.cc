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
// GSPMPhysicsList - instantiates all particles

#include "globals.hh"
#include "ExamplePhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
//#include <iomanip.h>                


ExamplePhysicsList::ExamplePhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

ExamplePhysicsList::~ExamplePhysicsList()
{
}

void ExamplePhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructIons();
}

void ExamplePhysicsList::ConstructBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExamplePhysicsList::ConstructLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExamplePhysicsList::ConstructMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExamplePhysicsList::ConstructBaryons()
{
  //  Construct all baryons
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExamplePhysicsList::ConstructIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}


void ExamplePhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
void ExamplePhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
     pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
     
     pmanager->AddProcess(new G4eIonisation(),-1,2,2);
     pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
     pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
     pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    } else { 
      if ((particle->GetPDGCharge() != 0.0) && 
          (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
       pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
     }
    }
  }
}

void ExamplePhysicsList::ConstructHad()
{;}

void ExamplePhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void ExamplePhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager->AddProcess(theDecayProcess, INT_MAX, -1, INT_MAX); 
   }
  }
}

void ExamplePhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "ExamplePhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  

  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
 


