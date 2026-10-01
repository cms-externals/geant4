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
/// \file QuasiOpticalPhysics.cc
/// \brief QuasiOpticalPhysics constructor
///
/// This class is equivalent to G4OpticalPhysics, with extensions for
/// quasi-optical processes to support offloading of optical photon
/// generation and transport.

#include "QuasiOpticalPhysics.hh"

#include "G4GeneralCerenkov.hh"
#include "G4LossTableManager.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4OpWLS.hh"
#include "G4OpWLS2.hh"
#include "G4OpticalParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4QuasiCerenkov.hh"
#include "G4QuasiOpticalPhoton.hh"
#include "G4QuasiScintillation.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(QuasiOpticalPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
QuasiOpticalPhysics::QuasiOpticalPhysics(G4int verbose, const G4String& name)
  : G4VPhysicsConstructor(name)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetVerboseLevel(verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QuasiOpticalPhysics::PrintStatistics() const
{
  G4OpticalParameters::Instance()->Dump();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QuasiOpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
  G4QuasiOpticalPhoton::QuasiOpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QuasiOpticalPhysics::ConstructProcess()
{
  if (verboseLevel > 0) 
  {
    G4cout << "QuasiOpticalPhysics: Add Optical Physics Processes" << G4endl;
  }
  auto params = G4OpticalParameters::Instance();

  // Add Optical Processes

  G4ProcessManager* pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  if (nullptr == pManager)
  {
    G4ExceptionDescription ed;
    ed << "Optical Photon without a Process Manager";
    G4Exception("QuasiOpticalPhysics::ConstructProcess()", "",
                FatalException, ed);
    return;
  }

  if (params->GetProcessActivation("OpAbsorption"))
  {
    auto absorption = new G4OpAbsorption();
    pManager->AddDiscreteProcess(absorption);
  }

  if (params->GetProcessActivation("OpRayleigh"))
  {
    auto rayleigh = new G4OpRayleigh();
    pManager->AddDiscreteProcess(rayleigh);
  }

  if (params->GetProcessActivation("OpMieHG"))
  {
    auto mie = new G4OpMieHG();
    pManager->AddDiscreteProcess(mie);
  }

  if (params->GetProcessActivation("OpBoundary"))
  {
    auto boundary = new G4OpBoundaryProcess();
    pManager->AddDiscreteProcess(boundary);
  }

  if (params->GetProcessActivation("OpWLS"))
  {
    auto wls = new G4OpWLS();
    pManager->AddDiscreteProcess(wls);
  }

  if (params->GetProcessActivation("OpWLS2"))
  {
    auto wls2 = new G4OpWLS2();
    pManager->AddDiscreteProcess(wls2);
  }

  // Processes for generating quasi-optical photons

  G4VProcess* theCerenkov{nullptr};
  if (params->GetProcessActivation("QuasiCerenkov"))
  {
    if(params->GetCerenkovOffloadPhotons() &&
       params->GetCerenkovStackPhotons())
    {
      G4Exception("QuasiOpticalPhysics::ConstructProcess()", 
                  "InvalidQuasiCerenkovConfiguration", JustWarning,
                  "Cannot stack Cerenkov photons during offloading: \n"
                  "Cerenkov photon stacking will be disabled.");
    }
    theCerenkov = new G4QuasiCerenkov();
  }

  G4VProcess* theScint{nullptr};
  if (params->GetProcessActivation("QuasiScintillation"))
  {
    if(params->GetScintOffloadPhotons() && params->GetScintStackPhotons())
    {
      G4Exception("QuasiOpticalPhysics::ConstructProcess()",
                  "InvalidQuasiScintillationConfiguration", JustWarning,
                  "Cannot stack scintillation photons during offloading: \n"
                  "Scintillation photon stacking will be disabled.");
    }
    auto scint = new G4QuasiScintillation();
    scint->AddSaturation(G4LossTableManager::Instance()->EmSaturation());
    theScint = scint;
  }
 
  auto myParticleIterator = GetParticleIterator();
  myParticleIterator->reset();

  while ((*myParticleIterator)())
  {
    auto particle = myParticleIterator->value();
    if (particle->IsShortLived())
    {
      continue;
    }

    pManager = particle->GetProcessManager();
    if (nullptr == pManager)
    {
      G4ExceptionDescription ed;
      ed << "Particle " << particle->GetParticleName() 
         << "without a Process Manager";
      G4Exception("QuasiOpticalPhysics::ConstructProcess()", "",
                  FatalException, ed);
      return; 
    }

    if (nullptr != theCerenkov && theCerenkov->IsApplicable(*particle))
    {
      pManager->AddDiscreteProcess(theCerenkov);
    }

    if (nullptr != theScint && theScint->IsApplicable(*particle))
    {
      pManager->AddProcess(theScint);
      pManager->SetProcessOrderingToLast(theScint, idxAtRest);
      pManager->SetProcessOrderingToLast(theScint, idxPostStep);
    }
  }

  if (verboseLevel > 0) 
  { 
    G4cout << "### " << namePhysics << " physics constructed." << G4endl;
  }
  if (verboseLevel > 1) 
  {
    PrintStatistics();
  }
}
