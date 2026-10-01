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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel and Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
// ********************************************************************
//
/// \file PhysicsConfigurator.cc
/// \brief Implementation of the CaTS::PhysicsConfigurator class

// Geant4 headers
#include "G4String.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PhysListFactoryAlt.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4PhysListRegistry.hh"
#include "G4OpticalParameters.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StepLimiter.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"
// project Headers
#include "PhysicsConfigurator.hh"
#include "ConfigurationManager.hh"
// c++ headers
#include <stdlib.h> 
PhysicsConfigurator* PhysicsConfigurator::instance = 0;

G4VModularPhysicsList* PhysicsConfigurator::Construct(G4String physName)
{
  //
  // Access to registries and factories
  //
  G4PhysicsConstructorRegistry* g4pcr =
    G4PhysicsConstructorRegistry::Instance();
  G4PhysListRegistry* g4plr = G4PhysListRegistry::Instance();
  G4bool verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  if(verbose)
  {
    G4cout << "Available Physics Constructors:  "
           << g4pcr->AvailablePhysicsConstructors().size() << G4endl;
    G4cout << "Available Physics Lists:         "
           << g4plr->AvailablePhysLists().size() << G4endl;
    G4cout << "Available Physics Extensions:    "
           << g4plr->AvailablePhysicsExtensions().size() << G4endl;
    G4cout << "Available Physics Lists Em:      "
           << g4plr->AvailablePhysListsEM().size() << G4endl;
    g4plr->SetVerbose(1);
  }
  else
  {
    g4plr->SetVerbose(0);
  }
  g4plr->AddPhysicsExtension("OPTICAL", "G4OpticalPhysics");
  g4plr->AddPhysicsExtension("QUASIOPTICAL", "QuasiOpticalPhysics");
  g4plr->AddPhysicsExtension("STEPLIMIT", "G4StepLimiterPhysics");
  g4plr->AddPhysicsExtension("NEUTRONLIMIT", "G4NeutronTrackingCut");

  if(verbose)
  {
    g4pcr->PrintAvailablePhysicsConstructors();
    g4plr->PrintAvailablePhysLists();
  }
  g4alt::G4PhysListFactory factory;
  G4VModularPhysicsList* phys = nullptr;
  if(verbose)
    G4cout << "Physics configuration: " << physName << G4endl;
  //
  // currently using the Constructor names doesn't work otherwise it would be:
  // G4String physName = "FTFP_BERT+G4OpticalPhysics+G4StepLimiterPhysics";
  // using the name doesn't work either
  // G4String physName = "FTFP_BERT+Optical+stepLimiter";
  // reference PhysicsList via its name
  //
  if(factory.IsReferencePhysList(physName))
  {
    phys = factory.GetReferencePhysList(physName);
  }
  else
  {
    G4cout << "Not a reference physics list" << G4endl;
    g4plr->PrintAvailablePhysLists();
    exit(EXIT_FAILURE);
  }
  if(verbose)
  {
    G4cout << phys->GetPhysicsTableDirectory() << G4endl;
  }

  // Default configuration; can be overridden via the optical messenger
  auto params = G4OpticalParameters::Instance();
  params->SetProcessActivation("Cerenkov", true);
  params->SetProcessActivation("Scintillation", true);
  params->SetProcessActivation("OpAbsorption", true);
  params->SetProcessActivation("OpRayleigh", true);
  params->SetProcessActivation("OpMieHG", false);
  params->SetProcessActivation("OpWLS", true);
  params->SetProcessActivation("OpWLS2", false);

  params->SetCerenkovStackPhotons(false);
  params->SetScintStackPhotons(false);

  // only relevant if we actually stack and trace the optical photons
  params->SetCerenkovTrackSecondariesFirst(true);
  params->SetScintTrackSecondariesFirst(true);

  params->SetCerenkovMaxPhotonsPerStep(100);
  params->SetCerenkovMaxBetaChange(10.0);
  if(verbose)
  {
    phys->DumpList();
  }
  return phys;
}

PhysicsConfigurator* PhysicsConfigurator::getInstance()
{
  if(instance == 0)
  {
    instance = new PhysicsConfigurator();
  }
  return instance;
}
