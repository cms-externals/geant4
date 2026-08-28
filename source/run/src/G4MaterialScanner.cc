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
// G4MaterialScanner implementation
//
// Author: M.Asai, May 2006
// --------------------------------------------------------------------

#include "G4MaterialScanner.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4GeometryManager.hh"
#include "G4MSSteppingAction.hh"
#include "G4MatScanMessenger.hh"
#include "G4NistManager.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RayShooter.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManagerKernel.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"

// --------------------------------------------------------------------
G4MaterialScanner::G4MaterialScanner()
{
  theRayShooter = new G4RayShooter();
  theMessenger = new G4MatScanMessenger(this);
  theEventManager = G4EventManager::GetEventManager();

  eyePosition = G4ThreeVector(0., 0., 0.);
  nTheta = 81;
  thetaMin = -80.0 * deg;
  thetaSpan = 160.0 * deg;
  nPhi = 36;
  phiMin = 0.0 * deg;
  phiSpan = 350.0 * deg;

  ost = &G4cout;
}

// --------------------------------------------------------------------
G4MaterialScanner::~G4MaterialScanner()
{
  delete theRayShooter;
  delete theMatScannerSteppingAction;
  delete theMessenger;
  if (ost != &G4cout)
  {
    delete ost;
  }
}

// --------------------------------------------------------------------
void G4MaterialScanner::Scan()
{
  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  G4ApplicationState currentState = theStateMan->GetCurrentState();
  if (currentState != G4State_Idle)
  {
    G4cerr << "Illegal application state - Scan() ignored." << G4endl;
    return;
  }

  if (theMatScannerSteppingAction == nullptr)
  {
    theMatScannerSteppingAction = new G4MSSteppingAction();
  }
  StoreUserActions();
  DoScan();
  RestoreUserActions();
}

// --------------------------------------------------------------------
void G4MaterialScanner::StoreUserActions()
{
  theUserEventAction = theEventManager->GetUserEventAction();
  theUserStackingAction = theEventManager->GetUserStackingAction();
  theUserTrackingAction = theEventManager->GetUserTrackingAction();
  theUserSteppingAction = theEventManager->GetUserSteppingAction();

  theEventManager->SetUserAction(theMatScannerEventAction);
  theEventManager->SetUserAction(theMatScannerStackingAction);
  theEventManager->SetUserAction(theMatScannerTrackingAction);
  theEventManager->SetUserAction(theMatScannerSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if (theSDMan != nullptr)
  {
    theSDMan->Activate("/", false);
  }
}

// --------------------------------------------------------------------
void G4MaterialScanner::RestoreUserActions()
{
  theEventManager->SetUserAction(theUserEventAction);
  theEventManager->SetUserAction(theUserStackingAction);
  theEventManager->SetUserAction(theUserTrackingAction);
  theEventManager->SetUserAction(theUserSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if (theSDMan != nullptr)
  {
    theSDMan->Activate("/", true);
  }
}

// --------------------------------------------------------------------
void G4MaterialScanner::DoScan()
{
  // Disable visualization
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/disable");
  // Confirm geometry is optimized and material table is updated
  G4UImanager::GetUIpointer()->ApplyCommand("/run/undertakeOptimisation");
  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  theStateMan->SetNewState(G4State_Init);
  G4RunManagerKernel::GetRunManagerKernel()->UpdateRegion();
  theStateMan->SetNewState(G4State_Idle);

  theMatScannerSteppingAction->SetThinMaterial(ignoreThinMaterial, thinDencity);

  G4ThreeVector center(0, 0, 0);
  G4Navigator* navigator =
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  navigator->LocateGlobalPointAndSetup(center, nullptr, false);

  theStateMan->SetNewState(G4State_GeomClosed);

  std::ios::fmtflags os_flags(G4cout.flags());
  // Event loop
  G4int iEvent = 0;
  for (G4int iTheta = 0; iTheta < nTheta; ++iTheta)
  {
    G4double theta = thetaMin;
    if (iTheta > 0) theta += G4double(iTheta) * thetaSpan / G4double(nTheta - 1);
    G4double aveLength = 0.;
    G4double aveMassThick = 0.;
    G4double aveX0 = 0.;
    G4double aveLambda = 0.;
    *(ost) << G4endl;
    *(ost) << "              Theta         Phi      Length  Mass length         x0     lambda0"
           << G4endl;
    *(ost) << "              (deg)       (deg)        (cm)      (g/cm2)" << G4endl;
    *(ost) << G4endl;
    for (G4int iPhi = 0; iPhi < nPhi; ++iPhi)
    {
      auto anEvent = new G4Event(iEvent++);
      G4double phi = phiMin;
      if (iPhi > 0) phi += G4double(iPhi) * phiSpan / G4double(nPhi - 1);
      eyeDirection = G4ThreeVector(std::cos(theta) * std::cos(phi), std::cos(theta) * std::sin(phi),
                                   std::sin(theta));
      theRayShooter->Shoot(anEvent, eyePosition, eyeDirection);
      theMatScannerSteppingAction->Initialize(regionSensitive, theRegion);
      theEventManager->ProcessOneEvent(anEvent);
      G4double length = theMatScannerSteppingAction->GetTotalStepLength();
      G4double massThick = theMatScannerSteppingAction->GetTotalMassThickness();
      G4double x0 = theMatScannerSteppingAction->GetX0();
      G4double lambda = theMatScannerSteppingAction->GetLambda0();

      *(ost) << "        " << std::setprecision(2) << std::setw(11) << std::scientific
             << theta / deg << " " << std::setw(11) << std::scientific << phi / deg << " "
             << std::setw(11) << std::scientific << length / cm << " " << std::setw(11)
             << std::scientific << massThick / (g / cm2) << " " << std::setw(11) << std::scientific
             << x0 << " " << std::setw(11) << std::scientific << lambda << G4endl;
      if (1 == verbosity)
      {
        theMatScannerSteppingAction->PrintIntegratedMaterialVerbose(*(ost));
      }
      else if (2 == verbosity)
      {
        theMatScannerSteppingAction->PrintEachMaterialVerbose(*(ost));
      }
      *(ost) << G4endl;
      delete anEvent;
      aveLength += length / cm;
      aveMassThick += massThick / (g / cm2);
      aveX0 += x0;
      aveLambda += lambda;
    }
    if (nPhi > 1)
    {
      *(ost) << G4endl;
      *(ost) << " ave. for theta = " << std::setprecision(2) << std::setw(11) << std::scientific
             << theta / deg << " : " << std::setw(11) << std::scientific << aveLength / nPhi << " "
             << std::setw(11) << std::scientific << aveMassThick / nPhi << " " << std::setw(11)
             << std::scientific << aveX0 / nPhi << " " << std::setw(11) << std::scientific
             << aveLambda / nPhi << G4endl;
    }
  }
  G4cout.flags(os_flags);

  theStateMan->SetNewState(G4State_Idle);
  // Enable visualization
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/enable");
  return;
}

// --------------------------------------------------------------------
G4bool G4MaterialScanner::SetRegionName(const G4String& val)
{
  G4Region* aRegion = G4RegionStore::GetInstance()->GetRegion(val);
  if (aRegion != nullptr)
  {
    theRegion = aRegion;
    regionName = val;
    return true;
  }
  return false;
}

// --------------------------------------------------------------------
G4bool G4MaterialScanner::SetThinMaterial(G4bool flg, G4String& mat)
{
  if (flg)
  {
    auto pMat = G4Material::GetMaterial(mat, false);
    if (!pMat) pMat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
    if (pMat)
    {
      ignoreThinMaterial = true;
      thinDencity = pMat->GetDensity();
      return true;
    }
    else
    {
      return false;
    }
  }

  ignoreThinMaterial = false;
  return true;
}

// --------------------------------------------------------------------
void G4MaterialScanner::SetDetDefaultParticle(const G4ParticleDefinition* ptcl, G4double eKin)
{
  if (theMatScannerSteppingAction == nullptr)
  {
    theMatScannerSteppingAction = new G4MSSteppingAction();
  }
  theMatScannerSteppingAction->SetDetDefaultParticle(ptcl, eKin);
}

// --------------------------------------------------------------------
void G4MaterialScanner::SetOutFile(G4String& fNm)
{
  if (fNm != fTarget)
  {
    if (fNm == "**Screen**")
    {
      ost = &G4cout;
    }
    else
    {
      if (ost != &G4cout) delete ost;
      ost = new std::ofstream(fNm);
    }
    fTarget = fNm;
  }
}
