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
// G4VolumeTreeWalker class Implementation
//
// Author: John Apostolakis (CERN), June 2026
// --------------------------------------------------------------------

#include "G4VolumeTreeWalker.hh"

#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VVolumeVisitor.hh"

#include "geomdefs.hh"

#include <cassert>

G4VolumeTreeWalker::G4VolumeTreeWalker(G4VVolumeVisitor& visitor) : fVisitor(visitor) {}

void G4VolumeTreeWalker::Walk(G4VPhysicalVolume* startVol)
{
  assert(startVol != nullptr);
  DescendPhysical(startVol, 0);
}

// ---------------------------------------------------------------------------

void G4VolumeTreeWalker::DescendPhysical(G4VPhysicalVolume* pv, G4int depth)
{
  // Dispatch to the typed visitor method.
  // kParameterised must be tested before kReplica because
  // G4PVParameterised inherits G4PVReplica and both share kReplica/kParam.
  //
  if (fVerbose)
  {
    G4cout << "Walker::DescendPhysical depth= " << depth << " - entering." << G4endl;
  }

  switch (pv->VolumeType())
  {
    case kNormal:
      fVisitor.VisitPlacement(static_cast<G4PVPlacement*>(pv), depth);
      break;

    case kParameterised:
      fVisitor.VisitParameterised(static_cast<G4PVParameterised*>(pv), depth);
      break;

    case kReplica:
      fVisitor.VisitReplica(static_cast<G4PVReplica*>(pv), depth);
      break;

    default:  // kExternal or any future type
      fVisitor.VisitUnknown(pv, depth);
      break;
  }

  // Descend into the logical volume regardless of the placement type.
  // For replicas/parameterisations this visits the shared logical volume once.
  DescendLogical(pv->GetLogicalVolume(), depth);

  if (fVerbose)
  {
    G4cout << "Visitor - Calling the Depart 'Physical' method " << G4endl;
  }

  switch (pv->VolumeType())
  {
    case kNormal:
      fVisitor.DepartPlacement(static_cast<G4PVPlacement*>(pv), depth);
      break;

    case kParameterised:
      fVisitor.DepartParameterised(static_cast<G4PVParameterised*>(pv), depth);
      break;

    case kReplica:
      fVisitor.DepartReplica(static_cast<G4PVReplica*>(pv), depth);
      break;

    default:  // kExternal or any future type
      fVisitor.DepartUnknown(pv, depth);
      break;
  }

  if (fVerbose)
  {
    G4cout << "Walker::DescendPhysical depth= " << depth << " - exiting." << G4endl;
  }
}

// ---------------------------------------------------------------------------

void G4VolumeTreeWalker::DescendLogical(G4LogicalVolume* lv, G4int depth)
{
  if (lv == nullptr)
  {
    return;
  }

  if (fVerbose)
  {
    G4cout << "Walker::DescendLogical depth= " << depth << " - entering." << G4endl;
  }

  fVisitor.EnterLogicalVolume(lv, depth);

  const std::size_t nDaughters = lv->GetNoDaughters();
  for (std::size_t i = 0; i < nDaughters; ++i)
  {
    DescendPhysical(lv->GetDaughter(i), depth + 1);
  }

  if (fVerbose)
  {
    G4cout << "Walker::DescendLogical depth= " << depth << " - exiting." << G4endl;
  }

  fVisitor.LeaveLogicalVolume(lv, depth);
}
