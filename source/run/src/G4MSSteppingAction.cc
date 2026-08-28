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
// G4MSSteppingAction implementation
//
// Author: M.Asai, 5 May 2006
// --------------------------------------------------------------------

#include "G4MSSteppingAction.hh"

#include "G4ElementVector.hh"
#include "G4HadronicProcessStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"

// --------------------------------------------------------------------
G4MSSteppingAction::G4MSSteppingAction()
{
  defParticle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  if (defParticle == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "Proton is not defined in the physics list. Material scanner "
       << "does not work without proton.";
    G4Exception("G4MSSteppingAction::G4MSSteppingAction()", "MatScan0001", FatalException, ed);
  }
  defEKin = 100.0 * CLHEP::GeV;
}

// --------------------------------------------------------------------
void G4MSSteppingAction::Initialize(G4bool rSens, G4Region* reg)
{
  regionSensitive = rSens;
  theRegion = reg;
  length = 0.;
  massThick = 0.;
  x0 = 0.;
  lambda = 0.;
  shape_mat_info_v.clear();
}

// --------------------------------------------------------------------
void G4MSSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4Region* region = preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetRegion();

  if (regionSensitive && (region != theRegion)) return;

  G4double stlen = aStep->GetStepLength();
  const G4Material* material = preStepPoint->GetMaterial();
  if (intLengthMap.find(material) == intLengthMap.end())
  {
    auto* store = G4HadronicProcessStore::Instance();
    auto intLen = store->GetInteractionLength(material, defParticle, defEKin);
    if (intLen == DBL_MAX)
    {
      G4ExceptionDescription ed;
      ed << " PANIC ! Failure to get the interaction length of " << defParticle->GetParticleName()
         << " for " << material->GetName()
         << ".\n It is likely that proper hadronic processes are not defined in the physics list.";
      G4Exception("G4MSSteppingAction::UpdateInteractionLength", "MatScan0010", FatalException, ed);
    }
    intLengthMap[material] = intLen;
#ifdef MatScannerCheck
    G4cout << "Updating interaction length for " << material->GetName() << " : original "
           << material->GetNuclearInterLength() / CLHEP::cm << " : updated " << intLen / CLHEP::cm
           << G4endl;
#endif
  }
  G4double nil = intLengthMap[material];
  length += stlen;
  if (!ignoreThinMaterial || material->GetDensity() > thinDencity)
  {
    massThick += material->GetDensity() * stlen;
    x0 += stlen / (material->GetRadlen());
    lambda += stlen / nil;
  }

  // store information per step (1 geantino step = 1 solid)
  {
    shape_mat_info_v.push_back({});

    shape_mat_info_t& thisMaterialInfo = shape_mat_info_v.back();

    thisMaterialInfo.density = material->GetDensity();
    thisMaterialInfo.radiation_length = material->GetRadlen();
    thisMaterialInfo.interaction_length = nil;
    thisMaterialInfo.thickness = aStep->GetStepLength();
    thisMaterialInfo.integrated_thickness = length;
    thisMaterialInfo.massThickness = material->GetDensity() * stlen;
    thisMaterialInfo.integrated_massThickness = massThick;
    thisMaterialInfo.lambda = stlen / nil;
    thisMaterialInfo.x0 = stlen / (material->GetRadlen());
    thisMaterialInfo.integrated_lambda = lambda;
    thisMaterialInfo.integrated_x0 = x0;
    thisMaterialInfo.entry_point = preStepPoint->GetPosition();
    thisMaterialInfo.exit_point = aStep->GetPostStepPoint()->GetPosition();
    thisMaterialInfo.material_name = material->GetName();
    thisMaterialInfo.physvol_name = preStepPoint->GetPhysicalVolume()->GetName();
  }
}

void G4MSSteppingAction::PrintEachMaterialVerbose(std::ostream& oss)
{
  G4int colwidth = 11;
  G4int matname_colwidth = 15;

  oss << G4endl << G4endl;
  oss
    << "Phys vol          Material          Density    Radiation Interaction Thickness  Mass       "
       "X0         Lambda                      Entry point                              Exit "
       "point\n"
    << "   name              name                       length     length               length     "
       "                                           (cm)                                   (cm)\n"
    << "                                    (g/cm3)     (cm)       (cm)       (cm)      (g/cm2)";
  oss << G4endl;

  std::ios::fmtflags os_flags(oss.flags());
  for (auto& matInfo : shape_mat_info_v)
  {
    oss << std::setw(matname_colwidth) << std::left << matInfo.GetPhysVolName() << "   ";
    oss << std::setw(matname_colwidth) << std::left << matInfo.GetName(matname_colwidth) << "   ";

    oss << std::setw(colwidth) << std::scientific << std::setprecision(2)
        << matInfo.density / (CLHEP::g / CLHEP::cm3);
    oss << std::setw(colwidth) << std::scientific << matInfo.radiation_length / CLHEP::cm;
    oss << std::setw(colwidth) << std::scientific << matInfo.interaction_length / CLHEP::cm;
    oss << std::setw(colwidth) << std::scientific << matInfo.thickness / CLHEP::cm;
    oss << std::setw(colwidth) << std::scientific
        << matInfo.massThickness / (CLHEP::g / CLHEP::cm2);
    oss << std::setw(colwidth) << std::scientific << matInfo.x0;
    oss << std::setw(colwidth) << std::scientific << matInfo.lambda;
    oss << std::setw(colwidth) << " ";
    oss << std::scientific << std::right << "(";
    oss << std::scientific << std::left << matInfo.entry_point.x() / CLHEP::cm;
    oss << ", " << matInfo.entry_point.y() / CLHEP::cm;
    oss << ", " << matInfo.entry_point.z() / CLHEP::cm << ")";
    oss << std::setw(colwidth) << " ";
    oss << std::scientific << std::right << "(";
    oss << std::scientific << std::left << matInfo.exit_point.x() / CLHEP::cm;
    oss << ", " << matInfo.exit_point.y() / CLHEP::cm;
    oss << ", " << matInfo.exit_point.z() / CLHEP::cm << ")";
    oss << G4endl;
  }
  oss.flags(os_flags);  // Restore original stream format
}

void G4MSSteppingAction::PrintIntegratedMaterialVerbose(std::ostream& oss)
{
  // create database (db) of material name (key) and information
  std::map<G4String, shape_mat_info_t> mat_db;
  // accumulate information for each material name into mat_db
  for (auto& matInfo : shape_mat_info_v)
  {
    G4String key = matInfo.material_name;
    if (0 == mat_db.count(key))
    {
      mat_db[key] = shape_mat_info_t{};
      mat_db[key].material_name = key;
      mat_db[key].density = matInfo.density;
      mat_db[key].radiation_length = matInfo.radiation_length;
      mat_db[key].interaction_length = matInfo.interaction_length;
    }

    mat_db[key].x0 += matInfo.x0;
    mat_db[key].thickness += matInfo.thickness;
    mat_db[key].lambda += matInfo.lambda;
    mat_db[key].massThickness += matInfo.massThickness;
  }

  std::ios::fmtflags os_flags(oss.flags());
  G4int colwidth = 11;
  G4int matname_colwidth = 15;

  oss << G4endl << G4endl;
  oss
    << "   Material       Density    Radiation  Interaction Total     Total mass Total      Total\n"
    << "     name                     length     length    thickness  length     X0         "
       "Lambda\n"
    << "                  (g/cm3)     (cm)       (cm)        (cm)      (g/cm2)";
  oss << G4endl;

  for (auto& [key, mat] : mat_db)
  {
    oss << std::setw(matname_colwidth) << std::left << mat.GetName(matname_colwidth) << "   ";

    oss << std::setw(colwidth) << std::scientific << std::setprecision(2)
        << mat.density / (CLHEP::g / CLHEP::cm3);
    oss << std::setw(colwidth) << std::scientific << mat.radiation_length / CLHEP::cm;
    oss << std::setw(colwidth) << std::scientific << mat.interaction_length / CLHEP::cm;
    oss << std::setw(colwidth) << std::scientific << mat.thickness / CLHEP::cm;
    oss << std::setw(colwidth) << std::scientific << mat.massThickness / (CLHEP::g / CLHEP::cm2);
    oss << std::setw(colwidth) << std::scientific << mat.x0;
    oss << std::setw(colwidth) << std::scientific << mat.lambda;
    oss << G4endl;
  }
  oss.flags(os_flags);  // Restore original stream format

  return;
}

void G4MSSteppingAction::SetDetDefaultParticle(const G4ParticleDefinition* ptcl, G4double eKin)
{
  defParticle = ptcl;
  defEKin = eKin;
  intLengthMap.clear();
}
