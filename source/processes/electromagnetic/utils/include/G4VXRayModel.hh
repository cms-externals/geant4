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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VXRayModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.04.2025
//
//
// Class Description:
//
// Abstract interface to a X-Ray production model

// -------------------------------------------------------------------
//

#ifndef G4VXRayModel_h
#define G4VXRayModel_h 1

#include "globals.hh"
#include <vector>

class G4LogicalVolume;
class G4ParticleDefinition;
class G4DynamicParticle;
class G4LossTableManager;

class G4VXRayModel
{

public:

  explicit G4VXRayModel(const G4String& nam);

  virtual ~G4VXRayModel();

  void Initialise(std::vector<const G4LogicalVolume*>*);

  virtual void InitialiseModel();

  virtual G4bool IsApplicable(const G4ParticleDefinition*,
		   	      const G4LogicalVolume*);

  // if secondary X-rays produced they are inside vector out
  // and "true" is returned
  virtual G4bool SampleXRays(std::vector<G4DynamicParticle*>* out,
			     const G4DynamicParticle* in,
                             const G4LogicalVolume*) = 0;

  // for automatic documentation
  virtual void ModelDescription(std::ostream& outFile) const;

  const G4String& GetName() const { return pName; };

  //  hide assignment operator
  G4VXRayModel& operator=(const G4VXRayModel& right) = delete;
  G4VXRayModel(const G4VXRayModel&) = delete;

protected:

  std::vector<const G4LogicalVolume*>* pLogicalVolumes{nullptr};
  G4LossTableManager* pEmManager;
  G4int verbose{1};
  G4bool isMaster;
  const G4String pName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
