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
// Geant4 header G4EmTableUtil
//
// Author V.Ivanchenko 14.03.2022
//
// Utilities used at initialisation of EM physics
//

#ifndef G4EmTableUtil_h
#define G4EmTableUtil_h 1

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4VMultipleScattering.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4MscStepLimitType.hh"

class G4EmTableUtil
{
public:

  static void PrepareMscProcess(G4VMultipleScattering* proc,
                                const G4ParticleDefinition& part,
			        G4EmModelManager* modelManager,
			        G4MscStepLimitType& stepLimit,
                                G4double& facrange,
			        G4bool& latDisplacement, G4bool& master,
			        G4bool& isIon, G4bool& baseMat);

  static void BuildMscProcess(G4VMultipleScattering* proc,
                              const G4VMultipleScattering* masterProc,
		              const G4ParticleDefinition& part,
		              const G4ParticleDefinition* firstPart,
		              G4int nModels, G4bool& master,
                              G4bool& baseMat);

  static G4bool StoreMscTable(G4VMultipleScattering* proc,
                              const G4ParticleDefinition* part,
                              const G4String& directory,
			      const G4int nModels, const G4int verb,
                              const G4bool ascii);

  static G4bool StoreTable(G4VProcess*, const G4ParticleDefinition*, 
                           G4PhysicsTable*, const G4String& dir,
                           const G4String& tname, G4int verb,
                           G4bool ascii);


};

#endif


