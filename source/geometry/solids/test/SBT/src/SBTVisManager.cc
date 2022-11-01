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
// SBTVisManager
//
// Implementation of visualization manager for SBT
//

#include "SBTVisManager.hh"
#include "G4ModelingParameters.hh"

// BuildFakeWorld
//
G4int SBTVisManager::BuildFakeWorld() const
{
        //
        // These are probably leaks...
        //
        G4ModelingParameters::DrawingStyle style = G4ModelingParameters::wf;
        G4ModelingParameters *model = new G4ModelingParameters
                                            (0,      // No default vis attributes.
                                             style,  // Wireframe
                                             true,   // Global culling.
                                             true,   // Cull invisible volumes.
                                             false,  // Density culling.
                                             0.,     // Density (not relevant if density culling false).
                                             true,   // Cull daughters of opaque mothers.
                                             24);    // No of sides (not relevant for this operation).
        SBTFakeModel *fakeModel = new SBTFakeModel(model);

        G4Scene *currentScene = GetCurrentScene();

        if (!currentScene) {
                G4cerr << "Please create a view first" << G4endl;
                return 1;
        }
        
        currentScene->AddRunDurationModel( (G4VModel *)fakeModel );
        return 0;
}
