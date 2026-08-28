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
//---------------------------------------------------------------------------
//
// ClassName:   G4NewExoticPhysics - test physics constructor (in namespace)
//
// Author: 2015-05-11 R. Hatcher
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4NewExoticPhysics_h
#  define G4NewExoticPhysics_h 1

#  include "G4Decay.hh"
#  include "G4VPhysicsConstructor.hh"
#  include "globals.hh"

namespace myns
{

class G4NewExoticPhysics : public G4VPhysicsConstructor
{
  public:

    G4NewExoticPhysics(G4int ver = 1);
    G4NewExoticPhysics(const G4String& name, G4int ver = 1);
    virtual ~G4NewExoticPhysics();

  public:

    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();

    virtual G4Decay* GetDecayProcess() { return fDecayProcess; }

  private:

    static G4ThreadLocal G4Decay* fDecayProcess;
    G4int verbose;
    static G4ThreadLocal G4bool wasActivated;
};

}  // namespace myns

#endif
