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
/// \file QuasiOpticalPhysics.hh
/// \brief QuasiOpticalPhysics constructor for quasi-optical processes
///
/// This class is equivalent to G4OpticalPhysics, with extensions for
/// quasi-optical processes to support offloading of optical photon
/// generation and transport.

#pragma once

#include "G4VPhysicsConstructor.hh"

class QuasiOpticalPhysics : public G4VPhysicsConstructor
{
  public:

    QuasiOpticalPhysics(G4int verb = 0, const G4String& name = "QuasiOptical");
    ~QuasiOpticalPhysics() override = default;
    void PrintStatistics() const;

    QuasiOpticalPhysics(const QuasiOpticalPhysics& right) = delete;
    QuasiOpticalPhysics& operator=(const QuasiOpticalPhysics& right) = delete;

  protected:

    // construct particle and physics
    void ConstructParticle() override;
    void ConstructProcess() override;

  private:

    void PrintWarning(G4ExceptionDescription&) const;
};
