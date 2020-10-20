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
// G4EqMagElectricField
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field.

// Created: V.Grichine, 10.11.1998
// -------------------------------------------------------------------
#ifndef G4EQMAGELECTRICFIELD_HH
#define G4EQMAGELECTRICFIELD_HH

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"
#include "G4ElectroMagneticField.hh"

class G4EqMagElectricField : public G4EquationOfMotion
{
  public:  // with description

    G4EqMagElectricField(G4ElectroMagneticField* emField );

    ~G4EqMagElectricField();

    void  SetChargeMomentumMass(G4ChargeState particleCharge, // in e+ units
                                G4double MomentumXc,
                                G4double mass);

    void EvaluateRhsGivenB(const G4double y[],
                           const G4double Field[],
                                 G4double dydx[] ) const;
      // Given the value of the electromagnetic field, this function 
      // calculates the value of the derivative dydx.

  private:

    G4double fElectroMagCof = 0.0;
    G4double fMassCof = 0.0;
};

#endif
