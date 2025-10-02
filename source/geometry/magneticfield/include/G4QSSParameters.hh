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
// G4QSSParameters
//
// Hold parameters for QSS Integrator driver -- used to create
//  all QSStepper objects.  ( directly or via IntegrationDriver )
// Checks consistency of values proposed.
//
// This design means that objects of only *one* order of QSS driver
//  can be created.  ( QSS2 or QSS3 must be used globally.  )
// 
// Author: John Apostolakis,  Aug 2025
// --------------------------------------------------------------------
#ifndef G4QSSParameters_H
#define G4QSSParameters_H

#include "G4Types.hh"

class G4QSSParameters 
{
  public:

    static G4QSSParameters* Instance();

    ~G4QSSParameters() = default;

    // Access methods
    G4int    GetQssOrder()    { return fQssOrder; }
    G4double Get_dQRel()      { return fdQRel; }
    G4double Get_dQMin()      { return fdQMin; }
    G4int    GetMaxSubsteps() { return fMaxSubsteps; }

    G4bool   SetQssOrder( G4int value,  G4bool onlyWarn= false );
    G4bool   Set_dQRel( G4double dQRel );
    G4bool   Set_dQMin( G4double dQMin );
    G4bool   SetMaxSubsteps( G4int maxSubsteps );

  private:

     G4QSSParameters() = default;

  private:

     G4int    fQssOrder = 2;
     G4double fdQMin    = 0.00001;
     G4double fdQRel    = 0.001;

     G4int    fMaxSubsteps = 5000;
};

#endif
