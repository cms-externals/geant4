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
// Alternative Stepping Verbose method that prints out:
//   - momentum direction
//   - safety  ( distance to nearest boundary at pre-step point )
// Additionally refined to check the energy conservation for this
//   test  ( printing E_0 - E_k  and   1-dir_y )
// ********************************************************************
// Code developed by:
//    J. Apostolakis,   April 2021, November 2022
// ********************************************************************

#ifndef G24SteppingVerboseWithDir_hh
#  define G24SteppingVerboseWithDir_hh 1

#  include "G4SteppingVerbose.hh"

class G24SteppingVerboseWithDir : public G4SteppingVerbose
{
  public:

    // Constructor/Destructor
    G24SteppingVerboseWithDir();
    ~G24SteppingVerboseWithDir();

    void StepInfo() override;
    void TrackingStarted() override;

    // implementation methods - simplify/reduce code
    void Banner();

  private:

    G4double fInitialTrackKinE = 0.0;  // Store initial kinE of each track (to check accuracy)
};

#endif
