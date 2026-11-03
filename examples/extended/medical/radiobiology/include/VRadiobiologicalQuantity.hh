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
/// \file VRadiobiologicalQuantity.hh
/// \brief Definition of the RadioBio::VRadiobiologicalQuantity class

#ifndef RADIOBIOLOGYVRADIOBIOLOGICALQUANTITY_H
#define RADIOBIOLOGYVRADIOBIOLOGICALQUANTITY_H

#include "VoxelizedSensitiveDetector.hh"

#include <valarray>

class G4VAccumulable;
class G4UImessenger;

namespace RadioBio
{

/// Defines the common interface for radiobiological scoring quantities.
class VRadiobiologicalQuantity
{
  public:
    VRadiobiologicalQuantity()
    {
      if (VoxelizedSensitiveDetector::GetInstance() == nullptr) {
        G4Exception("VRadiobiologicalQuantity::VRadiobiologicalQuantity", "NotVoxelized",
                    FatalException,
                    "Trying to create a radiobiological quantity before voxelizing the detector!");
      }
    }

    virtual ~VRadiobiologicalQuantity() = default;

    // If this is false (set with macro), nothing happens
    G4bool IsCalculationEnabled() const { return fCalculationEnabled; }

    // If this is true, calculation is over and values are updated
    G4bool HasBeenCalculated() const { return fCalculated; }

    // Parameter setting
    void SetCalculationEnabled(G4bool enabled) { fCalculationEnabled = enabled; }

    // Verbosity for output
    void SetVerboseLevel(G4int level) { fVerboseLevel = level; }
    G4int GetVerboseLevel() const { return fVerboseLevel; }

    // Set path for output file
    void SetPath(const G4String& path) { fPath = path; }

    // Alias for matrix type
    using array_type = std::valarray<G4double>;

    // Update data from the accumulable
    virtual void AddFromAccumulable(G4VAccumulable*) = 0;

    // Initialization
    virtual void Initialize() = 0;

    // Calculation
    virtual void Compute() = 0;

    // Clear accumulated data
    virtual void Reset() = 0;

    // Output to text files (called at the end of run)
    virtual void Store() = 0;

    // Some basic output to the screen
    void PrintVirtualParameters()
    {
      G4cout << "** Verbose level = " << fVerboseLevel << " **********************" << G4endl
             << "** Initialized = " << (fInitialized ? "true *" : "false ")
             << "********************" << G4endl << "** Saved = " << (fSaved ? "true *" : "false ")
             << "**************************" << G4endl
             << "** Calculation enabled = " << (fCalculationEnabled ? "true *" : "false ")
             << "************" << G4endl << "*******************************************" << G4endl;
    }
    virtual void PrintParameters() { PrintVirtualParameters(); }

  protected:
    // Proper messenger
    G4UImessenger* fMessenger = nullptr;

    // Verbose Level
    G4int fVerboseLevel = 1;

    // Bool to see if class has been initialized
    G4bool fInitialized = false;

    // Bool to see if output has already been saved
    G4bool fSaved = false;

    // Bool to see if quantity has been calculated
    G4bool fCalculated = false;

    // Output paths
    G4String fPath{};

    // Boolean for calculation
    G4bool fCalculationEnabled = false;
};

}  // namespace RadioBio

#endif  // RADIOBIOLOGYVRADIOBIOLOGICALQUANTITY_H
