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
/// \file Dose.hh
/// \brief Definition of the RadioBio::Dose class

#ifndef RADIOBIOLOGYDOSE_H
#define RADIOBIOLOGYDOSE_H

#include "globals.hh"

#include "VRadiobiologicalQuantity.hh"

#include <map>
#include <valarray>

namespace RadioBio
{

// Forward declaration of other radiobiology classes
class DetectorConstruction;
class DoseAccumulable;
class DoseMessenger;

/// Computes and stores voxel-resolved absorbed dose.
class Dose : public VRadiobiologicalQuantity
{
  public:
    Dose();
    ~Dose() override;

    // Virtual methods to override
    void AddFromAccumulable(G4VAccumulable*) override;
    void Initialize() override;
    void Compute() override;
    void Reset() override;
    void Store() override;
    void PrintParameters() override;

    // Get a copy of the total dose object
    const array_type& GetDose() const { return fDose; }

  private:
    using particle_map_type = std::map<G4String, array_type>;

    // Total, primary and all-secondary deposited energy and dose
    array_type fEnDep = {};
    array_type fPrimaryEnDep = {};
    array_type fSecondaryEnDep = {};
    array_type fDose = {};
    array_type fPrimaryDose = {};
    array_type fSecondaryDose = {};

    // Non-primary dose separated by particle species
    particle_map_type fSecondaryEnDepByParticle{};
    particle_map_type fSecondaryDoseByParticle{};

    // To be used for accumulation
    void AddEnergyDeposit(const array_type& dep) { fEnDep += dep; }
    void AddPrimaryEnergyDeposit(const array_type& dep) { fPrimaryEnDep += dep; }
    void AddSecondaryEnergyDeposit(const array_type& dep) { fSecondaryEnDep += dep; }
    void AddSecondaryEnergyDepositByParticle(const particle_map_type& depByParticle);
};

}  // namespace RadioBio

#endif  // RADIOBIOLOGYDOSE_H
