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
/// \file LET.hh
/// \brief Definition of the RadioBio::LET class

#ifndef RADIOBIOLOGYLET_H
#define RADIOBIOLOGYLET_H
#include "globals.hh"

#include "IonLet.hh"
#include "VRadiobiologicalQuantity.hh"

#include <vector>

namespace RadioBio
{

// Forward declaration of other radiobiology classes
class LETAccumulable;
class LETMessenger;

/// Computes and stores dose- and track-averaged LET.
class LET : public VRadiobiologicalQuantity
{
  public:
    LET();
    ~LET() override;

    // Virtual methods to override
    void AddFromAccumulable(G4VAccumulable*) override;
    void Initialize() override;
    void Compute() override;
    void Reset() override;
    void Store() override;
    void PrintParameters() override;

  private:
    // Numerator and Denominator of Track-averaged LET
    array_type fNTotalLETT = {};
    array_type fDTotalLETT = {};

    // Numerator and Denominator of Dose-averaged LET
    array_type fNTotalLETD = {};
    array_type fDTotalLETD = {};

    // Track-averaged and Dose-averaged LET
    array_type fTotalLETT = {};
    array_type fTotalLETD = {};

    std::vector<IonLet> fIonLetStore{};

    // To be used for accumulation
    void AddNTotalLETT(const array_type& value) { fNTotalLETT += value; }
    void AddNTotalLETD(const array_type& value) { fNTotalLETD += value; }
    void AddDTotalLETT(const array_type& value) { fDTotalLETT += value; }
    void AddDTotalLETD(const array_type& value) { fDTotalLETD += value; }
};

}  // namespace RadioBio

#endif  // RADIOBIOLOGYLET_H
