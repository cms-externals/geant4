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
/// \file G4MoleculeCounterScorer.hh
/// \brief Definition of the G4MoleculeCounterScorer class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025~26)
//

#ifndef G4MOLECULECOUNTERSCORER_HH
#define G4MOLECULECOUNTERSCORER_HH

#include "G4MoleculeCounter.hh"
#include "G4VPrimitiveMoleculeScorer.hh"

//------------------------------------------------------------------------------

class G4MoleculeCounterScorer
  : public G4VPrimitiveMoleculeScorer<G4MoleculeCounter, G4MoleculeCounterIndex, G4String>
{
  public:
    G4MoleculeCounterScorer(G4String, G4String, G4int = 0);
    ~G4MoleculeCounterScorer() override = default;

    std::tuple<G4String> ConvertCounterIndexToKey(const G4MoleculeCounterIndex&) const override;
    void WriteWithAnalysisManager() const override;
};

#endif
