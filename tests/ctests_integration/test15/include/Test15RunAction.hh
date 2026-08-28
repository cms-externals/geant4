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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo
//                    (27th November 2001)
//
// RunAction header
// --------------------------------------------------------------

#ifndef Test15RunAction_h
#  define Test15RunAction_h 1

#  include "G4AnalysisManager.hh"
#  include "G4UserRunAction.hh"
#  include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;
class Test15Run;

class Test15RunAction : public G4UserRunAction
{
  public:

    Test15RunAction();
    virtual ~Test15RunAction();

  public:

    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    void FillRadialExperimentalData();

    void createNeutronFluxHisto(G4int, const Test15Run*);
    void createRadialFluxHisto(G4int, const Test15Run*);
    void secondarySummary(G4int, const Test15Run*);

  private:

    G4double local_energy_integral[4];
    G4double copy_radial_fluence_1[76];
    G4double copy_radial_error_1[76];

    G4int number_shells;

    G4double outer_radius[26];
    G4double inner_radius[26];
    G4double shell_outer_radius;
    G4double shell_inner_radius;

    G4AnalysisManager* analysisManager;
};

#endif
