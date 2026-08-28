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
// This is a copy of the CMS hadronic configuration for Run-3
// Updated 05.11.2025
//

#include "CMSHadronPhysicsFTFP_BERT.hh"

#include "G4ios.hh"

#include <CLHEP/Units/SystemOfUnits.h>

CMSHadronPhysicsFTFP_BERT::CMSHadronPhysicsFTFP_BERT(G4int)
  : CMSHadronPhysicsFTFP_BERT(3. * CLHEP::GeV, 6. * CLHEP::GeV, 12 * CLHEP::GeV)
{}

CMSHadronPhysicsFTFP_BERT::CMSHadronPhysicsFTFP_BERT(G4double e1, G4double e2, G4double e3)
  : G4HadronPhysicsFTFP_BERT("hInelastic FTFP_BERT", false)
{
  minFTFP_pion = minFTFP_kaon = minFTFP_proton = minFTFP_neutron = e1;
  maxBERT_kaon = maxBERT_proton = maxBERT_neutron = e2;
  maxBERT_pion = e3;
}

void CMSHadronPhysicsFTFP_BERT::DumpBanner()
{
  G4cout << "### CMS version of FTFP_BERT : transition between BERT and FTFP is over the interval "
         << minFTFP_proton / CLHEP::GeV << " to " << maxBERT_proton / CLHEP::GeV << " GeV"
         << " GeV; for pions up to " << maxBERT_pion / CLHEP::GeV << " GeV" << G4endl;
}
