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
// V. Ivanchenko 23 September 2023 
//

#include "G4NeutronHPElasticXS.hh"
#include "G4Neutron.hh"
#include "G4ParticleHPManager.hh"
#include "G4SystemOfUnits.hh"

G4NeutronHPElasticXS::G4NeutronHPElasticXS()
  : G4CrossSectionHP(G4Neutron::Neutron(), "neutronElasticHP",
		     G4ParticleHPManager::GetInstance()->GetNeutronHPPath() + "/Elastic/CrossSection/",
                     20*CLHEP::MeV, 0, 100)
{
  SetMaxKinEnergy(20*CLHEP::MeV);
}

void G4NeutronHPElasticXS::CrossSectionDescription(std::ostream& outF) const
{
  outF << "High Precision cross data based on Evaluated Nuclear Data Files"
       << " (JEFF) for elastic scattering of neutrons below 20 MeV";
}
