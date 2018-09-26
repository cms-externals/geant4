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
// $Id: G4CrossSectionElastic.cc 111432 2018-08-10 07:32:25Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CrossSectionElastic
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 19.11.2010
// Modifications:
//
// Class Description:
//
// Wrapper for elastic cross section build from a component
//
// -------------------------------------------------------------------
//

#include "G4CrossSectionElastic.hh"
#include "G4VComponentCrossSection.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4HadronicParameters.hh"

G4CrossSectionElastic::G4CrossSectionElastic(G4VComponentCrossSection* c,
					     G4int zmin, G4int zmax, 
					     G4double Emin, G4double Emax)
  : G4VCrossSectionDataSet(c->GetName()), component(c),
    Zmin(zmin),Zmax(zmax)
{
  nist = G4NistManager::Instance();
  SetMinKinEnergy(Emin);
  SetMaxKinEnergy(Emax);
}

G4CrossSectionElastic::~G4CrossSectionElastic()
{
  G4CrossSectionDataSetRegistry::Instance()->DeleteComponent(component);
}
   
G4bool G4CrossSectionElastic::IsElementApplicable(const G4DynamicParticle* p, 
						  G4int Z, const G4Material*)
{
  G4double e = p->GetKineticEnergy();
  return 
    (Z >= Zmin && Z <= Zmax && e >= GetMinKinEnergy() && e <= GetMaxKinEnergy()); 
}

G4double 
G4CrossSectionElastic::GetElementCrossSection(const G4DynamicParticle* p, 
					      G4int Z, 
					      const G4Material*)
{
  return component->GetElasticElementCrossSection(p->GetDefinition(), 
						  p->GetKineticEnergy(), 
						  Z, nist->GetAtomicMassAmu(Z));
}

void G4CrossSectionElastic::BuildPhysicsTable(const G4ParticleDefinition& p)
{
  component->BuildPhysicsTable(p);
  SetMaxKinEnergy(G4HadronicParameters::Instance()->GetMaxEnergy());
}

void G4CrossSectionElastic::DumpPhysicsTable(const G4ParticleDefinition& p)
{
  component->DumpPhysicsTable(p);
}

void G4CrossSectionElastic::CrossSectionDescription(std::ostream& /*outFile*/) const
{
  component->Description();
}


