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
/// \file PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4ParticleTable.hh"
#include "G4SingleParticleSource.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* pPrimaryGenerator)
  : G4UImessenger(), fpPrimaryGenerator(pPrimaryGenerator)
{
  fparticle = std::make_unique<G4UIcmdWithAString>("/UHDR/source/particle", this);
  fparticle->SetGuidance("e-, proton, alpha");

  fenergy = std::make_unique<G4UIcmdWithADoubleAndUnit>("/UHDR/source/energy", this);
  fenergy->SetGuidance("Sets a monocromatic energy (same as  gps/energy)");
  fenergy->SetParameterName("monoenergy", false, false);
  fenergy->SetDefaultUnit("keV");

  //L. T. Anh: /UHDR/source/ion for defininng ion beam:
  fion = std::make_unique<G4UIcommand>("/UHDR/source/ion", this);
  fion->SetGuidance("Set properties of ion to be generated.");
  fion->SetGuidance("[usage] /UHDR/source/ion Z A Q E");
  fion->SetGuidance("        Z:(int) AtomicNumber");
  fion->SetGuidance("        A:(int) AtomicMass");
  fion->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
  fion->SetGuidance("        E:(double) Excitation energy (in keV)");
   
  auto* pZ = new G4UIparameter("Z", 'i', false);
  pZ->SetParameterRange("Z>0");
  pZ->SetDefaultValue(1);
  fion->SetParameter(pZ);
  auto* pA = new G4UIparameter("A", 'i', false);
  pA->SetParameterRange("A>0");
  pA->SetDefaultValue(1);
  fion->SetParameter(pA);
  auto* pQ = new G4UIparameter("Q", 'i', true);   
  pQ->SetDefaultValue(0);                      
  fion->SetParameter(pQ);
  auto* pE = new G4UIparameter("E", 'd', true);    
  pE->SetDefaultValue(0.);
  fion->SetParameter(pE);
  fion->AvailableForStates(G4State_Idle);          //after /run/initialize

  fenergyPerNucleon = std::make_unique<G4UIcmdWithADoubleAndUnit>(
                      "/UHDR/source/energyPerNucleon", this);
  fenergyPerNucleon->SetGuidance("Sets kinetic energy per nucleon (E_total = value * A).");
  fenergyPerNucleon->SetGuidance(
    "must be called after /UHDR/source/ion; A is taken from baryon number.");
  fenergyPerNucleon->SetParameterName("energyPerNucleon", false, false);
  fenergyPerNucleon->SetDefaultUnit("MeV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  const auto particleGun = fpPrimaryGenerator->GetSPGun();
  if (command == fparticle.get()) {
    G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->FindParticle(newValue);
    if (pd != nullptr) {
      particleGun->SetParticleDefinition(pd);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Particle [" << newValue << "] is not found.";
      command->CommandFailed(ed);
    }
  }
  else if (command == fenergy.get()) {
    const auto pEnergyDis = particleGun->GetEneDist();
    pEnergyDis->SetMonoEnergy(fenergy->GetNewDoubleValue(newValue));
  } else if (command == fion.get()) {
    G4Tokenizer next(newValue);          // splits on whitespace by default
    const G4int Z = StoI(next());
    const G4int A = StoI(next());
    G4int Q = Z;                         // default: fully stripped ion
    G4double Eexc = 0.;

    const G4String sQ = next();
    if (!sQ.empty() && sQ != "*") {
      Q = StoI(sQ);
      const G4String sE = next();
      if (!sE.empty()) Eexc = StoD(sE) * keV;
    }

    auto* ion = G4IonTable::GetIonTable()->GetIon(Z, A, Eexc);
    if (ion != nullptr) {
      particleGun->SetParticleDefinition(ion);
      particleGun->SetParticleCharge(Q * eplus);
    }
    else {
      G4ExceptionDescription ed;
      ed << "Ion Z=" << Z << " A=" << A << " cannot be created.";
      command->CommandFailed(ed);
    }
  } else if (command == fenergyPerNucleon.get()) {
    const auto* pd = particleGun->GetParticleDefinition();
    if (pd == nullptr) {
      G4ExceptionDescription ed;
      ed << "No particle/ion defined yet; call /UHDR/source/particle or /UHDR/source/ion first.";
      command->CommandFailed(ed);
      return;
    }

    const G4int A = pd->GetBaryonNumber();
    if (A <= 0) {
      G4ExceptionDescription ed;
      ed << "Particle [" << pd->GetParticleName() << "] has invalid baryon number A=" << A
        << "; energyPerNucleon không áp dụng được.";
      command->CommandFailed(ed);
      return;
    }

    const G4double ePerNucleon = fenergyPerNucleon->GetNewDoubleValue(newValue);
    const auto pEnergyDis = particleGun->GetEneDist();
    pEnergyDis->SetMonoEnergy(ePerNucleon * A);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......