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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <set>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction(),
  fEdepVectorCtor3(2U, 0.),
  fEdepVectorCtor3n("EdepVector3n", 2U, 0.),
  fEdepVectorCtor4(2U),
  fEdepVectorCtor4n("EdepVector4n", 2U),
  fEdepArrayCtor1n("EdepArray1n"),
  fProcCounterMapCtor1n("ProcCounterMap1n"),
  fProcCounterUMapCtor1n("ProcCounterUMap1n")
{
  // add new units for dose
  DefineUnits();

  // Get accumulable manager
  auto accumulableManager = G4AccumulableManager::Instance();

  // Values
  //

  // save all values in fDValues
  fDValues.push_back(&fEdep);               // accumulable_0
  fDValues.push_back(&fEdep2);              // Edep2
  fDValues.push_back(&fMaxEdep);            // accumulable_2
  fDValues.push_back(&fMinEdep);            // accumulable_3
  fDValues.push_back(&fEventCounter);       // accumulable_4
  fDValues.push_back(&fEventCounter2);      // accumulable_5
  fDValues.push_back(&fEventCounter2);      // refused as already registered
  // register all values
  for (auto& element : fDValues) {
    accumulableManager->Register(*element);
  }

  // Accumulables can be also created via accumulable manager
  // (they are registered automatically
  auto edepBis = accumulableManager->CreateAccValue<G4double>("EdepBis", 0.);
  auto edep2Bis = accumulableManager->CreateAccValue<G4double>(0.);
  fDValues.push_back(edepBis);              // EdepBis
  fDValues.push_back(edep2Bis);             // accumulable_7


  accumulableManager->Register(fNsec);      // accumulable_8
  accumulableManager->Register(fPassed);    // accumulable_9

  // User defined accumulable
  //

  fProcCounter = new ProcCounterAccumulable("ProcCounter");
  accumulableManager->Register(fProcCounter);     // ProcCounter
  // check attributing a name
  fProcCounterTest1 = new ProcCounterAccumulable();
  accumulableManager->Register(fProcCounterTest1); // accumulable_11
  G4cout << "fProcCounterTest1 name: " << fProcCounterTest1->GetName() << G4endl;

  // check that an accumulable with already existing name is not accepted
  ProcCounterAccumulable* procCounterTest2 = new ProcCounterAccumulable("ProcCounter");
  accumulableManager->Register(procCounterTest2);  // refused as already registered

  // Vectors
  //

  // initialize vectors, that were not initialized on creation
  fEdepVectorCtor1.push_back(0.);
  fEdepVectorCtor1.push_back(1.);
  fEdepVectorCtor4[0] = 0.;
  fEdepVectorCtor4[1] = 0.;

  // save all vectors in fVectors
  fVectors.push_back(&fEdepVectorCtor1);   // accumulable_12
  fVectors.push_back(&fEdepVectorCtor3);   // accumulable_13
  fVectors.push_back(&fEdepVectorCtor3n);  // EdepVector3n
  fVectors.push_back(&fEdepVectorCtor4);   // accumulable_15
  fVectors.push_back(&fEdepVectorCtor4n);  // EdepVector4n
  fVectors.push_back(&fEdepVectorCtor10);  // accumulable_17
  fVectors.push_back(&fEdepVectorCtor10n); // EdepVector10n

  // register all vectors
  for (auto& element : fVectors) {
    accumulableManager->Register(*element);
  }

  // Arrays
  //

  // save all arrays in fArrays
  fArrays.push_back(&fEdepArrayCtor1);    // accumulable_19
  fArrays.push_back(&fEdepArrayCtor1n);   // EdepArray1n
  fArrays.push_back(&fEdepArrayCtor2);    // accumulable_21
  fArrays.push_back(&fEdepArrayCtor3);    // EdepArray3

  // register all arrays
  for (auto& element : fArrays) {
    accumulableManager->Register(*element);
  }

  // Maps
  //

  // save all maps in fMaps
  fMaps.push_back(&fProcCounterMapCtor1);    // accumulable_23
  fMaps.push_back(&fProcCounterMapCtor1n);   // ProcCounterMap1n
  fMaps.push_back(&fProcCounterMapCtor10);   // accumulable_25
  fMaps.push_back(&fProcCounterMapCtor10n);  // ProcCounterMap10n

  // register all maps
  for (auto& element : fMaps) {
    accumulableManager->Register(*element);
  }

  // Unordered maps
  //

  // save all unordered maps in fUnorderedMaps
  fUnorderedMaps.push_back(&fProcCounterUMapCtor1);    // accumulable_27
  fUnorderedMaps.push_back(&fProcCounterUMapCtor1n);   // ProcCounterUMap1n
  fUnorderedMaps.push_back(&fProcCounterUMapCtor10);   // accumulable_29
  fUnorderedMaps.push_back(&fProcCounterUMapCtor10n);  // ProcCounterUMap10n

  // register all unordered maps
  for (auto& element : fUnorderedMaps) {
    accumulableManager->Register(*element);
  }

  // Print all defined accumulables names
  if (isMaster) {
    G4cout << "List of defined accumulables:" << G4endl;
    std::vector<G4VAccumulable*>::const_iterator it;
    for ( it = accumulableManager->BeginConst(); it != accumulableManager->EndConst(); it++) {
      G4cout << "Accumulable: " << (*it)->GetName() << G4endl;
    }
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fProcCounter;
  delete fProcCounterTest1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::DefineUnits()
{
  // add new units for dose
  //

  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//_____________________________________________________________________________
void RunAction::PrintAccRange(G4int startId, G4int count, const G4String& accType) const
{
  G4cout << "Accumulables of " << accType << " type:" << G4endl;
  auto accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Print(startId, count);
  G4cout << "-----" << G4endl;
}

//_____________________________________________________________________________
void RunAction::PrintAllAccumulables() const
{
  // Test printing by range
  PrintAccRange(fDValues[0]->GetId(), G4int(fDValues.size()), "Value<G4double>");
  PrintAccRange(fNsec.GetId(), 1, "Value<G4int>");
  PrintAccRange(fPassed.GetId(), 1, "Value<G4bool>");
  PrintAccRange(fProcCounter->GetId(), 1, "User");
  PrintAccRange(fVectors[0]->GetId(), G4int(fVectors.size()), "vector<G4double>");
  PrintAccRange(fArrays[0]->GetId(), G4int(fArrays.size()), "array<G4int, 2");
  PrintAccRange(fMaps[0]->GetId(), G4int(fMaps.size()), "map<G4String, G4int>");
  PrintAccRange(fUnorderedMaps[0]->GetId(), G4int(fUnorderedMaps.size()),
                    "unordered_map<G4String, G4int>");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Check access to accumulables by name
  for (auto it = accumulableManager->Begin(); it != accumulableManager->End(); ++it) {
    auto accumulable = *it;
    auto name = accumulable->GetName();
    auto findByName = accumulableManager->GetAccumulable(name);
    if (accumulable != findByName) {
      G4cerr << "Failed to find accumulable " << accumulable->GetName() << " by name." << G4endl;
    }
  }
  G4cout << "Check access to accumulables by name: done" << G4endl;
  G4cout << "-----" << G4endl;

  // Check access to accumulables by id via type specific access methods
  TestGetById();
  TestGetByName();

  // Print all accumulables (ourselves)
  if (isMaster) {
    G4cout << "Our printing" << G4endl;
    PrintAllAccumulables();
  }


  // Print all accumulables by accumulable manager
  if (isMaster) {
    G4cout << "G4AccumulableManager printing" << G4endl;
    accumulableManager->Print();
  }

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  const DetectorConstruction* detectorConstruction
   = static_cast<const DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Compute average number of secondary e-
  G4double nsec  = fNsec.GetValue();
  G4double nsecPerEvent = nsec/nofEvents;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const PrimaryGeneratorAction* generatorAction
   = static_cast<const PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : "
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << " Average number of e-  per event : " << nsecPerEvent
     << G4endl
     << " Max edep : " << G4BestUnit(fMaxEdep.GetValue(),"Energy")
     << G4endl
     << " Min edep : " << G4BestUnit(fMinEdep.GetValue(),"Energy")
     << G4endl
     << " Passed primary : " << std::boolalpha << fPassed.GetValue()
     << G4endl
     << " Event counter: " << fEventCounter.GetValue()
     << G4endl
     << " Event counter2: " << fEventCounter2.GetValue()
     << G4endl;

  //frequency of processes
  fProcCounter->Print();

  G4cout
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddNsec(G4int nsec)
{
  fNsec  += nsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;

  fMaxEdep = std::max(fMaxEdep.GetValue(), edep);
  if ( edep != 0. ) {
    fMinEdep = std::min (fMinEdep.GetValue(), edep);
  }

  // Accumulables defined via accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();

  // Get accumulables from manager
  auto edepBis = accumulableManager->GetAccValue<G4double>("EdepBis");
  auto edep2Bis = accumulableManager->GetAccValue<G4double>("accumulable_7");

  // Update accumulables from manager
  if (edepBis != nullptr) {
    (*edepBis)  += edep;
  }
  if (edep2Bis != nullptr) {
    (*edep2Bis) += edep*edep;
  }

  // Print error message if accumulables were not retrieved
  if (edepBis == nullptr) {
    G4cerr << "Failed to get EdepBis" << G4endl;
  }
  if (edep2Bis == nullptr) {
    G4cerr << "Failed to get Edep2Bis" << G4endl;
  }

  // Vectors
  for (auto& element : fVectors) {
    // G4cout << element->GetName() << ", vector size: " << element->size() << G4endl;
    (*element)[0] += edep;
    (*element)[1] += edep*edep;
  }

  // Arrays
  for (auto& element : fArrays) {
    // G4cout << element->GetName() << ", vector size: " << element->size() << G4endl;
    (*element)[0] += edep;
    (*element)[1] += edep*edep;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountEvent()
{
  if ( fEventCounter.GetValue() < 5 ) {
    G4cout << "Prefix counter: "  << (++fEventCounter).GetValue() << G4endl;
    G4cout << "Postfix counter: " << (fEventCounter2++).GetValue() << G4endl;
  } else {
    ++fEventCounter;
    fEventCounter2++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcess(G4String procName) 
{ 
  fProcCounter->CountProcesses(procName);

  // count process directly in accumulable maps
  for (auto procCounterMap : fMaps) {
    auto it = procCounterMap->find(procName);
    if ( it == procCounterMap->end()) {
      (*procCounterMap)[procName] = 1;
    }
    else {
      (*procCounterMap)[procName]++; 
    }
  }

  // count process directly in accumulable unordered maps
  for (auto procCounterUMap : fUnorderedMaps) {
    auto it = procCounterUMap->find(procName);
    if ( it == procCounterUMap->end()) {
      (*procCounterUMap)[procName] = 1;
    }
    else {
      (*procCounterUMap)[procName]++;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetPassed(G4bool passed)
{
  //fPassed = fPassed.GetValue() || passed;
  fPassed = passed;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

