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
// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestCalorimeterSD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestCalorimeterSD.hh"

#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Positron.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestCalorimeterSD::hTestCalorimeterSD(G4String name)
 :G4VSensitiveDetector(name),
  theHisto(hTestHisto::GetPointer()),
  evno(0),
  evnOld(-1),
  trIDold(-1),
  delta(1.0e-6*mm)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestCalorimeterSD::~hTestCalorimeterSD()
{
  energy.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  evno++;
  if(0 < theHisto->GetVerbose()) 
    G4cout << "hTestCalorimeterSD: Begin Of Event # " << evno << G4endl;

  numAbs = theHisto->GetNumberOfAbsorbers();
  energy.resize(numAbs);
  for(G4int i=0; i<numAbs; i++) { energy[i] = 0.0; }
  backEnergy = 0.0;
  leakEnergy = 0.0;
  G4double gap = theHisto->GetGap();
  zmax = (theHisto->GetAbsorberThickness() + gap)*numAbs - gap - delta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool hTestCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  //  if(0.0 == edep) return true;

  theHisto->AddTrackLength(aStep->GetStepLength());

  G4int j = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
   
  if(j < 0 || j >= numAbs) {
      G4cout << "Warning!!! hTestCalorimeterSD: cannot add " << edep/MeV
             << " MeV to the slice # " << j << G4endl;

  } else {
      energy[j] += edep;
  }

  if(1 < theHisto->GetVerbose()) {
      G4cout << "hTestCalorimeterSD: energy = " << edep/MeV
             << " MeV is deposited at " << j
             << "-th absorber slice " << G4endl;
  }

  const G4Track* track = aStep->GetTrack();
  G4int trIDnow  = track->GetTrackID();
  G4double tkin  = track->GetKineticEnergy(); 
  G4double theta = (track->GetMomentumDirection()).theta();
  G4double zend  = aStep->GetPostStepPoint()->GetPosition().z();
  //  G4double zstart= aStep->GetPreStepPoint()->GetPosition().z();

  G4bool stop = false;
  G4bool primary = false;
  G4bool outAbs = false;
  if(track->GetNextVolume()->GetName() == "World"
     && (zend <= delta || zend >= zmax-delta)  ) outAbs = true;

  if(tkin == 0.0) stop = true;
  if(0 == aStep->GetTrack()->GetParentID()) primary = true;

  // new particle
  if(trIDnow != trIDold || evno != evnOld) {
    trIDold = trIDnow;
    evnOld  = evno;
    part_is_out = true;
    tkinold = aStep->GetPreStepPoint()->GetKineticEnergy();
  }

  // Primary particle stop 
   
  if(primary && (stop || outAbs)) {

      G4double xend = aStep->GetPostStepPoint()->GetPosition().x();
      G4double yend = aStep->GetPostStepPoint()->GetPosition().y();
      theHisto->SaveToTuple("xend",xend/mm);      
      theHisto->SaveToTuple("yend",yend/mm);      
      theHisto->SaveToTuple("zend",zend/mm);      
      theHisto->SaveToTuple("ltpk",(theHisto->GetTrackLength())/mm);      
      theHisto->SaveToTuple("tend",tkin/MeV);
      theHisto->SaveToTuple("teta",theta);      
      theHisto->SaveToTuple("loss",(tkinold-tkin)/MeV);      
      theHisto->SaveToTuple("dedx",(tkinold-tkin)*mm/(zmax*MeV));      

    // exclude cases when track return back

      if(theHisto->GetTrackLength() < 2.0*zend) theHisto->AddEndPoint(zend);

  }

  // After step in absorber

  if(outAbs && part_is_out) {
    G4double e = tkin;
    if(track->GetDefinition() == G4Positron::Positron()) 
      e += 2.*electron_mass_c2;
    if(zend > zmax-delta)       leakEnergy += e;
    else if(zend < delta)       backEnergy += e;
    part_is_out = false;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int i, j;

  theHisto->SaveToTuple("back",backEnergy);      
  theHisto->SaveToTuple("leak",leakEnergy);      
  G4double etot = 0.0;

  // The histogramm on the energy deposition profile
  if(numAbs > 0) {
    G4double s = theHisto->GetAbsorberThickness();
    G4double z = -0.5 * s;
    for(i=0; i<numAbs; i++) {
      z += s; 
      etot += energy[i];
      theHisto->AddEnergy(energy[i], z);
    }
  }
  theHisto->SaveToTuple(G4String("edep"),etot/MeV);      
  //theHisto->SaveToTuple(G4String("Evt"),G4double(evno));      

  // Integrated energy deposition to nTuple
  G4int nMax = 60;
  G4double EE[60];
  G4String eSlice[60]={
      "S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", 
      "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", 
      "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", 
      "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", 
      "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", 
      "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59"};
  G4String eInteg[60]={
      "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", 
      "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", 
      "E20", "E21", "E22", "E23", "E24", "E25", "E26", "E27", "E28", "E29", 
      "E30", "E31", "E32", "E33", "E34", "E35", "E36", "E37", "E38", "E39", 
      "E40", "E41", "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49", 
      "E50", "E51", "E52", "E53", "E54", "E55", "E56", "E57", "E58", "E59"};

  G4int k = theHisto->GetNumAbsorbersSaved();
  if (nMax > k) nMax = k;

  if(nMax > 1) {
    for(i=0; i<nMax; i++){
      EE[i]=0.0;
      for(j=0; j<i+1; j++) {
        EE[i] += energy[j];
      }  
      //   theHisto->SaveToTuple(eSlice[i],energy[i]);      
      // theHisto->SaveToTuple(eInteg[i],EE[i]);      
    }
  }

  // Dump information about this event 
  theHisto->SaveEvent();

  if(0 < theHisto->GetVerbose()) { 
    G4cout << "hTestCalorimeterSD: Event #" << evno << " ended" << G4endl;
    G4cout << "Edep(MeV)= " << etot/MeV 
           << "; backE(MeV)= " << backEnergy/MeV
           << "; leakE(MeV)= " << leakEnergy/MeV
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestCalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void hTestCalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







