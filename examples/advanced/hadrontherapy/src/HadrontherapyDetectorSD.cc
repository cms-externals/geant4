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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy


#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "HadrontherapyAnalysis.hh"
#include "G4SDManager.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyDetectorHit.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyLet.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyRBE.hh"
#include "G4TrackVector.hh"
#include "HadrontherapySteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4UserEventAction.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "HadrontherapyRunAction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Run.hh"
#include "G4EmCalculator.hh"
#include <G4AccumulableManager.hh>
#include "HadrontherapyRBEAccumulable.hh"
//#include <TTree.h>

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorSD::HadrontherapyDetectorSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="HadrontherapyDetectorHitsCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorSD::~HadrontherapyDetectorSD()
{
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorSD::Initialize(G4HCofThisEvent*)
{
    HitsCollection = new HadrontherapyDetectorHitsCollection(sensitiveDetectorName, collectionName[0]);
}

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    
    if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "RODetectorZDivisionPhys") return false;
    
    
    G4Track* theTrack = aStep  ->  GetTrack();
    G4Step fStep = *theTrack -> GetStep();
    //const std::vector<const G4Track*>* secondary= fStep.GetSecondaryInCurrentStep();
    //size_t SecondarySize = (*secondary).size();
    
    G4double kineticEnergy =  theTrack -> GetKineticEnergy();
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    
    //Get particle name
    G4String particleName =  particleDef -> GetParticleName();
    G4int pdg = particleDef ->GetPDGEncoding();
    G4int Z = particleDef-> GetAtomicNumber();
    G4int A = particleDef-> GetAtomicMass();
    
    // Get unique track_id (in an event)
    G4int trackID = theTrack -> GetTrackID();
    G4double DX = aStep -> GetStepLength();
    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
    
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    
    // Read voxel indexes: i is the x index, k is the z index
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int k  = touchable->GetReplicaNumber(0);
    G4int i  = touchable->GetReplicaNumber(2);
    G4int j  = touchable->GetReplicaNumber(1);
    
    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    G4String namePre = volumePre->GetName();
    //G4double eKin = aStep -> GetPreStepPoint() -> GetKineticEnergy();
    G4double parentID =aStep->GetTrack()->GetParentID();
    
    
    //G4int eventNum = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
    //G4double PosX = aStep->GetTrack()->GetPosition().x();
    //G4double PosY = aStep->GetTrack()->GetPosition().y();
    //G4double PosZ = aStep->GetTrack()->GetPosition().z();
    
    const G4LogicalVolume* VolumeVertex = aStep -> GetTrack() -> GetLogicalVolumeAtVertex();
    G4String VolumeAtVertex = VolumeVertex -> GetName();
    
    //G4double vertexPosition_x = aStep -> GetTrack() -> GetVertexPosition().x();
    
    G4ThreeVector VertexPosition = theTrack->GetVertexPosition();
    //G4double VertexEnergy = theTrack->GetVertexKineticEnergy();
    
    // get Material
    G4String material= aStep -> GetTrack() -> GetMaterial() -> GetName();
    
    HadrontherapyMatrix* matrix = HadrontherapyMatrix::GetInstance();
    HadrontherapyLet* let = HadrontherapyLet::GetInstance();
    
    G4int* hitTrack = matrix -> GetHitTrack(i,j,k);
    G4int voxel = matrix -> Index(i,j,k);
    
    
    // The following section fill two histograms. One contain the total
    // energy deposited in the first slice (voxel) of the phantom and due to all
    // particles (primaries and secondary). The second histogram contains
    // the kinetic energy of the primary particle in the first slice (voxel).
    //
    auto analysisManager = G4AnalysisManager::Instance();
    if (voxel == 100){
        if (parentID == 0)
        {
            analysisManager -> FillH1(1, kineticEnergy *MeV);
            //G4cout << "$$$$$   " << kineticEnergy << "Energy Dep   " << TotalEnergyDeposit << G4endl;
        }
        analysisManager->FillH1(2, energyDeposit);
    }
    
    //instanciate EmCalculator
    G4EmCalculator emCal;
    
    // G4double dedx = emCal.ComputeTotalDEDX(kineticEnergy,particleName,material,0.1*mm);
    
    
    //  ******************** let ***************************
    if (let)
    {
        if ( !(Z==0 && A==1) ) // All but not neutrons
        {
            if( energyDeposit>0. && DX >0. )
            {
                if (pdg !=22) // not gamma
                {
                    let -> FillEnergySpectrum(trackID, particleDef,energyDeposit, DX, i, j, k);
                }
                else if (kineticEnergy > 50.*keV) // gamma cut
                {
                    let -> FillEnergySpectrum(trackID, particleDef,energyDeposit, DX, i , j, k);
                }
                
            }
        }
    }
    
    
    if (matrix)
    {
        
        // Increment Fluences & accumulate energy spectra
        // Hit voxels are marked with track_id throught hitTrack matrix
        //G4int* hitTrack = matrix -> GetHitTrack(i,j,k); // hitTrack MUST BE cleared at every eventAction!
        if ( *hitTrack != trackID )
        {
            *hitTrack = trackID;
            /*
             * Fill FLUENCE data for every single nuclide
             * Exclude e-, neutrons, gamma, ...
             */
            if ( Z >= 1) matrix -> Fill(trackID, particleDef, i, j, k, 0, true);
            
        }
        
        if(energyDeposit != 0)
        {
            /*
             *  This method will fill a dose matrix for every single nuclide.
             *  A method of the HadrontherapyMatrix class (StoreDoseFluenceAscii())
             *  is called automatically at the end of main (or via the macro command /analysis/writeDoseFile.
             *  It permits to store all dose/fluence data into a single plane ASCII file.
             */
            // if (A==1 && Z==1) // primary and sec. protons
            if ( Z>=1 )    //  exclude e-, neutrons, gamma, ...
                matrix -> Fill(trackID, particleDef, i, j, k, energyDeposit);
            /*
             * Create a hit with the information of position is in the detector
             */
            HadrontherapyDetectorHit* detectorHit = new HadrontherapyDetectorHit();
            detectorHit -> SetEdepAndPosition(i, j, k, energyDeposit);
            HitsCollection -> insert(detectorHit);
        }
    }
    
    auto rbe = HadrontherapyRBE::GetInstance();
    if (rbe->IsCalculationEnabled())
    {
        if (!fRBEAccumulable)
        {
            fRBEAccumulable = dynamic_cast<HadrontherapyRBEAccumulable*>(G4AccumulableManager::Instance()->GetAccumulable("RBE"));
            if (!fRBEAccumulable)
            {
                G4Exception("HadrontherapyDetectorSD::ProcessHits", "NoAccumulable", FatalException, "Accumulable RBE not found.");
            }
        }
        
        fRBEAccumulable->Accumulate(kineticEnergy / A, energyDeposit, DX, Z, i, j, k);
    }
    return true;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}
