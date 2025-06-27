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
// SteppingAction program
// --------------------------------------------------------------

#include "Test15SteppingAction.hh"

#include "Test15EventAction.hh"

// #include "Test15AnalysisManager.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "globals.hh"
#include "G4ios.hh"

#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test15SteppingAction::Test15SteppingAction(Test15EventAction* eventAction)
  : G4UserSteppingAction(),evtAction(eventAction)  {

  startEnergy = 0.;

  flag = false;

  number_generations = 0;

  // number_shells = 26;
  number_shells = 25;
  //  G4double shell_thickness = 10.0*mm;
  G4double shell_thickness = 2.0*mm;

  shell_outer_radius = 457.0*mm;
  shell_inner_radius = shell_outer_radius - shell_thickness;

  // G4double radii_start[] = {200.0*cm,190.0*cm,185.0*cm,175.0*cm,165.0*cm,150.0*cm,140.0*cm,130.0*cm,120.0*cm,110.0*cm,100.0*cm,90.0*cm,80.0*cm,70.0*cm,60.0*cm,50.0*cm,45.7*cm,40.0*cm,30.0*cm,25.0*cm,20.0*cm,15.0*cm,10.0*cm,8.0*cm,5.0*cm,3.0*cm};
  G4double radii_start[] = {200.0*cm,190.0*cm,185.0*cm,175.0*cm,165.0*cm,150.0*cm,140.0*cm,130.0*cm,120.0*cm,110.0*cm,100.0*cm,90.0*cm,80.0*cm,70.0*cm,60.0*cm,50.0*cm,40.0*cm,30.0*cm,25.0*cm,20.0*cm,15.0*cm,10.0*cm,8.0*cm,5.0*cm,3.0*cm};

  for(G4int i=0; i<number_shells; ++i) {
    outer_radius[i] = radii_start[i];
    inner_radius[i] = radii_start[i] - shell_thickness;
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test15SteppingAction::~Test15SteppingAction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test15SteppingAction::UserSteppingAction(const G4Step* fStep)
{

  // return;

  // G4cout << " Got here Step 1 " << G4endl;
  // still needed? 27/09/15
  /*
  Test15AnalysisManager* analysis =  Test15AnalysisManager::getInstance();
  */

  //xfull  if(analysis->GetNtupleFull()) G4cout << " ntuple full is changed in stepping action " << G4endl;

  // removed 28/11/01 - unnecessary unless program "freezes"
  // kill track if too many steps
  // NB: This is set to DBL_MAX - therefore may cause program to "hang"
  //  G4int MaxNoSteps = DBL_MAX;
  //  G4int StepNo = fStep->GetTrack()->GetCurrentStepNumber();
  //  if(StepNo >= MaxNoSteps) fStep->GetTrack()->SetTrackStatus(fStopAndKill);

  G4int StepNo = fStep->GetTrack()->GetCurrentStepNumber();
  G4double partEnergy = fStep->GetPreStepPoint()->GetKineticEnergy();
  G4double partTime = fStep->GetPreStepPoint()->GetGlobalTime();
  //  if(partTime > 100*ns) fStep->GetTrack()->SetTrackStatus(fStopAndKill);
  G4double partMomentum = fStep->GetPreStepPoint()->GetMomentum().mag();
  G4double zMomentum = fStep->GetPreStepPoint()->GetMomentum().z();
  G4double angle = fStep->GetPreStepPoint()->GetMomentum().angle(fStep->GetPreStepPoint()->GetPosition());
  G4double cos_angle = std::abs(cos(angle));
  
  G4ParticleDefinition* particleType = fStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();
  G4double primaryEnergy = 0.0;

  //  if(StepNo == 0) 

  if(StepNo == 1 && particleName == "neutron") 
    { 
      // G4cout << " NEUTRON Got here +++++++++++++++++ " << G4endl;
      //x2018 evtAction->analyseParticleSource(partEnergy, particleName, partTime, partMomentum, zMomentum);

      startEnergy = partEnergy;
      startTime = partTime;
      flag = true;
    }

  // still needed? 27/09/15
  /*
  if(StepNo == 1 && particleName == "proton") 
    { 
      analysis->analyseProtonSource(partEnergy, particleName, partTime, partMomentum);
    }
  */
  G4int TrackID = fStep->GetTrack()->GetTrackID();
  G4int ParentTrackID = fStep->GetTrack()->GetParentID();

  if(StepNo == 1 && TrackID == 1) {
    parent_energy.clear();
    parent_particle.clear();
    parent_particleID.clear();
    number_generations = 0;
  }

  parent_energy[TrackID] = partEnergy;
  parent_particle[TrackID] = particleName;
  parent_particleID[TrackID] = ParentTrackID;

  G4bool reduced_tally = false;
  
  //xfull  if(TrackID == 5243 && StepNo > 1990) G4cout << " step: " << StepNo << " ntuple_full: " << analysis->GetNtupleFull() << G4endl;

  if(TrackID == 1 && StepNo == 1) primaryEnergy = partEnergy;
  if(StepNo == 1 && TrackID != 1) 
    { 
      // find number of generations - hopefully:
      number_generations = 1;
      G4int temp_id = TrackID;
      while(parent_particleID[temp_id] != 1) {
	//	G4cout << " TrackID: " << TrackID << " parentID: " << parent_particleID[temp_id] << " TempID: " << temp_id << " parent:" << ParentTrackID << " number of generations: " << number_generations << G4endl;
	temp_id = parent_particleID[temp_id];
	number_generations++;
      }

      if(parent_particle[ParentTrackID] == "neutron" && particleName == "neutron") {
	reduced_tally = true;
	parent_particle.erase(ParentTrackID);
	//x	number_generations++;
	//	G4cout << " erasing parentID: " << ParentTrackID << " and trackID: " << TrackID << G4endl;
      }

      //xfull      if(analysis->GetNtupleFull()) G4cout << " BEFORE SECONDARIES ntuple full is changed in stepping action " << G4endl;

//       G4cout << " calling analyse secondaries with ntuple full: " << analysis->GetNtupleFull() << G4endl;

      evtAction->analyseSecondaries(partEnergy, particleName, partTime, partMomentum, ParentTrackID, primaryEnergy,parent_energy[ParentTrackID],parent_particle[ParentTrackID],reduced_tally,number_generations);

      //      fStep->GetTrack()->SetTrackStatus(fStopAndKill);
      //      if(particleName == "neutron") fStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }
  if(particleName == "neutron") {
  // if(particleName == "neutron" && flag) {
    //x  Test15AnalysisManager* analysis =  Test15AnalysisManager::getInstance();
    evtAction->NeutronEnergyTime(partEnergy,partTime,startEnergy);
  } else {
    // G4cout << " particleName is: " << particleName << G4endl;
    if(particleName == "Pb207" || particleName == "Pb208") evtAction->OtherEnergyTime(partEnergy,partTime,startEnergy);
  }

  G4double radiusPre  = fStep->GetPreStepPoint()->GetPosition().mag();
  G4double radiusPost = fStep->GetPostStepPoint()->GetPosition().mag();
  G4double zPos  = fStep->GetPreStepPoint()->GetPosition().z();
  G4double StepLength  = fStep->GetStepLength();
  G4String Volume = fStep->GetTrack()->GetVolume()->GetName();
  
  G4bool entering_sph = false;
  G4bool exiting_sph = false;
  G4bool entering_sph_front = false;
  G4bool exiting_sph_front = false;
  G4bool entering_cyl = false;
  G4bool exiting_cyl = false;
  
  //  G4cout << " GOT HERE 1" << G4endl;
  G4TouchableHistory* thePreTouchable = 
    (G4TouchableHistory*)(fStep->GetPreStepPoint()->GetTouchable());
  G4TouchableHistory* thePostTouchable = 
    (G4TouchableHistory*)(fStep->GetPostStepPoint()->GetTouchable());
  G4int PreReplica = thePreTouchable->GetReplicaNumber();
  G4int PostReplica = thePostTouchable->GetReplicaNumber();
  //  G4cout << " GOT HERE 2" << G4endl;

  //xxx  if(particleName == "neutron" || particleName == "proton") {
  if(particleName == "neutron") {

    G4String PreVolGrichine = thePreTouchable->GetVolume()->GetName();
    if(PreVolGrichine=="check_phys" && fStep->GetTrack()->GetTrackID() > 1 ) {
      G4bool exitingGrichine = true;
      evtAction->exitingGrichineTally(exitingGrichine);
    }

    if(fStep->GetTrack()->GetNextVolume()) {
      // step points
      G4String PreVol = thePreTouchable->GetVolume()->GetName();
      G4String PostVol = thePostTouchable->GetVolume()->GetName();
      
      // exiting lead volume:
      G4bool exiting = false;
      G4bool exiting_check = false;
      
      if(PostVol == "sample_phys" && PostReplica == 0) entering_sph = true;
      if(PreVol == "sample_phys" && PreReplica == 0) exiting_sph = true;
      //xtemp      if(PostVol == "sample_phys2" && PostReplica == 0 && PostParentReplica == 51) entering_sph_front = true;
      //xtemp      if(PreVol == "sample_phys2" && PreReplica == 0 && PreParentReplica == 51) exiting_sph_front = true;
      if(PostVol == "sampleTube_phys") entering_cyl = true;
      if(PreVol == "sampleTube_phys") exiting_cyl = true;
      
//       if(entering_sph) G4cout << " entered sphere " << G4endl;
//       if(entering_cyl) G4cout << " entered cylinder " << G4endl;
//       if(exiting_sph) G4cout << " exiting sphere " << G4endl;
//       if(exiting_cyl) G4cout << " exiting cylinder " << G4endl;

      if(PreVol=="lab_phys" && PostVol=="world_phys") {
	//	G4cout << " exiting stepping action, prevol: " << PreVol << " postvol: " << PostVol << G4endl;
	exiting = true;
	evtAction->exitingTally(exiting,partEnergy);
      }

      if(PostVol=="world_phys") {
	//	G4cout << " exiting stepping action, prevol: " << PreVol << " postvol: " << PostVol << G4endl;
	exiting_check = true;
	evtAction->exitingTallyCheck(exiting_check);
      }

    }
    
    //    G4double radius = 1.0*(pow((pow(start.x(),2.)+pow(start.y(),2.)+pow(start.z(),2.)),0.5));
    //    if(abs((int)radius-456) < 10) {
    //xxx    if(abs((int)radius-228) < 1) {
    G4double radius = -9999.;

    //    if(((radiusPre && radiusPost) > 44.6*cm) && ((radiusPre && radiusPost) < 46.6*cm) && particleName == "neutron")

  // still needed? 27/09/15
    // G4double outer_radius = analysis->GetShellOuterRadius();
    // G4double inner_radius = analysis->GetShellInnerRadius();

    //    G4cout << " outer radius: " << G4BestUnit(outer_radius,"Length") << " inner radius: " << G4BestUnit(inner_radius,"Length") << G4endl;

    //debug    G4cout << " outer_radius: " << outer_radius << " inner_radius: " << inner_radius << G4endl;
//     if(((radiusPre > inner_radius ) && (radiusPost > inner_radius)) && ((radiusPre < outer_radius) && (radiusPost < outer_radius)) && particleName == "neutron")
    G4double my_tolerance = 1e-9*mm;
    G4bool pre_inside = false;
    G4bool post_inside = false;

    if((radiusPre <= (shell_outer_radius+my_tolerance) ) && (radiusPre >= (shell_inner_radius-my_tolerance)) && particleName == "neutron") {
      //      G4cout << " pre-step is inside shell " << radiusPre << G4endl;
      pre_inside = true;
    }
    if((radiusPost <= (shell_outer_radius+my_tolerance) ) && (radiusPost >= (shell_inner_radius-my_tolerance)) && particleName == "neutron") {
      //      G4cout << " post-step is inside shell " << radiusPost << G4endl;
      post_inside = true;
    }

//     if(pre_inside && radiusPost <= outer_radius) post_inside = true;
//     if(pre_inside && radiusPost >= inner_radius) post_inside = true;
//     if(post_inside && radiusPre <= outer_radius) pre_inside = true;
//     if(post_inside && radiusPre >= inner_radius) pre_inside = true;

//     if(pre_inside && post_inside != true) G4cout << " PROBLEM pre? radiusPost= " << radiusPost << " radiusPre= " << radiusPre << G4endl;
//     if(post_inside && pre_inside != true) G4cout << " PROBLEM post? radiusPost= " << radiusPost << " radiusPre= " << radiusPre << G4endl;

//     if((((radiusPre < (outer_radius-my_tolerance) ) && (radiusPre > (inner_radius+my_tolerance))) || ((radiusPost < (outer_radius-my_tolerance)) && (radiusPost > (inner_radius-my_tolerance)))) && particleName == "neutron") {

    //27/09/15 still needed...
    // G4cout << " Test15SteppingAction:: got here 1 " << " post_inside: " << post_inside << " pre_inside: " << pre_inside << G4endl;
    if(post_inside && pre_inside) {
      // G4cout << " Test15SteppingAction:: got here 2 " << G4endl;
      // getchar();
//       if(radiusPre < inner_radius) G4cout << " pre inner failed " << radiusPre << " difference: " << radiusPre-inner_radius << " post_inside is: " << post_inside << " pre_inside is: " << pre_inside << G4endl;
//       if(radiusPost < inner_radius) G4cout << " post inner failed " << radiusPost << " difference: " << radiusPost-inner_radius << " post_inside is: " << post_inside << " pre_inside is: " << pre_inside << G4endl;
//       if(radiusPre > outer_radius) G4cout << " pre outer failed " << radiusPre << " difference: " << radiusPre-outer_radius << " post_inside is: " << post_inside << " pre_inside is: " << pre_inside << G4endl;
//       if(radiusPost > outer_radius) G4cout << " post outer failed " << radiusPost << " difference: " << radiusPost-outer_radius << " post_inside is: " << post_inside << " pre_inside is: " << pre_inside << G4endl;
//       if((radiusPre < inner_radius) || (radiusPost < inner_radius) || (radiusPre > outer_radius) || (radiusPost > outer_radius)) G4cout << " STEP IS STRADDLING " << " outer: " << outer_radius << " and post-step: " << radiusPost << " inner: " << inner_radius << " and pre-step: " << radiusPre << G4endl; 
//      G4cout << " calling neutron shell " << G4endl;
//xbugcheck      G4cout << " radiusPre: " << G4BestUnit(radiusPre,"Length") << " radiusPost: " << G4BestUnit(radiusPost,"Length") << G4endl;
      evtAction->analyseNeutronShellFluence(partEnergy,partTime,startEnergy,TrackID,ParentTrackID,zMomentum,startTime,radius,zPos,parent_energy[ParentTrackID],parent_particle[ParentTrackID],StepLength,entering_sph,entering_cyl,exiting_sph,exiting_cyl,Volume,entering_sph_front,exiting_sph_front,0,0,PreReplica,PostReplica);
      //xRadial      G4cout << " called neutron shell fluence with step: " << StepLength << G4endl;
    }
    
    // problem? 100907
    // still needed: 27/09/15:
    // for(G4int i=0; i<analysis->GetNumberShells(); ++i) {
    // for(G4int i=0; i<26; ++i) { 
    for(G4int i=0; i<25; ++i) {
    // still needed: 27/09/15:
      // G4double radial_outer_radius = analysis->GetRadialOuterRadius(i);
      // G4double radial_inner_radius = analysis->GetRadialInnerRadius(i);
      G4double radial_outer_radius = outer_radius[i];
      G4double radial_inner_radius = inner_radius[i];
      G4double radial_tolerance = 1e-6*mm;
      G4bool pre_inside_radial = false;
      G4bool post_inside_radial = false;

      if((radiusPre <= (radial_outer_radius+radial_tolerance) ) && (radiusPre >= (radial_inner_radius-radial_tolerance)) && particleName == "neutron") {
	//      if((radiusPre <= radial_outer_radius) && (radiusPre >= radial_inner_radius) && particleName == "neutron") {
	//      G4cout << " pre-step is inside shell " << radiusPre << G4endl;
	pre_inside_radial = true;
      }
      if((radiusPost <= (radial_outer_radius+radial_tolerance) ) && (radiusPost >= (radial_inner_radius-radial_tolerance)) && particleName == "neutron") {
	//      if((radiusPost <= radial_outer_radius) && (radiusPost >= radial_inner_radius) && particleName == "neutron") {
	//      G4cout << " post-step is inside shell " << radiusPost << G4endl;
	post_inside_radial = true;
      }

//xbugcheck       if(pre_inside_radial && radiusPost <= radial_outer_radius) post_inside_radial = true;
//xbugcheck       if(pre_inside_radial && radiusPost >= radial_inner_radius) post_inside_radial = true;
//xbugcheck       if(post_inside_radial && radiusPre <= radial_outer_radius) pre_inside_radial = true;
//xbugcheck       if(post_inside_radial && radiusPre >= radial_inner_radius) pre_inside_radial = true;

//       if(pre_inside_radial && post_inside_radial != true) G4cout << " PROBLEM pre? radiusPost= " << radiusPost << " radiusPre= " << radiusPre << " whilst radial_outer_radius = " << radial_outer_radius << " and inner = " << radial_inner_radius << G4endl;
//       if(post_inside_radial && pre_inside_radial != true) G4cout << " PROBLEM post? radiusPost= " << radiusPost << " radiusPre= " << radiusPre << G4endl;
      
      //xbugcheck      if(post_inside_radial || pre_inside_radial)
      if(post_inside_radial && pre_inside_radial)
	{
	  //	  G4cout << " SHELLS: radiusPre: " << G4BestUnit(radiusPre,"Length") << " radiusPost: " << G4BestUnit(radiusPost,"Length") << G4endl;
	  // analysis->analyseNeutronRadialFluence(partEnergy,partTime,StepLength,i);
	  evtAction->analyseNeutronRadialFluence(partEnergy,partTime,StepLength,i);
//xRadial 	  if(i==17) G4cout << " called neutron radial fluence with step: " << StepLength << G4endl;
	}
    }
    // problem? 100907

    if(Volume == "sample_phys" || Volume == "sampleTube_phys" || Volume == "sample_phys2") {
      //       G4cout << " GOT HERE VOLUME: " << Volume << G4endl;
//       if(entering_sph) G4cout << " entering sphere flag true " << G4endl;
//       if(entering_cyl) G4cout << " entering cylinder flag true " << G4endl;
//       if(exiting_sph) G4cout << " exiting sphere flag true " << G4endl;
//       if(exiting_cyl) G4cout << " exiting cylinder flag true " << G4endl;
      //      if(entering) G4cout << " entering " << G4endl;
      radius = 456.0;
      // G4int PreParentReplica = thePreTouchable->GetReplicaNumber(1);
      // G4int PostParentReplica = thePostTouchable->GetReplicaNumber(1);
      //      G4cout << " parentReplica: " << PreParentReplica << " and: " << PostParentReplica << " and volume: " << Volume << G4endl;
      //      G4cout << " GOT HERE 3" << G4endl;

  // still needed? 27/09/15
  /*
      if(particleName == "neutron") {
	analysis->analyseNeutronFluence(partEnergy,partTime,startEnergy,TrackID,ParentTrackID,zMomentum,startTime,radius,zPos,parent_energy[ParentTrackID],parent_particle[ParentTrackID],StepLength,entering_sph,entering_cyl,exiting_sph,exiting_cyl,Volume,entering_sph_front,exiting_sph_front,PreParentReplica,PostParentReplica,PreReplica,PostReplica);
      }

  */
      //      G4cout << " GOT HERE 4 " << G4endl;
    }


    //    static const 
    G4double radii[] = {16.8*cm,40.4*cm,45.6*cm,69.1*cm,81.1*cm,98.6*cm,105.3*cm,113.5*cm,124.8*cm,153.9*cm};
    
    for(G4int i=0; i<10; ++i){
      if((radiusPre<radii[i] && radiusPost>radii[i]) || (radiusPre>radii[i] && radiusPost<radii[i])) {
	radius = radii[i];
	//      G4cout << " calling analyse neutron flux, step: " <<  StepNo << " trackID: " << TrackID << " parent: " << ParentTrackID << G4endl;
	//    radius = 45.6*cm;
	//	if((radiusPre<radius && radiusPost>radius) || (radiusPre>radius && radiusPost<radius)) {

	//xfull	if(TrackID == 5243 && StepNo > 1990) G4cout << " step: " << StepNo << " ntuple_full: " << analysis->GetNtupleFull() << " and calling analyseNeutronFlux " << G4endl;
  // still needed? 27/09/15
	// metal band at Hallenstadion 26/09/15
      // G4cout << " Test15SteppingAction:: GOT HERE 1!!!! " << " for i: " << i << G4endl
      // 	     << " RadiusPre: " << radiusPre << " radiusPost: " << radiusPost << " vs. " << radii[i] << G4endl;
	evtAction->analyseNeutronFlux(partEnergy,partTime,startEnergy,TrackID,ParentTrackID,zMomentum,startTime,radius,zPos,parent_energy[ParentTrackID],parent_particle[ParentTrackID],cos_angle,number_generations,particleName,reduced_tally);
	//      analysis->NeutronEnergyTime(partEnergy,partTime,startEnergy);
	//      fStep->GetTrack()->SetTrackStatus(fStopAndKill);
	//      flag = false;
      }
    }
  }

}


