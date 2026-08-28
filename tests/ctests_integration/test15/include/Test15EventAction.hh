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
// EventAction header
// --------------------------------------------------------------

#ifndef Test15EventAction_h
#  define Test15EventAction_h 1

#  include "G4AnalysisManager.hh"
#  include "G4RunManager.hh"
#  include "G4SystemOfUnits.hh"
#  include "G4THitsMap.hh"
#  include "G4UnitsTable.hh"
#  include "G4UserEventAction.hh"
#  include "G4ios.hh"
#  include "globals.hh"

#  include "Test15Run.hh"

// #include "Test15LeadHit.hh"
// #include "Test15SampleHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class Test15EventAction : public G4UserEventAction
{
  public:

    Test15EventAction();
    virtual ~Test15EventAction();
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  public:

    void SetPrintModulo(G4int val) { printModulo = val; };

    // fill histograms with data from Test15ParticleSource / secondary history
    // x2018 void analyseParticleSource(G4double, G4String, G4double, G4double, G4double);
    // fill histograms with data from Test15ParticleSource / secondary history
    void analyseSecondaries(G4double, G4String, G4double, G4double, G4int, G4double, G4double,
                            G4String, G4bool, G4int);

    // fill historgram and tuple with neutron time vs. energy
    void NeutronEnergyTime(G4double, G4double, G4double);
    // fill historgram and tuple with other particle time vs. energy
    void OtherEnergyTime(G4double, G4double, G4double);

    // number of neutrons leaving volume:
    void exitingTally(G4bool, G4double);

    // number of neutrons leaving volume:
    void exitingGrichineTally(G4bool);

    // number of neutrons on stack:
    void AddToNeutronStack();

    // number of neutrons leaving volume:
    void exitingTallyCheck(G4bool);

    // fill histograms with data from Test15NeutronFlux / secondary history
    void analyseNeutronFlux(G4double, G4double, G4double, G4int, G4int, G4double, G4double,
                            G4double, G4double, G4double, G4String, G4double, G4int, G4String,
                            G4bool);

    // fill histograms with data from Test15NeutronFlux / secondary history
    // x2018 void analyseNeutronFluence(G4double, G4double, G4double, G4int, G4int, G4double,
    // G4double, G4double, G4double, G4double, G4String, G4double, G4bool, G4bool, G4bool, G4bool,
    // G4String,G4bool,G4bool,G4int,G4int,G4int,G4int);

    void analyseNeutronShellFluence(G4double, G4double, G4double, G4int, G4int, G4double, G4double,
                                    G4double, G4double, G4double, G4String, G4double, G4bool,
                                    G4bool, G4bool, G4bool, G4String, G4bool, G4bool, G4int, G4int,
                                    G4int, G4int);

    void analyseNeutronRadialFluence(G4double, G4double, G4double, G4int);

  private:

    // methods
    // x G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
    // x                                         const G4Event* event) const;

    // // methods
    // Test15LeadHitsCollection* GetLeadHitsCollection(G4int hcID,
    // 						const G4Event* event) const;

    G4int event_id;

    G4int printModulo;

    G4bool ntuple_full;

    G4int fNeutronStack;
};

// inline functions

// x2018 inline void Test15EventAction::analyseParticleSource(G4double source_energy, G4String
// source_name, G4double source_time, G4double source_momentum, G4double source_zMomentum)
//  {
/*
// get analysis manager
G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


G4int Iparticle=-9;
G4double temp_time = time/microsecond;
G4double temp_energy = energy/eV;
G4double temp_momentum = momentum;
//  G4double temp_momentum = momentum2;
// to access 1/p distribution will have to use uwfunc to weight bins

if(name == "gamma") {
  analysisManager->FillH1(1, energy/keV);
  // hGammaEdep->fill(energy/keV);
  Iparticle = 1;
}
if(name == "neutron") {
  // G4cout << " filling neutron lethargy with: energy = " << energy/eV << " and momentum = " <<
1/momentum << G4endl; analysisManager->FillH1(2, energy/eV,1/momentum);
  // hNeutronEdep->fill(energy/eV,1/momentum);  // fill(x,weight)
  Iparticle = 2;
  if(energy/MeV < 2.0) {
    neutflux[0]++;
    enflux[0]+=energy;
  }
  else if(energy/MeV > 2.0 && energy/MeV < 20.0) {
    neutflux[1]++;
    enflux[1]+=energy;
  }
  else if(energy/MeV > 20.0) {
    neutflux[2]++;
    enflux[2]+=energy;
  }
  if(energy/MeV > 1000.0) {
    neutflux[3]++;
    enflux[3]+=energy;
  }
}
if(name == "e-") {
  analysisManager->FillH1(3, energy/keV);
  // hElectronEdep->fill(energy/keV);  // fill(x,weight)
  Iparticle = 3;
}
if(name == "e+") {
  analysisManager->FillH1(4, energy/keV);
  // hPositronEdep->fill(energy/keV);  // fill(x,weight)
  Iparticle = 4;
}
if(name == "other") {
  analysisManager->FillH1(5, energy/keV);
  // hOtherEdep->fill(energy/keV);  // fill(x,weight)
  Iparticle = 5;
}


// ntupleSource->fill( ntupleSource->findColumn( "energy" ), (G4float) temp_energy);
// ntupleSource->fill( ntupleSource->findColumn( "time" ), (G4float) temp_time);
// ntupleSource->fill( ntupleSource->findColumn( "particle" ), (G4float) Iparticle);
// ntupleSource->fill( ntupleSource->findColumn( "momentum" ), (G4float) temp_momentum);
// ntupleSource->fill( ntupleSource->findColumn( "zmom" ), (G4float) zMomentum);
// ntupleSource->addRow();

//  G4cout << " Filled Source Ntuple " << G4endl;

*/
// x2018 }

inline void Test15EventAction::exitingTally(G4bool exiting_flag, G4double energy)
{
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (exiting_flag)
  {
    Test15Run* run =
      static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    // G4double test = 999.;
    run->AddExitingFlux(energy);
    // exiting_flux++;
    // 27/09/15:
    analysisManager->FillNtupleDColumn(2, 0, energy);
    analysisManager->AddNtupleRow(2);
    // ntupleExiting->fill( ntupleExiting->findColumn( "energy" ), (G4float) energy);
    // ntupleExiting->addRow();
  }
}

inline void Test15EventAction::exitingGrichineTally(G4bool exiting_flag)
{
  // get analysis manager
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (exiting_flag)
  {
    Test15Run* run =
      static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    // G4double test = 999.;
    run->AddExitingGrichineFlux();
    // // exiting_flux++;
    // //27/09/15:
    // analysisManager->FillNtupleDColumn(2,0, energy);
    // analysisManager->AddNtupleRow(2);
    // // ntupleExiting->fill( ntupleExiting->findColumn( "energy" ), (G4float) energy);
    // // ntupleExiting->addRow();
  }
}

inline void Test15EventAction::AddToNeutronStack()
{
  fNeutronStack++;
}

inline void Test15EventAction::exitingTallyCheck(G4bool exiting_flag_check)
{
  if (exiting_flag_check)
  {
    Test15Run* run =
      static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    // G4double test = 999.;
    run->AddExitingCheckFlux();
    // exiting_check_flux++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Test15EventAction::analyseSecondaries(G4double energy, G4String name, G4double time,
                                                  G4double momentum, G4int ParentID,
                                                  G4double primaryEnergy, G4double parentEnergy,
                                                  G4String parentParticle, G4bool reduced_flux,
                                                  G4int number_generations)
{
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // if(ntuple_full) {
  //   G4cout << " 1) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
  //   G4double break_it = std::sqrt(-100.0);
  // }

  G4int Iparticle = -9;
  G4double temp_time = time / microsecond;
  G4double temp_energy = energy / eV;
  G4double temp_momentum = momentum;
  //  G4double temp_momentum = momentum2;

  Test15Run* run = static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  if (name == "gamma")
  {
    run->AddFlux(name);
    // gamma_flux++;
    Iparticle = 1;
    //    return;
  }
  else if (name == "neutron")
  {
    if (!reduced_flux) run->AddFlux("neutron_check");  // neutron_check++;
    run->AddFlux(name);
    // neutron_flux++;
    Iparticle = 2;
    if (reduced_flux) Iparticle = -2;
    if (energy > 0.1 * eV && energy < 10.0 * keV)
      run->AddFlux("neutron_fluence");  // neutron_fluence++;
  }
  else if (name == "e-")
  {
    run->AddFlux(name);
    // electron_flux++;
    Iparticle = 3;
    //    return;
  }
  else if (name == "pi-")
  {
    run->AddFlux(name);
    // piminus_flux++;
    Iparticle = 4;
  }
  else if (name == "pi+")
  {
    run->AddFlux(name);
    // piplus_flux++;
    Iparticle = 5;
  }
  else if (name == "pi0")
  {
    run->AddFlux(name);
    // pizero_flux++;
    Iparticle = 6;
  }
  else if (name == "e+")
  {
    run->AddFlux(name);
    // positron_flux++;
    Iparticle = 7;
  }
  else if (name == "proton")
  {
    run->AddFlux(name);
    // proton_flux++;
    Iparticle = 8;
  }
  else if (name == "mu-")
  {
    run->AddFlux(name);
    // proton_flux++;
    Iparticle = 9;
  }
  else if (name == "mu+")
  {
    run->AddFlux(name);
    // proton_flux++;
    Iparticle = 10;
  }
  else
  {
    run->AddFlux("other");
    // other_flux++;
    Iparticle = 99;
    return;
  }

  G4int iParent = 0;
  if (parentParticle == "gamma")
    iParent = 1;
  else if (parentParticle == "neutron")
    iParent = 2;
  else if (reduced_flux)
    iParent = -2;
  else if (parentParticle == "e-")
    iParent = 3;
  else if (parentParticle == "pi-")
    iParent = 4;
  else if (parentParticle == "pi+")
    iParent = 5;
  else if (parentParticle == "pi0")
    iParent = 6;
  else if (parentParticle == "e+")
    iParent = 7;
  else if (parentParticle == "proton")
    iParent = 8;
  else if (parentParticle == "proton")
    iParent = 8;
  else if (parentParticle == "mu-")
    iParent = 9;
  else if (parentParticle == "mu+")
    iParent = 10;

  // if(ntuple_full) {
  //   G4cout << " 2) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
  //   G4double break_it = std::sqrt(-100.0);
  // }

  if (ntuple_full)
  {
    // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleDColumn(0, 0, temp_energy);
    analysisManager->FillNtupleDColumn(0, 1, temp_time);
    analysisManager->FillNtupleIColumn(0, 2, Iparticle);
    analysisManager->FillNtupleDColumn(0, 3, temp_momentum);
    analysisManager->FillNtupleIColumn(0, 4, ParentID);
    analysisManager->FillNtupleDColumn(0, 5, primaryEnergy);
    analysisManager->FillNtupleIColumn(0, 6, iParent);
    analysisManager->FillNtupleDColumn(0, 7, parentEnergy);
    analysisManager->FillNtupleIColumn(0, 8, number_generations);
    analysisManager->FillNtupleIColumn(0, 9, event_id);
    analysisManager->AddNtupleRow();

    // ntupleSecondary->fill( ntupleSecondary->findColumn( "energy" ), (G4float) temp_energy);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "time" ), (G4float) temp_time);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "particle" ), (G4float) Iparticle);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "momentum" ), (G4float) temp_momentum);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "parentid" ), (G4float) ParentID);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "e_prim" ), (G4float) primaryEnergy);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "parent" ), (G4float) iParent);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "e_parent" ), (G4float) parentEnergy);
    // ntupleSecondary->fill( ntupleSecondary->findColumn( "numgen" ), (G4float)
    // number_generations); ntupleSecondary->fill( ntupleSecondary->findColumn( "event" ), (G4float)
    // event_id); ntupleSecondary->addRow();
  }

  //  G4cout << " Filled Source Ntuple " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// x2018 inline void Test15EventAction::analyseNeutronFluence(G4double energy, G4double time,
// G4double startEnergy, G4int TrackID, G4int ParentID, G4double zMomentum, G4double startTime,
// G4double radius, G4double zPos, G4double parentEnergy, G4String parentParticle, G4double
// steplength, G4bool enter_sph, G4bool enter_cyl, G4bool exit_sph, G4bool exit_cyl, G4String
// Volume, G4bool enter_sph_front, G4bool exit_sph_front, G4int preParentReplica, G4int
// postParentReplica, G4int preReplica, G4int postReplica)
// {
// get analysis manager
// x2018 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

// 27/09/15:
/*
  if(enter_sph) {
    if(TrackID == oldTrackID) duplicate_neutron2++;
    else duplicate_neutron2 = 0;
    oldTrackID2 = TrackID;
    //x    G4cout << " duplicate_neutron2: " << duplicate_neutron2 << G4endl;
  }

  G4double temp_energy = energy/eV;
//   if(enter_sph) G4cout << " ENTERING SPHERE " << G4endl;
//   if(enter_cyl) G4cout << " ENTERING CYLINDER " << G4endl;
//   if(exit_sph) G4cout << " EXITING SPHERE " << G4endl;
//   if(exit_cyl) G4cout << " EXITING CYLINDER " << G4endl;

  G4double temp_time = time/microsecond;
  G4double temp_startenergy = startEnergy/eV;

  G4int iParent = 0;
  if (parentParticle == "gamma")   iParent = 1;
  if (parentParticle == "neutron") iParent = 2;
  if (parentParticle == "e-")      iParent = 3;
  if (parentParticle == "pi-")     iParent = 4;
  if (parentParticle == "pi+")     iParent = 5;
  if (parentParticle == "pi0")     iParent = 6;
  if (parentParticle == "e+")      iParent = 7;
  if (parentParticle == "proton")  iParent = 8;

  if(ntuple_full) {
    G4cout << " 2) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
    G4double break_it = std::sqrt(-100.0);
  }

  if(ntuple_full) {
    ntupleFluence->fill( ntupleFluence->findColumn( "energy" ), (G4float) temp_energy);
    ntupleFluence->fill( ntupleFluence->findColumn( "time" ), (G4float) temp_time);
    ntupleFluence->fill( ntupleFluence->findColumn( "starte" ), (G4float) temp_startenergy);
    ntupleFluence->fill( ntupleFluence->findColumn( "trackid" ), (G4float) TrackID);
    ntupleFluence->fill( ntupleFluence->findColumn( "parentid" ), (G4float) ParentID);
    ntupleFluence->fill( ntupleFluence->findColumn( "zmom" ), (G4float) zMomentum);
    ntupleFluence->fill( ntupleFluence->findColumn( "startt" ), (G4float) startTime);
    ntupleFluence->fill( ntupleFluence->findColumn( "radius" ), (G4float) radius);
    ntupleFluence->fill( ntupleFluence->findColumn( "e_parent" ), (G4float) parentEnergy);
    ntupleFluence->fill( ntupleFluence->findColumn( "parent" ), (G4float) iParent);
    ntupleFluence->fill( ntupleFluence->findColumn( "step" ), (G4float) steplength);
    ntupleFluence->fill( ntupleFluence->findColumn( "dupli" ), (G4float) duplicate_neutron2);

    ntupleFluence->addRow();
  }

    //  if(std::abs(radius-456.0) > 0.1 && std::abs(radius-50.0) > 0.1 ) return;

  G4double temp_flux_energy = -999.0;
  G4int flux_idx = -999;
  G4double temp_low_energy = -999.0;
  G4double temp_lithium_energy = -999.0;
  G4int flux_idx_low = -999;

  if(temp_energy<lithium_energy[100]) {
    for(G4int i = 0; i<100; ++i) {
      if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1]) {
  if(enter_sph) lithium_fluence[i]++;
  if(Volume == "sample_phys" && (preReplica == 0 || postReplica == 0)) lithium_fluence_step[i] +=
steplength/mm; if(Volume == "sample_phys2" && (preReplica == 0 || postReplica == 0) &&
(preParentReplica == 51 || postParentReplica == 51)) lithium_fluence_front_step[i] += steplength/mm;
  if(enter_cyl) lithium_fluence_cyl[i]++;
  if(Volume == "sampleTube_phys") lithium_fluence_step_cyl[i] += steplength/mm;
      }
    }
  }

    //  if(temp_energy<flux_energy[0]) {
  if(temp_energy<low_energy[100]) {
    for(G4int i = 0; i<100; ++i) {
      if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1]) {
  if(enter_sph) low_fluence[i]++;
  if(Volume == "sample_phys" && (preReplica == 0 || postReplica == 0)) low_fluence_step[i] +=
steplength/mm; if(Volume == "sample_phys2" && (preReplica == 0 || postReplica == 0) &&
(preParentReplica == 51 || postParentReplica == 51)) low_fluence_front_step[i] += steplength/mm;
  temp_low_energy = std::exp(0.5*(log(low_energy[i+1])+log(low_energy[i])));
  flux_idx_low = i;
  if(enter_cyl) low_fluence_cyl[i]++;
  if(Volume == "sampleTube_phys") low_fluence_step_cyl[i] += steplength/mm;

      }
    }
  }
  else if(temp_energy>flux_energy[0]) {
    for(G4int i = 0; i<21; ++i) {
      //     if(energy > *(flux_energy+i) && energy < *(flux_energy+i+1)) {
      //    if(energy/eV > flux_energy[i] && energy/eV < flux_energy[i+1]) {
      //    if(energy > flux_energy[i]/1000000.0 && energy < flux_energy[i+1]/1000000.0) {
      if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1]) {
  if(enter_sph) fluence[i]++;
  if(Volume == "sample_phys" && (preReplica == 0 || postReplica == 0)) fluence_step[i] +=
steplength/mm; if(Volume == "sample_phys2" && (preReplica == 0 || postReplica == 0) &&
(preParentReplica == 51 || postParentReplica == 51)) fluence_front_step[i] += steplength/mm;
    //xx	G4cout << " i: " << i << " temp_energy: " << temp_energy << " flux_energy[i]:" <<
flux_energy[i] << " flux_energy[i+1]: " << flux_energy[i+1] << " flux: " << flux[i] << G4endl;
  flux_idx = i;
  if(enter_cyl)fluence_cyl[i]++;
  if(Volume == "sampleTube_phys") fluence_step_cyl[i] += steplength/mm;
      }
    }
  }
  //     G4cout << " energy: " << energy/eV << " flux energy0: " << flux_energy[i] << " flux
energy1: " << flux_energy[i+1] << G4endl;

  if(ntuple_full) {
    G4cout << " 3) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
    G4double break_it = std::sqrt(-100.0);
  }

  if(ntuple_full) {
    nFluence3->fill( nFluence3->findColumn( "fluxe" ), (G4float) temp_flux_energy);
    nFluence3->fill( nFluence3->findColumn( "fluxidx" ), (G4float) flux_idx);
    nFluence3->fill( nFluence3->findColumn( "lowe" ), (G4float) temp_low_energy);
    nFluence3->fill( nFluence3->findColumn( "lowidx" ), (G4float) flux_idx_low);
    nFluence3->fill( nFluence3->findColumn( "step" ), (G4float) steplength);
    nFluence3->addRow();
  }

//   // integrate step lengths:
//     if(temp_energy<lithium_energy[100]) {
//       for(G4int i = 0; i<100; ++i) {
// 	if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1]) lithium_fluence_step[i]
+= steplength/mm;
//       }

//       if(temp_energy<low_energy[100]) {
// 	for(G4int i = 0; i<100; ++i) {
// 	  if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1]) low_fluence_step[i] +=
steplength/mm;
// 	}
//       }
//       else if(temp_energy>flux_energy[0]) {
// 	for(G4int i = 0; i<21; ++i) {
// 	  if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1]) fluence_step[i] +=
steplength/mm;
// 	}
//       }

//     }
*/
// x2018 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Test15EventAction::analyseNeutronRadialFluence(G4double energy, G4double time,
                                                           G4double steplength, G4int radius_index)
{
  // get analysis manager
  //  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  Test15Run* run = static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->analyseNeutronRadialFluence(energy, time, steplength, radius_index);

  // 27/09/15:
  /*
if(radius_index < 0 || radius_index > 26) G4cout << " WARNING radius index is wrong!!!!!! " <<
radius_index << G4endl;

G4double temp_energy = energy/eV;
G4double temp_time = time/microsecond;

//xRadial  if(radius_index == 17) G4cout << " temp_energy: " << temp_energy << "
lithium_radial_energy_upper[10]: " << lithium_radial_energy_upper[9] << G4endl;

if(temp_energy < lithium_radial_energy_upper[9] && temp_energy > lithium_radial_energy_lower[0]) {
  for (G4int i=0; i<10; ++i) {
    if(temp_energy > lithium_radial_energy_lower[i] && temp_energy < lithium_radial_energy_upper[i])
//xRadial{
radial_fluence_step[radius_index][i] += steplength/mm;
//xRadial	G4cout << " FILLING radial array " << G4endl;
//xRadial      } else {
//xRadial	G4cout << " SKPPED radial array " << temp_energy << G4endl;
//xRadial      }

  }
}
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Test15EventAction::NeutronEnergyTime(G4double energy, G4double time,
                                                 G4double startEnergy)
{
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (time > 0. && energy > 0.)
    analysisManager->FillH2(1, log10(time / microsecond), log10(energy / eV), 1.0);

  G4double temp_time = time / microsecond;
  G4double temp_energy = energy / eV;
  G4double temp_startenergy = startEnergy / eV;

  // ntupleEnergyTime->fill( ntupleEnergyTime->findColumn( "energy" ), (G4float) temp_energy);
  // ntupleEnergyTime->fill( ntupleEnergyTime->findColumn( "time" ), (G4float) temp_time);
  // ntupleEnergyTime->fill( ntupleEnergyTime->findColumn( "primary" ), (G4float) temp_startenergy);
  // ntupleEnergyTime->addRow();

  if (ntuple_full)
  {
    analysisManager->FillNtupleDColumn(1, 0, temp_energy);
    analysisManager->FillNtupleDColumn(1, 1, temp_time);
    analysisManager->FillNtupleDColumn(1, 2, temp_startenergy);
    analysisManager->AddNtupleRow(1);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Test15EventAction::OtherEnergyTime(G4double energy, G4double time, G4double startEnergy)
{
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (time > 0 && energy > 0)
  {
    analysisManager->FillH2(2, log10(time / microsecond), log10(energy / eV), 1.0);
  }

  G4double temp_time = time / microsecond;
  G4double temp_energy = energy / eV;
  G4double temp_startenergy = startEnergy / eV;

  // ntupleEnergyTime->fill( ntupleEnergyTime->findColumn( "energy" ), (G4float) temp_energy);
  // ntupleEnergyTime->fill( ntupleEnergyTime->findColumn( "time" ), (G4float) temp_time);
  // ntupleEnergyTime->fill( ntupleEnergyTime->findColumn( "primary" ), (G4float) temp_startenergy);
  // ntupleEnergyTime->addRow();

  if (ntuple_full)
  {
    analysisManager->FillNtupleDColumn(14, 0, temp_energy);
    analysisManager->FillNtupleDColumn(14, 1, temp_time);
    analysisManager->FillNtupleDColumn(14, 2, temp_startenergy);
    analysisManager->AddNtupleRow(14);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void Test15EventAction::analyseNeutronShellFluence(
  G4double energy, G4double time, G4double startEnergy, G4int TrackID, G4int ParentID,
  G4double zMomentum, G4double startTime, G4double radius, G4double zPos, G4double parentEnergy,
  G4String parentParticle, G4double steplength, G4bool enter_sph, G4bool enter_cyl, G4bool exit_sph,
  G4bool exit_cyl, G4String Volume, G4bool enter_sph_front, G4bool exit_sph_front,
  G4int preParentReplica, G4int postParentReplica, G4int preReplica, G4int postReplica)
{
  // G4cout << " HEADER Test15EventAction::analyseNeutronShellFluence got here 0 " << G4endl;
  // getchar();

  Test15Run* run = static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->analyseNeutronShellFluence(
    energy, time, startEnergy, TrackID, ParentID, zMomentum, startTime, radius, zPos, parentEnergy,
    parentParticle, steplength, enter_sph, enter_cyl, exit_sph, exit_cyl, Volume, enter_sph_front,
    exit_sph_front, preParentReplica, postParentReplica, preReplica, postReplica);

  // // get analysis manager
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // 27/09/15:
  /*
    number_shell_steps++;

    G4double temp_energy = energy/eV;
  //   if(enter_sph) G4cout << " ENTERING SPHERE " << G4endl;
  //   if(enter_cyl) G4cout << " ENTERING CYLINDER " << G4endl;
  //   if(exit_sph) G4cout << " EXITING SPHERE " << G4endl;
  //   if(exit_cyl) G4cout << " EXITING CYLINDER " << G4endl;

    G4double temp_time = time/microsecond;
    G4double temp_startenergy = startEnergy/eV;

    G4int iParent = 0;
    if (parentParticle == "gamma")   iParent = 1;
    if (parentParticle == "neutron") iParent = 2;
    if (parentParticle == "e-")      iParent = 3;
    if (parentParticle == "pi-")     iParent = 4;
    if (parentParticle == "pi+")     iParent = 5;
    if (parentParticle == "pi0")     iParent = 6;
    if (parentParticle == "e+")      iParent = 7;
    if (parentParticle == "proton")  iParent = 8;

    if(ntuple_full) {
      G4cout << " 4) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
      G4double break_it = std::sqrt(-100.0);
    }

    if(ntuple_full) {
      ntupleShellFluence->fill( ntupleShellFluence->findColumn( "energy" ), (G4float) temp_energy);
      ntupleShellFluence->fill( ntupleShellFluence->findColumn( "time" ), (G4float) temp_time);
      ntupleShellFluence->fill( ntupleShellFluence->findColumn( "starte" ), (G4float)
  temp_startenergy); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "trackid" ),
  (G4float) TrackID); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "parentid" ),
  (G4float) ParentID); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "zmom" ), (G4float)
  zMomentum); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "startt" ), (G4float)
  startTime); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "radius" ), (G4float)
  radius); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "e_parent" ), (G4float)
  parentEnergy); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "parent" ), (G4float)
  iParent); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "step" ), (G4float)
  steplength); ntupleShellFluence->fill( ntupleShellFluence->findColumn( "dupli" ), (G4float)
  duplicate_neutron2);

      ntupleShellFluence->addRow();
    }

      //  if(std::abs(radius-456.0) > 0.1 && std::abs(radius-50.0) > 0.1 ) return;

    G4double temp_flux_energy = -999.0;
    G4int flux_idx = -999;
    G4double temp_low_energy = -999.0;
    G4double temp_lithium_energy = -999.0;
    G4int flux_idx_low = -999;

    if(temp_energy<lithium_energy[100]) {
      for(G4int i = 0; i<100; ++i) {
        if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1])
    lithium_fluence_step_shell[i] += steplength/mm;
      }
    }

      //  if(temp_energy<flux_energy[0]) {
    if(temp_energy<low_energy[100]) {
      for(G4int i = 0; i<100; ++i) {
        if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1])
    low_fluence_step_shell[i] += steplength/mm;
      }
    }
    if(temp_energy>flux_energy[0]) {
      for(G4int i = 0; i<21; ++i) {
        //     if(energy > *(flux_energy+i) && energy < *(flux_energy+i+1)) {
        //    if(energy/eV > flux_energy[i] && energy/eV < flux_energy[i+1]) {
        //    if(energy > flux_energy[i]/1000000.0 && energy < flux_energy[i+1]/1000000.0) {
        if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1])
    fluence_step_shell[i] += steplength/mm;
      }
    }
    //     G4cout << " energy: " << energy/eV << " flux energy0: " << flux_energy[i] << " flux
  energy1: " << flux_energy[i+1] << G4endl;

    if(ntuple_full) {
      G4cout << " 5) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
      G4double break_it = std::sqrt(-100.0);
    }

    if(ntuple_full) {
      nShellFluence3->fill( nShellFluence3->findColumn( "fluxe" ), (G4float) temp_flux_energy);
      nShellFluence3->fill( nShellFluence3->findColumn( "fluxidx" ), (G4float) flux_idx);
      nShellFluence3->fill( nShellFluence3->findColumn( "lowe" ), (G4float) temp_low_energy);
      nShellFluence3->fill( nShellFluence3->findColumn( "lowidx" ), (G4float) flux_idx_low);
      nShellFluence3->fill( nShellFluence3->findColumn( "step" ), (G4float) steplength);
      nShellFluence3->addRow();
    }

  //   // integrate step lengths:
  //     if(temp_energy<lithium_energy[100]) {
  //       for(G4int i = 0; i<100; ++i) {
  // 	if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1])
  lithium_fluence_step[i] += steplength/mm;
  //       }

  //       if(temp_energy<low_energy[100]) {
  // 	for(G4int i = 0; i<100; ++i) {
  // 	  if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1]) low_fluence_step[i] +=
  steplength/mm;
  // 	}
  //       }
  //       else if(temp_energy>flux_energy[0]) {
  // 	for(G4int i = 0; i<21; ++i) {
  // 	  if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1]) fluence_step[i] +=
  steplength/mm;
  // 	}
  //       }

  //     }

  */
}

inline void Test15EventAction::analyseNeutronFlux(
  G4double energy, G4double time, G4double startEnergy, G4int TrackID, G4int ParentID,
  G4double zMomentum, G4double startTime, G4double radius, G4double zPos, G4double parentEnergy,
  G4String parentParticle, G4double cos_angle, G4int number_generations, G4String Particle,
  G4bool reduced_tally)
{
  // G4cout << " EventAction.hh, AnalyseNeutronFlux:: Got here, now do something!!! " << G4endl;
  // getchar();

  Test15Run* run = static_cast<Test15Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->analyseNeutronFlux(energy, time, startEnergy, TrackID, ParentID, zMomentum, startTime,
                          radius, zPos, parentEnergy, parentParticle, cos_angle, number_generations,
                          Particle, reduced_tally);

  // 27/09/15:
  /*
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int total_flux;
  G4double fluence_spectrum[1000];
  G4int n_max;
  G4double radii[10];
  G4double radii_energies[10];
  G4double flux_radius[10][10];
  G4double flux_energy[22];
  G4double fine_energy[43];
  G4double flux_data[21];
  G4double eflux_data[21];
  G4double flux_stat_error[21];
  G4double flux_syst_error[21];
  G4double flux[21];
  G4double cos_flux[21];
  G4double fluence[21];
  G4double fluence_step[21];
  G4double fluence_front_step[21];
  G4double fluence_cyl[21];
  G4double fluence_step_cyl[21];
  G4double fluence_step_shell[21];
  G4double eflux[21];
  G4double fine_eflux[42];
  G4double energy_integral[4];
  G4double enflux[4];
  G4int neutflux[4];

  G4int integral_scintillation;
  G4double integral_scintillation_E;
  G4int integral_lithium;
  G4double integral_lithium_E;
  G4int integral_helium;
  G4double integral_helium_E;

  G4int duplicate_neutron;
  G4int oldTrackID;
  G4int duplicate_neutron2;
  G4int oldTrackID2;

  G4double fractional_bin_width;


  //xfull  if(ntuple_full) G4cout << " entering analyseNeutronFlux " << ntuple_full << G4endl;
  //xfull  G4cout << " entering analyseNeutronFlux " << ntuple_full << G4endl;


  if(Particle == "neutron") {
    if(TrackID == oldTrackID && std::abs(radius - 456.0*mm)<0.1) {
      duplicate_neutron++;
      //x    G4cout << " duplicate_neutron " << duplicate_neutron << " trackID: " << TrackID <<
  G4endl; }  else { duplicate_neutron = 0;
    }
    oldTrackID = TrackID;
  }
  //xfull  G4cout << " inside analyseNeutronFlux 1 " << ntuple_full << G4endl;

  fractional_bin_width = 0.2;

  G4double temp_time = time/microsecond;
  G4double temp_energy = energy/eV;
  G4double temp_startenergy = startEnergy/eV;

  //    static const
  G4double radii_temp[] =
  {16.8*cm,40.4*cm,45.6*cm,69.1*cm,81.1*cm,98.6*cm,105.3*cm,113.5*cm,124.8*cm,153.9*cm}; G4double
  radii_energies_temp[] = {0.1,1.5,5.0,10.0,18.0,100.0,480.0,1000.0,10000.0,50000.0};

  for(G4int i=0; i<10; ++i) {
    radii[i] = radii_temp[i];
    radii_energies[i] = radii_energies_temp[i];
    for(G4int j=0; j<10; ++j) flux_radius[i][j] = 0;
  }

  for(G4int i=0; i<1000; ++i) fluence_spectrum[i] = 0.0;
  n_max = 0;

  if(Particle == "neutron") {
    //    G4cout << " radius is: " << radius << G4endl;

    for (G4int i=0; i<10; ++i) {
      //xfull      G4cout << " inside analyseNeutronFlux 1a, index: " << i << " ntuple_full: " <<
  ntuple_full << G4endl; if( std::abs(radius-radii[i])<0.1) { for (G4int j=0; j<10; ++j) { if (
  std::abs(energy-radii_energies[j])<fractional_bin_width*radii_energies[j]) flux_radius[i][j]
  += 1.0/std::abs(cos_angle);
  }
      }
    }

    //xfull    G4cout << " inside analyseNeutronFlux 1b " << ntuple_full << G4endl;

    if( std::abs(radius-45.6*cm)<0.1) {
      //      if(temp_energy > 0.345*eV && temp_energy < 1.e5*eV) {
      //xfull      G4cout << " inside analyseNeutronFlux 1b1 " << ntuple_full << G4endl;
      if(energy > 0.345*eV && energy < 1.e5*eV) {
  integral_scintillation += 1.0/std::abs(cos_angle);
  integral_scintillation_E += energy/eV*1.0/std::abs(cos_angle);
      }
      //xfull      G4cout << " inside analyseNeutronFlux 1b2 " << ntuple_full << G4endl;

      G4int n = (int) ((2.0 + std::log10(energy/eV))/0.09);
      //BUG BUG BUG BUG BUG!!!!
      if(n<0) n = 0;
      //xfull      G4cout << " inside analyseNeutronFlux 1b3 " << ntuple_full << G4endl;
      if(n < 100) {
  //xfull	G4cout << " inside analyseNeutronFlux 1b4 " << ntuple_full << " angle: " << cos_angle << "
  n: " << n << " fluence: " << fluence_spectrum[n] << " ntuple_full again: " << ntuple_full <<
  G4endl; fluence_spectrum[n] += 1.0/std::abs(cos_angle);
  //xfull	G4cout << " inside analyseNeutronFlux 1b5 " << ntuple_full << G4endl;
  if(n > n_max) n_max = n;
  //xfull	G4cout << " inside analyseNeutronFlux 1b6 " << ntuple_full << G4endl;
      }

      //xfull      G4cout << " inside analyseNeutronFlux 1c " << ntuple_full << G4endl;

      if(temp_energy > 0.0194 && temp_energy < 1.e5) {
  integral_lithium += 1.0/std::abs(cos_angle);
  integral_lithium_E += temp_energy*1.0/std::abs(cos_angle);
      }

      //xfull      G4cout << " inside analyseNeutronFlux 1d " << ntuple_full << G4endl;

      if(temp_energy > 59500 && temp_energy < 1825092) {
  integral_helium += 1.0/std::abs(cos_angle);
  integral_helium_E += temp_energy*1.0/std::abs(cos_angle);
      }
      //xfull      G4cout << " inside analyseNeutronFlux 1e " << ntuple_full << G4endl;
    }
  }

  //xfull  G4cout << " inside analyseNeutronFlux 2 " << ntuple_full << G4endl;

  G4int iParent = 0;
  if (!reduced_tally) iParent = -2;
  if (parentParticle == "gamma")   iParent = 1;
  if (parentParticle == "neutron") iParent = 2;
  if (parentParticle == "e-")      iParent = 3;
  if (parentParticle == "pi-")     iParent = 4;
  if (parentParticle == "pi+")     iParent = 5;
  if (parentParticle == "pi0")     iParent = 6;
  if (parentParticle == "e+")      iParent = 7;
  if (parentParticle == "proton")  iParent = 8;

  G4int iParticle = 0;
  if (Particle == "gamma")   iParticle = 1;
  if (Particle == "neutron") iParticle = 2;
  if (Particle == "e-")      iParticle = 3;
  if (Particle == "pi-")     iParticle = 4;
  if (Particle == "pi+")     iParticle = 5;
  if (Particle == "pi0")     iParticle = 6;
  if (Particle == "e+")      iParticle = 7;
  if (Particle == "proton")  iParticle = 8;

  //xfull  G4cout << " inside analyseNeutronFlux 3 " << ntuple_full << G4endl;

  if(ntuple_full) {
    G4cout << " 6) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
    G4double break_it = std::sqrt(-100.0);
  }
  if(ntuple_full) {
    ntupleFlux->fill( ntupleFlux->findColumn( "energy" ), (G4float) temp_energy);
    ntupleFlux->fill( ntupleFlux->findColumn( "time" ), (G4float) temp_time);
    ntupleFlux->fill( ntupleFlux->findColumn( "starte" ), (G4float) temp_startenergy);
    ntupleFlux->fill( ntupleFlux->findColumn( "trackid" ), (G4float) TrackID);
    ntupleFlux->fill( ntupleFlux->findColumn( "parentid" ), (G4float) ParentID);
    ntupleFlux->fill( ntupleFlux->findColumn( "zmom" ), (G4float) zMomentum);
    ntupleFlux->fill( ntupleFlux->findColumn( "startt" ), (G4float) startTime);
    ntupleFlux->fill( ntupleFlux->findColumn( "radius" ), (G4float) radius/mm);
    ntupleFlux->fill( ntupleFlux->findColumn( "e_parent" ), (G4float) parentEnergy);
    ntupleFlux->fill( ntupleFlux->findColumn( "parent" ), (G4float) iParent);
    ntupleFlux->fill( ntupleFlux->findColumn( "c_angle" ), (G4float) cos_angle);
    ntupleFlux->fill( ntupleFlux->findColumn( "numgen" ), (G4float) number_generations);
    ntupleFlux->fill( ntupleFlux->findColumn( "particle" ), (G4float) iParticle);
    ntupleFlux->fill( ntupleFlux->findColumn( "dupli" ), (G4float) duplicate_neutron);

    ntupleFlux->addRow();
  }

  //  if(std::abs(radius-456.0) > 0.1 && std::abs(radius-50.0) > 0.1 ) return;

  if(Particle == "neutron") {

    G4double temp_flux_energy = -999.0;
    G4int flux_idx = -999;
    G4double temp_low_energy = -999.0;
    G4double temp_lithium_energy = -999.0;
    G4int flux_idx_low = -999;

    if(std::abs(radius-50.0*mm) < 0.1) integral_flux_5cm++;
    if(std::abs(radius-100.0*mm) < 0.1) integral_flux_10cm++;
    if(std::abs(radius-456.0*mm) < 0.1) {
      integral_flux_46cm++;
      integral_Eflux_46cm += energy;
      if(energy > 0.1*eV && energy < 10.0*keV) neutron_fluence_46cm++;
      if(std::abs(zPos-75.0) < 15.0) integral_Zflux_46cm++;
      //     std::ofstream hitsfile("energyintegral.out", std::ios::app);
      //     hitsfile << " Hit Energy: " << G4BestUnit(energy,"Energy") << G4endl;
    }
    if(std::abs(radius-700.0*mm) < 0.1) integral_flux_70cm++;
    if(std::abs(radius-1000.0*mm) < 0.1) integral_flux_100cm++;
    if(std::abs(radius-1200.0*mm) < 0.1) integral_flux_120cm++;

    if(std::abs(radius-456.0*mm) < 0.1) {
      total_flux++;

      if(energy/eV < flux_energy[0]) energy_integral[0]+=energy/eV;

      if(temp_energy<lithium_energy[100]) {
  for(G4int i = 0; i<100; ++i) {
    if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1]) {
      lithium_flux[i]++;
      cos_lithium_flux[i] += 1.0/std::abs(cos_angle);
      integral_Eflux_46cm_restricted += energy;
      if(std::abs(zPos-75.0) < 15.0) lithium_Zflux[i]++;
      temp_lithium_energy = std::exp(0.5*(log(lithium_energy[i+1])+log(lithium_energy[i])));
    }
  }
      }

      //  if(temp_energy<flux_energy[0]) {
      if(temp_energy<low_energy[100]) {
  for(G4int i = 0; i<100; ++i) {
    if(temp_energy > low_energy[i] && temp_energy < low_energy[i+1]) {
      low_flux[i]++;
      cos_low_flux[i] += 1.0/std::abs(cos_angle);
      temp_low_energy = std::exp(0.5*(log(low_energy[i+1])+log(low_energy[i])));
      flux_idx_low = i;
    }
  }
      }
      else if(temp_energy>flux_energy[0]) {
  for(G4int i = 0; i<21; ++i) {
    //     if(energy > *(flux_energy+i) && energy < *(flux_energy+i+1)) {
    //    if(energy/eV > flux_energy[i] && energy/eV < flux_energy[i+1]) {
    //    if(energy > flux_energy[i]/1000000.0 && energy < flux_energy[i+1]/1000000.0) {
    if(temp_energy > flux_energy[i] && temp_energy < flux_energy[i+1]) {
      energy_integral[1]+=energy/eV;
      flux[i]++;
      if(cos_angle != 0.0) cos_flux[i] += 1.0/std::abs(cos_angle);
      //xx	G4cout << " i: " << i << " temp_energy: " << temp_energy << " flux_energy[i]:" <<
  flux_energy[i] << " flux_energy[i+1]: " << flux_energy[i+1] << " flux: " << flux[i] << G4endl;
      eflux[i] += energy/eV;
      temp_flux_energy = (flux_energy[i]+flux_energy[i+1])/2.;
      flux_idx = i;
    }
  }
  //     G4cout << " energy: " << energy/eV << " flux energy0: " << flux_energy[i] << " flux
  energy1: " << flux_energy[i+1] << G4endl;
      }
      for(G4int j = 0; j<42; ++j) {
  //     if(energy > *(flux_energy+i) && energy < *(flux_energy+i+1)) {
  //    if(energy/eV > flux_energy[i] && energy/eV < flux_energy[i+1]) {

  if(energy > fine_energy[j]/1000000.0 && energy < fine_energy[j+1]/1000000.0) fine_eflux[j] +=
  energy/eV;

  //     G4cout << " energy: " << energy/eV << " flux energy0: " << flux_energy[i] << " flux
  energy1: " << flux_energy[i+1] << G4endl;
      }
      //  G4cout << " Filled Source Ntuple " << G4endl;

      if(ntuple_full) {
  G4cout << " 7) WHY IS HISTOGRAM FLAG SET TO TRUE? " << G4endl;
  G4double break_it = std::sqrt(-100.0);
      }

      if(ntuple_full) {
  nFlux3->fill( nFlux3->findColumn( "fluxe" ), (G4float) temp_flux_energy);
  nFlux3->fill( nFlux3->findColumn( "fluxidx" ), (G4float) flux_idx);
  nFlux3->fill( nFlux3->findColumn( "lowe" ), (G4float) temp_low_energy);
  nFlux3->fill( nFlux3->findColumn( "lowidx" ), (G4float) flux_idx_low);
  nFlux3->addRow();
      }
      }


    if(std::abs(radius-50.0*mm) < 0.1) {

      if(temp_energy<lithium_energy[100]) {
  for(G4int i = 0; i<100; ++i) {
    if(temp_energy > lithium_energy[i] && temp_energy < lithium_energy[i+1]) lithium_flux_5cm[i]++;
  }
      }

    }

  }

*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
