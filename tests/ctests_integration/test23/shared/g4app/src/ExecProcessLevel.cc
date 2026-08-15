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

#include "ExecProcessLevel.hh"

#include "G4DynamicParticle.hh"
#include "G4ElementVector.hh"
#include "G4IsotopeVector.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "globals.hh"

// --> obsolete --> #include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"
// --> obsolete --> #include "G4HadronInelasticDataSet.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4Box.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4IonConstructor.hh"
#include "G4IonTable.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleTable.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "G4Region.hh"
#include "G4Step.hh"
#include "G4TrackingManager.hh"

#include "ProcessWrapper.hh"
#include "TstDiscreteProcessReader.hh"
#include "TstTarget.hh"

#ifdef G4_USE_FLUKA
// interface to FLUKA.CERN, if specified
#  include "FLUKAInelasticScatteringXS.hh"
#endif

ExecProcessLevel::ExecProcessLevel()
  : ExecBase(),
    fProcWrapper(0),
    fTarget(0),
    fRegion(0),
    fProcManager(0),
    fBeam(0),
    fTrack(0),
    fStep(0),
    fPartChange(0)
{}

ExecProcessLevel::~ExecProcessLevel()
{
  if (fTarget)
  {
    delete fTarget;
    fTarget = 0;
  }
  if (fRegion)
  {
    delete fRegion;
    fRegion = 0;
  }
  if (fProcWrapper)
  {
    delete fProcWrapper;
    fProcWrapper = 0;
  }
  if (fBeam)
  {
    delete fBeam;
    fBeam = 0;
  }
  if (fTrack)
  {
    delete fTrack;
    fTrack = 0;
  }
  if (fStep) delete fStep;
}

void ExecProcessLevel::InitSetup(const TstReader* pset)
{
  fTarget = new TstTarget();

  // in principle, the return type here is G4Material*
  // but the compiler was complaining about unsed variable
  //
  // G4Material* mat = fTarget->ResetMaterial( pset->GetTargetMaterial() );
  fTarget->ResetMaterial(pset->GetTargetMaterial());

  G4ThreeVector targetSize = pset->GetTargetSize();
  fTarget->SetDimentions(targetSize.x(), targetSize.y(), targetSize.z());

  fTarget->ResetGeom();
  G4VPhysicalVolume* targetPhys = fTarget->ConstructTarget();

  fRegion = new G4Region("Region");  // needed by tracking manager
  targetPhys->GetLogicalVolume()->SetRegion(fRegion);
  fRegion->AddRootLogicalVolume(targetPhys->GetLogicalVolume());

  return;
}

void ExecProcessLevel::InitBeam(const TstReader* pset)
{
  assert(fProcWrapper);

  assert(fTarget->GetCurrentMaterial());

  // not a very good design...
  // need to see if it can be refactored into subclasses
  //
  if (pset->GetBeamParticle() == "proton")
  {
    fProcManager = new G4ProcessManager(G4Proton::Proton());
  }
  else if (pset->GetBeamParticle() == "pi-")
  {
    fProcManager = new G4ProcessManager(G4PionMinus::PionMinus());
  }
  else if (pset->GetBeamParticle() == "pi+")
  {
    fProcManager = new G4ProcessManager(G4PionPlus::PionPlus());
  }
  else if (pset->GetBeamParticle() == "gamma")
  {
    fProcManager = new G4ProcessManager(G4Gamma::Definition());
  }
  else if (pset->GetBeamParticle() == "kaon-")
  {
    fProcManager = new G4ProcessManager(G4KaonMinus::KaonMinus());
  }
  else if (pset->GetBeamParticle() == "kaon+")
  {
    fProcManager = new G4ProcessManager(G4KaonPlus::KaonPlus());
  }
  // need to add also omega, sigma, mu-, pbar !!!!!

  assert(fProcManager);

  if (pset->IsDiscreteProcess())
  {
    fProcManager->AddDiscreteProcess(fProcWrapper);
  }

  G4ParticleDefinition* partDef =
    (G4ParticleTable::GetParticleTable())->FindParticle(pset->GetBeamParticle());
  G4double partMass = partDef->GetPDGMass();
  G4double partMom = pset->GetBeamMomentum();
  G4double partEnergy = std::sqrt(partMom * partMom + partMass * partMass);  // total energy

  if (!fBeam) fBeam = new Beam();
  fBeam->SetBeam(pset->GetBeamParticle(), partMass, partEnergy);

  // new addition for restructuring test47
  //
  fBeam->SetBeamPosition(pset->GetPosition());
  fBeam->SetBeamMomentum(partMom * pset->GetDirection());
  if (pset->IsBeamKineByEvt())
  {
    fBeam->SetAngleRange(pset->GetAngleX(), pset->GetdAngleX(), pset->GetAngleY(),
                         pset->GetdAngleY());
    fBeam->SetdpByp(pset->GetdpByp());
    fBeam->SetPosSigma(pset->GetBeamSigmaX(), pset->GetBeamSigmaY());
    fBeam->BeamKineByEvtOn();
  }
  // end new stuff for restructuring test47

  const G4Element* elm = fTarget->GetCurrentMaterial()->GetElement(0);
  G4int A = (G4int)(elm->GetN() + 0.5);
  G4int Z = (G4int)(elm->GetZ() + 0.5);
  G4double amass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);

  G4DynamicParticle dParticle(partDef, pset->GetDirection(),
                              partEnergy - partMass);  // this has to be EKin, not total

  fBeam->SetLabV(G4LorentzVector(0., 0.,
                                 std::sqrt(partEnergy * (partEnergy + 2.0 * partMass)) / GeV,
                                 (partEnergy + partMass + amass) / GeV));

  // if under an agnle, then like this:
  // labv = G4LorentzVector(mom.x()/CLHEP::GeV, mom.y()/CLHEP::GeV,
  //			      mom.z()/CLHEP::GeV, (e0+mass+amass)/CLHEP::GeV);

  fBeam->SetLabP(G4LorentzVector(0., 0.,
                                 std::sqrt(partEnergy * (partEnergy + 2.0 * partMass)) / GeV,
                                 (partEnergy + partMass + G4Proton::Proton()->GetPDGMass()) / GeV));

  // Track
  fTrack = new G4Track(
    new G4DynamicParticle(partDef, pset->GetDirection(),
                          partEnergy - partMass),  // the 3rd arg of G4DynPart has to be EKin
    pset->GetTime(), pset->GetPosition());
  fTrack->SetTouchableHandle(G4TouchableHandle(new G4TouchableHistory()));

  // Step
  fStep = new G4Step();
  fStep->SetTrack(fTrack);
  fTrack->SetStep(fStep);

  // step points
  //
  // pre-step
  //
  fStep->SetPreStepPoint(new G4StepPoint());
  fStep->GetPreStepPoint()->SetPosition(pset->GetPosition());
  fStep->GetPreStepPoint()->SetMaterial(fTarget->GetCurrentMaterial());
  fStep->GetPreStepPoint()->SetSafety(10000. * CLHEP::cm);
  //
  // post-step
  //
  // NOTE-1 (JVY): this (commented) fragment below has been copied long ago
  //               from even older code (perhaps originates in test35);
  //               I do NOT understand why step-points were done this way
  //               instead of making them 2 separate instances (objects)...
  //               (as it is now)
  //    G4StepPoint* bPoint;
  //    bPoint = fStepPoint; // ...=aPoint, as it was in the original code
  //    G4ThreeVector bPosition = pset->GetDirection() * pset->GetStep();
  //    bPosition += pset->GetPosition();
  //    bPoint->SetPosition(bPosition);
  //    fStep->SetPostStepPoint( bPoint );
  //
  G4ThreeVector bPosition = pset->GetDirection() * pset->GetStep();
  bPosition += pset->GetPosition();
  fStep->SetPostStepPoint(new G4StepPoint());
  fStep->GetPostStepPoint()->SetPosition(bPosition);
  fStep->GetPostStepPoint()->SetMaterial(fTarget->GetCurrentMaterial());
  fStep->GetPostStepPoint()->SetSafety(10000. * CLHEP::cm);
  //
  fStep->SetStepLength(pset->GetStep());

  // actually, it needs to be a bit more extensive...
  //

  G4VCrossSectionDataSet* cs = 0;

  if (pset->GetBeamParticle() == "gamma")
  {
    // special case:
    // will calculate gamma-N xsec later in the process,
    // via G4PhotoNuclearCrossSection::GetElementCrossSection(...)
    return;
  }

  std::string physics_lc = G4StrUtil::to_lower_copy(pset->GetPhysics());
  std::cout << " GetPhysics = " << pset->GetPhysics() << " ---> physics_lc = " << physics_lc
            << std::endl;
  if (physics_lc.find("fluka") == std::string::npos)
  {
    if (partDef == G4Proton::Definition() || partDef == G4Neutron::Definition())
    {
      cs = new G4BGGNucleonInelasticXS(partDef);
    }
    else if (partDef == G4PionPlus::Definition() || partDef == G4PionMinus::Definition())
    {
      cs = new G4BGGPionInelasticXS(partDef);
    }
    else if (partDef->GetBaryonNumber() > 1)
    {  // Ions
      cs = new G4CrossSectionInelastic(new G4ComponentGGNuclNuclXsc);
    }
    else if (partDef == G4AntiProton::Definition() || partDef == G4AntiNeutron::Definition()
             || partDef == G4AntiDeuteron::Definition() || partDef == G4AntiTriton::Definition()
             || partDef == G4AntiHe3::Definition() || partDef == G4AntiAlpha::Definition())
    {
      cs = new G4CrossSectionInelastic(new G4ComponentAntiNuclNuclearXS);
    }
    else
    {
      cs = new G4CrossSectionInelastic(new G4ComponentGGHadronNucleusXsc);
    }
  }
#ifdef G4_USE_FLUKA
  else
  {
    cs = new FLUKAInelasticScatteringXS();
  }
#endif

  if (cs)
  {
    cs->BuildPhysicsTable(*partDef);
    fXSecOnTarget = cs->GetCrossSection(&dParticle, elm);
    fProcWrapper->AddDataSet(cs);
    fProcWrapper->PreparePhysicsTable(*partDef);
  }

  return;
}

G4VParticleChange* ExecProcessLevel::DoEvent()
{
  // related to restructuring test47
  //
  if (fBeam->IsBeamKineByEvt()) GenerateBeamKine();

  fTrack->SetKineticEnergy(fBeam->GetBeamEnergy() - fBeam->GetBeamPartMass());

  fPartChange = fProcWrapper->PostStepDoIt(*fTrack, *fStep);

  return fPartChange;
}
// new addition - mainly for test47 restructuring
//
void ExecProcessLevel::GenerateBeamKine()
{
  G4double angleX = G4RandGauss::shoot(fBeam->GetAngleX(), fBeam->GetdAngleX());
  G4double angleY = G4RandGauss::shoot(fBeam->GetAngleY(), fBeam->GetdAngleY());
  G4double zv = fTarget->GetZLength() * (G4UniformRand() - 0.5);
  G4double xv = fTarget->GetRMax();
  G4double yv = xv;
  if (fTarget->GetRMax() > 0)
  {
    while (std::sqrt(xv * xv + yv * yv) > fTarget->GetRMax())
    {
      xv = G4RandGauss::shoot(fBeam->GetBeamPosition().x(), fBeam->GetBeamSigmaX());
      yv = G4RandGauss::shoot(fBeam->GetBeamPosition().y(), fBeam->GetBeamSigmaY());
    }
  }

  G4double pmom = (fBeam->GetBeamMomentum().mag()) * G4RandGauss::shoot(1.0, fBeam->GetdpByp());
  G4ThreeVector newPosition = fBeam->GetBeamPosition() + G4ThreeVector(xv, yv, zv);
  G4ThreeVector aDirection = G4ThreeVector(-std::cos(angleX) * std::sin(angleY), std::sin(angleX),
                                           std::cos(angleX) * std::cos(angleY));

  G4double mass = fBeam->GetBeamPartMass();
  G4double energy = sqrt(pmom * pmom + mass * mass);
  fBeam->ResetBeamEnergy(energy);
  const G4ThreeVector& mom = pmom * aDirection;

  fStep->GetPreStepPoint()->SetPosition(newPosition);
  fStep->GetPostStepPoint()->SetPosition((newPosition + aDirection * fStep->GetStepLength()));

  // FIXME !!!
  // I hate thicks like that !!!
  // Need to find a better way... than casting away const...
  //
  ((G4DynamicParticle*)(fTrack->GetDynamicParticle()))->SetMomentum(mom);

  fBeam->SetLabV(G4LorentzVector(mom.x() / GeV, mom.y() / GeV, mom.z() / GeV,
                                 (energy + mass + fTarget->GetAMass()) / GeV));  // need amass !!!

  fBeam->SetLabP(G4LorentzVector(mom.x() / GeV, mom.y() / GeV, mom.z() / GeV,
                                 (energy + mass + G4Proton::Proton()->GetPDGMass()) / GeV));

  // FIXME !!!
  // Technically speaking, XSec on target also needs to be redefined because pmom changes.
  // But it's likely to be a minor thing, so let's leave it as-is for now...

  return;
}
