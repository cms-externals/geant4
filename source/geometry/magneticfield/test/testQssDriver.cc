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
// $Id: $
//
// Started from testTemplatedStepper, an adapted version of testPropagateMagField
//  First version:  Josh Xie (Google Summer of Code 2014)
//  Supervisors:    Sandro Wenzel, John Apostolakis
//
//  Refinements/improvements:  John Apostolakis  October 2021

#include <cstdio>
#include <limits>

// #include "ApproxEqual.hh"

#include "G4FieldTrack.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

// Global defs
// #include "globals.hh"
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ios.hh"

#include <cassert>
#include <ctime>
#include <iomanip>

// Sample Parameterisation
class G4LinScale : public G4VPVParameterisation
{
    virtual void ComputeTransformation(const G4int n, G4VPhysicalVolume* pRep) const
    {
      pRep->SetTranslation(G4ThreeVector(0, (n - 1) * 15, 0));
    }

    virtual void ComputeDimensions(G4Box& pBox, const G4int n, const G4VPhysicalVolume*) const
    {
      pBox.SetXHalfLength(10);
      pBox.SetYHalfLength(5 + n);
      pBox.SetZHalfLength(5 + n);
    }

    virtual void ComputeDimensions(G4Tubs&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Trd&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Cons&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Trap&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Hype&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Orb&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Sphere&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Torus&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Para&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Polycone&, const G4int, const G4VPhysicalVolume*) const {}
    virtual void ComputeDimensions(G4Polyhedra&, const G4int, const G4VPhysicalVolume*) const {}
};
G4LinScale myParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{
  G4Box* myHugeBox = new G4Box("huge box", 15 * m, 15 * m, 25 * m);
  G4Box* myBigBox = new G4Box("big cube", 10 * m, 10 * m, 10 * m);
  G4Box* mySmallBox = new G4Box("smaller cube", 2.5 * m, 2.5 * m, 2.5 * m);
  G4Box* myTinyBox = new G4Box("tiny  cube", .25 * m, .25 * m, .25 * m);

  // G4Box *myVariableBox=
  new G4Box("Variable Box", 10, 5, 5);

  //  World Volume
  //
  G4LogicalVolume* worldLog = new G4LogicalVolume(myHugeBox, 0, "World", 0, 0, 0);
  // Logical with no material,field,
  // sensitive detector or user limits

  G4PVPlacement* worldPhys =
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "World", worldLog, 0, false, 0);
  // Note: no mother pointer set

  //  Create the logical Volumes
  //
  //  G4LogicalVolume(*pSolid, *pMaterial, Name, *pField, *pSDetector, *pULimits)
  //
  G4LogicalVolume* BigBoxLog = new G4LogicalVolume(myBigBox, 0, "Crystal Box (large)", 0, 0, 0);
  G4LogicalVolume* smallBoxLog = new G4LogicalVolume(mySmallBox, 0, "Crystal Box (small)");
  G4LogicalVolume* tinyBoxLog = new G4LogicalVolume(myTinyBox, 0, "Crystal Box (tiny)");

  //  Place them.
  //
  //  1) Two big boxes in the world volume
  //
  // G4PVPlacement *BigTg1Phys=
  new G4PVPlacement(0, G4ThreeVector(0, 0, -15 * m), "Big Target 1", BigBoxLog, worldPhys, false,
                    0);
  // G4PVPlacement *BigTg2Phys=
  new G4PVPlacement(0, G4ThreeVector(0, 0, 15 * m), "Big Target 2", BigBoxLog, worldPhys, false, 0);

  //  2) Four (medium) boxes in X & Y near the origin of the world volume
  //
  // G4PVPlacement *MedTg3a_Phys=
  new G4PVPlacement(0, G4ThreeVector(0, 7.5 * m, 0), "Target 3a", smallBoxLog, worldPhys, false, 0);
  // G4PVPlacement *MedTg3b_Phys=
  new G4PVPlacement(0, G4ThreeVector(0, -7.5 * m, 0), "Target 3b", smallBoxLog, worldPhys, false,
                    0);
  // G4PVPlacement *MedTg3c_Phys=
  new G4PVPlacement(0, G4ThreeVector(-7.5 * m, 0, 0), "Target 3c", smallBoxLog, worldPhys, false,
                    0);
  // G4PVPlacement *MedTg3d_Phys=
  new G4PVPlacement(0, G4ThreeVector(7.5 * m, 0, 0), "Target 3d", smallBoxLog, worldPhys, false, 0);

  //  3) Eight small boxes around the origin of the world volume
  //        (in +-X, +-Y & +-Z)
  //
  // G4PVPlacement *SmTg4a_Phys=
  new G4PVPlacement(0, G4ThreeVector(0.3 * m, 0.3 * m, 0.3 * m), "Target 4a", tinyBoxLog, worldPhys,
                    false, 0);
  // G4PVPlacement *SmTg4b_Phys=
  new G4PVPlacement(0, G4ThreeVector(0.3 * m, -0.3 * m, 0.3 * m), "Target 4b", tinyBoxLog,
                    worldPhys, false, 0);
  // G4PVPlacement *SmTg4c_Phys=
  new G4PVPlacement(0, G4ThreeVector(-0.3 * m, -0.3 * m, 0.3 * m), "Target 4c", tinyBoxLog,
                    worldPhys, false, 0);
  // G4PVPlacement *SmTg4d_Phys=
  new G4PVPlacement(0, G4ThreeVector(-0.3 * m, 0.3 * m, 0.3 * m), "Target 4d", tinyBoxLog,
                    worldPhys, false, 0);

  // G4PVPlacement *SmTg4e_Phys=
  new G4PVPlacement(0, G4ThreeVector(0.3 * m, 0.3 * m, -0.3 * m), "Target 4e", tinyBoxLog,
                    worldPhys, false, 0);
  // G4PVPlacement *SmTg4f_Phys=
  new G4PVPlacement(0, G4ThreeVector(0.3 * m, -0.3 * m, -0.3 * m), "Target 4f", tinyBoxLog,
                    worldPhys, false, 0);
  // G4PVPlacement *SmTg4g_Phys=
  new G4PVPlacement(0, G4ThreeVector(-0.3 * m, -0.3 * m, -0.3 * m), "Target 4g", tinyBoxLog,
                    worldPhys, false, 0);
  // G4PVPlacement *SmTg4h_Phys=
  new G4PVPlacement(0, G4ThreeVector(-0.3 * m, 0.3 * m, -0.3 * m), "Target 4h", tinyBoxLog,
                    worldPhys, false, 0);

  return worldPhys;
}

#include "G4CachedMagneticField.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagneticField.hh"
#include "G4PropagatorInField.hh"
#include "G4QSSDriver.hh"
#include "G4QSStepper.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4TCachedMagneticField.hh"
#include "G4TMagFieldEquation.hh"
#include "G4TQuadrupoleMagField.hh"
#include "G4TransportationManager.hh"
#include "G4UniformMagField.hh"

// typedef G4CachedMagneticField Field_t;
// typedef G4TCachedMagneticField<G4QuadrupoleMagField> Field_t;
using Field_t = G4TCachedMagneticField<G4TQuadrupoleMagField>;
using Equation_t = G4TMagFieldEquation<Field_t>;
// Experimenting to see whether G4TCashKarpRKF45 can be used to implement the old G4CashKarpRKF45

G4TQuadrupoleMagField tQuadrupoleMagField(10. * tesla / (50. * cm));
// G4QuadrupoleMagField   tQuadrupoleMagField( 10.*tesla/(50.*cm) );
Field_t tMagField(&tQuadrupoleMagField, 1.0 * cm);
//===============================================

// G4UniformMagField      uniformMagField(10.*tesla, 0., 0.);
//  G4CachedMagneticField  myMagField( &uniformMagField, 1.0 * cm);
//  G4String   fieldName("Uniform 10Tesla");

G4double distConstField =  // 1.0 * CLHEP::cm;
  1.0e-2 * CLHEP::millimeter;

G4QuadrupoleMagField quadrupoleMagField(10. * tesla / (50. * cm));
G4CachedMagneticField myMagField(&quadrupoleMagField, distConstField);
G4String fieldName("Cached Quadropole field, 20T/meter, cache_dist="
                   + std::to_string(distConstField / CLHEP::mm) + " mm");

G4bool verbose = 1;

G4FieldManager* SetupField(G4int qssOrder)
{
  G4FieldManager* pFieldMgr;
  G4ChordFinder* pChordFinder;
  G4Mag_UsualEqRhs* pEquation = new G4Mag_UsualEqRhs(&tMagField);
  //=============test template mode================
  // Equation_t *templEquation = new Equation_t(&tMagField);
  //===============================================

  G4cout << " Setting up field (name) = " << fieldName << "  Qss integration Order = " << qssOrder
         << G4endl;

  // G4bool ok=
  G4QSSMessenger::instance()->SetQssOrder(qssOrder);

  auto qssStepper = new G4QSStepper(pEquation, 6, qssOrder);
  auto qss_driver = new G4QSSDriver<G4QSStepper>(qssStepper);

  // QSS-3 is available as part of QSS-v2 (June 2025)
  // auto qssStepper3 = G4QSSDriverCreator::CreateQss3Stepper(pEquation);

  pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  pFieldMgr->SetDetectorField(&tMagField);

  tMagField.ClearCounts();

  pChordFinder = new G4ChordFinder(qss_driver);
  pChordFinder->SetVerbose(0);

  pFieldMgr->SetChordFinder(pChordFinder);

  return pFieldMgr;
}

#include "G4BrentLocator.hh"
#include "G4MultiLevelLocator.hh"
#include "G4SimpleLocator.hh"

G4PropagatorInField* SetupPropagator(G4int type)
{
  // G4FieldManager* fieldMgr=
  SetupField(type);

  // G4ChordFinder  theChordFinder( &MagField, 0.05*mm ); // Default stepper

  G4PropagatorInField* thePropagator =
    G4TransportationManager::GetTransportationManager()->GetPropagatorInField();

  // Let us test the new Minimum Epsilon Step functionality
  // thePropagator -> SetMinimumEpsilonStep( 1.0e-3 ) ;
  // thePropagator -> SetMaximumEpsilonStep( 1.0e-5 ) ;

  G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // Test the options for Locator
  G4VIntersectionLocator* pLocator = 0;
  // G4cout << "Over-riding  PropagatorInField to use ";
  pLocator = new G4MultiLevelLocator(theNavigator);  // G4cout << "Multi"; // default
  // pLocator= new G4SimpleLocator(theNavigator); //G4cout << "Simple";
  // pLocator= new G4BrentLocator(theNavigator); //G4cout << " Brent ";
  // G4cout << " Locator. ( In the unit test code. ) " << G4endl;

  thePropagator->SetIntersectionLocator(pLocator);

  return thePropagator;
}

G4PropagatorInField* pMagFieldPropagator = 0;

#define CHECK_RESULTS 1

//
// Test Stepping
//
G4bool testG4PropagatorInField(G4VPhysicalVolume*,  // *pTopNode,
                               G4int type)
{
  void report_endPV(G4ThreeVector Position, G4ThreeVector UnitVelocity, G4double step_len,
                    G4double physStep, G4double safety, G4ThreeVector EndPosition,
                    G4ThreeVector EndUnitVelocity, G4int Step, G4VPhysicalVolume * startVolume);

  G4UniformMagField MagField(10. * tesla, 0., 0.);
  G4Navigator* pNavig =
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  pMagFieldPropagator = SetupPropagator(type);

  G4double particleCharge = +1.0;  // in e+ units
  G4double spin = 0.0;  // ignore the spin
  G4double magneticMoment = 0.0;  // ignore the magnetic moment

  G4ChargeState chargeState(particleCharge,  // The charge can change (dynamic)
                            spin = 0.0, magneticMoment = 0.0);

  G4ChordFinder* pChordFndr = pMagFieldPropagator->GetChordFinder();

  G4EquationOfMotion* equationOfMotion = pChordFndr->GetIntegrationDriver()->GetEquationOfMotion();

  equationOfMotion->SetChargeMomentumMass(chargeState,
                                          0.5 * proton_mass_c2,  // Momentum in Mev/c
                                          proton_mass_c2);
  // pNavig->SetWorldVolume(pTopNode);

  G4VPhysicalVolume* located;
  G4double step_len, physStep, safety;
  G4ThreeVector xHat(1, 0, 0), yHat(0, 1, 0), zHat(0, 0, 1);
  G4ThreeVector mxHat(-1, 0, 0), myHat(0, -1, 0), mzHat(0, 0, -1);

  // physStep=kInfinity;
  G4ThreeVector Position(0., 0., 0.);
  G4ThreeVector UnitMomentum(0., 0.6, 0.8);
  G4ThreeVector EndPosition, EndUnitMomentum, EndMomentum;

  //
  // Test location & Step computation
  //
  /* assert(located->GetName()=="World"); */
  if (std::fabs(UnitMomentum.mag() - 1.0) > 1.e-8)
  {
    G4cerr << "UnitMomentum.mag() - 1.0 = " << UnitMomentum.mag() - 1.0 << G4endl;
  }

  // How to ensure that integrDriver->OnStartTracking(); is called ??  2022.10.20
  // pMagFieldpropagator->ClearPropagatorState();  // ???

  // G4cout << G4endl;

  for (int iparticle = 0; iparticle < 1; iparticle++)
  {
    physStep = 2.5 * mm;  // millimeters
    Position = G4ThreeVector(0., 0., 0.) + iparticle * G4ThreeVector(0.2, 0.3, 0.4);
    UnitMomentum =
      (G4ThreeVector(0., 0.6, 0.8) + (float)iparticle * G4ThreeVector(0.1, 0.2, 0.3)).unit();

    G4double momentum = (0.5 + iparticle * 10.0) * proton_mass_c2;

    G4double kineticEnergy =
      momentum * momentum
      / (std::sqrt(momentum * momentum + proton_mass_c2 * proton_mass_c2) + proton_mass_c2);
    G4double velocity = momentum / (proton_mass_c2 + kineticEnergy);
    G4double labTof = 10.0 * ns, properTof = 0.1 * ns;
    G4ThreeVector Spin(1.0, 0.0, 0.0);
    // Momentum in Mev/c ?
    // pMagFieldPropagator
    G4ChargeState charge(+1.0, 0.0, 0.0);  // electrical charge in e+ units, mdm, pdgSpin
    equationOfMotion->SetChargeMomentumMass(chargeState, momentum, proton_mass_c2);
    if (verbose)
    {
      G4cout << G4endl;
      G4cout << "Test PropagateMagField: ***********************" << G4endl
             << " New Particle: position= " << Position << " direction= " << UnitMomentum
             << " Momentum " << momentum / GeV << " GeV/c " << G4endl;
    }
    clock_t total = 0;
    for (int istep = 0; istep < 14; istep++)
    {
      // G4cerr << "UnitMomentum Magnitude is " << UnitMomentum.mag() << G4endl;
      located = pNavig->LocateGlobalPointAndSetup(Position);
      // G4cerr << "Starting Step " << istep << " in volume "
      // << located->GetName() << G4endl;

      G4FieldTrack initTrack(Position, UnitMomentum,
                             0.0,  // starting S curve len
                             kineticEnergy, proton_mass_c2, velocity, labTof, properTof,
                             0  // or &Spin
      );
      initTrack.SetChargeAndMoments(particleCharge);  //  Critical !! Else not set

      clock_t t;
      t = clock();

      step_len = pMagFieldPropagator->ComputeStep(initTrack, physStep, safety, located);
      total += clock() - t;
      //       --------------------
      EndPosition = pMagFieldPropagator->EndPosition();
      EndUnitMomentum = pMagFieldPropagator->EndMomentumDir();
      EndMomentum = pMagFieldPropagator->GetEndState().GetMomentum();
      //       --------

      if (std::fabs(EndUnitMomentum.mag2() - 1.0) > 1.e-8)
        G4cerr << "EndUnitMomentum.mag2() - 1.0 = " << EndUnitMomentum.mag2() - 1.0 << G4endl;

      G4ThreeVector MoveVec = EndPosition - Position;
      // assert( MoveVec.mag() < physStep*(1.+1.e-9) );
#ifdef CHECK_RESULTS
      const double maxMoveRatio = 1.0e-9;
      double epsMomentum = 1.0e-9;

      double momentumMagStart = initTrack.GetMomentum().mag();
      double momentumMagEnd = EndMomentum.mag();
      double diffMomMag = momentumMagEnd - momentumMagStart;
      double averMomMag = 0.5 * (momentumMagEnd - momentumMagStart);

      double move = MoveVec.mag();
      double displaceRatio = 0.0;
      if (physStep > 0.)
      {
        displaceRatio = move / physStep - 1.0;
      }
      if (displaceRatio > maxMoveRatio)
      {  // Was 1.0e-9
        int eprec = G4cerr.precision(10);
        G4ThreeVector diffUnitMom = EndUnitMomentum - UnitMomentum;
        G4double dotDir = UnitMomentum.dot(EndUnitMomentum);
        G4cerr << " Move displaced " << std::setw(13) << move
               << " - further than step= " << std::setw(8) << physStep << " (achieved= " << step_len
               << " ) "
               << "  by  " << std::setw(9) << std::setprecision(4) << move - physStep << ", "
               << " a fraction of " << displaceRatio << "  "
               << " |D dir|= " << diffUnitMom.mag() << " "
               << " p-dir: "
               << " dot-1 = " << dotDir - 1.0 << " angle= "
               << std::acos(dotDir)
               // << " start= " << std::setw(12) << std::setprecision(9) << UnitMomentum << "  "
               // << " end = "  << std::setw(12) << std::setprecision(9) << EndUnitMomentum << " "
               << " diff= " << diffUnitMom << G4endl;
        G4cerr.precision(eprec);
      }
      G4double maxEpsilon = pMagFieldPropagator->GetMaximumEpsilonStep();
      assert(displaceRatio < maxEpsilon);  // Revised check condition
      if (displaceRatio >= maxEpsilon)
      {
        G4cerr << "ERROR> Problem - displacement ratio " << displaceRatio
               << " exceeds maximum allowed value = " << maxEpsilon << G4endl;
        exit(1);
      }

      if (std::fabs(diffMomMag) > epsMomentum * averMomMag)
      {
        int eprec = G4cerr.precision(10);
        G4cerr << " Momentum magnitude change:  relative: " << std::setw(9) << std::setprecision(5)
               << diffMomMag / averMomMag << " absolute: " << std::setw(12) << std::setprecision(9)
               << diffMomMag << " versus start " << std::setw(12) << std::setprecision(9)
               << momentumMagStart << " and end "
               << " versus end " << std::setw(12) << std::setprecision(9) << momentumMagEnd
               << G4endl;
        G4cerr.precision(eprec);
      }

      if (verbose)
      {
        // G4cout << " testPropagatorInField: After stepI " << istep  << " : " << G4endl;
        report_endPV(Position, UnitMomentum, step_len, physStep, safety, EndPosition,
                     EndUnitMomentum, istep, located);
      }
#endif
      assert(safety >= 0);
      pNavig->SetGeometricallyLimitedStep();
      // pMagFieldPropagator->SetGeometricallyLimitedStep();

      Position = EndPosition;
      UnitMomentum = EndUnitMomentum;

      if (step_len == physStep)
      {
        physStep *= 2.;
      }
    }  // ...........................  end for ( istep )
    G4cout << "============== time= " << total << " =================" << G4endl;
    tMagField.ReportStatistics();

  }  // ..............................  end for ( iparticle )

  return (1);
}

void report_endPV(G4ThreeVector Position, G4ThreeVector InitialUnitVelocity, G4double step_len,
                  G4double physStep, G4double safety, G4ThreeVector EndPosition,
                  G4ThreeVector EndUnitVelocity, G4int Step, G4VPhysicalVolume* startVolume)
//   G4VPhysicalVolume* endVolume)
{
  const G4int verboseLevel = 1;

  if (Step == 0 && verboseLevel <= 3)
  {
    G4cout.precision(6);
    // G4cout.setf(ios_base::fixed,ios_base::floatfield);
    G4cout << std::setw(3) << "Step " << std::setw(9) << "X(mm)" << " " << std::setw(9) << "Y(mm)"
           << " " << std::setw(9) << "Z(mm)" << " " << std::setw(9) << " N_x " << " "
           << std::setw(9) << " N_y " << " " << std::setw(9) << " N_z " << " " << std::setw(9)
           << " Delta|N|" << " " << std::setw(9) << " Delta(N_z) "
           << " "
           // << std::setw( 9) << "KinE(MeV)" << " "
           << std::setw(7) << "StepLen" << " " << std::setw(8) << "PhsStep" << " " << std::setw(8)
           << "Safety" << " " << std::setw(11) << "Field:Calls" << " " << std::setw(7) << " Evals"
           << " " << std::setw(18) << "NextVolume" << " " << G4endl;
  }
  //
  //
  if (verboseLevel > 3)
  {
    G4cout << "End  Position is " << EndPosition << G4endl << " and UnitVelocity is "
           << EndUnitVelocity << G4endl;
    G4cout << "Step taken was " << step_len << " out of PhysicalStep= " << physStep << G4endl;
    G4cout << "Final safety is: " << safety << G4endl;

    G4cout << "Chord length = " << (EndPosition - Position).mag() << G4endl;
    G4cout << G4endl;
  }
  else  // if( verboseLevel > 0 )
  {
    G4cout.precision(6);
    G4cout << std::setw(3) << Step << " " << std::setw(9) << Position.x() << " " << std::setw(9)
           << Position.y() << " " << std::setw(9) << Position.z() << " " << std::setw(9)
           << EndUnitVelocity.x() << " " << std::setw(9) << EndUnitVelocity.y() << " "
           << std::setw(9) << EndUnitVelocity.z() << " ";
    G4cout.precision(2);
    G4cout << std::setw(9) << EndUnitVelocity.mag() - InitialUnitVelocity.mag() << " "
           << std::setw(9) << EndUnitVelocity.z() - InitialUnitVelocity.z() << " ";
    //    << std::setw( 9) << KineticEnergy << " "
    //    << std::setw( 9) << EnergyDifference << " "
    G4cout.precision(6);
    G4cout << std::setw(9) << step_len << " " << std::setw(9) << physStep << " " << std::setw(9)
           << safety << " ";
    G4cout << std::setw(6) << tMagField.GetCountCalls() << " " << std::setw(5)
           << tMagField /*myMagField*/.GetCountEvaluations() << " ";

    if (startVolume != 0)
    {
      G4cout << std::setw(12) << startVolume->GetName() << " ";
    }
    else
    {
      G4cout << std::setw(12) << "OutOfWorld" << " ";
    }
#if 0
        if( endVolume != 0) 
            G4cout << std::setw(12) << endVolume()->GetName() << " ";
        else 
            G4cout << std::setw(12) << "OutOfWorld" << " ";
#endif
    G4cout << G4endl;
  }
}
// Main program
// -------------------------------
int main(int argc, char** argv)
{
  G4VPhysicalVolume* myTopNode;
  G4int qssOrder, optim, optimSaf;
  G4bool optimiseVoxels = true;
  G4bool optimisePiFwithSafety = true;

  // G4cout << " Arguments:  stepper-no  optimise-Voxels optimise-PiF-with-safety" << G4endl;

  if (argc >= 2)
  {
    qssOrder = atoi(argv[1]);
  }
  else
  {
    qssOrder = 2;  // QSS-2 only one currently available
  }

  if (argc >= 3)
  {
    optim = atoi(argv[2]);
    if (optim == 0)
    {
      optimiseVoxels = false;
    }
  }

  if (argc >= 4)
  {
    optimSaf = atoi(argv[3]);
    if (optimSaf == 0)
    {
      optimisePiFwithSafety = false;
    }
  }

  int len = 1;
  // G4cout << " Testing QSS with order " << qssOrder << G4endl;

  for (int k = 0; k < len; k++)
  {
    // G4cout << "             " ;
    // G4cout << " voxel optimisation      " << (optimiseVoxels ? "On" : "Off" )
    // G4cout << G4endl << "             " ;
    // G4cout << " Propagator safety optim " ;
    //  const char* OnOff= (optimisePiFwithSafety ? "on" : "off") ;
    //  //G4cout << OnOff << G4endl;
    // G4cout << (optimisePiFwithSafety ? "On" : "Off")  << G4endl;

    // Create the geometry & field
    myTopNode = BuildGeometry();  // Build the geometry

    G4Navigator* pNavig =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    pNavig->SetWorldVolume(myTopNode);

    G4GeometryManager::GetInstance()->CloseGeometry(false);

    // Setup the propagator (will be overwritten by testG4Propagator ...)
    pMagFieldPropagator = SetupPropagator(qssOrder);
    // G4cout << " Using default values for "
    //    << " Min Eps = "  <<   pMagFieldPropagator->GetMinimumEpsilonStep()
    //      << " and "
    //<< " MaxEps = " <<  pMagFieldPropagator->GetMaximumEpsilonStep()
    //<< G4endl;

    pMagFieldPropagator->SetUseSafetyForOptimization(optimisePiFwithSafety);
    // Do the tests without voxels
    // G4cout << " Test with no voxels" << G4endl;
    testG4PropagatorInField(myTopNode, qssOrder);

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels);
    pMagFieldPropagator->SetVerboseLevel(0);

    // Repeat tests but with full voxels
    // G4cout << " Test with full voxels" << G4endl;

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    testG4PropagatorInField(myTopNode, qssOrder);

    G4GeometryManager::GetInstance()->OpenGeometry();

    // G4cout << G4endl
    //    << "----------------------------------------------------------"
    //  << G4endl;

    // Repeat tests with full voxels and modified parameters
    // G4cout << "Test with more accurate parameters " << G4endl;

    G4double maxEpsStep = 0.001;
    G4double minEpsStep = 2.5e-8;
    // G4cout << " Setting values for Min Eps = " << minEpsStep
    //      << " and MaxEps = " << maxEpsStep << G4endl;

    pMagFieldPropagator->SetMaximumEpsilonStep(maxEpsStep);
    pMagFieldPropagator->SetMinimumEpsilonStep(minEpsStep);

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    testG4PropagatorInField(myTopNode, qssOrder);

    G4GeometryManager::GetInstance()->OpenGeometry();

    optimiseVoxels = !optimiseVoxels;
    // Repeat tests but with the opposite optimisation choice
    // G4cout << " Now test with optimisation " ;
    if (optimiseVoxels)
      G4cout << "on";
    else
      G4cout << "off";
    G4cout << G4endl;

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels);
    testG4PropagatorInField(myTopNode, qssOrder);

    G4GeometryManager::GetInstance()->OpenGeometry();
  }

  return 0;
}
