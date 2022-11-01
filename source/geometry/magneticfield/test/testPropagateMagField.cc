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
// Test all methods of integration of tracks in pure magnetic fields:
//   - all general RK steppers
//   - all specific RK steppers ( G4NystromRK4 )
//   - specialised methods ( helix-based methods such as G4HelixHeum )
//   - non-RK methods, such as multi-step (including Bulirsch-Stoer).
//
// Uses full stack of G4 classes for field propagation including
//    field, equation, stepper, driver(s), chord-finder, propagator-in-field
//
// Started from testG4Navigator1.cc,
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

#include <assert.h>
#include <iomanip>

// #include "ApproxEqual.hh"

// Global defs
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4ios.hh"

// Sample Parameterisation
class G4LinScale : public G4VPVParameterisation
{
  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0,(n-1)*15,0));
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int n,
				 const G4VPhysicalVolume* ) const
  {
    pBox.SetXHalfLength(10);
    pBox.SetYHalfLength(5+n);
    pBox.SetZHalfLength(5+n);
  }

  virtual void ComputeDimensions(G4Tubs &,
				 const G4int ,
                                 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trd &, 
				 const G4int,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Cons &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trap &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Hype &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Orb &,
		                 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Sphere &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Torus &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Para &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polycone &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polyhedra &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
};
G4LinScale myParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myHugeBox=  new G4Box("huge box",15*m,15*m,25*m);
    G4Box *myBigBox=   new G4Box("big cube",10*m,10*m,10*m);
    G4Box *mySmallBox= new G4Box("smaller cube",2.5*m,2.5*m,2.5*m);
    G4Box *myTinyBox=  new G4Box("tiny  cube",.25*m,.25*m,.25*m);

    // G4Box *myVariableBox=
    new G4Box("Variable Box",10,5,5);

    //  World Volume
    //
    G4LogicalVolume *worldLog=new G4LogicalVolume(myHugeBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new 
         G4PVPlacement(0,G4ThreeVector(0,0,0), "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set

//  Create the logical Volumes
//
//  G4LogicalVolume(*pSolid, *pMaterial, Name, *pField, *pSDetector, *pULimits)
//
    G4LogicalVolume *BigBoxLog=new G4LogicalVolume(myBigBox,0,
						"Crystal Box (large)",0,0,0);
    G4LogicalVolume *smallBoxLog=new G4LogicalVolume(mySmallBox,0,
						 "Crystal Box (small)");
    G4LogicalVolume *tinyBoxLog=new G4LogicalVolume(myTinyBox,0,
						 "Crystal Box (tiny)");


//  Place them.
//
//  1) Two big boxes in the world volume
//
    // G4PVPlacement *BigTg1Phys=
    new G4PVPlacement(0,G4ThreeVector(0,0,-15*m),
						"Big Target 1",BigBoxLog,
						worldPhys,false,0);
    // G4PVPlacement *BigTg2Phys=
    new G4PVPlacement(0,G4ThreeVector(0,0, 15*m),
						"Big Target 2",BigBoxLog,
						worldPhys,false,0);

//  2) Four (medium) boxes in X & Y near the origin of the world volume
//
    // G4PVPlacement *MedTg3a_Phys=
    new G4PVPlacement(0,G4ThreeVector(0, 7.5*m,0),
					      "Target 3a",smallBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *MedTg3b_Phys=
    new G4PVPlacement(0,G4ThreeVector(0,-7.5*m,0),
					      "Target 3b",smallBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *MedTg3c_Phys=
    new G4PVPlacement(0,G4ThreeVector(-7.5*m,0,0),
					      "Target 3c",smallBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *MedTg3d_Phys=
    new G4PVPlacement(0,G4ThreeVector( 7.5*m,0,0),
					      "Target 3d",smallBoxLog,
					      worldPhys,false,0);


//  3) Eight small boxes around the origin of the world volume 
//        (in +-X, +-Y & +-Z)
//
    // G4PVPlacement *SmTg4a_Phys=
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m, 0.3*m,0.3*m), "Target 4a",tinyBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *SmTg4b_Phys=
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m,-0.3*m,0.3*m), "Target 4b",tinyBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *SmTg4c_Phys=
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m,-0.3*m,0.3*m), "Target 4c",tinyBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *SmTg4d_Phys=
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m, 0.3*m,0.3*m), "Target 4d",tinyBoxLog,
					      worldPhys,false,0);

    // G4PVPlacement *SmTg4e_Phys=
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m, 0.3*m,-0.3*m), "Target 4e",tinyBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *SmTg4f_Phys=
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m,-0.3*m,-0.3*m), "Target 4f",tinyBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *SmTg4g_Phys=
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m,-0.3*m,-0.3*m), "Target 4g",tinyBoxLog,
					      worldPhys,false,0);
    // G4PVPlacement *SmTg4h_Phys=
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m, 0.3*m,-0.3*m), "Target 4h",tinyBoxLog,
					      worldPhys,false,0);

    return worldPhys;
}

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4CachedMagneticField.hh"

#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"

#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4DoLoMcPriRK34.hh"
#include "G4TsitourasRK45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ConstRK4.hh"
#include "G4NystromRK4.hh"

#include "G4HelixMixedStepper.hh"

// Templated equation -- make call to field non-virtual
#include "G4TMagFieldEquation.hh"

// Templated steppers -- make call to equation non-virtual
#include "G4TExplicitEuler.hh"
#include "G4TSimpleRunge.hh"
#include "G4TSimpleHeum.hh"
#include "G4TClassicalRK4.hh"
#include "G4TCashKarpRKF45.hh"
#include "G4TDormandPrince45.hh"

#include "G4RK547FEq1.hh"
#include "G4RK547FEq2.hh"
#include "G4RK547FEq3.hh"

// For new BulirsrchStoer method - Feb 2018
#include "G4BulirschStoer.hh"
#include "G4BulirschStoerDriver.hh"

#include "G4IntegrationDriver.hh"

#include "G4BFieldIntegrationDriver.hh"
#include "G4InterpolationDriver.hh"
#include "G4MagIntegratorDriver.hh"

#include "globals.hh"

using CLHEP::millimeter;

G4double  gDistConstField = 0.05 * millimeter;

G4double globDeltaChord= 0.010 * millimeter;         //  New values -- check old ones
G4double globDeltaIntersection= 0.005 * millimeter;  // 

G4UniformMagField      uniformMagField(10.*tesla, 0., 0.); 
G4CachedMagneticField  myCachedUniformField( &uniformMagField, gDistConstField );
G4String   fieldNameUniform("Uniform 10 Tesla"); 

G4RotationMatrix       rotM( 0.1, 0.2, 0.05 );  // RotationMatrix( phi, theta, psi );
G4ThreeVector          origin( 5.0 * millimeter, 2.5 * millimeter, 1.25 * millimeter );
G4QuadrupoleMagField   simpleQuadrupoleMagField( 10.*tesla/(50.*cm));
G4QuadrupoleMagField   quadrupoleMagField( 10.*tesla/(50.*cm), origin, &rotM ); 
// G4CachedMagneticField  myCachedQuadField( &quadrupoleMagField, gDistConstField );
G4CachedMagneticField  myCachedQuadField( &simpleQuadrupoleMagField, gDistConstField );
G4String   fieldNameQuad("Cached Quadropole field, 20T/meter, cache=1cm"); 

bool checkResults(
   G4ThreeVector InitPosition, G4ThreeVector InitMomentum, G4ThreeVector InitMomDir,
   G4ThreeVector EndPosition,  G4ThreeVector EndMomentum,  G4ThreeVector EndMomDir,
   double physStep,  double maxEpsilon );

static int FieldChoice= 1;

void SetFieldType( int ft )
{
   FieldChoice = ft;
   G4cout << " SetFieldType called: " << G4endl;
   if( ft == 0 ) {
      G4cout << "  Chosen Uniform field type: " << fieldNameUniform << " ." << G4endl;
   } else {
      G4cout << "  Chosen Quad    field type: " << fieldNameQuad << " ." << G4endl;
   }
} 

G4VIntegrationDriver* CodeForDocumentation( int useInterpolation )
{
   // const int numVar= 6;
   double stepMin= 0.1 * CLHEP::millimeter;
   
    // Testing code for Documentation  26.11.2020
   using Field_t = G4QuadrupoleMagField;
   using Equation_t = G4TMagFieldEquation<Field_t>;

   Field_t*      quadMagField = new G4QuadrupoleMagField(1.*tesla/(1.*meter));       
   Equation_t*  equation= new Equation_t(quadMagField);
   
   auto  equation3= new G4TMagFieldEquation<G4QuadrupoleMagField>(quadMagField);
   
   using TemplatedStepper_t = G4TDormandPrince45<Equation_t>;
   
   TemplatedStepper_t* dopri5_stepper=
      new G4TDormandPrince45<Equation_t>( equation );
   
   // Using it with templated drivers
   auto integrationDrv =
      new G4IntegrationDriver<TemplatedStepper_t>(stepMin, dopri5_stepper); // , numVar);
   
   auto interpolatingDrv =
      new G4InterpolationDriver<TemplatedStepper_t>(stepMin, dopri5_stepper); // , numVar);

   // delete equation2;
   delete equation3;
   // delete quadMagField;

   G4VIntegrationDriver* drv = nullptr;
   if( useInterpolation )
      drv = interpolatingDrv;
   else
      drv = integrationDrv;
   
   return drv;
   // return std::pair(integrationDriver,interpolationDriver);
}

G4CachedMagneticField *pMyMagField= nullptr;
using Equation_t = G4TMagFieldEquation<G4CachedMagneticField>;
Equation_t            *fEquation= nullptr;

G4FieldManager* SetupFieldManager(G4int stepperType)
{
    // using Equation_t = G4Mag_UsualEqRhs;

    
    G4FieldManager   * pFieldMgr= nullptr;
    G4ChordFinder    * pChordFinder= nullptr;
    G4MagIntegratorStepper  * pStepper= nullptr;
    G4VIntegrationDriver    * pDriver = nullptr;
    G4BulirschStoer         * pBSstepper = nullptr;
    const int nVar= 6;
    
    // Parameters
    G4double  stepMinimum = 1.0e-2 * millimeter;
    G4double epsilon = 1.0e-5;    

    // -- Parameters used for 'tight muon stepper field mgr' - ATLAS since ca 2017
    // GetChordFinder()->SetDeltaChord(0.00000002);
    // SetDeltaOneStep(0.000001);
    // SetDeltaIntersection(0.00000002);
    // SetMinimumEpsilonStep(0.0000009);
    // SetMaximumEpsilonStep(0.000001);
    
    G4String fieldName;
    if( FieldChoice == 0 ) { 
       pMyMagField = &myCachedUniformField;
       fieldName = fieldNameUniform;
    } else {
       pMyMagField = &myCachedQuadField;
       fieldName = fieldNameQuad;       
    }

    if( ! fEquation ) {
      assert(pMyMagField);       
      // fEquation= new G4Mag_UsualEqRhs(pMyMagField);
      fEquation= new Equation_t(pMyMagField);
      G4cout << " Created Equation for Magnetic field." << G4endl;
    }
    G4cout << " Setting up field of type: " << fieldName << G4endl;

    // For use in 'B-field driver' options;
    int  numVar = 6;
    using LargeStepDriver = G4IntegrationDriver<G4HelixHeum>;
    // auto longStepper = std::unique_ptr<G4HelixHeum>(new G4HelixHeum(fEquation));
    auto longStepperPtr = new G4HelixHeum(fEquation);


    auto rk54eq1= std::unique_ptr<G4RK547FEq1>(new G4RK547FEq1(fEquation) );
    auto shortStepper2= /* std::unique_ptr<G4RK547FEq2>(*/ new G4RK547FEq2(fEquation); // );
#if 0            
    auto shortStepper3= /*std::unique_ptr<G4RK547FEq3>(*/  new G4RK547FEq3(fEquation); // );
#endif
    auto shortStepper5= /* std::unique_ptr<G4DormandPrince745>( */ new G4DormandPrince745(fEquation); // );    
    auto shortStepper6= /* std::unique_ptr<G4DormandPrinceRK56>( */new G4DormandPrinceRK56(fEquation); // );

    switch ( stepperType ) 
    {
      case 0: pStepper = new G4ExplicitEuler( fEquation ); break;
      case 1: pStepper = new G4ImplicitEuler( fEquation ); break;
      case 2: pStepper = new G4SimpleRunge( fEquation ); break;
      case 3: pStepper = new G4SimpleHeum( fEquation ); break;
      case 4: pStepper = new G4ClassicalRK4( fEquation ); break;
      case 5: pStepper = new G4HelixExplicitEuler( fEquation ); break;
      case 6: pStepper = new G4HelixImplicitEuler( fEquation ); break;
      case 7: pStepper = new G4HelixSimpleRunge( fEquation ); break;
      case 17: pStepper = new G4HelixHeum( fEquation ); break;
      case 8: pStepper = new G4CashKarpRKF45( fEquation );    break;
      case 9: pStepper = new G4ExactHelixStepper( fEquation );   break;
      case 10: pStepper = new G4RKG3_Stepper( fEquation );       break;
      case 11: pStepper = new G4HelixMixedStepper( fEquation );  break;
      case 12: pStepper = new G4ConstRK4( fEquation ); break;
      case 13: pStepper = new G4NystromRK4( fEquation ); break;
      case 14: pStepper = new G4RK547FEq1( fEquation ); break;
      case 15: pStepper = new G4RK547FEq2( fEquation ); break;
      case 16: pStepper = new G4RK547FEq3( fEquation ); break;
      case 19: pStepper = new G4NystromRK4( fEquation, gDistConstField ); break;         
         
      case 23: pStepper = new G4BogackiShampine23( fEquation ); break;
      case 34: pStepper = new G4DoLoMcPriRK34( fEquation ); break;         
      case 45: pStepper = new G4BogackiShampine45( fEquation ); break;
      case 145: pStepper = new    G4TsitourasRK45( fEquation ); break;
      case 745: pStepper = new G4DormandPrince745( fEquation ); break;
      case 56: pStepper = new G4DormandPrinceRK56( fEquation ); break;
      case 78: pStepper = new G4DormandPrinceRK78( fEquation ); break;

      case 101: pStepper = new G4TExplicitEuler<Equation_t,6>( fEquation ); break;
      case 102: pStepper = new G4TSimpleRunge<Equation_t,6>( fEquation ); break;
      case 103: pStepper = new G4TSimpleHeum<Equation_t,6>( fEquation ); break;
      case 104: pStepper = new G4TClassicalRK4<Equation_t,6>( fEquation ); break;
      case 108: pStepper = new G4TCashKarpRKF45<Equation_t,6>( fEquation );    break;
      case 345: pStepper = new G4TDormandPrince45<Equation_t,6>( fEquation );  break;

      case 99: pStepper = nullptr;
         pBSstepper = new G4BulirschStoer( fEquation, nVar, epsilon );
         pDriver = new G4IntegrationDriver<G4BulirschStoer>( stepMinimum, pBSstepper, nVar );
         G4cout << "Using Bulirsch Stoer method (and driver) - Not Runge Kutta" << G4endl;
         break;

      case  205:
      case 1000:
          pDriver = new G4InterpolationDriver<G4DormandPrince745>(
                                               stepMinimum, new G4DormandPrince745(fEquation));
          break;

      case 245:
           pDriver = new G4InterpolationDriver<G4TDormandPrince45<Equation_t,6>>(
                                               stepMinimum,
                                               new G4TDormandPrince45<Equation_t,6>(fEquation) );                                    break;
          
      case 100:         
         pDriver = new G4IntegrationDriver<G4RK547FEq1>( stepMinimum,
                                                         rk54eq1.release()
                                                         , numVar );
         G4cout << "Using Integration Driver with RK547FEq1 Runge Kutta." << G4endl;
         break;

      case 201:
         pDriver = new G4BFieldIntegrationDriver(
            std::unique_ptr<G4IntegrationDriver<G4RK547FEq1>>(
               new G4IntegrationDriver<G4RK547FEq1>( stepMinimum, rk54eq1.release()
                                                       , numVar ) ),
            std::unique_ptr<LargeStepDriver>(new LargeStepDriver(stepMinimum,
                                                                 longStepperPtr, numVar) )
            );
         G4cout << "Using B-Field Driver with RK547FEq1 Runge Kutta (short) and Helix Heum (for long steps)" << G4endl;
         break;
      case 202:  

         pDriver = new G4BFieldIntegrationDriver(
            std::unique_ptr<G4IntegrationDriver<G4RK547FEq2>>(
               new G4IntegrationDriver<G4RK547FEq2>( stepMinimum, shortStepper2,  numVar ) ),
            std::unique_ptr<LargeStepDriver>(new LargeStepDriver(stepMinimum,
                                                                 longStepperPtr, numVar) )
            );
         G4cout << "Using B-Field Driver with RK547FEq2 Runge Kutta (short) and Helix Heum (for long steps)" << G4endl;
         break;
    /******
      case 203:
         pDriver = new G4BFieldIntegrationDriver(
            std::unique_ptr<G4InterpolationDriver<G4RK547FEq2>>(
               new G4InterpolationDriver<G4RK547FEq2>( stepMinimum, shortStepper3, numVar ) ),
            std::unique_ptr<LargeStepDriver>(new LargeStepDriver(stepMinimum,
                                                                 longStepperPtr, numVar) )
            );         
         G4cout << "Using B-Field Driver with RK547FEq3 Runge Kutta (short) and Helix Heum (for long steps)" << G4endl;
         break;
    ***/
      case 405:
         pDriver = new G4BFieldIntegrationDriver(
            std::unique_ptr<G4InterpolationDriver<G4DormandPrince745>>(
               new G4InterpolationDriver<G4DormandPrince745>( stepMinimum, shortStepper5, numVar ) ),
            std::unique_ptr<LargeStepDriver>(new LargeStepDriver(stepMinimum,
                                                                 longStepperPtr, // longStepper.get(),
                                                                 numVar) )
            );         
         G4cout << "Using B-Field Driver with Dormand-Prince 5/4 (short) and Helix Heum (for long steps)" << G4endl;
         break;

      case 106:  
         pDriver = new G4BFieldIntegrationDriver(
               std::unique_ptr<G4IntegrationDriver<G4DormandPrinceRK56>>(
                           new G4IntegrationDriver<G4DormandPrinceRK56>( stepMinimum,
                                                                         shortStepper6 // .get()
                                                                         , numVar )
                  ),
               std::unique_ptr<LargeStepDriver>(new LargeStepDriver(stepMinimum,
                                                                 longStepperPtr, numVar) )
            );         
         G4cout << "Using B-Field Driver with Dormand-Prince 6/5 (short) and Helix Heum (for long steps)" << G4endl;
         break;

      //  case 108: // pDriver = new G4BFieldIntegrationDriver<G4DormandPrinceRK78>(
              //                 stepMinimum, new G4DormandPrinceRK78(fEquation)); break;

      default: 
          pStepper = nullptr;   // Can use default= new G4ClassicalRK4( fEquation );
          G4ExceptionDescription ErrorMsg;
          ErrorMsg << " Incorrect Stepper type requested. Value was id= " 
                   << stepperType << G4endl;
          ErrorMsg << " NO replacement stepper chosen! " << G4endl;
          G4Exception("application::SetupFieldManager",
                      "Runtime Error",
                      FatalErrorInArgument,       //  use JustWarning,
                      " Invalid value of stepper type" );
          break; 
    }

    pFieldMgr= G4TransportationManager::GetTransportationManager()->
       GetFieldManager();

    pFieldMgr->SetDetectorField( pMyMagField ); // ( &myMagField );

    if( pStepper )
    {
      pChordFinder = new G4ChordFinder( pMyMagField,
                                        stepMinimum,
                                        pStepper);

    }
    else if ( pDriver )
    {
       pChordFinder = new G4ChordFinder( pDriver );
    }
    else
    {
       G4cerr << "ERROR> Neither valid stepper nor driver for Bulirsch Stoer method configured."
              << G4endl << " FATAL - exiting "          
              << G4endl;
    }
    // pFieldMgr->SetDeltaIntersection( deltaIntersection );
    
    assert( pChordFinder != 0 );
    pChordFinder->SetDeltaChord( globDeltaChord );
    // pChordFinder->SetDeltaIntersection( globDeltaIntersection );
    
    pChordFinder->SetVerbose(1);

    pFieldMgr->SetChordFinder( pChordFinder );

    return    pFieldMgr;
}

#include "G4SimpleLocator.hh"
#include "G4BrentLocator.hh"
#include "G4MultiLevelLocator.hh"

G4PropagatorInField*  SetupPropagator( G4int StepperType, G4int LocatorType = 1 )
{
    G4FieldManager* fieldMgr= 
       SetupFieldManager( StepperType ) ;

    // G4ChordFinder  theChordFinder( &MagField, 0.05*mm ); // Default stepper
 
    G4PropagatorInField *thePropagator = 
      G4TransportationManager::GetTransportationManager()->
       GetPropagatorInField ();

    G4bool useTightParameters = false;
    if( useTightParameters )     
    {
       // Parameters - old values
       G4double epsilon = 1.0e-7;        // new: 1.0e-6
       G4double reductionFactor = 1.0 ;  // new: 0.9
       G4double maxEpsilon = epsilon;
       G4double minEpsilon = reductionFactor * epsilon;  // was = epsilon

       //  ATLAS tight muon values
       // GetChordFinder()->SetDeltaChord(0.00000002);
       // SetDeltaOneStep(0.000001);
       // SetDeltaIntersection(0.00000002);
       // SetMinimumEpsilonStep(0.0000009);
       // SetMaximumEpsilonStep(0.000001);
       
       // Let us test the Minimum Epsilon Step functionality
       /* thePropagator*/ fieldMgr -> SetMinimumEpsilonStep( minEpsilon ) ; 
       /* thePropagator*/ fieldMgr -> SetMaximumEpsilonStep( maxEpsilon ) ;

       // fieldMgr->SetDeltaOneStep(0.000001);  --> Object not available
       globDeltaChord= 0.00000002;

       // GetChordFinder()->SetDeltaIntersection(0.00000002);
       globDeltaIntersection= 0.00000002;
      
       G4cout << " Installed values for Min Eps = " // " Found values for Min Eps = "
              << thePropagator->GetMinimumEpsilonStep()
              << " and Max Eps = " << thePropagator->GetMaximumEpsilonStep()
              << G4endl;
       G4cout << "  -- storing requested values : " << G4endl
              << "DeltaChord = " << globDeltaChord   << G4endl
              <<  "DeltaIntersection = " << globDeltaIntersection << G4endl;
    }
    else
    {
       G4cout << " Found values for Min Eps = "
              << thePropagator->GetMinimumEpsilonStep()
              << " and Max Eps = " << thePropagator->GetMaximumEpsilonStep()
              << G4endl; 
    }
    G4Navigator *theNavigator= G4TransportationManager::GetTransportationManager()->
       GetNavigatorForTracking();
    // Test the options for Locator
    G4VIntersectionLocator *pLocator=0;
    G4cout << "Over-riding  PropagatorInField to use ";
    // pLocator= new G4MultiLevelLocator(theNavigator); G4cout << "Multi"; // default
    // pLocator= new G4SimpleLocator(theNavigator); G4cout << "Simple";
    // pLocator= new G4BrentLocator(theNavigator); G4cout << " Brent "; 

    switch ( LocatorType ) {
    case 0: 
       pLocator= new G4SimpleLocator(theNavigator); G4cout << "Simple";
       break;       

    case 2:
      pLocator= new G4BrentLocator(theNavigator); G4cout << " Brent ";
      break;      
       
    case 1: 
    default : 
       pLocator= new G4MultiLevelLocator(theNavigator); G4cout << "Multi";
       break;
    };
    G4cout << " Locator. ( In the unit test code. ) " << G4endl;
    
    thePropagator->SetIntersectionLocator(pLocator);

    return thePropagator;
}

G4PropagatorInField *pMagFieldPropagator=0; 
//
// Test Stepping
//
G4bool testG4PropagatorInField(G4VPhysicalVolume*,     // *pTopNode, 
			       G4int               stepperType,
                               G4int               locatorType = 1)
{
    void report_endPV(G4ThreeVector    Position, 
                  G4ThreeVector UnitVelocity,
		  G4double step_len, 
                  G4double physStep, 
                  G4double safety,
		  G4ThreeVector EndPosition, 
                  G4ThreeVector EndUnitVelocity,
                  G4int             Step, 
                  G4VPhysicalVolume* startVolume);

    G4UniformMagField MagField(10.*tesla, 0., 0.);
    G4Navigator   *pNavig= G4TransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    
    // Force the Navigator to do a full search from the top of the tree
    G4bool relativeSearch =false;
    pNavig->LocateGlobalPointAndSetup(G4ThreeVector(2000., 20000., 2000.),
                                      nullptr, relativeSearch, true );
  
    pMagFieldPropagator= SetupPropagator(stepperType, locatorType );

    G4double particleCharge= +1.0;  // in e+ units
    G4double spin=0.0;              // ignore the spin
    G4double magneticMoment= 0.0;   // ignore the magnetic moment

    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0); 

    G4ChordFinder *pChordFndr= pMagFieldPropagator->GetChordFinder();
    
    G4EquationOfMotion* equationOfMotion =
        pChordFndr->GetIntegrationDriver()->GetEquationOfMotion();
    
    equationOfMotion->SetChargeMomentumMass( chargeState, 
			            0.5 * proton_mass_c2, // Momentum in Mev/c
					 proton_mass_c2 );
    // pNavig->SetWorldVolume(pTopNode);

    G4VPhysicalVolume *located;
    G4double step_len, physStep, safety;
    G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    // physStep=kInfinity;
    G4ThreeVector Position(0.,0.,0.); 
    G4ThreeVector UnitMomentum(0.,0.6,0.8);  
    G4ThreeVector EndPosition, EndUnitMomentum, EndMomentum;

//
// Test location & Step computation
//  
    /* assert(located->GetName()=="World"); */
    if( std::fabs(UnitMomentum.mag() - 1.0) > 1.e-8 ) 
    {
      G4cerr << "UnitMomentum.mag() - 1.0 = " << UnitMomentum.mag() - 1.0 <<
	G4endl;
    }

    G4cout << G4endl; 

    for( int iparticle=0; iparticle < 2; iparticle++ )
    { 
       physStep=  2.5 * mm ;  // millimeters 
       Position = G4ThreeVector(0.,0.,0.) 
	        + iparticle * G4ThreeVector(0.2, 0.3, 0.4); 
       UnitMomentum = (G4ThreeVector(0.,0.6,0.8) 
		    + (float)iparticle * G4ThreeVector(0.1, 0.2, 0.3)).unit();

       G4double momentum = (0.5+iparticle*10.0) * proton_mass_c2; 

       G4double kineticEnergy =  momentum*momentum /
                  ( std::sqrt( momentum*momentum + proton_mass_c2 * proton_mass_c2 ) 
		    + proton_mass_c2 );
       G4double velocity = momentum / ( proton_mass_c2 + kineticEnergy );
       G4double labTof= 10.0*ns, properTof= 0.1*ns;
       G4ThreeVector Spin(1.0, 0.0, 0.0);
                                                   // Momentum in Mev/c ?

       G4ChargeState chargeSt(1.0, 0.0, 0.5 ); 
       // pMagFieldPropagator
       equationOfMotion->SetChargeMomentumMass(
		      chargeSt,                    // charge in e+ units
		      momentum, 
		      proton_mass_c2); 
       G4cout << G4endl;
       G4cout << "Test PropagateMagField: ***********************" << G4endl
            << " Starting New Particle with Position " << Position << G4endl 
	    << " and UnitVelocity " << UnitMomentum << G4endl;
       G4cout << " Momentum in GeV/c is " << momentum / GeV
	      << " = " << (0.5+iparticle*10.0)*proton_mass_c2 / MeV << " MeV"
              << G4endl;

       G4bool canRelaxDeltaChord= true;
       pMagFieldPropagator->SetIterationsToIncreaseChordDistance(5); // 10);

       for( int istep=0; istep < 14; istep++ ){ 
          // G4cerr << "UnitMomentum Magnitude is " << UnitMomentum.mag() << G4endl;
	  located= pNavig->LocateGlobalPointAndSetup(Position);
	  // G4cerr << "Starting Step " << istep << " in volume " 
	       // << located->GetName() << G4endl;

          G4FieldTrack  initTrack( Position, 
				   UnitMomentum,
				   0.0,            // starting S curve len
				   kineticEnergy,
				   proton_mass_c2,
				   velocity,
				   labTof, 
				   properTof,
				   0              // or &Spin
				   ); 
         
	  step_len=pMagFieldPropagator->ComputeStep( initTrack, 
						     physStep, 
						     safety,
						     located,
                                                     canRelaxDeltaChord
                                                    );
	  //       --------------------
	  EndPosition=     pMagFieldPropagator->EndPosition();
	  EndUnitMomentum= pMagFieldPropagator->EndMomentumDir();
	  EndMomentum= pMagFieldPropagator->GetEndState().GetMomentum();
	  //       --------

          G4double maxEpsilon = pMagFieldPropagator -> GetMaximumEpsilonStep();
          G4bool ok = checkResults( Position, initTrack.GetMomentum(), initTrack.GetMomentumDir(),
                                    EndPosition, EndMomentum, EndUnitMomentum,
                                    physStep,    maxEpsilon );
          
	  // G4cout << " testPropagatorInField: After stepI " << istep  << " : " << G4endl;
          // if ( !ok && !verbose )
	  // report_endPV(Position, UnitMomentum, -1.0, -1.0, previousSafety,
	  //       Position, UnitMomentum, istep-1, located );
          G4bool verbose = true;
          if ( verbose || !ok )
             report_endPV(Position, UnitMomentum, step_len, physStep, safety,
                          EndPosition, EndUnitMomentum, istep, located );

	  assert(safety>=0);
	  pNavig->SetGeometricallyLimitedStep();
	  // pMagFieldPropagator->SetGeometricallyLimitedStep();

	  Position= EndPosition;
	  UnitMomentum= EndUnitMomentum;
	  physStep *= 2.; 
       } // ...........................  end for ( istep )

       // myMagField.
       G4cout << "Accumulated statistics: " << G4endl;
       pMyMagField->ReportStatistics(); 
       // pMyMagField->ClearCounts(); // --> don't accumulate

    }    // ..............................  end for ( iparticle )

    return(1);
}

void reportPositionDirection(G4ThreeVector    Position, 
                             G4ThreeVector    UnitVelocity)
{
   G4cout.precision(6);
   G4cout << std::setw( 9) << Position.x() << " "
          << std::setw( 9) << Position.y() << " "
          << std::setw( 9) << Position.z() << " "
          << std::setw( 9) << UnitVelocity.x() << " "
          << std::setw( 9) << UnitVelocity.y() << " "
          << std::setw( 9) << UnitVelocity.z() << " ";
}

void report_endPV(G4ThreeVector    Position, 
                  G4ThreeVector    InitialUnitVelocity,
		  G4double step_len, 
                  G4double physStep, 
                  G4double safety,
		  G4ThreeVector EndPosition, 
                  G4ThreeVector EndUnitVelocity,
                  G4int             Step, 
                  G4VPhysicalVolume* startVolume)
		  //   G4VPhysicalVolume* endVolume)
{
    const G4int verboseLevel=1;

    void reportPositionDirection(G4ThreeVector , G4ThreeVector);
    
    if( Step == 0 && verboseLevel <= 3 )
    {
       G4cout.precision(6);
       // G4cout.setf(ios_base::fixed,ios_base::floatfield);
       G4cout << std::setw( 5) << "Step#" << " "
            << std::setw( 9) << "X(mm)" << " "
            << std::setw( 9) << "Y(mm)" << " "  
            << std::setw( 9) << "Z(mm)" << " "
            << std::setw( 9) << " N_x " << " "
            << std::setw( 9) << " N_y " << " "
            << std::setw( 9) << " N_z " << " "
            << std::setw( 9) << " Delta|N|" << " "
            << std::setw( 9) << " Delta(N_z) " << " "
	   // << std::setw( 9) << "KinE(MeV)" << " "
	   // << std::setw( 9) << "dE(MeV)" << " "  
            << std::setw( 9) << "StepLen" << " "  
            << std::setw( 9) << "PhsStep" << " "  
            << std::setw( 9) << "Safety" << " "
            << std::setw(10) << "FieldEvals" << " "
            << std::setw( 9) << "NextVolume" << " "
            << G4endl;

       G4cout.precision(6);
       G4cout << std::setw( 5) << Step << " ";
       reportPositionDirection( Position, InitialUnitVelocity );   // Report starting values
       G4cout << std::setw( 9) << " " << " "
              << std::setw( 9) << " " << " ";
       G4cout 
	    << std::setw( 9) << 0.0 << " "
	    << std::setw( 9) << 0.0 << " ";
       //   << std::setw( 9) << safety << " ";
       G4cout << G4endl;
    }
    //
    //
    if( verboseLevel > 3 )
    {
       G4cout << "End  Position is " << EndPosition << G4endl 
	    << " and UnitVelocity is " << EndUnitVelocity << G4endl;
       G4cout << "Step taken was " << step_len  
	    << " out of PhysicalStep= " <<  physStep << G4endl;
       G4cout << "Final safety is: " << safety << G4endl;

       G4cout << "Chord length = " << (EndPosition-Position).mag() << G4endl;
       G4cout << G4endl; 
    }
    else // if( verboseLevel > 0 )
    {
       G4cout.precision(6);
       G4cout << std::setw( 5) << Step << " ";
       // reportPositionDirection( Position, EndUnitVelocity );  // Old code
       reportPositionDirection( EndPosition, EndUnitVelocity );  // Consistent
       G4cout.precision(2); 
       G4cout
	    << std::setw( 9) << EndUnitVelocity.mag()-InitialUnitVelocity.mag() << " "
	    << std::setw( 9) << EndUnitVelocity.z() - InitialUnitVelocity.z() << " ";
	 //    << std::setw( 9) << KineticEnergy << " "
	 //    << std::setw( 9) << EnergyDifference << " "
       G4cout.precision(6);
       G4cout 
	    << std::setw( 9) << step_len << " "
	    << std::setw( 9) << physStep << " "
	    << std::setw( 9) << safety << " ";

       // Report number of field evaluations
       G4cout // << std::setw(5) << pMyMagField->GetCountCalls()
              << std::setw(8) << pMyMagField->GetCountEvaluations();
       G4cout << " ";
       if( startVolume != 0) {
	 G4cout << std::setw(11) << startVolume->GetName() << " ";
       } else {
	 G4cout << std::setw(11) << "OutOfWorld" << " ";
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

bool checkResults( G4ThreeVector InitialPosition, G4ThreeVector InitialMomentum, G4ThreeVector InitialMomentumDir,
                   G4ThreeVector EndPosition,     G4ThreeVector EndMomentum,  G4ThreeVector EndMomentumDir,
                   double physStep,    // double stepTaken,
                   double maxEpsilon
   )
{
   G4bool failed= false;
   
   if( std::fabs(EndMomentumDir.mag2() - 1.0) > 1.e-8 )
      G4cerr << "EndMomentumDir.mag2() - 1.0 = " <<
         EndMomentumDir.mag2() - 1.0 << G4endl;

   G4ThreeVector MoveVec = EndPosition - InitialPosition;

   // assert( MoveVec.mag() < physStep*(1.+1.e-9) );

   const  double maxMoveRatio = 1.0e-9; 
   const  double epsMomentum = 1.0e-9;
   
   double momentumMagStart = InitialMomentum.mag();
   double momentumMagEnd   = EndMomentum.mag();
   double diffMomMag = momentumMagEnd - momentumMagStart;
   double averMomMag = 0.5 * ( momentumMagEnd - momentumMagStart );
   
   double move= MoveVec.mag();
   double displaceRatio = 0.0;
   if( physStep > 0. ) {
      displaceRatio = move / physStep - 1.0;
   }
   if( displaceRatio > maxMoveRatio ) {  // Was 1.0e-9
      int eprec= G4cerr.precision(10);
      G4ThreeVector diffUnitMom = EndMomentumDir - InitialMomentumDir;
      G4double  dotDir = InitialMomentumDir.dot( EndMomentumDir );
      G4cerr << " Move displaced " << std::setw(13)
             << move << " - further than step= "
             << std::setw(6) << physStep << "  by  "
             << std::setw(9) << std::setprecision(4) << move - physStep << ", "
             << " a fraction of " << displaceRatio << "  "
             << " |D dir|= " << diffUnitMom.mag() << " "
             << " p-dir: "
             << " dot-1 = " << dotDir - 1.0 << " angle= " << std::acos( dotDir )
         // << " start= " << std::setw(12) << std::setprecision(9) << InitialMomentumDir << "  "
         // << " end = "  << std::setw(12) << std::setprecision(9) << EndMomentumDir << " "
             << " diff= " << diffUnitMom
             << G4endl;
      G4cerr.precision(eprec);
      failed= true;
   }

   assert( displaceRatio < maxEpsilon ); // Revised check condition 
   if( displaceRatio >= maxEpsilon ){
      G4cerr << "ERROR> Problem - displacement ratio " << displaceRatio
             << " exceeds maximum allowed value = " << maxEpsilon
             << G4endl;
      exit(1);
   }
          
   if( std::fabs(diffMomMag) > epsMomentum * averMomMag )
   {
      int eprec= G4cerr.precision(10);
      G4cerr << " Momentum magnitude change:  relative: "
             << std::setw(9) << std::setprecision(5)
             << diffMomMag / averMomMag
             << " absolute: " << std::setw(12) << std::setprecision(9) << diffMomMag
             << " versus start " << std::setw(12) << std::setprecision(9) << momentumMagStart << " and end " 
             << " versus end " << std::setw(12) << std::setprecision(9) << momentumMagEnd
             << G4endl;
      G4cerr.precision(eprec);
      failed= true;      
      
   }
   
   return ! failed;
}

void repeatWith( double minEpsStep, double maxEpsStep, G4VPhysicalVolume *topNode, int stepperType )
{
    G4cout << "Test with more accurate parameters " << G4endl; 
    pMagFieldPropagator->SetMaximumEpsilonStep(maxEpsStep);
    pMagFieldPropagator->SetMinimumEpsilonStep(minEpsStep);

    G4cout << " Setting values for Min Eps = "
           << pMagFieldPropagator->GetMinimumEpsilonStep()
           << " ( intended = " << minEpsStep << " ) "
           << " and MaxEps = " << pMagFieldPropagator->GetMaximumEpsilonStep()
           << " ( intended = " << maxEpsStep << " ). "
           << G4endl; 

    testG4PropagatorInField(topNode, stepperType);   
}


// Main program
// -------------------------------
int main(int argc, char **argv)
{
    G4VPhysicalVolume *myTopNode;
    G4int  stepperType;
    G4int  optim, optimSaf;
    G4bool optimiseVoxels=true;
    G4bool optimisePiFwithSafety=true;

    stepperType = 8 ;  // Default stepper - Cash Karp RKF 45
    int    locatorType= 1;
    int    fieldType= 1;  // Default field type:  Quadrupole (1)  or Constant (0)

    G4cout << " Arguments:  stepper-no  optimise-Voxels optimise-PiF-with-safety" << G4endl;

    if( argc >= 2 ){
       stepperType = atoi(argv[1]);
    } 

    if( argc >=3 ){
       locatorType = atoi(argv[2]);
    }
    
    if( argc >=4 ){
      optim= atoi(argv[3]);
      if( optim == 0 ) { optimiseVoxels = false; }
    }

    if( argc >=5 ){
      optimSaf= atoi(argv[4]);
      if( optimSaf == 0 ) { optimisePiFwithSafety= false; }
    }

    if( argc >=6 ){
      G4int fieldTypeIn= atoi(argv[5]);
      if( fieldTypeIn < 0 || fieldTypeIn > 1 ){
         G4cerr << "Warning: Invalid field type input= " << fieldTypeIn << G4endl;
         G4cerr << "          Keeping default value =  " << fieldType << G4endl;
      }else{
         fieldType = fieldTypeIn;
      }
      SetFieldType(fieldType);      
    }
    
    G4cout << " Testing with stepper number   = " << stepperType << G4endl;
    G4cout << "       Using 'locator'-type    = " << locatorType; 
    G4cout << "        ( Key: 0 = simple, 1 = Multi, 2= Brent ) " << G4endl;
    G4cout << "       voxel optimisation      : "
           << (optimiseVoxels ? "On" : "Off")  << G4endl;
    G4cout << "       Propagator safety optim : " 
           << (optimisePiFwithSafety ? "On" : "Off")  << G4endl;

    // Create the geometry & field 
    myTopNode=BuildGeometry();	// Build the geometry
 
    G4Navigator *pNavig= G4TransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    pNavig->SetWorldVolume(myTopNode);

    G4GeometryManager::GetInstance()->CloseGeometry(false);
    
    // pMagFieldPropagator->SetMinimumEpsilonStep();
    // pMagFieldPropagator->SetMaximumEpsilonStep();
       
    // Setup the propagator (will be overwritten by testG4Propagator ...)
    pMagFieldPropagator= SetupPropagator(stepperType);
    G4cout << " Using default values for " 
	   << " Min Eps = "  <<   pMagFieldPropagator->GetMinimumEpsilonStep()
           << " and "
	   << " MaxEps = " <<  pMagFieldPropagator->GetMaximumEpsilonStep()
	   << G4endl; 

    pMagFieldPropagator->SetUseSafetyForOptimization(optimisePiFwithSafety); 

// Do the tests without voxels
    G4cout << " Test with no voxels" << G4endl; 
    testG4PropagatorInField(myTopNode, stepperType);

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels); 
    pMagFieldPropagator->SetVerboseLevel( 1 ); 

// Repeat tests but with full voxels
    G4cout << " Test with full voxels" << G4endl; 

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    testG4PropagatorInField(myTopNode, stepperType);

    G4GeometryManager::GetInstance()->OpenGeometry();

    G4cout << G4endl
	   << "----------------------------------------------------------"
	   << G4endl; 

// Repeat tests with full voxels and modified parameters
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 

    G4double  maxEps= 1.0e-5;
    G4double  minEps= 1.0e-5;
    G4cout << " Setting values for Min Eps = " << minEps 
           << " and MaxEps = " << maxEps << G4endl; 

    pMagFieldPropagator->SetMaximumEpsilonStep(maxEps);
    pMagFieldPropagator->SetMinimumEpsilonStep(minEps);

    testG4PropagatorInField(myTopNode, stepperType);

    
// Repeat tests with full voxels and modified parameters
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-6;
    minEps= 1.0e-6;
    repeatWith( minEps, maxEps, myTopNode, stepperType );
    
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-7;
    minEps= 1.0e-7;
    repeatWith( minEps, maxEps, myTopNode, stepperType );
    
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-8;
    minEps= 1.0e-8;
    repeatWith( minEps, maxEps, myTopNode, stepperType );
    
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-9;
    minEps= 1.0e-9;
    repeatWith( minEps, maxEps, myTopNode, stepperType );
    
    G4GeometryManager::GetInstance()->OpenGeometry();
    optimiseVoxels = ! optimiseVoxels;
// Repeat tests but with the opposite optimisation choice
    G4cout << " Now test with optimisation " ; 
    if (optimiseVoxels)   G4cout << "on"; 
    else            G4cout << "off"; 
    G4cout << G4endl;

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels); 
    testG4PropagatorInField(myTopNode, stepperType);

    G4GeometryManager::GetInstance()->OpenGeometry();

    // Cannot delete G4TransportationManager::GetInstance(); 
    return 0;
}


  
