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
//
//  
//
// Started from testG4Navigator1.cc,v 1.7 1996/08/29 15:42 pkent 
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

#include <assert.h>
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
#include <iomanip>

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
#include "G4HelixImplicitEuler.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "DormandPrince745.hh"
#include "DormandPrinceRK56.hh"
#include "DormandPrinceRK78.hh"
#include "DoLoMcPriRK34.hh"
#include "TsitourasRK45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ConstRK4.hh"
#include "G4NystromRK4.hh"
#include "G4HelixMixedStepper.hh"

// FSAL-interface versions of steppers
#include "FBogackiShampine45.hh"
#include "fDormandPrince745.hh"

#include "globals.hh"

G4UniformMagField      uniformMagField(10.*tesla, 0., 0.); 
G4CachedMagneticField  myCachedUniformField( &uniformMagField, 1.0 * cm); 
G4String   fieldNameUniform("Uniform 10 Tesla"); 

G4QuadrupoleMagField   quadrupoleMagField( 10.*tesla/(50.*cm) ); 
G4CachedMagneticField  myCachedQuadField( &quadrupoleMagField, 1.0 * cm); 
G4String   fieldNameQuad("Cached Quadropole field, 20T/meter, cache=1cm"); 

G4CachedMagneticField *pMyMagField;

static int FieldChoice= 0;

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

G4Field* SetupField(G4int type)
{
    G4FieldManager   *pFieldMgr;
    G4ChordFinder    *pChordFinder;

    G4MagIntegratorStepper   *pStepper= 0; // nullptr;
    FSALMagIntegratorStepper *pFSALStepper= 0; // nullptr;    
    G4String fieldName;

    // fEquation= new G4Mag_UsualEqRhs(&myMagField); 
    if( FieldChoice == 0 ) { 
       pMyMagField = &myCachedUniformField;
       fieldName = fieldNameUniform;
    } else {
       pMyMagField = &myCachedQuadField;
       fieldName = fieldNameQuad;       
    }
    
    G4cout << " Set  up field of type: " << fieldName << G4endl;    
}

// static G4Mag_UsualEqRhs *fEquation= nullptr;

G4EquationOfMotion *SetupEquation( G4MagneticField* magField )
{
    if( ! fEquation ) {
      fEquation= new G4Mag_UsualEqRhs(pMyMagField);
      G4cout << " Created Equation for Magnetic field." << G4endl;
    }
}

G4MagIntegratorStepper* CreateStepper( int StepperType, G4EquationOfMotion* equation )
{
    G4MagIntegratorStepper* pStepper;
    switch ( type ) 
    {
      case 0: pStepper = new G4ExplicitEuler( equation ); break;
      case 1: pStepper = new G4ImplicitEuler( equation ); break;
      case 2: pStepper = new G4SimpleRunge( equation ); break;
      case 3: pStepper = new G4SimpleHeum( equation ); break;
      case 4: pStepper = new G4ClassicalRK4( equation ); break;
      case 5: pStepper = new G4HelixExplicitEuler( equation ); break;
      case 6: pStepper = new G4HelixImplicitEuler( equation ); break;
      case 7: pStepper = new G4HelixSimpleRunge( equation ); break;
      case 8: pStepper = new G4CashKarpRKF45( equation );    break;
      case 9: pStepper = new G4ExactHelixStepper( equation );   break;
      case 10: pStepper = new G4RKG3_Stepper( equation );       break;
      case 11: pStepper = new G4HelixMixedStepper( equation );  break;
      case 12: pStepper = new G4ConstRK4( equation ); break;
      case 13: pStepper = new G4NystromRK4( equation ); break;
      case 23: pStepper = new G4BogackiShampine23( equation ); break;
      case 34: pStepper = new DoLoMcPriRK34( equation ); break;         
      case 45: pStepper = new G4BogackiShampine45( equation ); break;
      case 145: pStepper = new      TsitourasRK45( equation ); break;
      case 745: pStepper = new DormandPrince745( equation ); break;
      case 56: pStepper = new DormandPrinceRK56( equation ); break;
      case 78: pStepper = new DormandPrinceRK78( equation ); break;
      default: 
          pStepper = 0;   // Can use default= new G4ClassicalRK4( equation );
          G4ExceptionDescription ErrorMsg;
          ErrorMsg << " Incorrect Stepper type requested. Value was id= " 
                   << type << G4endl;
          ErrorMsg << " NO replacement stepper chosen! " << G4endl;
          G4Exception("application::CreateStepper()", "Runtime Error",
                      FatalErrorInArgument,       //  use JustWarning,
                      " Invalid value of stepper type" );
          break; 
    }
    return pStepper;
}

G4VFSALIntegratorStepper* CreateFSALStepper( int StepperType, G4EquationOfMotion* equation )
{
    G4VFSALIntegratorStepper* pFSALStepper= nullptr;
    switch ( StepperType )
    {
         // FSAL Steppers
      case 245: pFSALStepper = new G4FSALBogackiShampine45( equation ); break;
      case 345: pFSALStepper = new G4FSALDormandPrince745( equation ); break;
      default: 
          pStepper = 0;   // Could use default= new 
          G4ExceptionDescription ErrorMsg;
          ErrorMsg << " Incorrect Stepper type requested. Value was id= " 
                   << type << G4endl;
          ErrorMsg << " NO replacement stepper chosen! " << G4endl;
          G4Exception("application::CreateFSALStepper",   "Runtime Error", 
                      FatalErrorInArgument, " Invalid value of stepper type" );
          break; 
    }
    return pFSALStepper;
}

G4EquationOfMotion *SetupEquation( G4MagneticField* magField );
G4MagIntegratorStepper* CreateStepper( int StepperType, G4EquationOfMotion* equation );
G4VFSALIntegratorStepper* CreateFSALStepper( int StepperType, G4EquationOfMotion* equation );


#include "G4SimpleLocator.hh"
#include "G4BrentLocator.hh"
#include "G4MultiLevelLocator.hh"

// G4PropagatorInField*  SetupPropagator( G4int type)

G4NewChordFinder*  SetupNewChordFinder( G4int stepperType, G4MagIntegratorStepper )
{
    // G4FieldManager* fieldMgr= 
    SetupField( stepperType) ;

    G4NewChordFinder  *pChordFinder= nullptr; // new G4NewChordFinder( &MagField, 0.05*mm ); // Default stepper

    G4double  stepMinimum = 1.0e-2 * millimeter;    
    if( pStepper )      
       pChordFinder = new G4NewChordFinder( pMyMagField, // &myMagField,
                                            stepMinimum, // 1.0e-2 * mm,
                                            pStepper);
    else {
       pChordFinder = new G4NewChordFinder( pMyMagField,
                                            pFSALStepper,
                                            stepMinimum);
    }
    pChordFinder->SetVerbose(1);  // ity();
    
    // G4PropagatorInField *thePropagator = 
    //   G4TransportationManager::GetTransportationManager()->
    //    GetPropagatorInField ();

    // G4double epsilon = 1.0e-4;
    // G4double minEpsilon = epsilon;
    // G4double maxEpsilon = epsilon;
    
    // Let us test the new Minimum Epsilon Step functionality
    // thePropagator -> SetMinimumEpsilonStep( minEpsilon ) ; 
    // thePropagator -> SetMaximumEpsilonStep( maxEpsilon ) ;

    // G4Navigator *theNavigator= G4TransportationManager::GetTransportationManager()->
    //      GetNavigatorForTracking();
    // Test the options for Locator

    return pChordFinder; // thePropagator;
}


G4FieldManager* SetupFieldManager( )
{    
    pFieldMgr= G4TransportationManager::GetTransportationManager()->
       GetFieldManager();

    pFieldMgr->SetDetectorField( pMyMagField ); // ( &myMagField );



    pFieldMgr->SetChordFinder( pChordFinder );

    return    pFieldMgr;
}


G4PropagatorInField *pMagFieldPropagator=0; 
//
// Test Stepping
//
G4bool testG4PropagatorInField(G4VPhysicalVolume*,     // *pTopNode, 
			       G4int             type)
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
    
    pMagFieldPropagator= SetupPropagator(type);

    G4double particleCharge= +1.0;  // in e+ units
    G4double spin=0.0;              // ignore the spin
    G4double magneticMoment= 0.0;   // ignore the magnetic moment

    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0); 

    G4EquationOfMotion* equationOfMotion = 
        ( pMagFieldPropagator->GetChordFinder()
                             ->GetIntegrationDriver()
                             ->GetEquationOfMotion()  );
    
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
    G4ThreeVector EndPosition, EndUnitMomentum;

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
						     located);
	  //       --------------------
	  EndPosition=     pMagFieldPropagator->EndPosition();
	  EndUnitMomentum= pMagFieldPropagator->EndMomentumDir();
	  //       --------
	  
	  if( std::fabs(EndUnitMomentum.mag2() - 1.0) > 1.e-8 )
	    G4cerr << "EndUnitMomentum.mag2() - 1.0 = " <<
	      EndUnitMomentum.mag2() - 1.0 << G4endl;

	  G4ThreeVector MoveVec = EndPosition - Position;

	  // assert( MoveVec.mag() < physStep*(1.+1.e-9) );          

          double move= MoveVec.mag();
          double displaceRatio = 0.0;
          if( physStep > 0. ) {
             displaceRatio = move / physStep - 1.0;
          }
          if( displaceRatio > 1.0e-9 ) {
             G4cerr.precision(8); 
             G4cerr << " Move displaced " << move << " - further than step= "
                    << physStep << " by " << move - physStep << ", "
                    << " a fraction of " << displaceRatio << G4endl;
          }
          G4double maxEpsilon = pMagFieldPropagator -> GetMaximumEpsilonStep();
          assert( displaceRatio < maxEpsilon ); // Revised check condition 
             
	  // G4cout << " testPropagatorInField: After stepI " << istep  << " : " << G4endl;
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
       pMyMagField->          
          ReportStatistics(); 

    }    // ..............................  end for ( iparticle )

    return(1);
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
            << std::setw(18) << "NextVolume" << " "
            << G4endl;
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
       G4cout << std::setw( 5) << Step << " "
	    << std::setw( 9) << Position.x() << " "
	    << std::setw( 9) << Position.y() << " "
	    << std::setw( 9) << Position.z() << " "
	    << std::setw( 9) << EndUnitVelocity.x() << " "
	    << std::setw( 9) << EndUnitVelocity.y() << " "
	      << std::setw( 9) << EndUnitVelocity.z() << " ";
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
       if( startVolume != 0) {
	 G4cout << std::setw(12) << startVolume->GetName() << " ";
       } else {
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
    G4int type, optim, optimSaf;
    G4bool optimiseVoxels=true;
    G4bool optimisePiFwithSafety=true;

    type = 8 ;
    G4cout << " Arguments:  stepper-no  optimise-Voxels optimise-PiF-with-safety" << G4endl;

    if( argc >= 2 ){
       type = atoi(argv[1]);
    } 

    if( argc >=3 ){
      optim= atoi(argv[2]);
      if( optim == 0 ) { optimiseVoxels = false; }
    }

    if( argc >=4 ){
      optimSaf= atoi(argv[3]);
      if( optimSaf == 0 ) { optimisePiFwithSafety= false; }
    }

    G4cout << " Testing with stepper number    " << type << G4endl; 
    G4cout << "             " ; 
    G4cout << " voxel optimisation      " ; 
    // if (optimiseVoxels)   G4cout << "On"; 
    // else                  G4cout << "Off"; 
    G4cout << (optimiseVoxels ? "On" : "Off")  << G4endl;
    G4cout << "             " ; 
    G4cout << " Propagator safety optim " ; 
    // const char* OnOff= (optimisePiFwithSafety ? "on" : "off") ; 
    // G4cout << OnOff << G4endl;
    G4cout << (optimisePiFwithSafety ? "On" : "Off")  << G4endl;

    // Create the geometry & field 
    myTopNode=BuildGeometry();	// Build the geometry
 
    G4Navigator *pNavig= G4TransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    pNavig->SetWorldVolume(myTopNode);

    G4GeometryManager::GetInstance()->CloseGeometry(false);
    
    // pMagFieldPropagator->SetMinimumEpsilonStep();
    // pMagFieldPropagator->SetMaximumEpsilonStep();
       
    // Setup the propagator (will be overwritten by testG4Propagator ...)
    pMagFieldPropagator= SetupPropagator(type);
    G4cout << " Using default values for " 
	   << " Min Eps = "  <<   pMagFieldPropagator->GetMinimumEpsilonStep()
           << " and "
	   << " MaxEps = " <<  pMagFieldPropagator->GetMaximumEpsilonStep()
	   << G4endl; 

    pMagFieldPropagator->SetUseSafetyForOptimization(optimisePiFwithSafety); 

// Do the tests without voxels
    G4cout << " Test with no voxels" << G4endl; 
    testG4PropagatorInField(myTopNode, type);

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels); 
    pMagFieldPropagator->SetVerboseLevel( 1 ); 

// Repeat tests but with full voxels
    G4cout << " Test with full voxels" << G4endl; 

    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    testG4PropagatorInField(myTopNode, type);

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

    testG4PropagatorInField(myTopNode, type);

    
// Repeat tests with full voxels and modified parameters
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-6;
    minEps= 1.0e-6;
    repeatWith( minEps, maxEps, myTopNode, type );
    
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-7;
    minEps= 1.0e-7;
    repeatWith( minEps, maxEps, myTopNode, type );
    
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-8;
    minEps= 1.0e-8;
    repeatWith( minEps, maxEps, myTopNode, type );
    
    G4cout << "Test with more accurate parameters ( voxels on )" << G4endl; 
    maxEps= 1.0e-9;
    minEps= 1.0e-9;
    repeatWith( minEps, maxEps, myTopNode, type );
    
    G4GeometryManager::GetInstance()->OpenGeometry();
    optimiseVoxels = ! optimiseVoxels;
// Repeat tests but with the opposite optimisation choice
    G4cout << " Now test with optimisation " ; 
    if (optimiseVoxels)   G4cout << "on"; 
    else            G4cout << "off"; 
    G4cout << G4endl;

    pMagFieldPropagator->SetUseSafetyForOptimization(optimiseVoxels); 
    testG4PropagatorInField(myTopNode, type);

    G4GeometryManager::GetInstance()->OpenGeometry();

    // Cannot delete G4TransportationManager::GetInstance(); 
    return 0;
}


  
