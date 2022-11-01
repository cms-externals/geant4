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

#include "G4UniformMagField.hh"

#include "G4ChordFinder.hh"

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

#include "G4Mag_SpinEqRhs.hh"

#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"

#include "G4CashKarpRKF45.hh"
// #include "G4RKG3_Stepper.hh"

#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4DoLoMcPriRK34.hh"
#include "G4TsitourasRK45.hh"

G4UniformMagField myMagField(10.*tesla, 0., 0.); 
G4Mag_SpinEqRhs *gEquation = new G4Mag_SpinEqRhs(&myMagField);

G4FieldManager* SetupField(G4int type)
{
    G4FieldManager   *pFieldMgr;
    G4ChordFinder    *pChordFinder;

    auto fEquation = new G4Mag_SpinEqRhs(&myMagField);
    gEquation = fEquation;
    
    G4MagIntegratorStepper *pStepper;

    const int ncompspin=12;

    switch ( type ) 
    {
      case 0: pStepper = new G4ExplicitEuler( fEquation, ncompspin ); break;
      case 1: pStepper = new G4ImplicitEuler( fEquation, ncompspin ); break;
      case 2: pStepper = new G4SimpleRunge( fEquation, ncompspin ); break;
      case 3: pStepper = new G4SimpleHeum( fEquation, ncompspin ); break;
      case 4: pStepper = new G4ClassicalRK4( fEquation, ncompspin ); break;
      case 8: pStepper = new G4CashKarpRKF45( fEquation, ncompspin ); break;

      case 23: pStepper = new G4BogackiShampine23( fEquation, ncompspin ); break;
      case 34: pStepper = new G4DoLoMcPriRK34( fEquation, ncompspin ); break;         
      case 45: pStepper = new G4BogackiShampine45( fEquation, ncompspin ); break;
      case 145: pStepper = new    G4TsitourasRK45( fEquation, ncompspin ); break;
      case 745: pStepper = new G4DormandPrince745( fEquation, ncompspin ); break;
      case 56: pStepper = new G4DormandPrinceRK56( fEquation, ncompspin ); break;
      case 78: pStepper = new G4DormandPrinceRK78( fEquation, ncompspin ); break;
         
      default: pStepper = new G4ClassicalRK4( fEquation, ncompspin ); break;
    }
    
    pFieldMgr= G4TransportationManager::GetTransportationManager()->
       GetFieldManager();

    pFieldMgr->SetDetectorField( &myMagField );

    pChordFinder = new G4ChordFinder( &myMagField,
				      1.0e-2 * mm,
				      pStepper);

    pFieldMgr->SetChordFinder( pChordFinder );

    return    pFieldMgr;
}

G4PropagatorInField*  SetupPropagator( G4int type)
{
    // G4FieldManager* fieldMgr= 
    SetupField( type) ;

    // G4ChordFinder  theChordFinder( &MagField, 0.05*mm ); // Default stepper
 
    G4PropagatorInField *thePropagator = 
      G4TransportationManager::GetTransportationManager()->
       GetPropagatorInField ();

    // Let us test the new Minimum Epsilon Step functionality
    thePropagator -> SetMinimumEpsilonStep( 1.0e-8 ) ; 
    thePropagator -> SetMaximumEpsilonStep( 1.0e-8 ) ; 

    return thePropagator;
}

// //  This is Done only for this test program ... the transportation does it.
// //  The method is now obsolete -- as propagator in Field has this method,
// //    in order to message the correct field manager's chord finder.
// //
// void  SetChargeMomentumMass( G4ChordFinder* pChordFinder, 
//                              G4double charge, G4double MomentumXc, G4double Mass)
// {

//     // G4ChordFinder*  pChordFinder= 0; 
//     G4EquationOfMotion* equation= 0;

//     //  pChordFinder= G4TransportationManager::GetTransportationManager()->
//     //		   GetFieldManager()->GetChordFinder();

//     G4MagInt_Driver* driver= pChordFinder->GetIntegrationDriver();
//     equation = driver->GetEquationOfMotion(); 

// // equation GetEquationOfMotion(); 

//     G4ChargeState chargeSt(charge, 0.0, 0.5 ); // Assume M-dipole=0, spin=0.5
//     equation->SetChargeMomentumMass(
//                       chargeSt, 
// 		      MomentumXc,   // Momentum in Mev/c ?
//                       Mass );
// }

//
// Test Stepping
//
G4bool testG4PropagatorInField(G4VPhysicalVolume *pTopNode, G4int ) // type
{
    void report_endPV(G4ThreeVector    Position, 
                  G4ThreeVector UnitVelocity,
                  G4ThreeVector Spin,
		  G4double step_len, 
                  G4double physStep, 
                  G4double safety,
		  G4ThreeVector EndPosition, 
                  G4ThreeVector EndUnitVelocity,
                  G4ThreeVector EndSpin,
                  G4int             Step, 
                  G4VPhysicalVolume* startVolume);

    G4UniformMagField MagField(10.*tesla, 0., 0.);  // Tesla Defined ? 
    G4TransportationManager* transportMgr = G4TransportationManager::GetTransportationManager();
    G4Navigator* pNavig=  transportMgr-> GetNavigatorForTracking();
    // G4Navigator   *pNavig= G4TransportationManager::
    //               GetTransportationManager()-> GetNavigatorForTracking();

    // G4PropagatorInField *pMagFieldPropagator= SetupPropagator(type);
    G4PropagatorInField *pMagFieldPropagator= 
                            transportMgr->GetPropagatorInField(); 
      
    // pMagFieldPropagator->SetChargeMomentumMass(
    G4ChargeState chargeSt( +1.,         // charge in e+ units
                            0.0,         // M-dipole=0
                            0.5 );       // Polarisation    
    gEquation->SetChargeMomentumMass(  
                            chargeSt, 
			    0.1*GeV,                // Momentum in Mev/c ?
			    0.105658387*GeV );
    pNavig->SetWorldVolume(pTopNode);

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

    const G4double threshold= 1.e-6; 

    if( std::fabs(UnitMomentum.mag() - 1.0) > threshold ) 
    {
      G4cout << "UnitMomentum.mag() - 1.0 = " << UnitMomentum.mag() - 1.0 <<
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

       G4double momentum_val= (0.5+iparticle*1.0) * 0.1*GeV;  // As energy/c
       G4double rest_mass   = 0.105658387*GeV;                // A muon
                                             
       G4double kineticEnergy =  momentum_val*momentum_val /
                  ( std::sqrt( momentum_val*momentum_val + rest_mass * rest_mass ) 
		    + rest_mass );
       G4double labTof= 10.0*ns, properTof= 0.1*ns;
       // pMagFieldPropagator->SetChargeMomentumMass(
       G4ChargeState chargeState( +1.,         // charge in e+ units
                                  0.0,         // M-dipole=0
                                  0.5 );       // Polarisation           
       gEquation->SetChargeMomentumMass( chargeState,
                                         momentum_val, 
                                         rest_mass);

       G4double  beta = momentum_val / std::sqrt( rest_mass*rest_mass + momentum_val*momentum_val );
       G4double       velocity_magnitude = beta * c_light;
       // G4ThreeVector  Velocity = UnitMomentum * velocity_magnitude;

       G4cout << G4endl;
       G4cout << "Test PropagateMagField: ***********************" << G4endl
            << " Starting New Particle with Position " << Position << G4endl 
	    << " and UnitVelocity " << UnitMomentum << G4endl;
       G4cout << " Momentum in MeV/c is "<< (0.5+iparticle*1.0)*0.1*GeV/MeV;
       G4cout << G4endl;

   //  G4ThreeVector initialSpin = UnitMomentum; 
   //  G4ThreeVector initialSpin(UnitMomentum.y(),-UnitMomentum.x(),UnitMomentum.z()); 
       G4ThreeVector initialSpin = (UnitMomentum.cross(G4ThreeVector(1.0,0.,0.))).unit(); 

       G4cout << " Initial spin is set to = " << initialSpin << G4endl;
       G4cout << " Initial dot product of spin and momentum = " 
              << initialSpin.dot(UnitMomentum) << G4endl;

       for( int istep=0; istep < 14; istep++ ){ 
   //        // G4cout << "UnitMomentum Magnitude is " << UnitMomentum.mag() << G4endl;
	  located= pNavig->LocateGlobalPointAndSetup(Position);
	  //  Is the following better ?? It would need "changes"
	  // located= pMagFieldPropagator->LocateGlobalPointAndSetup(Position);
	  // G4cout << "Starting Step " << istep << " in volume " 
	       // << located->GetName() << G4endl;

          G4ThreeVector spinValue( initialSpin ); 

          G4FieldTrack  initTrack( Position, 
				   UnitMomentum,
				   0.0,            // starting S curve len
				   kineticEnergy,
				   rest_mass,
				   velocity_magnitude,
				   labTof, 
				   properTof,
				   &spinValue
				   ); 

	  step_len=pMagFieldPropagator->ComputeStep( initTrack, 
						     physStep, 
						     safety
						     ,located);
	  //       --------------------
	  EndPosition=     pMagFieldPropagator->EndPosition();
	  EndUnitMomentum= pMagFieldPropagator->EndMomentumDir();
	  //       --------
	  G4FieldTrack  EndFieldTrack= pMagFieldPropagator->GetEndState();
	  G4ThreeVector EndSpin=         EndFieldTrack.GetSpin();
          EndUnitMomentum  = EndFieldTrack.GetMomentumDir();

//          G4cout << " EndPosition " << EndPosition << G4endl;
//          G4cout << " EndUnitMomentum " << EndUnitMomentum << G4endl;
//          G4cout << " initialSpin " << initialSpin.mag() << G4endl;
//          G4cout << " EndSpin     " << EndSpin.mag()     << G4endl;

	  if( std::fabs(EndUnitMomentum.mag2() - 1.0) > threshold )
	    G4cout << "EndUnitMomentum.mag2() - 1.0 = " <<
	      EndUnitMomentum.mag2() - 1.0 << G4endl;

          // In this case spin should be parallel (equal) to momentum
          G4double endDot= EndSpin.dot(EndUnitMomentum) ;
          G4cout << " dot product of spin and momentum = " << endDot << G4endl;
	  if( std::fabs(endDot) > threshold ){
	     G4cout << " $$$$$$$$$$$$ Spin dot Momentum is above threshold of "
		    << threshold << G4endl ;
	     G4cout << " Spin dot UnitMomentum= " << endDot << " ";
	     G4cout << " Spin magnitude= " << EndSpin.mag() << " "
		    << "(-1=" << (EndSpin.mag()-1.0) << ") ";
	     G4cout << " UnitMom mag= " << EndUnitMomentum.mag() << " "
		    << "(-1=" << (EndUnitMomentum.mag()-1.0) << ") ";
	     G4cout << G4endl;
	  }
	  G4ThreeVector MoveVec = EndPosition - Position;
	  assert( MoveVec.mag() < physStep*(1.+1.e-9) );

	  // G4cout << " testPropagatorInField: After stepI " << istep  << " : " << G4endl;
	  report_endPV(Position, UnitMomentum, initialSpin, step_len, physStep, safety,
		       EndPosition, EndUnitMomentum, EndSpin, istep, located );

	  assert(safety>=0);
	  pNavig->SetGeometricallyLimitedStep();
	  // pMagFieldPropagator->SetGeometricallyLimitedStep();

	  Position=     EndPosition;
          // Velocity=     EndVelocity;
	  UnitMomentum= EndUnitMomentum;
	  initialSpin = EndSpin;
          
	  physStep *= 2.; 
       } // ...........................  end for ( istep )
    }    // ..............................  end for ( iparticle )

    return(1);
}

// int main(int argc, char** argv)
int main(int argc, char **argv)
{
    G4VPhysicalVolume *myTopNode;
    G4int type;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);

    type = 8 ;

    if( argc == 2 )
      type = atoi(argv[1]);

    // Setup the Propagator and register it with the Transportation Manager
    G4PropagatorInField *pMagFieldPropagator= SetupPropagator(type);

    G4cout << " Using the following values for " 
	   << " Min Eps = " <<  pMagFieldPropagator->GetMinimumEpsilonStep()
           << " and "
	   << " Max Eps = " <<  pMagFieldPropagator->GetMaximumEpsilonStep()
	   << G4endl; 

// Do the tests without voxels
    G4cout << " Test with no voxels" << G4endl; 
    testG4PropagatorInField(myTopNode, type);

// Repeat tests but with full voxels
    G4cout << " Test with full voxels" << G4endl; 
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    testG4PropagatorInField(myTopNode, type);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}


void report_endPV(G4ThreeVector    Position, 
                  G4ThreeVector , // UnitVelocity,
                  G4ThreeVector , // Spin,
		  G4double step_len, 
                  G4double physStep, 
                  G4double safety,
		  G4ThreeVector EndPosition, 
                  G4ThreeVector EndUnitVelocity,
                  G4ThreeVector EndSpin,
                  G4int             Step, 
                  G4VPhysicalVolume* startVolume)
		  //   G4VPhysicalVolume* endVolume)
{
    const G4int verboseLevel=1;
    
    if( Step == 0 && verboseLevel <= 3 )
    {
       G4cout.precision(3);
       // G4cout.setf(ios_base::fixed,ios_base::floatfield);
       G4cout << std::setw( 5) << "Step#" << " "
            << std::setw( 9) << "X(mm)" << " "
            << std::setw( 9) << "Y(mm)" << " "  
            << std::setw( 9) << "Z(mm)" << " "
            << std::setw( 7) << " N_x " << " "
            << std::setw( 7) << " N_y " << " "
            << std::setw( 7) << " N_z " << " "
            << std::setw( 7) << " S_x " << " "
            << std::setw( 7) << " S_y " << " "
            << std::setw( 7) << " S_z " << " "
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
       G4cout.precision(3);
       G4cout << std::setw( 5) << Step << " "
	    << std::setw( 9) << Position.x() << " "
	    << std::setw( 9) << Position.y() << " "
	    << std::setw( 9) << Position.z() << " "
	    << std::setw( 7) << EndUnitVelocity.x() << " "
	    << std::setw( 7) << EndUnitVelocity.y() << " "
	    << std::setw( 7) << EndUnitVelocity.z() << " "
	    << std::setw( 7) << EndSpin.x() << " "
	    << std::setw( 7) << EndSpin.y() << " "
	    << std::setw( 7) << EndSpin.z() << " "
	 //    << std::setw( 9) << KineticEnergy << " "
	 //    << std::setw( 9) << EnergyDifference << " "
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
       {
	 G4cout << std::setw(12) << endVolume()->GetName() << " ";
       } 
       else 
       {
	 G4cout << std::setw(12) << "OutOfWorld" << " ";
       }
#endif
       G4cout << G4endl;
    }
}

int readin_particle( )
{
 static const
 double pmass[5] = {
                    0.00051099906 ,         //  electron
                    0.105658389   ,         //  muon
                    0.13956995    ,         //  pion
                    0.493677      ,         //  kaon
                    0.93827231              //  proton
                   } ;
 int pCharge, i ;
 double pMomentum, pTeta, pPhi, h ;
 G4cout<<"Enter particle type: 0 - electron, 1 - muon, 2 - pion, \n"
     <<"3 - kaon, 4 - proton "<< G4endl ;
 G4cin>>i ;
 double pMass = pmass[i] ;
 G4cout<<"Enter particle charge in units of the positron charge "<< G4endl ;
 G4cin>>pCharge ;
 G4cout<<"Enter particle momentum in GeV/c"<<G4endl ;
 G4cin>>pMomentum ;
 G4cout<<"Enter particle teta & phi in degrees"<<G4endl ;
 G4cin>>pTeta ;
 G4cin>>pPhi ;
 G4cout<<"Enter particle Step in centimeters"<<G4endl ;
 G4cin>>h ;

 h *=  10.; // G4 units are in millimeters.

 double betaGamma = pMomentum/pMass ;
 double pSpeed = betaGamma*c_light/std::sqrt(1 + betaGamma*betaGamma) ;
 double pEnergy = pMomentum*c_light/pSpeed ;
        pEnergy *= 1.60217733e-10  ; // energy in J (SI units)
 pTeta *= pi/180 ;
 pPhi  *= pi/180 ;

#if 0
 for(i=0;i<3;i++) ystart[i] = 0 ;            // initial coordinates
 ystart[3] = pSpeed*std::sin(pTeta)*std::cos(pPhi) ;   // and speeds
 ystart[4] = pSpeed*std::sin(pTeta)*std::sin(pPhi) ;
 ystart[5] = pSpeed*std::cos(pTeta) ;
#endif

 return 1;
}

  
