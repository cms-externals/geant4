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

// #include "G4UniformMagField.hh"
#include "G4UniformElectricField.hh"

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
				 const G4VPhysicalVolume*) const
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
  //  virtual void ComputeDimensions(G4Polycone &, const G4int , const G4VPhysicalVolume*) const {}
  //  virtual void ComputeDimensions(G4Polyhedra &, const G4int , const G4VPhysicalVolume*) const {}
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
    
    G4PVPlacement *worldPhys=
         new G4PVPlacement(0,G4ThreeVector(0,0,0), "World",worldLog,
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
    new G4PVPlacement(0,G4ThreeVector( 7.5*m,0,0),
					      "Target 3d",smallBoxLog,
					      worldPhys,false,0);


//  3) Eight small boxes around the origin of the world volume 
//        (in +-X, +-Y & +-Z)
//
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m, 0.3*m,0.3*m), "Target 4a",tinyBoxLog,
					      worldPhys,false,0);
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m,-0.3*m,0.3*m), "Target 4b",tinyBoxLog,
					      worldPhys,false,0);
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m,-0.3*m,0.3*m), "Target 4c",tinyBoxLog,
					      worldPhys,false,0);
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m, 0.3*m,0.3*m), "Target 4d",tinyBoxLog,
					      worldPhys,false,0);

    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m, 0.3*m,-0.3*m), "Target 4e",tinyBoxLog,
					      worldPhys,false,0);
    new G4PVPlacement
          (0,G4ThreeVector( 0.3*m,-0.3*m,-0.3*m), "Target 4f",tinyBoxLog,
					      worldPhys,false,0);
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m,-0.3*m,-0.3*m), "Target 4g",tinyBoxLog,
					      worldPhys,false,0);
    new G4PVPlacement
          (0,G4ThreeVector(-0.3*m, 0.3*m,-0.3*m), "Target 4h",tinyBoxLog,
					      worldPhys,false,0);

    return worldPhys;
}

//  Equation of motion 
// #include "G4Mag_UsualEqRhs.hh"
#include "G4EqMagElectricField.hh"

#include "G4ChordFinder.hh"

#include "G4MagIntegratorDriver.hh"
#include "G4IntegrationDriver.hh"
#include "G4FSALIntegrationDriver.hh"

#include "G4PropagatorInField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"

#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4DoLoMcPriRK34.hh"
#include "G4TsitourasRK45.hh"

#include "G4RK547FEq1.hh"
#include "G4RK547FEq2.hh"
#include "G4RK547FEq3.hh"

#include "G4TCashKarpRKF45.hh"
#include "G4TDormandPrince45.hh"

// G4UniformMagField myMagField(10.*tesla, 0., 0.); 

G4UniformElectricField myElectricField(10.*kilovolt/cm, 0., 0.); 

G4EqMagElectricField *fEquation = 0;

bool gNewTypeDriver= true;   //  Use new templated (RK) integration Driver

constexpr G4int nvar = 8;  // Use 8 variable in order to integrate Time!!

G4FieldManager* SetupField(G4int type)
{
    G4FieldManager   *pFieldMgr;
    G4ChordFinder    *pChordFinder;

    G4double stepMinimum=  1.0e-2 * mm;  // hmin
    G4VIntegrationDriver* pDriver= nullptr;
    
    fEquation = new G4EqMagElectricField(&myElectricField);
 
    pFieldMgr= G4TransportationManager::GetTransportationManager()->
       GetFieldManager();

    pFieldMgr->SetDetectorField( &myElectricField );

    if( !gNewTypeDriver )
    {
       G4MagIntegratorStepper *pStepper;
       G4cout << "Using old type of Integration Driver: G4MagInt_Driver (does not know type of stepper.)"
              << G4endl;
       
       switch ( type ) 
       {
       case 0: G4cout << "ExplicitEuler         1st order stepper used (with old driver) ." << G4endl;
          pStepper = new G4ExplicitEuler( fEquation, nvar ); break;
       case 1: G4cout << "ImplicitEuler         1st order stepper used (with old driver) ." << G4endl;
          pStepper = new G4ImplicitEuler( fEquation, nvar  ); break;
       case 2: G4cout << "SimpleRunge          2nd order stepper used (with old driver) ." << G4endl;
          pStepper = new G4SimpleRunge( fEquation, nvar  ); break;
       case 3: G4cout << "SimpleHeum           3rd order stepper used (with old driver) ." << G4endl;
          pStepper = new G4SimpleHeum( fEquation, nvar  ); break;
       case 4: G4cout << "ClassicalRK4         4th order stepper used (with old driver) ." << G4endl;
          pStepper = new G4ClassicalRK4( fEquation, nvar  ); break;
       case 8: G4cout << "CaskKarpRKF45        5th order stepper used (with old driver) ." << G4endl;
          pStepper = new G4CashKarpRKF45( fEquation, nvar  ); break;
       case 23: G4cout << "BogackiShampine23   3rd order stepper used (with old driver) ." << G4endl;          
          pStepper = new G4BogackiShampine23( fEquation, nvar ); break;
       case 34: G4cout << "DoLoMcPriRK34       4th order stepper used (with old driver) ." << G4endl;                    
          pStepper = new G4DoLoMcPriRK34( fEquation, nvar ); break;         
       case 45:  G4cout << "BogackiShampine45  5th order stepper used (with old driver) ." << G4endl;
          pStepper = new G4BogackiShampine45( fEquation, nvar ); break;
       case 145: G4cout << "TsitourasRK45      5th order stepper used (with old driver) ." << G4endl;          
          pStepper = new    G4TsitourasRK45( fEquation, nvar ); break;
       case 245:
       case 745: G4cout << "DormandPrince745   5th order 7-stage stepper used (with old driver) ." << G4endl;       
          pStepper = new G4DormandPrince745( fEquation, nvar ); break;
       case 56: G4cout << "DormandPrinceRK56   6th order stepper used (with old driver) ." << G4endl;          
          pStepper = new G4DormandPrinceRK56( fEquation, nvar ); break;
       case 78: G4cout << "DormandPrinceRK78   8th order stepper used (with old driver) ." << G4endl;
          pStepper = new G4DormandPrinceRK78( fEquation, nvar ); break;

       case 241: G4cout << "RKF745Eq1  stepper used (with old driver) ." << G4endl;
          pStepper = new G4RK547FEq1(fEquation, nvar ); break;
       case 242: G4cout << "RKF745Eq2  stepper used (with old driver) ." << G4endl;
          pStepper = new G4RK547FEq2(fEquation, nvar ); break;          
       case 243: G4cout << "RKF745Eq3  stepper used (with old driver) ." << G4endl;
          pStepper = new G4RK547FEq3(fEquation, nvar ); break;
          
       case 108: G4cout << "TCaskKarpRKF45    (templated) 5th order stepper used (with old driver) ." << G4endl;
          pStepper = new G4TCashKarpRKF45<G4EqMagElectricField,nvar>( fEquation ); break;
       case 345: G4cout << "TDormandPrince45  (templated) 5th order stepper used (with old driver) ." << G4endl;
          pStepper = new G4TDormandPrince45<G4EqMagElectricField,nvar>( fEquation ); break;
          
          // --- case 9: pStepper = new G4RKG3_Stepper( fEquation, nvar  );    break;
       default: pStepper = 0;
          G4cout << "Chosen stepper " << type << " does not exist. " << G4endl;
          G4Exception("SetupField()", "InvalidArgument", FatalException,
                      "SetupField: incorrect argument for type"); 
       }
       pDriver = new G4MagInt_Driver(stepMinimum, pStepper, 
                             pStepper->GetNumberOfVariables() );
    }
    else
    {
       auto pStepper0 = new G4ExplicitEuler( fEquation, nvar );
       auto pStepper1 = new G4ImplicitEuler( fEquation, nvar  );
       auto pStepper2 = new G4SimpleRunge( fEquation, nvar  );
       auto pStepper3 = new G4SimpleHeum( fEquation, nvar  );
       auto pStepper4 = new G4ClassicalRK4( fEquation, nvar  );
       auto pStepper8 = new G4CashKarpRKF45( fEquation, nvar  );
       auto pStepper23 = new G4BogackiShampine23( fEquation, nvar );
       auto pStepper34 = new G4DoLoMcPriRK34( fEquation, nvar );        
       auto pStepperBS45 = new G4BogackiShampine45( fEquation, nvar );
       auto pStepper145 = new    G4TsitourasRK45( fEquation, nvar );
       auto pStepperDP45 = new G4DormandPrince745( fEquation, nvar );
       auto pStepperDP56 = new G4DormandPrinceRK56( fEquation, nvar );       
       auto pStepperDP78 = new G4DormandPrinceRK78( fEquation, nvar );
       
       auto pStepperEq1  = new G4RK547FEq1( fEquation, nvar );
       
       using TemplatedDoPri5_t = G4TDormandPrince45<G4EqMagElectricField,nvar>;
       TemplatedDoPri5_t*     pStepperTDP45= nullptr;
       using TemplatedCashKarp45_t = G4TCashKarpRKF45<G4EqMagElectricField,nvar>;
       TemplatedCashKarp45_t* pStepperTCashKarp45= nullptr;
       // auto pStepperTDP45 = new TemplatedDoPri5_t( fEquation );
                         // new G4TDormandPrince45<G4EqMagElectricField,nvar>(fEquation);         

       G4cout << "Using new type of Integration Driver: templated G4IntegrationDriver<StepperType> "
              << " or G4FSALIntegrationDriver<FsalStepperType> ." << G4endl;
       // using StepperType= G4CashKarpRKF45;
       // auto pDriver = new G4IntegrationDriver<StepperType>(stepMinimum, pStepper,
       //                                pStepper->GetNumberOfVariables() );
       switch ( type ) 
       {
       case 0:  G4cout << "ExplicitEuler     stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4ExplicitEuler>(stepMinimum, pStepper0, nvar );
          pStepper0 = nullptr;
          break;
       case 1:  G4cout << "ImplicitEuler     stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4ImplicitEuler>(stepMinimum, pStepper1, nvar  );
          pStepper1 = nullptr;
          break;
       case 2:  G4cout << "SimpleRunge       stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4SimpleRunge>(stepMinimum, pStepper2, nvar  );
          pStepper2 = nullptr;
          break;
       case 3:  G4cout << "SimpleHeum        stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4SimpleHeum>(stepMinimum, pStepper3, nvar  );
          pStepper3 = nullptr; 
          break;
       case 4:  G4cout << "ClassicalRK4      stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4ClassicalRK4>(stepMinimum, pStepper4, nvar  );
          pStepper4 = nullptr;
          break;

       case 8:  G4cout << "CashKarpRKF45     stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4CashKarpRKF45>(stepMinimum, pStepper8, nvar  );
          pStepper8 = nullptr;          
          break;
       case 23:  G4cout << "BogackiShampine23 stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4BogackiShampine23>(stepMinimum, pStepper23, nvar );
          pStepper23 = nullptr;          
          break;
       case 34:  G4cout << "DoLoMcPriRK34     stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4DoLoMcPriRK34>(stepMinimum, pStepper34, nvar );
          pStepper34 = nullptr;
          break;
       case 45:  G4cout << "BogackiShampine45 stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4BogackiShampine45>(stepMinimum, pStepperBS45, nvar );
          pStepperBS45 = nullptr;
          break;
       case 145: G4cout << "TsitourasRK45     stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver<    G4TsitourasRK45>(stepMinimum, pStepper145, nvar );
          pStepper145 = nullptr;          
          break;
       case 745: G4cout << "DormandPrince745  stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4DormandPrince745>(stepMinimum, pStepperDP45, nvar );
          pStepperDP45 = nullptr;          
          break;
       case 56:  G4cout << "DormandPrinceRK56 stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4DormandPrinceRK56>(stepMinimum, pStepperDP56, nvar );
          pStepperDP56 = nullptr;
          break;
       case 78:  G4cout << "DormandPrinceRK78 stepper used (with templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4DormandPrinceRK78>(stepMinimum, pStepperDP78, nvar );
          pStepperDP78 = nullptr;          
          break;
       case 44:  G4cout << "ClassicalRK4 stepper used (with initial templated driver) ." << G4endl;
          pDriver = new G4IntegrationDriver< G4ClassicalRK4>(stepMinimum, pStepper4, nvar  );
          pStepper4 = nullptr;          
          G4cout << " Renewing stepper with different instance of G4ClassicalRK4 " << G4endl;
          pDriver->RenewStepperAndAdjust( new G4ClassicalRK4( fEquation, nvar ) );
          G4cerr << " Renewing stepper with different instance of G4ImplicitEuler (MUST FAIL)." << G4endl;
          pDriver->RenewStepperAndAdjust( pStepper1 );  // Should FAIL          
          break;

       case 108:
          G4cout << "Cash-Karp 4th/5th order *templated* stepper used (with simple templated driver) ." << G4endl;
          pStepperTCashKarp45 = new TemplatedCashKarp45_t( fEquation );
          pDriver = new G4IntegrationDriver<TemplatedCashKarp45_t>(stepMinimum, pStepperTCashKarp45, nvar );
          break;

       case 241: G4cout << "RKF745Eq1  stepper used (with FSAL templated driver) ." << G4endl;
          pDriver = new G4FSALIntegrationDriver<G4RK547FEq1>(stepMinimum, pStepperEq1, nvar );
          pStepperEq1 = nullptr;
          break;          
       case 245: G4cout << "DormandPrince745  stepper used (with FSAL templated driver) ." << G4endl;
          pDriver = new G4FSALIntegrationDriver< G4DormandPrince745>(stepMinimum, pStepperDP45, nvar );
          pStepperDP45 = nullptr;       
          break;

       case 345: G4cout << "DormandPrince745  (templated) stepper used (with FSAL templated driver) ." << G4endl;
          pStepperTDP45 = new TemplatedDoPri5_t( fEquation );
          pDriver = new G4FSALIntegrationDriver<TemplatedDoPri5_t>(stepMinimum, pStepperTDP45, nvar );
          pStepperTDP45 = nullptr;  // Redundant if it's not deleted ... 
          break;

	// --- case 9: pDriver = new G4IntegrationDriver< G4RKG3_Stepper>(stepMinimum, pStepper, nvar  );    break;
      default: // pStepper = 0;
	G4cout << "Chosen stepper " << type << " does not exist. " << G4endl;
	G4Exception("SetupField()", "InvalidArgument", FatalException,
                    "SetupField: incorrect argument for stepper type"); 
       }

       // Delete unused steppers -- G4IntegrationDriver templated class must now own stepper
       delete pStepper0 ;
       delete pStepper1 ; 
       delete pStepper2 ; 
       delete pStepper3 ; 
       delete pStepper4 ; 
       delete pStepper8 ; 
       delete pStepper23 ;
       delete pStepper34 ; 
       delete pStepperBS45 ;
       delete pStepper145 ; 
       delete pStepperDP45 ;
       delete pStepperDP56 ;
       delete pStepperDP78 ;
       // delete pStepperTDP45 ;
       delete pStepperEq1 ;       
    }
    pChordFinder = new G4ChordFinder( pDriver ); 

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

    return thePropagator;
}

//  This is Done only for this test program ... the transportation does it.
//  The method is now obsolete -- as propagator in Field has this method,
//    in order to message the correct field manager's chord finder.
//
void  SetChargeMomentumMass(G4double charge, G4double MomentumXc, G4double Mass)
{
   G4FieldManager fieldMgr;

   // fieldMgr= G4TransportationManager::GetTransportationManager()->
   //          GetFieldManager(); 

   static G4ChargeState* pChargeState= new G4ChargeState ( -100, 0, 0 );  ;
   // delete pChargeState;
   // pChargeState= new G4ChargeState ( charge, 0, 0 );  // no Mag Dip Moment
   pChargeState->SetCharge(charge);

    // pMagFieldPropagator->set_magnetic_field();
   fEquation->SetChargeMomentumMass(
		      *pChargeState,  // charge in e+ units
		      MomentumXc,   // Momentum in Mev/c ?
                      Mass );
}

//
// Test Stepping
//
G4bool testG4PropagatorInField(G4VPhysicalVolume *pTopNode, G4int type)
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

    // G4UniformMagField MagField(10.*tesla, 0., 0.);  // Tesla Defined ? 
    G4Navigator   *pNavig= G4TransportationManager::
                    GetTransportationManager()-> GetNavigatorForTracking();
    G4PropagatorInField *pMagFieldPropagator= SetupPropagator(type);

    // pMagFieldPropagator->
    SetChargeMomentumMass(  
			   +1.,                    // charge in e+ units
			   0.5 * proton_mass_c2,    // Momentum in Mev/c ?
			    proton_mass_c2 );
    pNavig->SetWorldVolume(pTopNode);

    // Force the Navigator to do a full search from the top of the tree
    G4bool relativeSearch =false;
    pNavig->LocateGlobalPointAndSetup(G4ThreeVector(10000., 10000., 10000.),
                                      nullptr, relativeSearch, true );

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
       SetChargeMomentumMass(
		      +1,                    // charge in e+ units
		      momentum, 
		      proton_mass_c2); 
       G4cout << G4endl;
       G4cout << "***********************" << G4endl
            << " Starting New Particle with Position " << Position << G4endl 
	    << " and UnitVelocity " << UnitMomentum << G4endl;
       G4cout << " Momentum in GeV/c is "<< (0.5+iparticle*10.0)*proton_mass_c2;
       G4cout << G4endl;


       for( int istep=0; istep < 14; istep++ ){ 
   //        // G4cerr << "UnitMomentum Magnitude is " << UnitMomentum.mag() << G4endl;
	  located= pNavig->LocateGlobalPointAndSetup(Position);
	  //  Is the following better ?? It would need "changes"
	  // located= pMagFieldPropagator->LocateGlobalPointAndSetup(Position);
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
						     safety
#ifdef G4MAG_CHECK_VOLUME
						     ,located);
#else
	                                             );
#endif
	  //       --------------------
	  EndPosition=     pMagFieldPropagator->EndPosition();
	  EndUnitMomentum= pMagFieldPropagator->EndMomentumDir();
	  //       --------
	  
	  if( std::fabs(EndUnitMomentum.mag2() - 1.0) > 1.e-8 )
	    G4cerr << "EndUnitMomentum.mag2() - 1.0 = " <<
	      EndUnitMomentum.mag2() - 1.0 << G4endl;

	  G4ThreeVector MoveVec = EndPosition - Position;
	  assert( MoveVec.mag() < physStep*(1.+1.e-9) );

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
    }    // ..............................  end for ( iparticle )

    return(1);
}

int main(int argc, char **argv)
{
    G4VPhysicalVolume *myTopNode;
    G4int type;
    G4cout << "------------------ Test PropagateElectroMagField: ------------------" << G4endl;

    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);

    type = 4 ;

    if( argc == 2 )
      type = atoi(argv[1]);

    G4cout << "Pass 1: un- voxelised geometry" << G4endl;
    
    gNewTypeDriver= true;   //  Use new templated (RK) integration Driver
    testG4PropagatorInField(myTopNode, type);

    gNewTypeDriver= false;   //  Use old type of integration Driver - not templated
    testG4PropagatorInField(myTopNode, type);
    
// Repeat tests but with full voxels
    G4cout << "Pass 2:     voxelised geometry - i.e. full voxels" << G4endl;
    
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);

    gNewTypeDriver= true;   //  Use new templated (RK) integration Driver
    testG4PropagatorInField(myTopNode, type);
    
    gNewTypeDriver= false;   //  Use old type of integration Driver - not templated    
    testG4PropagatorInField(myTopNode, type);

    G4GeometryManager::GetInstance()->OpenGeometry();

    G4cout << G4endl;   // Add a final newline
    return 0;
}


void report_endPV(G4ThreeVector    Position, 
                  G4ThreeVector, // UnitVelocity,
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

    static double sumStepLen= 0.0;

    if( Step == 0 ) {
       sumStepLen= 0.0;
    }
    sumStepLen += step_len;
    
    if( Step == 0 && verboseLevel <= 3 )
    {
       // G4cout.precision(6);
       // G4cout.setf(ios_base::fixed,ios_base::floatfield);
       G4cout << std::setw( 5) << "Step#" << " "
            << std::setw( 10) << "AccumStep" << " "            
            << std::setw( 12) << "X(mm)" << " "
            << std::setw( 12) << "Y(mm)" << " "  
            << std::setw( 12) << "Z(mm)" << " "
            << std::setw( 10) << " N_x " << " "
            << std::setw( 10) << " N_y " << " "
            << std::setw( 10) << " N_z " << " "
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
       G4cout.precision(8);
       G4cout << std::setw( 3) << Step << " ";
       G4cout << std::setw( 9) << sumStepLen << " ";
       // G4cout.setf(ios_base::fixed,ios_base::floatfield);       
       G4cout.precision(10);       
       G4cout << std::setw( 12) << Position.x() << " "
              << std::setw( 12) << Position.y() << " "
              << std::setw( 12) << Position.z() << " ";
       G4cout.precision(8);
       G4cout << std::setw( 10) << EndUnitVelocity.x() << " "
              << std::setw( 10) << EndUnitVelocity.y() << " "
              << std::setw( 10) << EndUnitVelocity.z() << " ";
	 //    << std::setw( 9) << KineticEnergy << " "
	 //    << std::setw( 9) << EnergyDifference << " "
       G4cout.precision(8);
       G4cout << std::setw( 9) << step_len << " "
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

  
