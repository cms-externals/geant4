//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PropagatorInField.hh,v 1.28 2002/11/09 00:25:08 jacek Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// 
// class G4PropagatorInField 
//
// Class description:
// 
// This class performs the navigation/propagation of a particle/track 
// in a magnetic field. The field is in general non-uniform.
// For the calculation of the path, it relies on the class G4MagTr.
// To create an object, must have an object that calculates the Curved 
// paths and also must know the value of the maximum displacement allowed.
//
// Methods:
//        ComputeStep(..)
//        CalculateStepTimeAndAccuracy(..) 
//        LocateIntersectionPoint(..)

// History:
// -------
// 25.10.96 John Apostolakis,  design and implementation 
// 25.03.97 John Apostolakis,  adaptation for G4Transportation and cleanup
// ------------------------------------------------------------------------

#ifndef G4PropagatorInField_hh 
#define G4PropagatorInField_hh  1

#include "globals.hh"
#include "G4FieldTrack.hh"
#include "G4Navigator.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

class G4VCurvedTrajectoryFilter;

class G4PropagatorInField
{

 public:  // with description

   G4PropagatorInField( G4Navigator    *theNavigator, 
                        G4FieldManager *detectorFieldMgr );
  ~G4PropagatorInField();

   G4double ComputeStep( G4FieldTrack      &pFieldTrack,
                         G4double           pCurrentProposedStepLength,
                         G4double          &pNewSafety, 
                         G4VPhysicalVolume *pPhysVol=0 );
     // Compute the next geometric Step

   inline G4ThreeVector  EndPosition() const;       
   inline G4ThreeVector  EndMomentumDir() const;
   inline G4bool         IsParticleLooping() const;
     // Return the state after the Step

   inline G4double  GetEpsilonStep() const;
     // Relative accuracy for current Step (Calc.)
   inline void      SetEpsilonStep(G4double newEps);
     // The ratio DeltaOneStep()/h_current_step

   inline void SetChargeMomentumMass( G4double charge,     // in e+ units
                                      G4double momentum,   // in Geant4 units
                                      G4double pMass );  

   inline G4ChordFinder* GetChordFinder();

   inline G4int  SetVerboseLevel( G4int verbose );
   inline G4int  Verbose() const;

   inline G4double  GetDeltaIntersection() const;
     // Accuracy for boundary intersection.
   inline G4double  GetDeltaOneStep() const;
     // Accuracy for one tracking/physics step.
                                    
   inline void    SetAccuraciesWithDeltaOneStep( G4double deltaOneStep );  
     // Sets both accuracies for the Global (Detector) field, 
     // maintaining a particular ratio for accuracties 
     // of volume Intersection and Integration (in One Step).

   inline void    SetDeltaIntersection( G4double deltaIntersection );
     // Set accuracy of  intersection of a volume.  (only)
   inline void    SetDeltaOneStep( G4double deltaOneStep );  
     // Set accuracy for integration of one step.   (only)

   inline G4int   GetMaxLoopCount() const;
   inline void    SetMaxLoopCount( G4int new_max );
     // A maximum for the number of steps that a (looping) particle can take.

   void printStatus( const G4FieldTrack&        startFT,
                     const G4FieldTrack&        currentFT, 
                           G4double             requestStep, 
                           G4double             safety,
                           G4int                step, 
                           G4VPhysicalVolume*   startVolume);
     // Print Method - useful mostly for debugging.

   inline G4FieldTrack GetEndState() const;

   inline G4double  GetMinimumEpsilonStep() const;
   inline void      SetMinimumEpsilonStep( G4double newEpsMin );
     // Minimum for Relative accuracy of any Step 

   inline G4double  GetMaximumEpsilonStep() const;
   inline void      SetMaximumEpsilonStep( G4double newEpsMax );

   inline void      SetLargestAcceptableStep( G4double newBigDist );
   inline G4double  GetLargestAcceptableStep();

 public:  // without description

   inline G4FieldManager*  GetCurrentFieldManager();

 public:  // no description

   inline void SetThresholdNoZeroStep( G4int noAct,
                                       G4int noHarsh,
                                       G4int noAbandon );
   inline G4int GetThresholdNoZeroSteps( G4int i ); 

 protected:  // with description

   G4bool LocateIntersectionPoint( 
        const  G4FieldTrack&       curveStartPointTangent,  //  A
        const  G4FieldTrack&       curveEndPointTangent,    //  B
        const  G4ThreeVector&      trialPoint,              //  E
               G4FieldTrack&       intersectPointTangent);  // Output
     // If such an intersection exists, this function 
     // calculate the intersection point of the true path of the particle 
     // with the surface of the current volume (or of one of its daughters). 
     // (Should use lateral displacement as measure of convergence). 

   void PrintStepLengthDiagnostic( G4double      currentProposedStepLength,
                                   G4double      decreaseFactor,
                                   G4double      stepTrial,
                             const G4FieldTrack& aFieldTrack);

 private:

  // ----------------------------------------------------------------------
  //  DATA Members
  // ----------------------------------------------------------------------

   G4FieldManager *fDetectorFieldMgr; 
     // The  Field Manager of the whole Detector.  (default)

   G4FieldManager *fCurrentFieldMgr;
     // The  Field Manager of the current volume (may be the one above.)

   G4Navigator   *fNavigator;

   //  STATE information
   //  -----------------

   G4double    fEpsilonStep;
     // Relative accuracy for current Step (Calc.)

   G4FieldTrack    End_PointAndTangent;
     // End point storage

   G4bool      fParticleIsLooping;

   G4int  fVerboseLevel;
     // For debuging purposes

   //  Values for the small possible relative accuracy of a step
   //  (corresponding to the greatest possible integration accuracy)

   G4double  fEpsilonMinDefault;   // Can be 1.0e-5 to 1.0e-10 ...
   G4double  fEpsilonMaxDefault;   // Can be 1.0e-1 to 1.0e-3 ...
   G4double  fEpsilonMin; 
   G4double  fEpsilonMax;
     // Limits for the Relative accuracy of any Step 

   G4int  fmax_loop_count;

   //  Variables to keep track of "abnormal" case - which causes loop
   //
   G4int     fNoZeroStep;                        //  Counter of zeroStep
   G4int     fActionThreshold_NoZeroSteps;       //  Threshold: above this - act
   G4int     fSevereActionThreshold_NoZeroSteps; //  Threshold to act harshly
   G4int     fAbandonThreshold_NoZeroSteps;      //  Threshold to abandon

   G4double  fFull_CurveLen_of_LastAttempt; 
   G4double  fLast_ProposedStepLength; 
   G4double  fLargestAcceptableStep;

   G4double  fCharge, fInitialMomentumModulus, fMass;

  // Introducing smooth curved trajectory display (jacek 30/10/2002)
public:
  // 
  void SetTrajectoryFilter(G4VCurvedTrajectoryFilter* filter);

  // Access the points which have passed through the filter. The
  // points are stored as ThreeVectors for the initial impelmentation
  // only (jacek 30/10/2002)
  // Responsibility for deleting the points lies with
  // SmoothTrajectoryPoint, which is the points' final
  // destination. The points pointer is set to NULL, to ensure that
  // the points are not re-used in subsequent steps, therefore THIS
  // METHOD MUST BE CALLED EXACTLY ONCE PER STEP. (jacek 08/11/2002)
  G4std::vector<G4ThreeVector>* GimmeTrajectoryVectorAndForgetIt() const;
private:
  // The filter encapsulates the algorithm which selects which
  // intermediate points should be stored in a trajectory. If the
  // pointer is not NULL, then PIF should submit all the intermediate
  // points it calculates, to the filter. The pointer should be set to
  // NULL when no intermediate points need to be stored. This
  // (jacek 04/11/2002)
  G4VCurvedTrajectoryFilter* fpTrajectoryFilter;
};

// ********************************************************************
// Inline methods.
// ********************************************************************

#include "G4PropagatorInField.icc"

#endif 
