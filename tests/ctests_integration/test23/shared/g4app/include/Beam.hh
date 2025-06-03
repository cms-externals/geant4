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
#ifndef Beam_H
#define Beam_H 1

#include "G4ThreeVector.hh" // new - to put in changes for refactoring test47

#include "G4LorentzVector.hh"
#include "G4String.hh"

class Beam
{

   public:
   
      Beam() : fBeamPartName(""), fBeamPartMass(0.), fBeamEnergy(0.),
               fBeamMomentum(G4ThreeVector()), fBeamPosition(G4ThreeVector()), // new for test47
                fLabV(G4LorentzVector()), fLabP(G4LorentzVector()),
		fBeamKineByEvt(false), fdpByp(0.) // new - for test47 
		{ // new - for test47
		  fAngleX[0]=fAngleX[1]=0.;
	          fAngleY[0]=fAngleY[1]=0.;
		  fPosSigma[0]=fPosSigma[1]=0.;
		}
      ~Beam() {}
      
      void SetBeam( G4String name, G4double mass, G4double energy ) { fBeamPartName=name;
                                                                      fBeamPartMass=mass;
						                      fBeamEnergy=energy;
						                      return; }
						 
      void SetLabV( G4LorentzVector  lv ) { fLabV=lv; return; }
      void SetLabP( G4LorentzVector  lp ) { fLabP=lp; return; }
      
// new - for test47
//
      void SetAngleRange( G4double x, G4double dx, G4double y, G4double dy ) { fAngleX[0]=x; 
                                                                               fAngleX[1]=dx; 
									       fAngleY[0]=y;
									       fAngleY[1]=dy;
									       return; }
      void SetBeamMomentum( const G4ThreeVector& bmom ) { fBeamMomentum=bmom; return; } 
      void SetdpByp( G4double dp )       { fdpByp = dp; return; }
      void SetBeamPosition( const G4ThreeVector& pos ) { fBeamPosition=pos; return; }
      void SetPosSigma( G4double dx, G4double dy )  { fPosSigma[0]=dx;
						      fPosSigma[1]=dy;
						      return; }
      
      void ResetBeamEnergy( G4double e ) { fBeamEnergy=e; return; }
      
      void BeamKineByEvtOn() { fBeamKineByEvt=true; return; }
//

      G4String         GetBeamPartName()   const { return fBeamPartName; }
      G4double         GetBeamPartMass()   const { return fBeamPartMass; }
      G4double         GetBeamEnergy()     const { return fBeamEnergy; }
      const G4LorentzVector&  GetLabV()    const { return fLabV; }
      const G4LorentzVector&  GetLabP()    const { return fLabP; }

// new - for test47
//
      const G4ThreeVector&    GetBeamMomentum() const { return fBeamMomentum; } 
      G4bool                  IsBeamKineByEvt() const { return fBeamKineByEvt; }
      G4double                GetAngleX()       const { return fAngleX[0]; }
      G4double                GetdAngleX()      const { return fAngleX[1]; }
      G4double                GetAngleY()       const { return fAngleY[0]; }
      G4double                GetdAngleY()      const { return fAngleY[1]; }
      G4double                GetdpByp()        const { return fdpByp; }
      const G4ThreeVector&    GetBeamPosition() const { return fBeamPosition; }
      G4double                GetBeamSigmaX()   const { return fPosSigma[0]; }
      G4double                GetBeamSigmaY()   const { return fPosSigma[1]; }
//


   private:
   
      G4String        fBeamPartName;
      G4double        fBeamPartMass;
      G4double        fBeamEnergy;
      
      // new addition - related to refactoring test47
      //
      G4ThreeVector   fBeamMomentum;
      G4ThreeVector   fBeamPosition;      
      
      G4LorentzVector fLabV;
      G4LorentzVector fLabP;

      // new addition - related to refactoring test47
      //
      G4bool          fBeamKineByEvt;
      G4double        fAngleX[2];
      G4double        fAngleY[2];
      G4double        fdpByp;
      G4double        fPosSigma[2];
      
};

#endif
