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
#ifndef TestHARPHisto_h
#define TestHARPHisto_h

#include "TstHistoSet.hh"

#include <vector>

#include "TH1F.h"

class TestHARPHisto : public TstHistoSet
{

public:

   TestHARPHisto( G4String );
   ~TestHARPHisto();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );

private:
   
   TH1D*              fHistoNSec;
   
   std::vector<TH1D*> fHistoSecPiMinusFW; 
   std::vector<TH1D*> fHistoSecPiPlusFW; 
   std::vector<TH1D*> fHistoSecPiMinusLA; 
   std::vector<TH1D*> fHistoSecPiPlusLA;    
   std::vector<TH1D*> fHistoSecProtonFW;
   std::vector<TH1D*> fHistoSecProtonLA;
      
//   std::vector<TH1D*> fHistoSecPiMinusTot; 
//   std::vector<TH1D*> fHistoSecPiPlusTot; 
   
   // mainly for Mu2e - inspired by G4 review in Aug.2015
   //
   std::vector<TH1D*> fHistoSecPiMinusMomBins;
   std::vector<TH1D*> fHistoSecPiPlusMomBins;

   G4int                fNThetaBinsFW;
   G4double             fThetaMinFW;
   G4double             fDeltaThetaFW;   
   G4int                fNThetaBinsLA;
   G4double             fThetaMinLA;
   G4double             fDeltaThetaLA;   

   // mainly for Mu2e - inspired by G4 review in Aug.2015
   //
   G4double*            fMomBinsLA;
   G4int                fNMomBinsLA;
   G4int                fNTheta1;
   G4double             fThetaMin1;
   G4double             fDeltaTheta1;

};

#endif
