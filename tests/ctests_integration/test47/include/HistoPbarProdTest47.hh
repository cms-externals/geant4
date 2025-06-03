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
#ifndef HistoPbarProdTest47_H
#define HistoPbarProdTest47_H

#include "TstHistoSet.hh"
#include "TH1F.h"
#include <string>

// fwd declaration
// 
class G4VParticleChange;

class HistoPbarProdTest47 : public TstHistoSet {

public:

  HistoPbarProdTest47(G4String);
  virtual ~HistoPbarProdTest47();

  virtual void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );
  virtual void Write( G4int, G4double );

private:

   std::vector<G4double> fThetaSecOut;
   std::vector<G4double> fCosThetaMinSecOut, fCosThetaMaxSecOut, fDeltaCosThetaSecOut;
   TH1F*                 fMultEvt;
   TH1F*                 fMultEvtWithPbar;
   std::vector<TH1F*>    fInvXSecVsMomPbar;
   std::vector<TH1F*>    fInvXSecTst;


};

#endif
