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
#ifndef TestNA61PiplusBeamHisto_h
#define TestNA61PiplusBeamHisto_h

#include "TstHistoSet.hh"

#include <vector>

#include "TH1F.h"

class TestNA61PiplusBeamHisto : public TstHistoSet
{

public:

   TestNA61PiplusBeamHisto( G4String );
   ~TestNA61PiplusBeamHisto();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );

private:

   int      fNPiplusBinsP2019[11];
   double*  fPiplusBinsP2019[11];
   int      fNPiminusBinsP2019[11];
   double*  fPiminusBinsP2019[11];
   int      fNKplusBinsP2019[8];
   double*  fKplusBinsP2019[8];
   int      fNKminusBinsP2019[8];
   double*  fKminusBinsP2019[8];
      
   std::vector<TH1F*> fHistoProtonP2019; 
   std::vector<TH1F*> fHistoPiMinusP2019; 
   std::vector<TH1F*> fHistoPiPlusP2019;
   std::vector<TH1F*> fHistoKPlusP2019;
   std::vector<TH1F*> fHistoKMinusP2019; 
   std::vector<TH1F*> fHistoK0sP2019; 
   std::vector<TH1F*> fHistoLambdaP2019; 

};

#endif
