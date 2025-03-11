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
#ifndef Tst23NA61Histo_h
#define Tst23NA61Histo_h

#include "TstHistoSet.hh"
// #include "TstHistoSetForNu.hh"

#include <vector>
#include <map>

//#include "TProfile.h"
#include "TH1.h"

class Tst23NA61Histo : public TstHistoSet
{

public:

   Tst23NA61Histo( G4String );
   ~Tst23NA61Histo();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );


private:
      
   TH1D*              fHistoNSec;
   
   //
   // FIXME !!!
   // 
   // Causing compilation warning as is not used for the moment;
   // histos are booked with regular binnig, but it's likely to
   // change in the near future
   // int*     fNSecProtonThetaBins;
   // double** fSecProtonBins;
   //
   int*     fNSecPiplusThetaBins;
   double** fSecPiplusBins;
   int*     fNSecPiminusThetaBins;
   double** fSecPiminusBins;

   int      fNProtonBinsP2015[10];
   double*  fProtonBinsP2015[10];
   int      fNPiplusBinsP2015[11];
   double*  fPiplusBinsP2015[11];
   int      fNPiminusBinsP2015[11];
   double*  fPiminusBinsP2015[11];
   int      fNKplusBinsP2015[8];
   double*  fKplusBinsP2015[8];
   int      fNKminusBinsP2015[8];
   double*  fKminusBinsP2015[8];
   int      fNK0sBinsP2015[7];
   double*  fK0sBinsP2015[7];
   int      fNLambdaBinsP2015[7];
   double*  fLambdaBinsP2015[7];

   std::vector<TH1F*> fHistoSecProton; 
   std::vector<TH1F*> fHistoSecPiMinus; 
   std::vector<TH1F*> fHistoSecPiPlus;
   std::vector<TH1F*> fHistoSecKPlus;
   std::vector<TH1F*> fHistoSecPiPlus2; 

   std::vector<TH1F*> fHistoProtonP2015; 
   std::vector<TH1F*> fHistoPiMinusP2015; 
   std::vector<TH1F*> fHistoPiPlusP2015;
   std::vector<TH1F*> fHistoKPlusP2015;
   std::vector<TH1F*> fHistoKMinusP2015; 
   std::vector<TH1F*> fHistoK0sP2015; 
   std::vector<TH1F*> fHistoLambdaP2015; 
   
};

#endif
