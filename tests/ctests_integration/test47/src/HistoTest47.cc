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

#include "TSystem.h"

#include "HistoTest47.hh"
#include "TstReader.hh"

#include "HistoITEPTest47.hh"
#include "HistoBNLTest47.hh"
#include "HistoMIPSTest47.hh"
#include "HistoPbarProdTest47.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


HistoTest47::HistoTest47( const TstReader* pset )
   : TstHisto(pset)
{

   TH1::SetDefaultSumw2();
   
   // ugly trick to deal with old naming conventions
   //
   fBeamNameTmp = fBeam;
   if ( fBeam == "pi+" )
   {
      fBeamNameTmp = "piplus";
   }
   if ( fBeam == "pi-" )
   {
      fBeamNameTmp = "piminus";
   }

   // another ugly trick to deal with old (original) naming convention:
   // override fHistoTitle
   //
   // fisrt of all, pad fBeamMomentum (with extra zero's), to match old naming style
   //
   std::string::size_type cnt1 = fBeamMomentum.find("GeV");
   G4String tmp = fBeamMomentum.substr( 0, cnt1 );
   G4int counter = 4 - tmp.length();
   for ( G4int i=0; i<counter; ++i )
   {
      tmp += "0";
   }
   fBeamMomentum = tmp + "GeV";
   fHistoTitle = fBeamNameTmp + "+" + fTarget + fModel + "-" + fBeamMomentum;
   
   std::cout << " fHistoTitle = " << fHistoTitle << std::endl;
   
   if ( pset->GetExpDataSet() == "ITEP" )
   {
      fHistoSet = new HistoITEPTest47( fHistoTitle );
   }
   else if ( pset->GetExpDataSet() == "BNL" )
   {
      fHistoSet = new HistoBNLTest47( fHistoTitle );
   }
   else if ( pset->GetExpDataSet() == "MIPS" )
   {
      fHistoSet = new HistoMIPSTest47( fHistoTitle );
   }
   else if ( pset->GetExpDataSet() == "PbarProd" )
   {
      fHistoTitle = fBeamMomentum + " " + fBeam + "+" + fTarget;
      fHistoSet = new HistoPbarProdTest47( fHistoTitle );
   }
   
   // fHistoSet->SetDoResDecay(fDoResDecay);   

}

TFile* HistoTest47::OpenHistoFile()
{

   G4String fname = "./" + fBeamNameTmp + fTarget + fModel + fBeamMomentum; 
   
   if ( fJobID > -1 )
   {
      fname += "-" + std::to_string( fJobID );
   }  
//   if ( fDoResDecay )
//   {
//      fname += "-with-decays";
//   }
   fname += ".root";

   return new TFile( fname.c_str(), "recreate" );

}
