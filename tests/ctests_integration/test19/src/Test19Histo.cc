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

#include "Test19Histo.hh"
#include "TstReader.hh"

#include "TestNA49Histo.hh"
#include "TestNA61ProtonBeamHisto.hh"
#include "TestNA61PiplusBeamHisto.hh"
#include "TestHARPHisto.hh"
#include "TestSASM6EHisto.hh"
#include "TestIAEAHisto.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


Test19Histo::Test19Histo( const TstReader* pset )
   : TstHisto(pset)
{

   TH1::SetDefaultSumw2();

   if ( pset->GetExpDataSet() == "NA61" )
   {
      if (pset->GetBeamParticle() == "proton")
      { 
         std::cout << "Initialize NA61 histo for proton beam" << std::endl;
	 fHistoSet = new TestNA61ProtonBeamHisto( fHistoTitle );
      }
      else if ( pset->GetBeamParticle() == "pi+" )
      {
         std::cout << "Initialize NA61 histo for pi+ beam" << std::endl;
         fHistoSet = new TestNA61PiplusBeamHisto( fHistoTitle );         
      }
      else
      {
         std::cout << "NA61 BENCHMARK : UNKNOWN BEAM TYPE "  
	           << pset->GetBeamParticle() << std::endl;
         exit(1);
      }
      fHistoDirName  = "na61-histo";
   }
   else if ( pset->GetExpDataSet() == "NA49" )
   {
      fHistoSet = new TestNA49Histo( fHistoTitle );
      // fHistoSet = new Tst23NA49Histo( fHistoTitle );
      fHistoDirName  = "na49-histo";
   }
   else if ( pset->GetExpDataSet() == "HARP" )
   {
      fHistoSet = new TestHARPHisto( fHistoTitle );
      fHistoDirName = "harp-histo";
   }
   else if ( pset->GetExpDataSet() == "SASM6E" )
   {
      fHistoSet = new TestSASM6EHisto( fHistoTitle );
      fHistoDirName = "sasm6e-histo";
   }
   else if ( pset->GetExpDataSet() == "IAEA" )
   {
      fHistoSet = new TestIAEAHisto( fHistoTitle );
      fHistoDirName = "iaea-histo";
   }   
   // fHistoSet->SetDoResDecay(pset->ForseResDecay());
   fHistoSet->SetDoResDecay(fDoResDecay);   

}

TFile* Test19Histo::OpenHistoFile()
{

   // AccessPathName actually returns o (false) is directory ** exists **
   // and 1 (true) if it does NOT
   //
   if ( gSystem->AccessPathName( fHistoDirName.c_str() ) )
   {   
      gSystem->mkdir( fHistoDirName.c_str() );
   }
   
   G4String tmp_bname;
   if ( fBeam == "pi+" )
   {
      tmp_bname = "piplus";
   }
   else if ( fBeam == "pi-" )
   {
      tmp_bname = "piminus";
   }
   else if ( fBeam == "kaon+" )
   {
      tmp_bname = "kplus";
   }
   else if ( fBeam == "kaon-" )
   {
      tmp_bname = "kminus";
   }
   else
   {
      tmp_bname = fBeam;
   }
//.   G4String fname = fHistoDirName + "/" + fBeam + fTarget + fBeamMomentum + fModel;
   G4String fname = fHistoDirName + "/" + tmp_bname + fTarget + fBeamMomentum + fModel;
   
   if ( fJobID > -1 )
   {
      fname += "-" + std::to_string( fJobID );
   }  
   
   if ( fDoResDecay )
   {
      fname += "-with-decays";
   }
   fname += ".root";

   return new TFile( fname.c_str(), "recreate" );

}
