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

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"

#include "G4Timer.hh"

//#include "TstDiscreteProcessReader.hh"
#include "Tst47Reader.hh"

#include "Tst47ExecProcessLevel.hh"

#include "G4VParticleChange.hh"

#include "HistoTest47.hh"

#include "G4SystemOfUnits.hh"

using namespace std;

int main(int argc, char** argv) 
{

//   TstDiscreteProcessReader* theConfigReader = new TstDiscreteProcessReader();
   TstDiscreteProcessReader* theConfigReader = new Tst47Reader();

   // Choose the Random engine
   //
//   CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
//   CLHEP::HepRandom::setTheSeed( 1234 ); // just something... to be on the safe side
   
   // Control on input
   if (argc < 2) 
   {
      G4cout << "Input file is not specified! Exit" << G4endl;
      exit(1);
   }
  
   theConfigReader->OpenAppConfig( argv[1] );
   theConfigReader->Help();
  
   do 
   {

      theConfigReader->ProcessConfig();
      if ( theConfigReader->IsDone() ) break; // double protection because 
                                              // in the previous cycle it stopped at #run
           
//      CLHEP::HepRandom::setTheSeed( theConfigReader->GetRndmSeed() ); // reset the seed

      Tst47ExecProcessLevel* exec = new Tst47ExecProcessLevel( theConfigReader );
            
      TstHisto* histo = 0;
      histo = new HistoTest47( theConfigReader );
      if ( !histo )
      {
         G4cout << " Histograms NOT booked" << G4endl;
         exit(1);
      }

      G4VParticleChange* aChange = 0;
      
      G4int NEvts = theConfigReader->GetNEvents();
      
      G4Timer timer;
      timer.Start();
      
      for (G4int iter=0; iter<NEvts; ++iter) 
      {

         if ( theConfigReader->GetVerbosity() > -1 && iter%10000 == 0 )
	 {
	    G4cout << " Processing event # " << iter << G4endl;
	 }
	 
	 aChange = exec->DoEvent();
	 	 
	 histo->FillEvt( aChange, exec->GetBeam()->GetLabV(), exec->GetBeam()->GetLabP() );
	 
         G4int nsec = aChange->GetNumberOfSecondaries();
         for (G4int i=0; i<nsec; ++i) 
         {   
            delete aChange->GetSecondary(i);
         } 
         aChange->Clear();

      } // end loop over events
      
      timer.Stop();
      G4cout << " CPU = " << timer.GetUserElapsed() << G4endl;
      G4cout << " Real Time = " << timer.GetRealElapsed() << G4endl;
      
      histo->Write( theConfigReader->GetNEvents(), (exec->GetXSecOnTarget())/millibarn  );
      
      if ( histo ) delete histo; // delete histograms after writing them out
      
      delete exec;

   } while(!theConfigReader->IsDone()); // end loop over "runs"

   G4cout << "###### End of test #####" << G4endl;
  
   delete theConfigReader; 
   
   return 0;

}
