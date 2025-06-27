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

#include "TestSASM6EHisto.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

TestSASM6EHisto::TestSASM6EHisto( G4String htitle )
   : TstHistoSet(htitle)
{

  G4String title;
  G4String outcome;
  
  // at most, 7 bins, with bin center being 30, 40, 50, 60, 70, 80, 88
  //
  fNBins = 7;
  
  fBins = new double[fNBins+1];
  fBins[0] = 25.;
  for ( int i=1; i<7; ++i )
  {
     fBins[i] = fBins[0] + 10.*((double)i);
  }
  fBins[7] = 91.; // (this will make center of bin at 88, i.e. 85-91) 
  
  fDeltaPt = 0.1;
  
  outcome = " -> pi+ + X ";
  title = htitle + outcome + "(pt=0.18GeV/c)";
  fHistoSecPiPlus.push_back( new TH1F( "piplus_pt180", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.3GeV/c)";
  fHistoSecPiPlus.push_back( new TH1F( "piplus_pt300", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.5GeV/c)";
  fHistoSecPiPlus.push_back( new TH1F( "piplus_pt500", title.c_str(), fNBins, fBins ) );
  
  outcome = " -> pi- + X ";
  title = htitle + outcome + "(pt=0.18GeV/c)";
  fHistoSecPiMinus.push_back( new TH1F( "piminus_pt180", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.3GeV/c)";
  fHistoSecPiMinus.push_back( new TH1F( "piminus_pt300", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.5GeV/c)";
  fHistoSecPiMinus.push_back( new TH1F( "piminus_pt500", title.c_str(), fNBins, fBins ) );
  
  outcome = " -> K+ + X ";
  title = htitle + outcome + "(pt=0.18GeV/c)";
  fHistoSecKPlus.push_back( new TH1F( "kplus_pt180", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.3GeV/c)";
  fHistoSecKPlus.push_back( new TH1F( "kplus_pt300", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.5GeV/c)";
  fHistoSecKPlus.push_back( new TH1F( "kplus_pt500", title.c_str(), fNBins, fBins ) );
  
  outcome = " -> K- + X ";
  title = htitle + outcome + "(pt=0.18GeV/c)";
  fHistoSecKMinus.push_back( new TH1F( "kminus_pt180", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.3GeV/c)";
  fHistoSecKMinus.push_back( new TH1F( "kminus_pt300", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.5GeV/c)";
  fHistoSecKMinus.push_back( new TH1F( "kminus_pt500", title.c_str(), fNBins, fBins ) );

  outcome = " -> proton + X ";
  title = htitle + outcome + "(pt=0.18GeV/c)";
  fHistoSecProton.push_back( new TH1F( "proton_pt180", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.3GeV/c)";
  fHistoSecProton.push_back( new TH1F( "proton_pt300", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.5GeV/c)";
  fHistoSecProton.push_back( new TH1F( "proton_pt500", title.c_str(), fNBins, fBins ) );

  outcome = " -> antiproton + X ";
  title = htitle + outcome + "(pt=0.18GeV/c)";
  fHistoSecAntiProton.push_back( new TH1F( "antiproton_pt180", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.3GeV/c)";
  fHistoSecAntiProton.push_back( new TH1F( "antiproton_pt300", title.c_str(), fNBins, fBins ) );
  title = htitle + outcome + "(pt=0.5GeV/c)";
  fHistoSecAntiProton.push_back( new TH1F( "antiproton_pt500", title.c_str(), fNBins, fBins ) );
  
}

TestSASM6EHisto::~TestSASM6EHisto()
{

   for (size_t i=0; i<fHistoSecPiPlus.size(); i++ )
   {
      delete fHistoSecPiPlus[i];
   }
   for (size_t i=0; i<fHistoSecPiMinus.size(); i++ )
   {
      delete fHistoSecPiMinus[i];
   }
   for ( size_t i=0; i<fHistoSecKPlus.size(); i++ )
   {
      delete fHistoSecKPlus[i];
   }
   for ( size_t i=0; i<fHistoSecKMinus.size(); i++ )
   {
      delete fHistoSecKMinus[i];
   }
   for (size_t i=0; i<fHistoSecProton.size(); i++ )
   {
      delete fHistoSecProton[i];
   }
   for (size_t i=0; i<fHistoSecAntiProton.size(); i++ )
   {
      delete fHistoSecAntiProton[i];
   }

}

void TestSASM6EHisto::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& ) 
{

   G4int NSec = aChange->GetNumberOfSecondaries();
           
   for (G4int i=0; i<NSec; i++) 
   {

      G4Track* strk =  aChange->GetSecondary(i);
      const G4DynamicParticle* sdp = strk->GetDynamicParticle();          
      const G4String& pname = sdp->GetDefinition()->GetParticleName();
					
      if (    pname == "pi+"    || pname == "pi-" 
           || pname == "proton" || pname == "anti_proton"
	   || pname == "kaon+"  || pname == "kaon-" )
      {
	
         G4ThreeVector mom = sdp->GetMomentum()/CLHEP::GeV;    
	 G4double pt = mom.perp();
	
	 int id = -1;
	
	 if ( fabs( pt - 0.18 ) < fDeltaPt/2. )
	 {
	    id = 0;
	 }
	 else if ( fabs( pt - 0.3 ) < fDeltaPt/2. )
	 {
	    id = 1;
	 }
	 else if ( fabs( pt - 0.5 ) < fDeltaPt/2. )
	 {
	    id = 2;
	 }
	 if ( id == -1 ) continue;

         G4double plab  = mom.mag();
	 
	 if ( plab < fBins[0] ) continue;
	 if ( plab > fBins[fNBins] ) continue;
	 
	 int ib = -1;
	 for ( int j=0; j<fNBins; ++j )
	 {
	    if ( plab > fBins[j] && plab <= fBins[j+1] )
	    {
	       ib = j;
	       break;
	    }
	 }
	 if ( ib == -1 ) continue; // nothing found
	 
         G4double etot = strk->GetTotalEnergy()/CLHEP::GeV;
	 
         // G4double wt   = etot/(plab*plab);	 
	 // or ???
	 // G4double wt = 1./plab;
	 //		 
	 // wt /= ( 2. * CLHEP::pi * (fDeltaPt/plab) );
	 // or ???
	 // wt /= (fDeltaPt/plab);
	 
	 // try the dp3 stuff
	 //
	 G4double ptmin = pt - fDeltaPt/2.;
	 G4double ptmax = pt + fDeltaPt/2.;
	 G4double dpt2 = ptmax*ptmax - ptmin*ptmin;
	 G4double pzmin = std::sqrt( fBins[ib]*fBins[ib] - ptmax*ptmax );
	 G4double pzmax = std::sqrt( fBins[ib+1]*fBins[ib+1] - ptmin*ptmin );
	 G4double dpz = std::fabs( pzmax - pzmin );
	 G4double dp3 = CLHEP::pi * dpz * dpt2 ;
	 
	 G4double wt = etot / dp3; 
	 
//	 std::cout << " etot = " << etot << " plab = " << plab << " pt = " << pt << std::endl;
//	 std::cout << " dpt2 = " << dpt2 << " dpz = " << dpz << " wt = " << wt << std::endl;
	 	 	 
	if ( pname == "pi+" )
	{
	   if ( id >= 0 && id < 3 ) fHistoSecPiPlus[id]->Fill( plab, wt );
	}
	if ( pname == "pi-" )
	{
	   if ( id >=0 && id < 3 ) fHistoSecPiMinus[id]->Fill( plab, wt );
	}
	if ( pname == "proton" )
	{	   
	   if ( id >=0 && id < 3 ) fHistoSecProton[id]->Fill( plab, wt );
	}
	if ( pname == "anti_proton" )
	{	   
	   if ( id >=0 && id < 3 ) fHistoSecAntiProton[id]->Fill( plab, wt );
	}
	if ( pname == "kaon+" )
	{	   
	   if ( id >=0 && id < 3 ) fHistoSecKPlus[id]->Fill( plab, wt );
	}
	if ( pname == "kaon-" )
	{	   
	   if ( id >=0 && id < 3 ) fHistoSecKMinus[id]->Fill( plab, wt );
	}
      }
               
   }
      
   return;
   
}

void TestSASM6EHisto::Write( G4int stat, G4double xsec )
{

   double scale = xsec / ( (double)stat ) ;
   
   std::cout << " xsec = " << xsec << std::endl;
   std::cout << " nevents = " << stat << std::endl;
   
   std::cout << " scale = " << scale << std::endl;
   
   // NOTE: no need to scale with the width because 
   //       it's already taken into account by the dp3 ("cell size")
      
   for ( size_t i=0; i<fHistoSecPiPlus.size(); ++i )
   {
      //fHistoSecPiPlus[i]->Scale(scale,"width");
      fHistoSecPiPlus[i]->Scale(scale);
      fHistoSecPiPlus[i]->Write();
   }
  
   for ( size_t i=0; i<fHistoSecPiMinus.size(); ++i )
   {
      //fHistoSecPiMinus[i]->Scale(scale,"width");
      fHistoSecPiMinus[i]->Scale(scale);
      fHistoSecPiMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecKPlus.size(); ++i )
   {
      //fHistoSecKPlus[i]->Scale(scale,"width");
      fHistoSecKPlus[i]->Scale(scale);
      fHistoSecKPlus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecKMinus.size(); ++i )
   {
      //fHistoSecKMinus[i]->Scale(scale,"width");
      fHistoSecKMinus[i]->Scale(scale);
      fHistoSecKMinus[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecProton.size(); ++i )
   {
      //fHistoSecProton[i]->Scale(scale,"width");
      fHistoSecProton[i]->Scale(scale);
      fHistoSecProton[i]->Write();
   }

   for ( size_t i=0; i<fHistoSecAntiProton.size(); ++i )
   {
      //fHistoSecAntiProton[i]->Scale(scale,"width");
      fHistoSecAntiProton[i]->Scale(scale);
      fHistoSecAntiProton[i]->Write();
   }
   
   return;

}
