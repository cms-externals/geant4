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
#include "G4ios.hh"

#include "HistoPbarProdTest47.hh"
#include "G4VParticleChange.hh"
#include "G4UnitsTable.hh"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <fstream>
#include <iomanip>

HistoPbarProdTest47::HistoPbarProdTest47( G4String htitle )
   : TstHistoSet( htitle )
{

   fThetaSecOut.push_back(   0.*CLHEP::deg  ); // leave it for later - need sensible calc. of dcth
   fThetaSecOut.push_back(   3.5*CLHEP::deg );
   fThetaSecOut.push_back(  10.5*CLHEP::deg );
   fThetaSecOut.push_back(  10.8*CLHEP::deg );
   fThetaSecOut.push_back(  59.0*CLHEP::deg );
   fThetaSecOut.push_back(  61.0*CLHEP::deg );
   fThetaSecOut.push_back(  90.0*CLHEP::deg );
   fThetaSecOut.push_back(  97.0*CLHEP::deg );
   fThetaSecOut.push_back( 119.0*CLHEP::deg ); 
   
   G4double dTheta = 4.0*CLHEP::deg;
//   G4double dTheta = 3.5*CLHEP::deg;
//   G4double dTheta = 0.15*CLHEP::deg;
   
   fMultEvt         = new TH1F( "MultEvt", "Interaction multiplicity", 50, 0., 50. );
   fMultEvtWithPbar = new TH1F( "MultEvtWithPbar", "Multiplicity of interactions that have pbar", 50, 0., 50. );
   
   for ( unsigned int i=0; i<fThetaSecOut.size(); ++i )
   {
      
      G4double cth1 = std::cos(std::min(fThetaSecOut[i]+dTheta,180.*CLHEP::deg));
      G4double cth2 = std::cos(std::max(fThetaSecOut[i]-dTheta,  0.*CLHEP::deg));
      fDeltaCosThetaSecOut.push_back(std::abs(cth1-cth2));
      fCosThetaMinSecOut.push_back(std::min(cth1,cth2));
      fCosThetaMaxSecOut.push_back(std::max(cth1,cth2));
      
      std::ostringstream tmp, tmp1;
      G4double theta = fThetaSecOut[i]/CLHEP::deg;
      tmp << theta;
      tmp1 << i;
      std::string name  = "XSecVsMomPbarAt" + tmp.str() + "deg";
      std::string title = htitle + " -> pbar + X, at " + tmp.str() + " deg";
      fInvXSecVsMomPbar.push_back( new TH1F( name.c_str(), title.c_str(), 80, 0., 8. ) );
      name = "XSecTst" + tmp1.str();
      fInvXSecTst.push_back( new TH1F( name.c_str(), title.c_str(), 80, 0., 8. ) );
      
   }
   
   // now book histos

}

HistoPbarProdTest47::~HistoPbarProdTest47()
{

   delete fMultEvt;
   delete fMultEvtWithPbar;
   
   for ( unsigned int i=0; i<fInvXSecVsMomPbar.size(); ++i )
   {
      delete fInvXSecVsMomPbar[i];
   }
   for ( unsigned int i=0; i<fInvXSecTst.size(); ++i )
   {
      delete fInvXSecTst[i];
   }

}

void HistoPbarProdTest47::FillEvt(G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& ) 
{

   G4int n = aChange->GetNumberOfSecondaries();
   
   fMultEvt->Fill( (float)n );
   
   bool haspbar = false;
   
   for (G4int j=0; j<n; j++) // loop over secondaries
   {

      G4Track* strk =  aChange->GetSecondary(j);
      const G4DynamicParticle* sdp = strk->GetDynamicParticle();          
      const G4String& pname = sdp->GetDefinition()->GetParticleName();
      
      if ( pname != "anti_proton" ) continue;
      
      haspbar = true;
      
      G4ThreeVector mom = sdp->GetMomentum()/CLHEP::GeV;
      
      G4double plab  = mom.mag();
//      G4double theta = mom.theta();
//      G4double cth   = std::cos(theta);
      G4double cth = strk->GetMomentumDirection().z();
      G4double etot = strk->GetTotalEnergy()/CLHEP::GeV;
      G4double wt   = etot/(plab*plab);
//      G4double wt = 1./plab;

//      G4cout << " found " << pname << " out of " << n << " secondaries" << G4endl;
      
      for ( unsigned int i=0; i<fThetaSecOut.size(); ++i )
      {
         if ( cth > fCosThetaMinSecOut[i] && cth < fCosThetaMaxSecOut[i] )
	 {
	    fInvXSecVsMomPbar[i]->Fill(  plab, wt ); 
	    fInvXSecTst[i]->Fill( plab, wt );
	    // break; // NOTE: do NOT break since theta-bins may overlap
	 }
      }

   } 
   
   if ( haspbar ) fMultEvtWithPbar->Fill( (float)n );  
   
   return;

}

void HistoPbarProdTest47::Write( G4int nevt /* stat */, G4double xsec ) 
{

   fMultEvt->Write();
   fMultEvtWithPbar->Write();
   

// it's passed in in millibarn (exec->GetXSecOnTarget())/millibarn)
// need to convert to microbarn as the 10GeV data are given in these units 
// also need to consider what to do with !GeV data as they're given in nanobarns
//
//   G4double xsec1 = xsec * CLHEP::millibarn;
//   xsec1 /= CLHEP::microbarn;
   
   for ( unsigned int i=0; i<fInvXSecVsMomPbar.size(); ++i )
   {
// the xsec is passed in millibarn (exec->GetXSecOnTarget())/millibarn)
// need to convert to microbarn as the 10GeV data are given in these units 
// also need to consider what to do with !GeV data as they're given in nanobarns
      G4double scale = (xsec*CLHEP::millibarn/CLHEP::microbarn) / ( (double)nevt *  2.*CLHEP::pi*fDeltaCosThetaSecOut[i] ) ;
      fInvXSecVsMomPbar[i]->Scale(scale,"width");
      fInvXSecVsMomPbar[i]->Write();
   }
   for ( unsigned int i=0; i<fInvXSecTst.size(); ++i )
   {
      G4double scale = (xsec*CLHEP::millibarn/CLHEP::microbarn) / ( (double)nevt *  2.*CLHEP::pi*fDeltaCosThetaSecOut[i] ) ;
      fInvXSecTst[i]->Scale(scale,"width");
      fInvXSecTst[i]->Write();
   }
   
   return;

}
