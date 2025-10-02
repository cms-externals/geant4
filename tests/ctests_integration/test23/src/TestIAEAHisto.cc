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
#include "TestIAEAHisto.hh"
#include "Tst23ParticleChange.hh"

#include "G4VParticleChange.hh"
#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

TestIAEAHisto::TestIAEAHisto( G4String htitle )
   : TstHistoSet(htitle)
{

   fHistoNSec = new TH1D( "NSec", htitle.c_str(), 100, 0., 100. );

   // NOTE (JVY): this binning will be goof for 1.5GeV & 3GeV/c proton beam
   //             for 0.8GeV beam the binning is different (much finer)
   
   
   fNEKinBins = 27;
   fEKinBins = new double[fNEKinBins+1];
   
   fEKinBins[0]  =    1.;
   fEKinBins[1]  =    1.5;
   fEKinBins[2]  =    2.;
   fEKinBins[3]  =    3.;
   fEKinBins[4]  =    4.;
   fEKinBins[5]  =    5.;
   fEKinBins[6]  =    6.;
   fEKinBins[7]  =    8.;
   fEKinBins[8]  =   10.;
   fEKinBins[9]  =   15.;
   fEKinBins[10] =   20.;
   fEKinBins[11] =   30.;
   fEKinBins[12] =   40.;
   fEKinBins[13] =   50.;
   fEKinBins[14] =   60.;
   fEKinBins[15] =   80.;
   fEKinBins[16] =  100.;
   fEKinBins[17] =  150.;
   fEKinBins[18] =  200.;
   fEKinBins[19] =  300.;
   fEKinBins[20] =  400.;
   fEKinBins[21] =  500.;
   fEKinBins[22] =  600.;
   fEKinBins[23] =  800.;
   fEKinBins[24] = 1000.;
   fEKinBins[25] = 1500.;
   fEKinBins[26] = 2000.;
   fEKinBins[27] = 3000.;
   
   // find smallest bin
   double binsmall = 1.e+10;
   for ( int i=0; i<fNEKinBins; ++i )
   {
      double deltax = fEKinBins[i+1] - fEKinBins[i];
      if ( deltax < binsmall ) binsmall = deltax;
   }
   
   int nbinstmp = 3000. / binsmall;
   
   
   fTheta.push_back( 15.);
   fTheta.push_back( 30.);
   fTheta.push_back( 60.);
   fTheta.push_back( 90.);
   fTheta.push_back(120.);
   fTheta.push_back(150.);
   fDeltaTheta = 5.;
   
   for ( size_t i=0; i<fTheta.size(); ++i )
   {

      double cth1 = std::cos( std::min( (fTheta[i]+fDeltaTheta)*CLHEP::deg, 180.*CLHEP::deg ) );
      double cth2 = std::cos( std::max( (fTheta[i]-fDeltaTheta)*CLHEP::deg,   0.*CLHEP::deg ) );
      fCosThMin.push_back( std::min(cth1,cth2) );
      fCosThMax.push_back( std::max(cth1,cth2) );
      fDeltaCosTh.push_back( std::abs(cth1-cth2) );

      int angle = fTheta[i];
      std::string hname  = "EKinNeutron" + std::to_string(angle);
      // std::string htit   = "theta=" + std::to_string( fTheta[i] ) + "deg";
      std::string htit   = "theta = " + std::to_string( angle ) + " deg";
      fHistoNeutron.push_back( new TH1D( hname.c_str(), htit.c_str(), nbinstmp, 0., 3000.) );
      // fHistoNeutron.push_back( new TH1D( hname.c_str(), htit.c_str(), fNEKinBins, fEKinBins ) );    
   
   }
      
}

TestIAEAHisto::~TestIAEAHisto()
{
   
   if ( fHistoNSec ) delete fHistoNSec;
    
   delete [] fEKinBins;
   
   for ( size_t i=0; i<fHistoNeutron.size(); ++i )
   {
      delete fHistoNeutron[i];
   }

}

void TestIAEAHisto::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector& ) 
{

   G4bool isFirst = (dynamic_cast<Tst23ParticleChange*>(aChange))->IsFisrtInteraction();
   
   int NSec = aChange->GetNumberOfSecondaries();
   
   if ( isFirst && NSec > 0 ) fHistoNSec->Fill( (double)NSec );

   const G4DynamicParticle* sec = 0;
   
   for ( int is=0; is<NSec; ++is )
   {
   
      sec = aChange->GetSecondary(is)->GetDynamicParticle();
      const G4String& pname = sec->GetDefinition()->GetParticleName();
      
      double ke = (sec->GetKineticEnergy())/CLHEP::MeV;
      // if ( ke < 0. ) ke = 0.;     
      
      G4ThreeVector mom = sec->GetMomentum()/CLHEP::MeV;      
      double        theta = mom.theta();
      double        cth = std::cos(theta);
      
      for ( size_t ith=0; ith<fTheta.size(); ++ith )
      {
         if ( cth > fCosThMin[ith] && cth <= fCosThMax[ith] )
	 {
	    if ( pname == "neutron" )
	    {
	       if ( isFirst ) fHistoNeutron[ith]->Fill(ke);
	    }
	 }
      }
   
   }
      
   return;

}

void TestIAEAHisto::Write( G4int, G4double xsec )
{
   
   int nevt = 1.;
   
   // NOTE: bear in mind that fHistoNSec->Integral() is is the number of entries WITHIN histo xmin-xmax
   //       while fHistoNSec->GetEntries() is the total number of entries, incl. under/overflow (i.e. the actual NEvents)
   //       however, AFTER the NORMALIZATION of fHistoSec, GetEntries() and Integral() give the same result
   //
   if ( fHistoNSec->GetEntries() > 0 ) nevt = fHistoNSec->GetEntries();
      
   double xbin = (G4double)fHistoNSec->GetBinWidth(1);
   fHistoNSec->Scale( (1./(xbin*nevt)) ) ;
   fHistoNSec->Write();

   double* err2  = new double[fNEKinBins];
   double* bcont = new double[fNEKinBins];
   for ( int ib=0; ib<fNEKinBins; ++ib )
   {
      err2[ib] = 0.;
      bcont[ib] = 0.;
   }
   
   for ( size_t i=0; i<fHistoNeutron.size(); ++i )
   {
      for ( int ib=0; ib<fNEKinBins; ++ib )
      {
         err2[ib] = 0.;
         bcont[ib] = 0.;
      }
      double scale = xsec / ( (double)nevt * 2.*CLHEP::pi*fDeltaCosTh[i] );
      int angle = fTheta[i];
      std::string hname  = "KENeutron" + std::to_string(angle);
      TH1D h( hname.c_str(), fHistoNeutron[i]->GetTitle(), fNEKinBins, fEKinBins );
      for ( int ib=1; ib<=fHistoNeutron[i]->GetNbinsX(); ++ib )
      {
         double xx  = fHistoNeutron[i]->GetBinCenter(ib);
	 double yy  = fHistoNeutron[i]->GetBinContent(ib);
	 double eyy = fHistoNeutron[i]->GetBinError(ib);
	 int    ibf = h.FindBin(xx);
	 if ( ibf > 0 && ibf <= fNEKinBins ) 
	 {
	    // h.AddBinContent( ibf, yy );
	    err2[ibf-1] += eyy*eyy;
	    bcont[ibf-1] += yy;
	 }
      }
      for ( int ib=1; ib<=fNEKinBins; ++ib )
      {
         h.SetBinContent( ib, bcont[ib-1] );
	 double err = std::sqrt( err2[ib-1] );
	 h.SetBinError( ib, err ); 
      }
      h.Scale( scale, "width" );
      h.Write();
      fHistoNeutron[i]->Scale( scale, "width" );
      fHistoNeutron[i]->Write();
   }
   
   return;

}
