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

#include "G4PhysicalConstants.hh"


#include "Tst23HARPHisto.hh"
#include "Tst23ParticleChange.hh"

#include "G4TrackVector.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <sstream>

Tst23HARPHisto::Tst23HARPHisto( G4String ht )
   : TstHistoSet( ht ), 
     fNThetaBinsFW(4), fThetaMinFW(0.05), fDeltaThetaFW(0.05),
     fNThetaBinsLA(9), fThetaMinLA(0.35), fDeltaThetaLA(0.2),
     fNTheta1(14), fThetaMin1(0.15), fDeltaTheta1(0.2)
{

   G4String htitle;
   G4String hname;
  
   fHistoNSec = new TH1D( "NSec", ht.c_str(), 100, 0., 100. );

//   std::ostringstream osTitle1(std::ios_base::out|std::ios_base::app);
//   std::ostringstream osTitle2(std::ios_base::out|std::ios_base::app);
//   std::ostringstream osTitle3(std::ios_base::out|std::ios_base::app);
       
   G4double thetaMin = 0.;
   G4double thetaMax = 0.;
   std::string theta_bin_fw;
   std::string theta_bin_la;

   double parbins_fw[] = { 0.5, 1.0, 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5, 8. };
   int    nparbins_fw = sizeof(parbins_fw) / sizeof(double) - 1;

   for ( G4int i=0; i<fNThetaBinsFW; i++ )
   {
      thetaMin = fThetaMinFW + fDeltaThetaFW*i;
      thetaMax = thetaMin + fDeltaThetaFW;
      
      std::ostringstream osTitle1;
      std::ostringstream osTitle2;
      std::ostringstream osTitle3;
      
      osTitle1.clear();
      osTitle1 << thetaMin;
      theta_bin_fw = osTitle1.str() + " < theta < ";
      osTitle2.clear();
      osTitle2 << thetaMax;
      theta_bin_fw += osTitle2.str();
      theta_bin_fw += "(rad)";
   
      osTitle3.clear();
      osTitle3 << i;
      
      htitle = ht + " -> X + pi-, " + theta_bin_fw;  
      hname = "piminus_FW_" + osTitle3.str();         
//      fHistoSecPiMinusFW.push_back( new TH1D( hname.c_str(), htitle.c_str(), 80, 0., 8.0 ) );
      fHistoSecPiMinusFW.push_back( new TH1D( hname.c_str(), htitle.c_str(), nparbins_fw, parbins_fw ) );

      htitle = ht + " -> X + pi+, " + theta_bin_fw;
      hname = "piplus_FW_" + osTitle3.str();
//      fHistoSecPiPlusFW.push_back( new TH1D( hname.c_str(), htitle.c_str(), 80, 0., 8.0 ) );
      fHistoSecPiPlusFW.push_back( new TH1D( hname.c_str(), htitle.c_str(), nparbins_fw, parbins_fw ) );

      htitle = ht + " -> X + p, " + theta_bin_fw; 
      hname = "proton_FW_" + osTitle3.str();
      fHistoSecProtonFW.push_back( new TH1D( hname.c_str(), htitle.c_str(),nparbins_fw, parbins_fw ) ); 

   }
   
   double parbins_la[] = { 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 };
   int    nparbins_la = sizeof(parbins_la) / sizeof(double) - 1;

   for ( G4int i=0; i<fNThetaBinsLA; i++ )
   {
     
      thetaMin = fThetaMinLA + fDeltaThetaLA*i;
      thetaMax = thetaMin + fDeltaThetaLA; 
     
      std::ostringstream osTitle1;
      std::ostringstream osTitle2;
      std::ostringstream osTitle3;

      osTitle1.clear();
      osTitle1 << thetaMin;
      theta_bin_la = osTitle1.str() + " < theta < ";
      osTitle2.clear();
      osTitle2 << thetaMax;
      theta_bin_la += osTitle2.str();
      theta_bin_la += "(rad)";
      
      osTitle3.clear();
      osTitle3 << i;
      
      htitle = ht + " -> X + pi-, " + theta_bin_la;  
      hname = "piminus_LA_" + osTitle3.str();         
//      fHistoSecPiMinusLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5 ) );
      fHistoSecPiMinusLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), nparbins_la, parbins_la ) );
//      hname = "pimimus_mom_total_" + osTitle3.str();
//      fHistoSecPiMinusTot.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5) );	

      htitle = ht + " -> X + pi+, " + theta_bin_la;
      hname = "piplus_LA_" + osTitle3.str();
//      fHistoSecPiPlusLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5 ) );
      fHistoSecPiPlusLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), nparbins_la, parbins_la ) );
//      hname = "piplus_mom_total_" + osTitle3.str();
//      fHistoSecPiPlusTot.push_back( new TH1D( hname.c_str(), htitle.c_str(), 50, 0., 2.5 ) );	

      htitle = ht + " -> X + p, " + theta_bin_la;
      hname = "proton_LA_" + osTitle3.str();
      fHistoSecProtonLA.push_back( new TH1D( hname.c_str(), htitle.c_str(), nparbins_la, parbins_la ) );

   }
   
//---------------
//
// Extras for benchmarking Shielding/ShieldingM vs the "Mu2e way"
// (Mu2e looks at backward pi- production in p+Ta at 8GeV/c, 
//  in a specific pi- momentum bin of 100-150GeV/c
//

   double* thetabins = new double[fNTheta1+1];
   for ( int ith=0; ith<=fNTheta1; ++ith )
   {
      thetabins[ith] = fThetaMin1 + fDeltaTheta1*ith;
   }
   
   std::string momtitle = "";
   fNMomBinsLA = nparbins_la+1;
   fMomBinsLA  = new G4double[fNMomBinsLA+1];
   
   fMomBinsLA[0] = 0.05;
   for ( int ip=0; ip<nparbins_la; ++ip )
   {
      fMomBinsLA[ip+1] = parbins_la[ip];
   }
   fMomBinsLA[fNMomBinsLA] = parbins_la[nparbins_la];
   
   for (int ip=0; ip<fNMomBinsLA; ++ip )
   {
      
      std::ostringstream osTitle1;
      std::ostringstream osTitle2;
      std::ostringstream osTitle3;
      
      osTitle1.clear();
      osTitle1 << fMomBinsLA[ip];
      momtitle = osTitle1.str() + " < P(pi-), GeV/c < ";
      osTitle2.clear();
      osTitle2 << fMomBinsLA[ip+1];
      momtitle += osTitle2.str();
      
      osTitle3.clear();
      osTitle3 << ip;
      
      //htitle = ht + " -> X + pi-, " + momtitle;  
      htitle = momtitle;
      hname = "piminus_mom_" + osTitle3.str();
      fHistoSecPiMinusMomBins.push_back( new TH1D( hname.c_str(), htitle.c_str(), fNTheta1, thetabins ) );
      
      momtitle = osTitle1.str() + " < P(pi+), GeV/c < ";
      momtitle += osTitle2.str();
      
      //htitle = ht + " -> X + pi+, " + momtitle;
      htitle = momtitle;
      hname = "piplus_mom_" + osTitle3.str();
      fHistoSecPiPlusMomBins.push_back( new TH1D( hname.c_str(), htitle.c_str(), fNTheta1, thetabins ) );
      
   }
   
   delete [] thetabins;

}

Tst23HARPHisto::~Tst23HARPHisto()
{

   if ( fHistoNSec ) delete fHistoNSec;
   
   for (size_t i=0; i<fHistoSecPiMinusFW.size(); i++ )
   {
      delete fHistoSecPiMinusFW[i];
   }
   for (size_t i=0; i<fHistoSecPiPlusFW.size(); i++ )
   {
      delete fHistoSecPiPlusFW[i];
   }
   for (size_t i=0; i<fHistoSecProtonFW.size(); i++ )
   {
      delete fHistoSecProtonFW[i];
   }

   for (size_t i=0; i<fHistoSecPiMinusLA.size(); i++ )
   {
      delete fHistoSecPiMinusLA[i];
   }
   for (size_t i=0; i<fHistoSecPiPlusLA.size(); i++ )
   {
      delete fHistoSecPiPlusLA[i];
   }
   for (size_t i=0; i<fHistoSecProtonLA.size(); i++ )
   {
      delete fHistoSecProtonLA[i];
   }
   
   for ( size_t i=0; i<fHistoSecPiMinusMomBins.size(); ++i )
   {
      delete fHistoSecPiMinusMomBins[i];
   }
   for ( size_t i=0; i<fHistoSecPiPlusMomBins.size(); ++i )
   {
      delete fHistoSecPiPlusMomBins[i];
   }
   
   delete [] fMomBinsLA;

}

void Tst23HARPHisto::FillEvt( G4VParticleChange* aChange, const G4LorentzVector&, const G4LorentzVector&  ) 
{

   
   G4bool isFirst = (dynamic_cast<Tst23ParticleChange*>(aChange))->IsFisrtInteraction();
      
   G4int NSec = aChange->GetNumberOfSecondaries();

   if ( isFirst && NSec > 0 ) fHistoNSec->Fill( (double)NSec );
     
   const G4DynamicParticle* sec = 0;
          
   for (G4int i=0; i<NSec; i++) 
   {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();			
	const G4String& pname = sec->GetDefinition()->GetParticleName();
		
	// G4ThreeVector mom = sec->GetMomentum() / GeV ;
	G4double pmom = sec->GetTotalMomentum() / GeV ;
	//double mass = sec->GetDefinition()->GetPDGMass() / GeV;
	//double ekin = sec->GetKineticEnergy() / GeV ;
	G4double theta = (sec->GetMomentum()).theta();
	
//-------------------
// this part is for Shielding/ShieldingM benchmarking the "Mu2e way"...

	int ipbin = -1; 
	for ( G4int ip=0; ip<fNMomBinsLA; ++ip )
	{
	   if ( pmom > fMomBinsLA[ip] && pmom <= fMomBinsLA[ip+1] )
	   {
	      ipbin = ip;
	      break;
	   }
	}
	if ( ipbin >= 0 && ipbin < fNMomBinsLA ) 
	{
	   if ( pname == "pi-" )
	   {
	      if ( isFirst ) fHistoSecPiMinusMomBins[ipbin]->Fill( theta );
	   }
	   else if ( pname == "pi+" )
	   {
	      if ( isFirst ) fHistoSecPiPlusMomBins[ipbin]->Fill( theta );
	   }
	}
//---------------------
// ...and from here on, it's the standard HARP binning, etc.

	if ( theta < fThetaMinFW ) continue;
	if ( theta < fThetaMinFW+fDeltaThetaFW*fNThetaBinsFW )
	{
	   G4int ith = ( theta - fThetaMinFW ) / fDeltaThetaFW;
	   if ( pname == "pi-" )
	   {
	      if ( isFirst ) fHistoSecPiMinusFW[ith]->Fill( pmom );
	   }
	   else if ( pname == "pi+" )
	   {
	      if ( isFirst ) fHistoSecPiPlusFW[ith]->Fill( pmom );
	   }
	   else if ( pname == "proton" )
	   {
	      if ( isFirst ) fHistoSecProtonFW[ith]->Fill( pmom );
	   }
	}
		
	if ( theta < fThetaMinLA ) continue;
	if ( theta > fThetaMinLA+fDeltaThetaLA*fNThetaBinsLA ) continue;
	G4int    itheta = ( theta - fThetaMinLA ) / fDeltaThetaLA;
	if ( itheta < 0 || itheta >= fNThetaBinsLA ) continue;
	        
	if ( pname == "pi-" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiMinusLA[itheta]->Fill( pmom );
	   }
	}
	else if ( pname == "pi+" )
	{
	   if ( isFirst )
	   {
	      fHistoSecPiPlusLA[itheta]->Fill( pmom );
	   }
	}
	else if ( pname == "proton" )
	{
	   if ( isFirst )
	   {
	      fHistoSecProtonLA[itheta]->Fill( pmom );
	   }
	}
	
   }
      
   return;
   
}

void Tst23HARPHisto::Write( G4int, G4double xsec )
{

   // fHistoFile->cd();

   G4double xbin = 1.;
   G4double norm = 1.;
   G4double scale = 1.;
   
   // NOTE: bear in mind that fHistoNSec->Integral() is is the number of entries WITHIN histo xmin-xmax
   //       while fHistoNSec->GetEntries() is the total number of entries, incl. under/overflow (i.e. the actual NEvents)
   //       however, AFTER the NORMALIZATION of fHistoSec, GetEntries() and Integral() give the same result
   //
   norm = (double)(fHistoNSec->GetEntries());
   
   std::cout << " xsec = " << xsec << std::endl;
   std::cout << " norm = " << norm << std::endl;
   std::cout << " integral = " << fHistoNSec->Integral() << std::endl;
   
   xbin = (G4double)fHistoNSec->GetBinWidth(1);
   scale = 1. / (xbin*norm);
   fHistoNSec->Scale( scale ) ;
   fHistoNSec->Write();
   
   // secondary pi-
   //
   for ( size_t i=0; i<fHistoSecPiMinusFW.size(); ++i )
   {
//      xbin = fHistoSecPiMinusFW[i]->GetBinWidth(1);
//      scale = xsec / (xbin*norm*fDeltaThetaFW);
      scale = xsec / (norm*fDeltaThetaFW);
      fHistoSecPiMinusFW[i]->Scale(scale,"width");
      fHistoSecPiMinusFW[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiMinusLA.size(); ++i )
   {
//      xbin = fHistoSecPiMinusLA[i]->GetBinWidth(1);
//      scale = xsec / (xbin*norm*fDeltaThetaLA);
      scale = xsec / (norm*fDeltaThetaLA);
      fHistoSecPiMinusLA[i]->Scale(scale,"width");
      fHistoSecPiMinusLA[i]->Write();
   }

   // secondary pi+
   //
   for ( size_t i=0; i<fHistoSecPiPlusFW.size(); ++i )
   {
//      xbin = fHistoSecPiPlusFW[i]->GetBinWidth(1);
//      scale = xsec / (xbin*norm*fDeltaThetaFW);
      scale = xsec / (norm*fDeltaThetaFW);
      fHistoSecPiPlusFW[i]->Scale(scale,"width");
      fHistoSecPiPlusFW[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiPlusLA.size(); ++i )
   {
//      xbin = fHistoSecPiPlusLA[i]->GetBinWidth(1);
//      scale = xsec / (xbin*norm*fDeltaThetaLA);
      scale = xsec / (norm*fDeltaThetaLA);
      fHistoSecPiPlusLA[i]->Scale(scale,"width");
      fHistoSecPiPlusLA[i]->Write();
   }

   // secondary proton
   //
   for ( size_t i=0; i<fHistoSecProtonFW.size(); ++i )
   {
      scale = xsec / (norm*fDeltaThetaFW);
      fHistoSecProtonFW[i]->Scale( scale, "width" );
      fHistoSecProtonFW[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecProtonLA.size(); ++i )
   {
      scale = xsec / (norm*fDeltaThetaLA);
      fHistoSecProtonLA[i]->Scale( scale, "width" );
      fHistoSecProtonLA[i]->Write();
   }
   
// ----------------
//  for Shielding/ShieldingM benchmarking the "Mu2e way"

   for ( size_t i=0; i<fHistoSecPiMinusMomBins.size(); ++i )
   {
      scale = xsec / (norm*(fMomBinsLA[i+1]-fMomBinsLA[i]));
      fHistoSecPiMinusMomBins[i]->Scale( scale, "width" );
      fHistoSecPiMinusMomBins[i]->Write();
   }
   for ( size_t i=0; i<fHistoSecPiPlusMomBins.size(); ++i )
   {
      scale = xsec / (norm*(fMomBinsLA[i+1]-fMomBinsLA[i]));
      fHistoSecPiPlusMomBins[i]->Scale( scale, "width" );
      fHistoSecPiPlusMomBins[i]->Write();
   }
   

   return;

}
