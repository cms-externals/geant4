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
//               FTF test: Pbar+P interaction channels
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Test30
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 March 2002
//
//      Modifications:
//      14.11.03 Renamed to cascade
//      09.05.06 Return back to test30
// -------------------------------------------------------------------
#include "globals.hh"
#include "G4Version.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "FTFtest1.icc"
#include "G4ChipsComponentXS.hh"                  // Uzhi 29.01.13

#include "UZHI_diffraction.hh"

int main(int argc, char** argv)
{
  CLHEP::RanluxEngine defaultEngine( 1234567, 4 );
  CLHEP::HepRandom::setTheEngine( &defaultEngine );

  G4cout << "========================================================" << G4endl;
  G4cout << "======              FTF Test Start              ========" << G4endl;
  G4cout << "========================================================" << G4endl;

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();
  G4String fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

//-----------------------------------------------------------------------
  #include "FTFtest2.icc"   // Initialization
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
G4double sigTot = 0;
//G4double sigEl  = 0;
G4double sigIn  = 0;

//-------------------------- Global histograms  -------------------------

  G4int Uzhi_run=0;

  G4double Xs[50][10];
  for(G4int ii=0; ii<50; ii++)
   {
    for(G4int jj=0; jj<10; jj++)
       {Xs[ii][jj]=0.;}
   };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;

  for(G4int run=0; run<100; run++) {
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//-------------------------- Current histograms -------------------------

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do {
    #include "FTFtest3.icc"  // -------- Read input file
    #include "FTFtest4.icc"  // -------- Start run processing
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    G4cout << "cross(mb)in= " << cross_sec*1000./barn << G4endl
           << "cross(mb)el= " << cross_secel*1000./barn<<G4endl<<G4endl;

    cross_inel=cross_sec-cross_secel; // +++++++++++++++++++++++++

    cross_sec/=millibarn;   // Inel Cross section in mb
    cross_secel/=millibarn; // Elas Cross section in mb
    cross_inel/=millibarn;  // Inel Cross section in mb

    G4cout<<"Element A Z N: "<<A<<" "<<Z<<" "<<A-Z<<G4endl;
    G4cout<<"Proposed Xs (mb): Tot El In: "
          <<cross_sec<<" "<<cross_secel<<" "<<cross_inel<<G4endl;

//---------------------------------------------------------------------------
// Kossov cross sections      ---------------------------
    G4double chipsTot, chipsEl, chipsIn;

    static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();
    G4ChipsComponentXS* CHIPSxsManager = _instance;

    G4bool CHIPapplic=false; //true;                          //false;   Uzhi 29.01.13
    if(CHIPapplic)
    {
     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);
     chipsTot/=millibarn; chipsEl/=millibarn; chipsIn/=millibarn;

     G4cout<<"CHIPS cross sections are used:----------------------"<<G4endl<<
             "Plab          Total        Elastic      Inelastic"   <<G4endl;
     G4cout<<" "<<Plab/GeV<<" "<< chipsTot<<" "<<chipsEl<<" "<<chipsIn <<G4endl<<G4endl;

     sigTot=chipsTot;
     //sigEl=chipsEl;
     sigIn=chipsIn;
    } else
    {
     sigTot = cross_sec;
     //sigEl  = cross_secel;
     sigIn  = cross_inel;

     G4cout<<"Proposed Xs (mb) are used: Tot El In: "
           <<cross_sec<<" "<<cross_secel<<" "<<cross_inel<<G4endl;
    }

//+++++++++++++++++++++++++++++++++ For each energy +++++++++++++++++++++
    Xs[Uzhi_run][0]=Plab/GeV;

    G4int Ntotal=nevt;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//-------------------------------------------------------

    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, theta;
    G4VParticleChange* aChange = 0;

//  G4double E=energy+part->GetPDGMass();                                  // Elab Proj
//  G4double SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E); // per  Proj+N
//  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                         //      Proj+N

    // -------- Event loop
    G4cout<<"Events start "<<nevt<<G4endl;
//  G4int Ninelast=nevt;                              // Uzhi
//=================================================================
    for (G4int iter=0; iter<nevt; ++iter) {
//=================================================================
      if(verbose > 0) G4cout<<"Start events loop***********************"<<G4endl;

      if(verbose>=1 || iter == modu*(iter/modu)) {
        G4cout << "### " << iter << "-th event start " <<Plab/GeV<<G4endl;
      }

      if(saverand) {defaultEngine.saveStatus("initial.conf");}

      G4double e0 = energy;
      do {
        if(sigmae > 0.0) e0 = G4RandGauss::shoot(energy,sigmae);
      } while (e0 < 0.0);

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      G4double amass = phys->GetNucleusMass();
      // note: check of 4-momentum balance for CHIPS is not guranteed due to
      // unknown isotope
      aChange = proc->PostStepDoIt(*gTrack,*step);

      G4double mass = part->GetPDGMass();

      if ( ionParticle )
      {
       e0/=ionA; G4double mass_N=938.*MeV;                              // Init 4-mom
       labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass_N)), //   NN
	 		      e0 + mass_N +  mass_N);
      } else
      {
       labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)),   //   hA
  		              e0 + mass + amass);
      }

      G4ThreeVector bst = labv.boostVector();          // To CMS NN in AA or hA
//------------
      G4LorentzVector labNN(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)),e0 + mass + amass);
      G4ThreeVector boostNN = labNN.boostVector();
//------------

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de);
      labv -= dee;

      G4int n = aChange->GetNumberOfSecondaries();     // Multiplicity of prod. part.

      if(verbose>=1) G4cout<<" Uzhi ------------ N prod. part "<<n<<G4endl;
//++++++++++++++++ Variables for each event +++++++++++++++++++++++++++++
      if((verbose > 0) && (n < 2)) {G4cout<<"Multiplicity of produced < 2!!!"<<G4endl;}
      if(n < 2) {Ntotal--;}

//    G4int nbar = 0;
      G4int Npim = 0;
      G4int Npip = 0;
      G4int Npi0 = 0;

      G4int Nrhop = 0;
      G4int Nrho0 = 0;
      G4int Nrhom= 0;

      G4int Nomega = 0;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      for(G4int i=0; i<n; ++i)              // Loop over produced particles
      {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        G4String pname=pd->GetParticleName();

        if(verbose>=1) G4cout<<" Part  "<<i<<" "<<pname
                             <<" "<<sec->Get4Momentum()/GeV
                             <<sec->Get4Momentum().mag()/GeV<<G4endl;

        fm  = sec->Get4Momentum();

        mom = sec->GetMomentum();
//      G4double mas = pd->GetPDGMass();
//	G4double p = mom.mag();

        labv -= fm;   // For checking energy-momentum conservation

	// electron can come only from internal conversion
	// its mass should be added to initial state
        if(pd == electron) {

	  labv += G4LorentzVector(0.0,0.0,0.0,electron_mass_c2);
	}

        e  = fm.e() - m;
        theta = mom.theta();

//        G4double CosTheta=std::cos(theta);

        theta=theta*180./pi;

        fm.boost(-bst);

//        G4double costcm = std::cos(fm.theta());
//+++++++++++++++++ For each particle in the event ++++++++++++++++++++++
	if      ( pname == "pi-" )    {Npim++;}
	else if ( pname == "pi+" )    {Npip++;}
	else if ( pname == "pi0" )    {Npi0++;            }
        else if ( pd->GetPDGEncoding() == 213 ) {Nrhop++;}
        else if ( pd->GetPDGEncoding() == 113 ) {Nrho0++;}
        else if ( pd->GetPDGEncoding() ==-213 ) {Nrhom++;}
        else if ( pd->GetPDGEncoding() == 223 ) {Nomega++;}
	else{	}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	de += e;

	//	delete sec;
        delete aChange->GetSecondary(i);

      } //     end of the loop on particles
//+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++
                 ;
      if((n == 2) && (Npip == 1)  && (Npim == 1)                 ) {Xs[Uzhi_run][1]++;} // Pi+ Pi-
      if((n == 2) && (Nrho0 == 1) && (Nomega == 1)               ) {Xs[Uzhi_run][2]++;} // Rho0 omega
      if((n == 2) && (Nrho0 == 1) && (Npi0 == 1)                 ) {Xs[Uzhi_run][3]++;} // Rho0 Pi0
      if((n == 2) && (Npim == 1)  && (Nrhop == 1)                ) {Xs[Uzhi_run][4]++;} // Pi- Rho+
      if((n == 2) && (Npip == 1)  && (Nrhom == 1)                ) {Xs[Uzhi_run][4]++;} // Pi+ Rho-
      if((n == 3) && (Npip == 1)  && (Npim == 1) && (Npi0 == 1)  ) {Xs[Uzhi_run][5]++;} // Pi+ Pi- Pi0
      if((n == 3) && (Npip == 1)  && (Npim == 1) && (Nrho0 == 1) ) {Xs[Uzhi_run][6]++;} // Pi+ Pi- Rho0
      if((n == 3) && (Npip == 1)  && (Npim == 1) && (Nomega == 1)) {Xs[Uzhi_run][7]++;} // Pi+ Pi- omega
      if((n == 3) && (Npi0 == 1)  && (Nrho0 == 2)                ) {Xs[Uzhi_run][8]++;} // Pi0 2 Rho0
      if((n == 3) && (Npim == 1)  && (Nrhop == 1) && (Nrho0 == 1)) {Xs[Uzhi_run][9]++;} // Pi- Rho+ Rho0
      if((n == 3) && (Npip == 1)  && (Nrhom == 1) && (Nrho0 == 1)) {Xs[Uzhi_run][9]++;} // Pi- Rho+ Rho0

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(verbose > 0)
        G4cout << "Energy/Momentum balance= " << labv << G4endl;

      aChange->Clear();

      if(verbose > 0)
      {
       G4cout << "End event =====================================" <<Plab<< G4endl; // Uzhi
       G4int Uzhi_i;                                                        // Uzhi
       G4cin >> Uzhi_i;                                                     // Uzhi
      }
//

    }   // End of the event loop ------------------------------------

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++
//sigTot=sigTot; sigEl=sigEl;
//sigTot=1.;
if(Ntotal != 0)
{
     for(G4int ii=1;ii<10; ii++) {Xs[Uzhi_run][ii]*=sigTot/Ntotal;}
}

//-------------------------------------------------------------------
G4cout<< "***********************************************************"<< G4endl;

G4cout<<"nevt Ninel "<<nevt<<" "<<Ntotal<<G4endl;
G4cout<<"Plab "<<Plab/GeV<<" SigIn "<<sigIn<<G4endl;
//-------------------------------------------------------------------


    if(verbose > 0) {
      G4cout << "###### End of run # " << run << "     ######" << G4endl;
    }
//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi
Uzhi_run++;                                                            // Uzhi
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi
//    G4cerr << "###### End of run # " << run << "     ######" << G4endl;

  } while(end);
  }  // End of job ----------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//----------------------------- Write distributions------------------

std::ofstream PbarPchan1("PbarPrho.dat",std::ios::out);
//PbarPchan1 ################### Channels with Rho-mesons in final sttates ###################################
PbarPchan1<<G4Version<<G4endl;
PbarPchan1<<" Plab    Pi+Pi- Rho0Om  Rho0Pi0  Pi-Rho+  Pi+Pi-Pi0  Pi+Pi-Rho0  Pi_Pi-Omega  Pi02Rho0  Pi-Rho+Rho0"<< G4endl;
//                      1       2       3        4          5         6           7           8            9
for(G4int ii=0;ii<Uzhi_run;ii++) //--------------------------------- Uzhi
  {
   PbarPchan1<<" "<<Xs[ii][0];
   for(G4int jj=1;jj<10;jj++) {PbarPchan1<<" "<<Xs[ii][jj];};

   PbarPchan1<<G4endl;
  };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  delete pFrame;
//  delete lFrame;
//  delete sFrame;
  delete mate;
  delete fin;
  delete phys;
  partTable->DeleteAllParticles();

  G4cout << "###### End of test #####" << G4endl;
}

#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"


// Decay all input tracks, put daughters onto end of list

G4DecayKineticTracks::G4DecayKineticTracks(G4KineticTrackVector *tracks) {
  if (tracks) Decay(tracks);
}

void G4DecayKineticTracks::Decay(G4KineticTrackVector *tracks) const {
  if (!tracks) return;
return;
  G4KineticTrackVector* daughters = 0;
  for (size_t i=0; i<tracks->size(); ++i) {
    G4KineticTrack* track = (*tracks)[i];
    if (!track) continue;

    // Select decay of current track, put daughters at end of vector
    daughters = track->GetDefinition()->IsShortLived() ? track->Decay() : 0;
//
if((!daughters) && ((track->GetDefinition()->GetPDGEncoding() == 221) ||    // Eta
                    (track->GetDefinition()->GetPDGEncoding() == 331)   ))  // Eta_prime
daughters = track->Decay();
//
    if (daughters) {
      tracks->insert(tracks->end(), daughters->begin(), daughters->end());
      delete track;		// Remove parent track
      delete daughters;
      (*tracks)[i] = NULL;	// Flag parent's slot for removal
    }
  }

  // Find and remove null pointers created by decays above
  for (int j=tracks->size()-1; j>=0; --j) {
    if (NULL == (*tracks)[j]) tracks->erase(tracks->begin()+j);
  }
}
