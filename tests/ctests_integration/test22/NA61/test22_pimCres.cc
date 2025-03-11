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
//               FTF test: Pi- + C12 interactions -> Rho, Omega, K*0(852)
// -------------------------------------------------------------------
#include "globals.hh"
#include "G4Version.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "FTFtest1.icc"
#include "G4ChipsComponentXS.hh"

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
//G4double sigTot = 0;
//G4double sigEl  = 0;
G4double sigIn  = 0;

//-------------------------- Global histograms  -------------------------
std::ofstream hAx("PimCna61.dat",std::ios::out);

G4double SqrtS;
//G4double Ybeam;

G4double XUzhi[40][3];
for(G4int ii=0; ii<40; ii++){for(G4int j=0;j<3;j++) XUzhi[ii][j]=0.;}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;

  for(G4int run=0; run<100; run++) {
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

    G4bool CHIPapplic=true;                          //false;   Uzhi 29.01.13
    if(CHIPapplic)
    {
     chipsTot=CHIPSxsManager->GetTotalElementCrossSection(part,energy,Z,A-Z);
     chipsEl =CHIPSxsManager->GetElasticElementCrossSection(part,energy,Z,A-Z);
     chipsIn =CHIPSxsManager->GetInelasticElementCrossSection(part,energy,Z,A-Z);
     chipsTot/=millibarn; chipsEl/=millibarn; chipsIn/=millibarn;

     G4cout<<"CHIPS cross sections are used:----------------------"<<G4endl<<
             "Plab          Total        Elastic      Inelastic"   <<G4endl;
     G4cout<<" "<<Plab/GeV<<" "<< chipsTot<<" "<<chipsEl<<" "<<chipsIn <<G4endl<<G4endl;

     //sigTot=chipsTot;
     //sigEl=chipsEl;
     sigIn=chipsIn;
    } else
    {
     //sigTot = cross_sec;
     //sigEl  = cross_secel;
     sigIn  = cross_inel;

     G4cout<<"Proposed Xs (mb) are used: Tot El In: "
           <<cross_sec<<" "<<cross_secel<<" "<<cross_inel<<G4endl;
    }

//+++++++++++++++++++++++++++++++++ For each energy +++++++++++++++++++++
   G4double E=energy+part->GetPDGMass();
   SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E);
   //Ybeam=0.5*std::log((E+Plab)/(E-Plab));

//  G4int Ntotal=nevt;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//-------------------------------------------------------

    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, pz;  //px, py, pt, pt2, theta;
    G4VParticleChange* aChange = 0;

    // -------- Event loop
    G4cout<<"Events start "<<nevt<<G4endl;
    G4int Ninelast=nevt;                              // Uzhi

//=================================================================
    for (G4int iter=0; iter<nevt; ++iter) {
//=================================================================
      if(verbose > 0) G4cout<<"Start events loop***********************"<<G4endl;

      if(verbose>=1 || iter == modu*(iter/modu)) {
        G4cout << "### " << iter << "-th event start " <<Plab/GeV<<" Nevents "<<nevt<<G4endl;
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
      G4LorentzVector labNN(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)),e0 + mass + 939.); // amass
      G4ThreeVector boostNN = labNN.boostVector();
//------------

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de);
      labv -= dee;

      G4int n = aChange->GetNumberOfSecondaries();     // Multiplicity of prod. part.

      if(verbose>=1) G4cout<<" Uzhi ------------ N prod. part "<<n<<G4endl;
      if((verbose > 0) && (n < 2)) {G4cout<<"Multiplicity of produced < 2!!!"<<G4endl;}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for(G4int i=0; i<n; ++i)              // Loop over produced particles
      {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        G4String pname=pd->GetParticleName();
        G4int partPDG = pd->GetPDGEncoding();

        if(verbose>=1) G4cout<<" Part  "<<i<<" "<<pname
                             <<" "<<sec->Get4Momentum()/GeV
                             <<sec->Get4Momentum().mag()/GeV<<G4endl;

        fm  = sec->Get4Momentum();
        labv -= fm;   // For checking energy-momentum conservation
	// electron can come only from internal conversion
	// its mass should be added to initial state
        if(pd == electron) {

	  labv += G4LorentzVector(0.0,0.0,0.0,electron_mass_c2);
	}


        fm.boost(-boostNN);

        mom = fm.vect();

//        px = mom.x();
//        py = mom.y();
        pz = mom.z();
//        G4double Pmod=mom.mag(); Pmod=Pmod;
//        pt = std::sqrt(px*px +py*py); pt2=sqr(pt/GeV); pt2=pt2;
        e  = fm.e();
//        theta = mom.theta();
//        theta=theta*180./pi;

        G4double xF=2.*pz/SqrtS;

        G4int NxUzhi=int(xF/0.05);
        if(NxUzhi > 19) NxUzhi=19;

//+++++++++++++++++ For each particle in the event ++++++++++++++++++++++

	 if  (partPDG == 113)     // Rho-0
         {
            if(NxUzhi >= 0) XUzhi[NxUzhi][0]+=xF;
         };

	 if  (partPDG == 223)     // omega
         {
            if(NxUzhi >= 0) XUzhi[NxUzhi][1]+=xF;
         };

	 if  (std::abs(partPDG) == 313)     // K*0(892)
         {
//G4cout<<"K*0 "<<partPDG<<G4endl;  G4int Uzhi; G4cin>>Uzhi;
            if(NxUzhi >= 0) XUzhi[NxUzhi][2]+=xF;
         };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	de += e;

	//	delete sec;
        delete aChange->GetSecondary(i);

      } //     end of the loop on particles

//+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++

      if(n == 0) {Ninelast--; G4cout<<"n=0 !!! "<<G4endl;}
      if(n == 1) {Ninelast--; G4cout<<"n=1 !!! "<<G4endl;}
      if(n == 2) {Ninelast--;} // G4cout<<"n=2 !!! "<<G4endl;}
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
G4cout<< "***********************************************************"<< G4endl;

G4cout<<"nevt Ninel "<<nevt<<" "<<Ninelast<<G4endl;
G4cout<<"Plab "<<Plab/GeV<<" SigIn "<<sigIn<<G4endl;

//-------------------------------------------------------------------


    if(verbose > 0) {
      G4cout << "###### End of run # " << run << "     ######" << G4endl;
    }
//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi


// ----------------------------- xF distributions----------------------
    G4cout<< "******** xF distr ******* at Plab "<<Plab<<" Xin " << sigIn<< G4endl;
    hAx<<G4Version<<G4endl;
    hAx<<" xF  Rho Omega K0892"<< G4endl;

    for(G4int ii=0; ii <20; ii++)
    {for(G4int jj=0;jj < 3; jj++){XUzhi[ii][jj]*= 1./Ninelast/0.05;}}  //*sigIn/12.;}}

    G4double xF=-0.025;
    for(G4int ii=0; ii <20; ii++)
    {
     xF+=0.05;
     hAx<<xF<<" ";
     for(G4int jj=0;jj < 3; jj++){hAx << XUzhi[ii][jj]<<" ";}
     hAx<<G4endl;
    }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi
//    G4cerr << "###### End of run # " << run << "     ######" << G4endl;

  } while(end);
  }  // End of job ----------------------------------------------
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

//=============================
#include "G4DecayKineticTracks.hh"
#include "G4KineticTrackVector.hh"
#include "G4KineticTrack.hh"


// Decay all input tracks, put daughters onto end of list

G4DecayKineticTracks::G4DecayKineticTracks(G4KineticTrackVector *tracks) {
  if (tracks) Decay(tracks);
}

void G4DecayKineticTracks::Decay(G4KineticTrackVector *tracks) const {

  if (!tracks) return;

  G4KineticTrackVector* daughters = 0;
  for (size_t i=0; i<tracks->size(); ++i) {
    G4KineticTrack* track = (*tracks)[i];
    if (!track) continue;

    // Select decay of current track, put daughters at end of vector
//
if( (track->GetDefinition()->GetPDGEncoding() == 113 ) ||
    (track->GetDefinition()->GetPDGEncoding() == 223 ) ||
    (track->GetDefinition()->GetPDGEncoding() ==-313 ) ||
    (track->GetDefinition()->GetPDGEncoding() == 313 )   ) {daughters =0;} else
{
//
    daughters = track->GetDefinition()->IsShortLived() ? track->Decay() : 0;
}
if(track->GetDefinition()->GetParticleName() == "eta"      ) daughters = track->Decay();
if(track->GetDefinition()->GetParticleName() == "eta_prime") daughters = track->Decay();
/*
if(track->GetDefinition()->GetPDGEncoding() == 111 ) daughters = track->Decay(); // Pi0

if(track->GetDefinition()->GetPDGEncoding() == 3122) daughters = track->Decay(); // Lambda
if(track->GetDefinition()->GetPDGEncoding() == 3222) daughters = track->Decay(); // Sigma+
if(track->GetDefinition()->GetPDGEncoding() == 3212) daughters = track->Decay(); // Sigma0
if(track->GetDefinition()->GetPDGEncoding() == 3112) daughters = track->Decay(); // Sigma-
*/
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
