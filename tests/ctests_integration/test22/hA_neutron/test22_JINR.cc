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
//
//               FTF test: p+C --> Pi+/Pi- ; NA49
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
std::ofstream hAn("jAn.dat",std::ios::out);

G4double Anles [7] = {30., 53., 60., 90., 100., 120., 150.};
G4double AnLow [7] = {25., 48., 55., 85.,  95., 115., 145.};
G4double AnHig [7] = {35., 58., 65., 95., 105., 125., 155.};

G4double dOmega [7];
for(G4int i=0; i < 7; i++)
{ Anles [i]*=pi/180.; AnLow [i]*=pi/180.; AnHig [i]*=pi/180.;
  dOmega[i]=2.*pi*(std::cos(AnLow[i])-std::cos(AnHig[i]));}
                                                                               // 10                                                                20                                                     30                                                          40
G4double Tkin [55] ={ 0.1 , 0.2  , 0.3  , 0.4  , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
                     1.250, 1.500, 1.750, 2.0  , 2.50, 3.00, 3.50, 4.00, 4.50, 5.00,
                       6.0,   7.0,   8.0,   9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 22.0,
                      25.0,  28.0,  31.0,  36.0, 41.0, 46.0, 51.0, 56.0, 61.0, 71.0,
                      81.0,  91.0, 101.0, 125.0,150.0,175.0,200.0,240.0,280.0,320.0,
                     360.0, 400.0, 440.0, 480.0,520.0};
G4double TkinL[55] ={ 0.05,  0.15,  0.25,  0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95,
                     1.125, 1.375, 1.625, 1.875, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75,
                       5.5,   6.5,   7.5,   8.5, 10.0, 12.0, 14.0, 16.0, 18.0, 21.0,
                      23.0,  26.0,  29.0,  33.0, 38.0, 43.0, 48.0, 53.0, 58.0, 66.0,
                      76.0,  86.0,  96.0, 110.0,135.0,160.0,180.0,220.0,260.0,300.0,
                     340.0, 380.0, 420.0, 460.0,500.0};
G4double TkinH[55] ={ 0.15,  0.25,  0.35,  0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05,
                     1.375, 1.625, 1.875, 2.125, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25,
                       6.5,   7.5,   8.5,   9.5, 12.0, 14.0, 16.0, 18.0, 20.0, 23.0,
                      27.0,  30.0,  33.0,  39.0, 44.0, 49.0, 54.0, 59.0, 64.0, 76.0,
                      86.0,  96.0, 106.0, 140.0,165.0,190.0,220.0,260.0,300.0,340.0,
                     380.0, 420.0, 460.0, 500.0,540.0};

//G4double SqrtS;

G4double XUzhiP[55][7], XUzhiN[55][7];
for(G4int i=0; i<55; i++){for(G4int j=0;j<7;j++) {XUzhiP[i][j]=0.; XUzhiN[i][j]=0.;}}

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
   //SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E);

//-------------------------------------------------------

    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
//    G4double e, px, py, pz, pt, pt2, theta;
    G4double e, theta;
    G4VParticleChange* aChange = 0;

//  G4double E=energy+part->GetPDGMass();                                  // Elab Proj
//  G4double SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E); // per  Proj+N
//  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                         //      Proj+N

    // -------- Event loop
    G4cout<<"Events start "<<nevt<<G4endl;
    G4int Ninelast=nevt;                              // Uzhi
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
      G4LorentzVector labNN(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)),e0 + mass + 938.); //amass);
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

      for(G4int i=0; i<n; ++i)              // Loop over produced particles
      {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        G4String pname=pd->GetParticleName();

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

        mom = fm.vect();

/*
        px = mom.x();
        py = mom.y();
        pz = mom.z();
        G4double Pmod=mom.mag();
        pt = std::sqrt(px*px +py*py); pt2=sqr(pt/GeV); pt2=pt2;
*/
        e  = fm.e();

        G4double Mhadron=pd->GetPDGMass();
        G4double T  = fm.e()- Mhadron;

        theta = mom.theta();
//      G4double CosTheta=std::cos(theta);

//        theta=theta*180./pi;
//G4cout<<"theta T "<<theta*180./pi<<" "<<T<<G4endl;
//      G4double costcm = std::cos(fm.theta());

//+++++++++++++++++ For each particle in the event ++++++++++++++++++++++

	 if  (pname == "proton")
         {
//
          for(G4int iT=0; iT <55; iT++)
          {
           if((TkinL[iT] <= T) && (T < TkinH[iT]))
           {
            for(G4int iTh=0; iTh < 7; iTh++)
            {
             if((AnLow[iTh] <= theta) && (theta < AnHig[iTh])) {XUzhiP[iT][iTh]+=1.;}
            }
           }
          }
//
         };

	 if ( pname == "neutron" )
         {
//
          for(G4int iT=0; iT <55; iT++)
          {
           if((TkinL[iT] <= T) && (T < TkinH[iT]))
           {
            for(G4int iTh=0; iTh < 7; iTh++)
            {
             if((AnLow[iTh] <= theta) && (theta < AnHig[iTh])) {XUzhiN[iT][iTh]+=1.;}
            }
           }
          }
//
         };

	 if ( pname == "pi+" )
         {
         };

	 if ( pname == "pi-" )
         {
         };

	 if ( pname == "kaon+" )
         {
         };

	 if ( pname == "kaon-" )
         {
         };

	 if ( pname == "anti_proton" )
         {
         };

	 if ( pname == "lambda" )
         {
         };

	 if ( pname == "sigma+" )
         {
         };

	 if ( pname == "sigma-" )
         {
         };
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	de += e;

	//	delete sec;
        delete aChange->GetSecondary(i);

      } //     end of the loop on particles
//}; //if(Uzhi_EQnex == 1)                           // Uzhi March 2015
//+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++

      if(n == 0) {Ninelast--; G4cout<<"n=0 !!! "<<G4endl;}// G4int Uzhi_i; G4cin >> Uzhi_i;}
      if(n == 1) {Ninelast--; G4cout<<"n=1 !!! "<<G4endl;}// G4int Uzhi_i; G4cin >> Uzhi_i;}
      if(n == 2) {Ninelast--; G4cout<<"n=2 !!! "<<G4endl;}// G4int Uzhi_i; G4cin >> Uzhi_i;}
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


// ----------------------------- Tkin distributions----------------------
    G4cout<< "******** Tkin distr ******* at Plab "<<Plab<<" Xin " << sigIn<< G4endl;

//    sigIn/=1000.;

//    hAp<<G4Version<<G4endl;
//    hAp<<" Tkin A30 A53 A60 A90 A100 A120 A150"<<G4endl;

    hAn<<G4Version<<G4endl;
    hAn<<" Tkin A30 A53 A60 A90 A100 A120 A150"<<G4endl;


    for(G4int ii=0; ii <55; ii++)
    {for(G4int jj=0;jj < 7; jj++)
     {
      XUzhiP[ii][jj]*= 1./Ninelast*sigIn/(TkinH[ii]-TkinL[ii])/dOmega[jj];
      XUzhiN[ii][jj]*= 1./Ninelast*sigIn/(TkinH[ii]-TkinL[ii])/dOmega[jj];
     }
    }

    for(G4int ii=0; ii <55; ii++)
    {
      hAn<<Tkin[ii];

      for(G4int jj= 0; jj <  7; jj++) { //hAp<<" "<<XUzhiP[ii][jj];
                                        hAn<<" "<<XUzhiN[ii][jj];}
      hAn<<G4endl;
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
