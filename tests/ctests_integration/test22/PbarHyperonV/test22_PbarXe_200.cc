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
//      edition  29.08.2014  A.Galoyan
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
#include "G4Timer.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "FTFtest1.icc"
#include "G4ChipsComponentXS.hh"                  // Uzhi 29.01.13
#include "UZHI_diffraction.hh"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>

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

//int npart;
/*
 //   Root initialization
 TFile f1("PbXe200.root","RECREATE");


TH1F *h1=new TH1F("YKs0","",50, -1.0, 6.0);
//h1->SetStats(kFALSE);
h1->GetYaxis()->SetTitle(" d#sigma/dY(mub)  ");
h1->GetXaxis()->SetTitle(" Y ");
h1->GetXaxis()->SetTitleSize(0.06);
h1->GetYaxis()->SetTitleSize(0.06);
h1->GetXaxis()->SetTitleColor(4);
h1->GetYaxis()->SetTitleColor(4);
h1->GetXaxis()->SetTitleOffset(0.7);
h1->GetYaxis()->SetTitleOffset(0.8);
h1->GetXaxis()->SetLabelSize(0.05);
h1->GetYaxis()->SetLabelSize(0.05);

TH1F *h2=new TH1F("YLambda","",50, -1.0, 6.0);
//h2->SetStats(kFALSE);
h2->GetYaxis()->SetTitle(" d#sigma/dY(mub)  ");
h2->GetXaxis()->SetTitle(" Y ");
h2->GetXaxis()->SetTitleSize(0.06);
h2->GetYaxis()->SetTitleSize(0.06);
h2->GetXaxis()->SetTitleColor(4);
h2->GetYaxis()->SetTitleColor(4);
h2->GetXaxis()->SetTitleOffset(0.7);
h2->GetYaxis()->SetTitleOffset(0.8);
h2->GetXaxis()->SetLabelSize(0.05);
h2->GetYaxis()->SetLabelSize(0.05);

TH1F *h3=new TH1F("YLambdabar","",50, -1.0, 6.0);
//h2->SetStats(kFALSE);
h3->GetYaxis()->SetTitle(" d#sigma/dY(mub)  ");
h3->GetXaxis()->SetTitle(" Y ");
h3->GetXaxis()->SetTitleSize(0.06);
h3->GetYaxis()->SetTitleSize(0.06);
h3->GetXaxis()->SetTitleColor(4);
h3->GetYaxis()->SetTitleColor(4);
h3->GetXaxis()->SetTitleOffset(0.7);
h3->GetYaxis()->SetTitleOffset(0.8);
h3->GetXaxis()->SetLabelSize(0.05);
h3->GetYaxis()->SetLabelSize(0.05);

TH1F *h4=new TH1F("PtKs0","",50, 0, 2.0);
h4->SetStats(kFALSE);
h4->GetYaxis()->SetTitle(" d#sigma/dPt2  ");
h4->GetXaxis()->SetTitle(" Pt2 ");

TH1F *h5=new TH1F("PtLambda","",50, 0.0, 2.0);
h5->SetStats(kFALSE);
h5->GetYaxis()->SetTitle(" d#sigma/dPt2  ");
h5->GetXaxis()->SetTitle(" Pt2 ");

TH1F *h6=new TH1F("PtLambdabar","",50, 0.0, 2.0);
h6->SetStats(kFALSE);
h6->GetYaxis()->SetTitle(" d#sigma/dPt2  ");
h6->GetXaxis()->SetTitle(" Pt2 ");
*/

//-------------------------- Current histograms -------------------------
G4double DistrY[50][3]; for(G4int i=0; i<50; i++){for(G4int j=0; j<3; j++){DistrY[i][j]=0.;}}
G4double Ynbin(50.0), Ylow(-1.0), Yhigh( 6.0);
G4double dY=(Yhigh - Ylow)/Ynbin;

G4double DistrP[50][3]; for(G4int i=0; i<50; i++){for(G4int j=0; j<3; j++){DistrP[i][j]=0.;}}
G4double Pnbin(50), Plow( 0.0), Phigh(2.0);
G4double dP=(Phigh - Plow)/Pnbin;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, px, py, pt, theta;  // pz,
    G4VParticleChange* aChange = 0;

//  G4double E=energy+part->GetPDGMass();                                  // Elab Proj
//  G4double SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E); // per  Proj+N
//  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                         //      Proj+N

    // -------- Event loop
    G4cout<<"Events start "<<nevt<<G4endl;
//  G4Timer timer;
//  timer.Start();

//  G4int Ninelast=nevt;                              // Uzhi

//    G4double weight=1./nevt/(7./50.);
//    G4double weight1=1./nevt/(2./50.);
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

      G4LorentzVector Proj4Mom = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)), e0 + mass);
//    G4cout<<Proj4Mom<<G4endl;
      Proj4Mom.boost(-boostNN);

//      G4cout<<Proj4Mom<<G4endl;
//      G4double Pmax=labNN.mag()/2.; //Proj4Mom.vect().mag();
//      G4double  weightX=sigIn/Pmax/nevt/(2./50.)/3.1416*1000.;
//      G4cout<<Pmax<<"  "<<weightX<<G4endl;

//------------

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de);
      labv -= dee;

      G4int n = aChange->GetNumberOfSecondaries();     // Multiplicity of prod. part.
//      npart = n;

      if(verbose>=1) G4cout<<" Uzhi ------------ N prod. part "<<n<<G4endl;
//++++++++++++++++ Variables for each event +++++++++++++++++++++++++++++
      if((verbose > 0) && (n < 2)) {G4cout<<"Multiplicity of produced < 2!!!"<<G4endl;}

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

        px = mom.x()/GeV;
        py = mom.y()/GeV;
//      pz = mom.z()/GeV;
        pt = std::sqrt(px*px +py*py);
        e  = fm.e()/GeV; // - m;
        theta = mom.theta();

        G4double rapidity=fm.rapidity();
        G4int id = pd->GetPDGEncoding();

        G4int RapBin=(rapidity - Ylow)/dY;
        G4int PtBin=pt/dP;


        if(id==310) //h4->Fill(pt, weight1/pt);
        {
         if((0 <= RapBin)&&(RapBin <= 49)) DistrY[RapBin][0] +=1.;
         if((0 <= PtBin )&&(PtBin  <= 49)) DistrP[PtBin][0] +=1./pt;
        }

        if(id==3122) //h5->Fill(pt, weight1/pt);
        {
         if((0 <= RapBin)&&(RapBin <= 49)) DistrY[RapBin][1] +=1.;
         if((0 <= PtBin )&&(PtBin  <= 49)) DistrP[PtBin][1] +=1./pt;
        }

        if(id==-3122) //h6->Fill(pt, weight1/pt);
        {
         if((0 <= RapBin)&&(RapBin <= 49)) DistrY[RapBin][2] +=1.;
         if((0 <= PtBin )&&(PtBin  <= 49)) DistrP[PtBin][2] +=1./pt;
        }

//         if(id==310)  h1->Fill(rapidity, weight);
//        if(id==3122) h2->Fill(rapidity, weight);
//        if(id==-3122) h3->Fill(rapidity, weight);

        theta=theta*180./pi;
//        fm.boost(-bst);       // Transformation to CM system

//        G4double Xf=fm.z()/Pmax;
//        G4cout<<Xf<<G4endl;


//        G4double costcm = std::cos(fm.theta());
	de += e;

	//	delete sec;
        delete aChange->GetSecondary(i);

      } //     end of the loop on particles


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

//    timer.Stop();
//    G4cout << "  "  << timer.GetUserElapsed() << G4endl;

/*
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    h5->Write();
    h6->Write();
*/

    if(verbose > 0) {
      G4cout << "###### End of run # " << run << "     ######" << G4endl;
    }

//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi

//sigTot=sigTot; sigEl=sigEl; sigIn=sigIn;
sigIn=1.0;
//-----------------------------------------------------------------------
//----------------------------- Rapidity distributions------------------// Uzhi ++++
std::ofstream Ydistr("PbarXeY.dat",std::ios::out);
    G4cout<< "******** Rapidity ******* at Plab "<<Plab<<" Xin " << sigIn<<" "<< G4endl;
    Ydistr<<G4Version<<G4endl;
    Ydistr<<"  Y    K0s Lambda  LambdaBar "<< G4endl;

    for(G4int i=0; i <50; i++)
    {
     Ydistr<<Ylow + i*dY +dY/2.<<" ";

     DistrY[i][0] *= sigIn/nevt/dY  *2.0;        // 2 K0s
     DistrY[i][1] *= sigIn/nevt/dY;
     DistrY[i][2] *= sigIn/nevt/dY;
     Ydistr<<DistrY[i][0]<<" "<<DistrY[i][1]<<" "<<DistrY[i][2]<<G4endl;
    }
//

//----------------------------- Pt -Y distributions----------------------
std::ofstream PTdistr("PbarXePt.dat",std::ios::out);

    G4cout<< "******** Pt ******* at Plab "<<Plab<<" Xin " << sigIn<< G4endl;
    PTdistr<<G4Version<<G4endl;
    PTdistr<<"  Pt2 K0s Lambda LambdaBar"<< G4endl;

    for(G4int i=0; i <50; i++)
    {
     PTdistr<<Plow + i*dP +dP/2.<<" ";
     DistrP[i][0] *= sigIn/nevt/dP  *2.0;        // 2 K0s
     DistrP[i][1] *= sigIn/nevt/dP;
     DistrP[i][2] *= sigIn/nevt/dP;
     PTdistr<<DistrP[i][0]<<" "<<DistrP[i][1]<<" "<<DistrP[i][2]<<G4endl;
    }
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi


  } while(end);
  }  // End of job ----------------------------------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//  delete pFrame;
//  delete lFrame;
//  delete sFrame;
G4cout<<G4Version<<G4endl;

  delete mate;
  delete fin;
  delete phys;
  partTable->DeleteAllParticles();
//  f1.Write();

  G4cout << "###### End of test #####" << G4endl;
}
