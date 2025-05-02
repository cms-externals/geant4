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
//               FTF test: Pbar+P interaction --> Lambda, LambdaBar, K0s
//
//      edition  25.11.2017  A.Galoyan
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
/*
 //   Root initialization
 TFile f1("Hyperon3_6N.root","RECREATE");

TH1F *h1=new TH1F("YKs0","",50, -2.0, 2.0);
TH1F *h2=new TH1F("YLambda","",50, -2.0, 2.0);
TH1F *h3=new TH1F("YLambdabar","",50, -2.0, 2.0);
TH1F *h4=new TH1F("Pt2Ks0","",50, 0.0, 1.0);
TH1F *h5=new TH1F("Pt2Lambda","",50, 0.0, 1.0);
TH1F *h6=new TH1F("Pt2Lambdabar","",50, 0.0, 1.0);
TH1F *h7=new TH1F("XKs0","",50, -1.0, 1.0);
TH1F *h8=new TH1F("XLam","",50, -1.0, 1.0);
TH1F *h9=new TH1F("XLam_bar","",50, -1.0, 1.0);
*/
//-------------------------- Current histograms -------------------------
G4double DistrY[50][3]; for(G4int i=0; i<50; i++){for(G4int j=0; j<3; j++){DistrY[i][j]=0.;}}
G4double Ynbin(50.0), Ylow(-5.0), Yhigh( 5.0);
G4double dY=(Yhigh - Ylow)/Ynbin;

G4double DistrP[50][3]; for(G4int i=0; i<50; i++){for(G4int j=0; j<3; j++){DistrP[i][j]=0.;}}
G4double Pnbin(50), Plow( 0.0), Phigh(2.0);
G4double dP=(Phigh - Plow)/Pnbin;

G4double DistrX[50][3]; for(G4int i=0; i<50; i++){for(G4int j=0; j<3; j++){DistrX[i][j]=0.;}}
G4double Xnbin(50.0), Xlow(-1.0), Xhigh( 1.0);
G4double dX=(Xhigh - Xlow)/Xnbin;

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
    G4int Ntotal=nevt;
    G4int Ninelast=nevt;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//-------------------------------------------------------

    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, px, py, pt, pt2;
    G4VParticleChange* aChange = 0;

    // -------- Event loop
    G4cout<<"Events start "<<nevt<<G4endl;

//    G4Timer timer;
//    timer.Start();

//    G4double weight=sigIn/3.1416/nevt/(4./50.);
//    G4double weight1=sigIn/nevt*50.;
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
      Proj4Mom.boost(-boostNN);

 G4double MKs=493.677;    //mass of Ks0
 G4double Mlam=1116.683;  //mass of Lambda

 G4double Pmax1=std::sqrt(labNN.mag2()/4.-sqr(MKs));
 G4double Pmax2=std::sqrt(labNN.mag2()/4.-sqr(Mlam));

// G4double  weightX1=sigIn/Pmax1/nevt/(2./50.)/3.1416*1000. ;
// G4double  weightX2=sigIn/Pmax2/nevt/(2./50.)/3.1416*1000. ;
//------------

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de);
      labv -= dee;

      G4int n = aChange->GetNumberOfSecondaries();     // Multiplicity of prod. part.

      if(verbose>=1) G4cout<<" Uzhi ------------ N prod. part "<<n<<G4endl;
//++++++++++++++++ Variables for each event +++++++++++++++++++++++++++++
      if((verbose > 1) && (n < 2)) {G4cout<<"Multiplicity of produced < 2!!!"<<G4endl;}
      if(n < 2) Ninelast--;

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

        labv -= fm;   // For checking energy-momentum conservation

	// electron can come only from internal conversion
	// its mass should be added to initial state
        if(pd == electron) {

	  labv += G4LorentzVector(0.0,0.0,0.0,electron_mass_c2);
	}

        px = mom.x()/GeV;
        py = mom.y()/GeV;
        pt = std::sqrt(px*px +py*py); pt2 = pt*pt;

        G4int Pt2Bin=pt2/dP; if(Pt2Bin < 0) Pt2Bin=0; if(Pt2Bin > 49) Pt2Bin=49;

        e  = fm.e()/GeV; // - m;

        G4double rapidity;
if(((11000.0 < Plab) && (Plab < 13000.0)))
{rapidity=fm.rapidity();                      // Lab rapidity

 fm.boost(-bst);}
else
{fm.boost(-bst);                              // Transformation to CM system
 rapidity=fm.rapidity();                      // CMS rapidity
}

        G4int RapBin=(rapidity - Ylow)/dY; // if(RapBin < 0.) RapBin=0; if(RapBin > 49) RapBin=49;

        G4int id = pd->GetPDGEncoding();

        if(id==310)  // K0s meson
        {
         if((0 <=RapBin)&&(RapBin <=49))  DistrY[RapBin][0] +=1.;
         DistrP[Pt2Bin][0] +=1.;

         G4double Xf = fm.z()/Pmax1;
         G4int XfBin=(Xf+1.0)/dX;    // if(XfBin < 0) XfBin=0; if(XfBin > 49) XfBin=49;
         if(Plab < 1000.0)
         {if((0 <= XfBin) && (XfBin <= 49)) DistrX[XfBin][0] +=1.;}
         else
         {if((0 <= XfBin) && (XfBin <= 49)) DistrX[XfBin][0] +=fm.e()/Pmax1;}
        }

        if(id==3122) // Lambda
        {
         DistrY[RapBin][1] +=1.;
         DistrP[Pt2Bin][1] +=1.;
         G4double Xf = fm.z()/Pmax2;
         G4int XfBin=(Xf+1.0)/dX;     if(XfBin < 0) XfBin=0; if(XfBin > 49) XfBin=49;
         if((0 <= XfBin) && (XfBin <= 49)) DistrX[XfBin][1] +=fm.e()/Pmax2;
        }

        if(id==-3122) // Anti-Lambda
        {
         DistrY[RapBin][2] +=1.;
         DistrP[Pt2Bin][2] +=1.;

         G4double Xf = fm.z()/Pmax2;
         G4int XfBin=(Xf+1.0)/dX;    if(XfBin < 0) XfBin=0; if(XfBin > 49) XfBin=49;

         if((0 <= XfBin) && (XfBin <= 49)) DistrX[XfBin][2] +=fm.e()/Pmax2;
        }

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
       G4int Uzhi;  G4cin>>Uzhi;                                                    // Uzhi
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

//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++
G4cout<< "***********************************************************"<< G4endl;

G4cout<<"nevt Ntotal Ninelast "<<nevt<<" "<<Ntotal<<" "<<Ninelast<<G4endl;
G4cout<<"Plab "<<Plab/GeV<<" SigIn "<<sigIn<<G4endl;

//-------------------------------------------------------------------


    if(verbose > 0) {
      G4cout << "###### End of run # " << run << "     ######" << G4endl;
    }
//++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi

//sigTot=sigTot; sigEl=sigEl; sigIn=sigIn;
//-----------------------------------------------------------------------

/*
G4String Version=G4Version; G4String Ins(nameGen+" $");
G4int Kins=0; do{Kins++;} while(Version[Kins]!='$');
Version.replace(Kins-1,Ins.size(),Ins);
*/

G4double FactorY=1.; if(( 3000.0 < Plab) && (Plab <  4000.0)) FactorY=1000./3.1416;
                     if((11000.0 < Plab) && (Plab < 13000.0)) FactorY=1000./3.1416;

G4double FactorX=1.; if((  700.0 < Plab) && (Plab <   800.0)) FactorX=3.1416;
                     if(( 3000.0 < Plab) && (Plab <  4000.0)) FactorX=1000.;
                     if((11000.0 < Plab) && (Plab < 13000.0)) FactorX=1000.;

G4double FactorP=1.; if((11000.0 < Plab) && (Plab < 13000.0)) FactorP=1000.;

//----------------------------- Rapidity distributions------------------// Uzhi ++++
std::ofstream Ydistr("Ydistr.dat",std::ios::out);
    G4cout<< "******** Rapidity ******* at Plab "<<Plab<<" Xin " << sigIn<<" "<< G4endl;
    Ydistr<<G4Version<<G4endl;
    Ydistr<<"  Y    K0s Lambda  LambdaBar "<< G4endl;

    for(G4int i=0; i <50; i++)
    {
     Ydistr<<Ylow + i*dY +dY/2.<<" ";

     DistrY[i][0] *= sigIn/Ninelast/dY*FactorY;
     DistrY[i][1] *= sigIn/Ninelast/dY*FactorY;
     DistrY[i][2] *= sigIn/Ninelast/dY*FactorY;
     Ydistr<<DistrY[i][0]<<" "<<DistrY[i][1]<<" "<<DistrY[i][2]<<G4endl;
    }
//

//----------------------------- Pt -Y distributions----------------------
std::ofstream PTdistr("PTdistr.dat",std::ios::out);

    G4cout<< "******** Pt-Y ******* at Plab "<<Plab<<" Xin " << sigIn<< G4endl;
    PTdistr<<G4Version<<G4endl;
    PTdistr<<"  Pt2 K0s Lambda LambdaBar"<< G4endl;

    for(G4int i=0; i <50; i++)
    {
     PTdistr<<Plow + i*dP +dP/2.<<" ";
     DistrP[i][0] *= sigIn/Ninelast/dP*FactorP;
     DistrP[i][1] *= sigIn/Ninelast/dP*FactorP;
     DistrP[i][2] *= sigIn/Ninelast/dP*FactorP;
     PTdistr<<DistrP[i][0]<<" "<<DistrP[i][1]<<" "<<DistrP[i][2]<<G4endl;
    }

// ----------------------------- Xf distributions----------------------// Uzhi ++++
std::ofstream XFdistr("XFdistr.dat",std::ios::out);
    G4cout<< "******** Xf distr ******* at Plab "<<Plab<<" Xin " << sigIn<< G4endl;

    XFdistr<<G4Version<<G4endl;
    XFdistr<<"  Xf K0s Lambda LambdaBar"<< G4endl;

    for(G4int i=0; i <50; i++)
    {
     XFdistr<<Xlow + i*dX +dX/2.<<" ";
     DistrX[i][0] *= sigIn/Ninelast/dX/pi*FactorX;
     DistrX[i][1] *= sigIn/Ninelast/dX/pi*FactorX;
     DistrX[i][2] *= sigIn/Ninelast/dX/pi*FactorX;
     XFdistr<<DistrX[i][0]<<" "<<DistrX[i][1]<<" "<<DistrX[i][2]<<G4endl;
    }

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi

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
//  f1.Write();

  G4cout << "###### End of test #####" << G4endl;
}
