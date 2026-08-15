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
#include "G4ChipsComponentXS.hh"  // Uzhi 29.01.13
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "FTFtest1.icc"
#include "UZHI_diffraction.hh"

#include <fstream>
#include <iomanip>

int main(int argc, char** argv)
{
  CLHEP::RanluxEngine defaultEngine(1234567, 4);
  CLHEP::HepRandom::setTheEngine(&defaultEngine);

  G4cout << "========================================================" << G4endl;
  G4cout << "======              FTF Test Start              ========" << G4endl;
  G4cout << "========================================================" << G4endl;

  // -------------------------------------------------------------------
  // Control on input

  if (argc < 2)
  {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();
  G4String fname = argv[1];
  fin->open(fname.c_str());
  if (!fin->is_open())
  {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  //-----------------------------------------------------------------------
#include "FTFtest2.icc"  // Initialization
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  G4double sigTot = 0;
  G4double sigEl = 0;
  G4double sigIn = 0;

  //-------------------------- Global histograms  -------------------------

  G4int Uzhi_run = 0;
  G4double TopoUzhi[50][16];
  for (G4int ii = 0; ii < 50; ii++)
  {
    for (G4int jj = 0; jj < 16; jj++)
    {
      TopoUzhi[ii][jj] = 0.;
    }
  };

  G4double Xs[50][34];
  for (G4int ii = 0; ii < 50; ii++)
  {
    for (G4int jj = 0; jj < 34; jj++)
    {
      Xs[ii][jj] = 0.;
    }
  };
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;

  for (G4int run = 0; run < 100; run++)
  {
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //-------------------------- Current histograms -------------------------

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do
    {
#include "FTFtest3.icc"  // -------- Read input file
#include "FTFtest4.icc"  // -------- Start run processing
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      G4cout << "cross(mb)in= " << cross_sec * 1000. / barn << G4endl
             << "cross(mb)el= " << cross_secel * 1000. / barn << G4endl << G4endl;

      cross_inel = cross_sec - cross_secel;  // +++++++++++++++++++++++++

      cross_sec /= millibarn;  // Inel Cross section in mb
      cross_secel /= millibarn;  // Elas Cross section in mb
      cross_inel /= millibarn;  // Inel Cross section in mb

      G4cout << "Element A Z N: " << A << " " << Z << " " << A - Z << G4endl;
      G4cout << "Proposed Xs (mb): Tot El In: " << cross_sec << " " << cross_secel << " "
             << cross_inel << G4endl;

      //---------------------------------------------------------------------------
      // Kossov cross sections      ---------------------------
      G4double chipsTot, chipsEl, chipsIn;

      static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();
      G4ChipsComponentXS* CHIPSxsManager = _instance;

      G4bool CHIPapplic = false;  // true;                          //false;   Uzhi 29.01.13
      if (CHIPapplic)
      {
        chipsTot = CHIPSxsManager->GetTotalElementCrossSection(part, energy, Z, A - Z);
        chipsEl = CHIPSxsManager->GetElasticElementCrossSection(part, energy, Z, A - Z);
        chipsIn = CHIPSxsManager->GetInelasticElementCrossSection(part, energy, Z, A - Z);
        chipsTot /= millibarn;
        chipsEl /= millibarn;
        chipsIn /= millibarn;

        G4cout << "CHIPS cross sections are used:----------------------" << G4endl
               << "Plab          Total        Elastic      Inelastic" << G4endl;
        G4cout << " " << Plab / GeV << " " << chipsTot << " " << chipsEl << " " << chipsIn << G4endl
               << G4endl;

        sigTot = chipsTot;
        sigEl = chipsEl;
        sigIn = chipsIn;
      }
      else
      {
        sigTot = cross_sec;
        sigEl = cross_secel;
        sigIn = cross_inel;

        G4cout << "Proposed Xs (mb) are used: Tot El In: " << cross_sec << " " << cross_secel << " "
               << cross_inel << G4endl;
      }

      //+++++++++++++++++++++++++++++++++ For each energy +++++++++++++++++++++
      TopoUzhi[Uzhi_run][0] = Plab / GeV;

      Xs[Uzhi_run][0] = Plab / GeV;
      Xs[Uzhi_run][1] = sigTot;
      Xs[Uzhi_run][2] = sigEl;
      Xs[Uzhi_run][3] = sigIn;

      G4int Ntotal = nevt;
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-------------------------------------------------------

      const G4DynamicParticle* sec = 0;
      G4ParticleDefinition* pd;
      G4ThreeVector mom;
      G4LorentzVector labv, fm;
      G4double e, theta;
      G4VParticleChange* aChange = 0;

      //  G4double E=energy+part->GetPDGMass();                                  // Elab Proj
      //  G4double SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E); // per  Proj+N
      //  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                         //      Proj+N

      // -------- Event loop
      G4cout << "Events start " << nevt << G4endl;
      //  G4int Ninelast=nevt;                              // Uzhi
      //=================================================================
      for (G4int iter = 0; iter < nevt; ++iter)
      {
        //=================================================================
        if (verbose > 0) G4cout << "Start events loop***********************" << G4endl;

        if (verbose >= 1 || iter == modu * (iter / modu))
        {
          G4cout << "### " << iter << "-th event start " << Plab / GeV << G4endl;
        }

        if (saverand)
        {
          defaultEngine.saveStatus("initial.conf");
        }

        G4double e0 = energy;
        do
        {
          if (sigmae > 0.0) e0 = G4RandGauss::shoot(energy, sigmae);
        } while (e0 < 0.0);

        dParticle.SetKineticEnergy(e0);

        gTrack->SetStep(step);
        gTrack->SetKineticEnergy(e0);
        G4double amass = phys->GetNucleusMass();
        // note: check of 4-momentum balance for CHIPS is not guranteed due to
        // unknown isotope
        aChange = proc->PostStepDoIt(*gTrack, *step);

        G4double mass = part->GetPDGMass();

        if (ionParticle)
        {
          e0 /= ionA;
          G4double mass_N = 938. * MeV;  // Init 4-mom
          labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0 * (e0 + 2. * mass_N)),  //   NN
                                 e0 + mass_N + mass_N);
        }
        else
        {
          labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0 * (e0 + 2. * mass)),  //   hA
                                 e0 + mass + amass);
        }

        G4ThreeVector bst = labv.boostVector();  // To CMS NN in AA or hA
        //------------
        G4LorentzVector labNN(0.0, 0.0, std::sqrt(e0 * (e0 + 2. * mass)), e0 + mass + amass);
        G4ThreeVector boostNN = labNN.boostVector();
        //------------

        // take into account local energy deposit
        G4double de = aChange->GetLocalEnergyDeposit();
        G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de);
        labv -= dee;

        G4int n = aChange->GetNumberOfSecondaries();  // Multiplicity of prod. part.

        if (verbose >= 1) G4cout << " Uzhi ------------ N prod. part " << n << G4endl;
        //++++++++++++++++ Variables for each event +++++++++++++++++++++++++++++
        if ((verbose > 0) && (n < 2))
        {
          G4cout << "Multiplicity of produced < 2!!!" << G4endl;
        }
        if (n < 2)
        {
          Ntotal--;
        }

        //    G4int nbar = 0;
        G4int Npim = 0;
        G4int Npip = 0;
        G4int Npi0 = 0;

        G4int NKm = 0;
        G4int NKp = 0;
        G4int NK0s = 0;
        G4int NK0l = 0;

        G4int NLambda = 0;
        G4int NLambdaBar = 0;

        G4int NSigma0 = 0;
        G4int NSigma0Bar = 0;

        G4int NSigma_p = 0;
        G4int NSigma_pBar = 0;

        G4int NSigma_m = 0;
        G4int NSigma_mBar = 0;

        G4int NXi0 = 0;
        G4int NXi0Bar = 0;

        G4int NXi_m = 0;
        G4int NXi_mBar = 0;

        G4int Nproton = 0;
        G4int NprotonBar = 0;

        G4int Nneutron = 0;
        G4int NneutronBar = 0;

        G4int Neta = 0;
        G4int Neta_prime = 0;

        G4int Ngamma = 0;

        G4int Ncharged = 0;  // Uzhi
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        for (G4int i = 0; i < n; ++i)  // Loop over produced particles
        {
          sec = aChange->GetSecondary(i)->GetDynamicParticle();
          pd = sec->GetDefinition();
          G4String pname = pd->GetParticleName();

          if (verbose >= 1)
            G4cout << " Part  " << i << " " << pname << " " << sec->Get4Momentum() / GeV
                   << sec->Get4Momentum().mag() / GeV << G4endl;

          fm = sec->Get4Momentum();

          mom = sec->GetMomentum();
          //      G4double mas = pd->GetPDGMass();
          //	G4double p = mom.mag();

          labv -= fm;  // For checking energy-momentum conservation

          // electron can come only from internal conversion
          // its mass should be added to initial state
          if (pd == electron)
          {
            labv += G4LorentzVector(0.0, 0.0, 0.0, electron_mass_c2);
          }

          e = fm.e() - m;
          theta = mom.theta();

          //        G4double CosTheta=std::cos(theta);

          theta = theta * 180. / pi;

          fm.boost(-bst);

          //        G4double costcm = std::cos(fm.theta());
          //+++++++++++++++++ For each particle in the event ++++++++++++++++++++++
          if (pname == "pi-")
          {
            Npim++;
            Ncharged++;
          }
          else if (pname == "pi+")
          {
            Npip++;
            Ncharged++;
          }
          else if (pname == "pi0")
          {
            Npi0++;
          }

          else if (pname == "kaon-")
          {
            NKm++;
            Ncharged++;
          }
          else if (pname == "kaon+")
          {
            NKp++;
            Ncharged++;
          }
          else if (pname == "kaon0S")
          {
            NK0s++;
          }
          else if (pname == "kaon0L")
          {
            NK0l++;
          }

          else if (pname == "lambda")
            NLambda++;
          else if (pname == "anti_lambda")
            NLambdaBar++;

          else if (pname == "sigma0")
            NSigma0++;
          else if (pname == "anti_sigma0")
            NSigma0Bar++;
          else if (pname == "sigma+")
          {
            NSigma_p++;
            Ncharged++;
          }
          else if (pname == "sigma-")
          {
            NSigma_m++;
            Ncharged++;
          }
          else if (pname == "anti_sigma+")
          {
            NSigma_pBar++;
            Ncharged++;
          }
          else if (pname == "anti_sigma-")
          {
            NSigma_mBar++;
            Ncharged++;
          }

          else if (pname == "xi0")
            NXi0++;
          else if (pname == "xi-")
          {
            NXi_m++;
            Ncharged++;
          }
          else if (pname == "anti_xi0")
            NXi0Bar++;
          else if (pname == "anti_xi-")
          {
            NXi_mBar++;
            Ncharged++;
          }

          else if (pname == "omega-")
          {
            Ncharged++;
          }
          else if (pname == "anti_omega-")
          {
            Ncharged++;
          }

          else if (pname == "proton")
          {
            Nproton++;
            Ncharged++;
          }
          else if (pname == "anti_proton")
          {
            NprotonBar++;
            Ncharged++;
          }

          else if (pname == "neutron")
            Nneutron++;
          else if (pname == "anti_neutron")
            NneutronBar++;
          else if (pname == "eta")
            Neta++;
          else if (pname == "eta_prime")
            Neta_prime++;
          else if (pname == "gamma")
            Ngamma++;
          else if (pd->GetParticleType() == "nucleus")
            ;
          else
          {
            G4cout << "****Found " << pname;
            if (pd->IsShortLived()) G4cout << "  is Shortlived";
            G4cout << G4endl << " .... width, Shortlived... " << pd->GetPDGWidth() << " "
                   << pd->IsShortLived() << G4endl;
            pd->DumpTable();
          }
          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          de += e;

          //	delete sec;
          delete aChange->GetSecondary(i);

        }  //     end of the loop on particles
        //+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++
        G4int JUzhi = (Ncharged + 2) / 2;  // Uzhi
        if (JUzhi < 16) TopoUzhi[Uzhi_run][JUzhi]++;  // Uzhi

        if (JUzhi == 2)
        {  // Uzhi Apr. 2016
          G4int FSbaryon = Nproton + Nneutron + NprotonBar + NneutronBar + NLambda + NLambdaBar
                           + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar
                           + NXi0 + NXi0Bar + NXi_m + NXi_mBar;
          if (FSbaryon != 0) TopoUzhi[Uzhi_run][JUzhi]--;
        }
        //    G4int SumQ=0;                                               // Uzhi
        /*
        //--------5 N Nbar
              Nother=     Nproton    + Nneutron    +
                          NprotonBar + NneutronBar +
                          Npim + Npip+ Npi0        +
                          NKm + NKp  + NK0s + NK0l +
                          NLambda    + NLambdaBar  +
                          NSigma0    + NSigma0Bar  +
                          NSigma_p   + NSigma_pBar +
                          NSigma_m   + NSigma_mBar +
                          NXi0       + NXi0Bar     +
                          NXi_m      + NXi_mBar    +
                          Neta       + Neta_prime  +
                          Ngamma                    ;

        if((NprotonBar == 1) && (Nproton == 1) && (Nother == 0))
        {Xs[Uzhi_run][]++;} // G4cout<<" 5 "<<G4endl;}

        */
        // if((n == 2) && (NSigma_m == 1)) {G4cout<<" NSigma_mBar "<<NSigma_mBar<<" NSigma_pBar
        // "<<NSigma_pBar<<G4endl; G4int Uzhi; G4cin>>Uzhi;}
        //--------4   Elastic Pbar P
        G4int Nother = Nneutron + NneutronBar + Npim + Npip + Npi0 + NKm + NKp + NK0s + NK0l
                       + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                       + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta
                       + Neta_prime + Ngamma;
        if ((NprotonBar == 1) && (Nproton == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][4]++;
        }  // G4cout<<" 4 "<<G4endl;} // TopoUzhi[Uzhi_run][JUzhi]--;

        //--------5 N Nbar

        Nother = Nproton + NprotonBar + Npim + Npip + Npi0 + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NneutronBar == 1) && (Nneutron == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][5]++;
        }  // G4cout<<" 5 "<<G4endl; G4int Uzhi; G4cin >> Uzhi;}

        //--------6 L Lbar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NKm + NKp
                 + NK0s + NK0l + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NLambdaBar == 1) && (NLambda == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][6]++;
        }  // G4cout<<" 6 L Lbar"<<G4endl; G4int Uzhi; G4cin>>Uzhi;}			//Vova

        //--------7 P Pi- Nbar
        Nother = Nneutron + NprotonBar + Npip + Npi0 + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Nproton == 1) && (NneutronBar == 1) && (Npim == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][7]++;
        }  // G4cout<<" 7 "<<G4endl;}

        //--------8 N Pi+ Pbar
        Nother = Nproton + NneutronBar + Npim + Npi0 + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NprotonBar == 1) && (Nneutron == 1) && (Npip == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][8]++;
        }  // G4cout<<" 8 "<<G4endl;}

        //--------9 P Pi0 Pbar
        Nother = Nneutron + NneutronBar + Npim + Npip + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NprotonBar == 1) && (Nproton == 1) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][9]++;
        }  // G4cout<<" 9 "<<G4endl;}

        //--------10 P Pi+ Pi- Pbar
        Nother = Nneutron + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l + NLambda + NLambdaBar
                 + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0
                 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NprotonBar == 1) && (Nproton == 1) && (Npip == 1) && (Npim == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][10]++;
        }  // G4cout<<" 10 "<<G4endl;}

        //--------11 P Pi+ Pi- Pi0 Pbar
        Nother = Nneutron + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda + NLambdaBar + NSigma0
                 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar
                 + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NprotonBar == 1) && (Nproton == 1) && (Npip == 1) && (Npim == 1) && (Npi0 == 1)
            && (Nother == 0))
        {
          Xs[Uzhi_run][11]++;
        }  // G4cout<<" 11 "<<G4endl;}

        //--------12 Pi+ Pi-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;
        if ((Npip == 1) && (Npim == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][12]++;
        }  //  G4cout<<" 12 Pi+ Pi-"<<Xs[Uzhi_run][12]<<G4endl; G4int Uzhi; G4cin>>Uzhi;}

        //--------13 K+ K-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NKm == 1) && (NKp == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][13]++;
        }  // G4cout<<" 13 K+ K-"<<G4endl;}

        //-------------14 Pi+ Pi- Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][14]++;
        }  // G4cout<<" 14 Pi+ Pi- Pi0"<<Xs[Uzhi_run][14]<<G4endl; G4int Uzhi; G4cin>>Uzhi;}

        //-----------15 2Pi+ 2Pi-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Nother == 0))
        {
          Xs[Uzhi_run][15]++;
        }  // G4cout<<" 15 2Pi+ 2Pi-"<<Xs[Uzhi_run][15]<<G4endl; G4int Uzhi; G4cin>>Uzhi;}

        //-----------16 2Pi+ 2Pi- Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][16]++;
        }  // G4cout<<" 16 2Pi+ 2Pi- Pi0 "<<Xs[Uzhi_run][16]<<G4endl; G4int Uzhi; G4cin>>Uzhi;}

        //-----------17 3Pi+ 3Pi-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 3) && (Npim == 3) && (Nother == 0))
        {
          Xs[Uzhi_run][17]++;
        }  // G4cout<<" 17 3Pi+ 3Pi-"<<G4endl;  G4int Uzhi; G4cin>>Uzhi;}

        //-----------18 3Pi+ 3Pi- Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 3) && (Npim == 3) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][18]++;
        }  // G4cout<<" 18 3Pi+ 3Pi- Pi0"<<G4endl;  G4int Uzhi; G4cin>>Uzhi;}

        //--------19 Lambda Pi0 LambdaBar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar
                 + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NLambda == 1) && (NLambdaBar == 1) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][19]++;
        }  // G4cout<<" 19 L Pi0 aL "<<G4endl; G4int Uzhi; G4cin>>Uzhi;}

        //--------20 Lambda Pi+ Pi- LambdaBar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0
                 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NLambda == 1) && (NLambdaBar == 1) && (Npip == 1) && (Npim == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][20]++;
        }  // G4cout<<" 20 L Pi+ Pi- aL "<<G4endl; G4int Uzhi; G4cin>>Uzhi;}

        //--------21 Sigma+ Sigma-Bar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NKm + NKp
                 + NK0s + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar +

                 NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((NSigma_p == 1) && (NSigma_pBar == 1) && (Nother == 0))  // NSigma_mBar -> NSigma_pBar
        {
          Xs[Uzhi_run][21]++;
        }  // G4cout<<" 21 "<<G4endl;}

        //--------22 Sigma- Sigma+Bar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NKm + NKp
                 + NK0s + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p
                 + NSigma_pBar +

                 NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NSigma_m == 1) && (NSigma_mBar == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][22]++;
          G4cout << " 22 " << G4endl;
          G4int Uzhi;
          G4cin >> Uzhi;
        }  //} // G4cout<<" 22 "<<G4endl;}

        //--------23 Sigma0 LambsaBar + Lambda Sigma0Bar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NKm + NKp
                 + NK0s + NK0l + NLambda + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NSigma0 == 1) && (NLambdaBar == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][23]++;
        }  // G4cout<<" 23  Sigma0 Lbar "<<G4endl; G4int Uzhi; G4cin>>Uzhi;}			//Vova

        //--------23 Sigma0 LambsaBar + Lambda Sigma0Bar
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NKm + NKp
                 + NK0s + NK0l + NLambdaBar + NSigma0 + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NLambda == 1) && (NSigma0Bar == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][23]++;
        }  // G4cout<<" 23 "<<G4endl;}

        //--------24 Sigma+Pi-LambdaBar Sigma-BarPi+Lambda
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npip + Npi0 + NKm + NKp + NK0s
                 + NK0l + NLambda + NSigma0 + NSigma0Bar + NSigma_pBar + NSigma_m + NSigma_mBar
                 + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NSigma_p == 1) && (Npim == 1) && (NLambdaBar == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][24]++;
        }  // G4cout<<" 24 "<<G4endl;}

        //--------24 Sigma+Pi-LambdaBar Sigma-BarPi+Lambda
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npi0 + NKm + NKp + NK0s
                 + NK0l + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NSigma_pBar == 1) && (Npip == 1) && (NLambda == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][24]++;
          G4cout << " 24 " << G4endl;
          G4int Uzhi;
          G4cin >> Uzhi;
        }  // G4cout<<" 24 "<<G4endl;}

        //--------25 Lambda K0 Nbar + N K0 Lambda0bar
        Nother = Nproton + Nneutron + NprotonBar + Npim + Npip + Npi0 + NKm + NKp + NLambdaBar
                 + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0
                 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((NLambda == 1) && (NneutronBar == 1) && (NK0s + NK0l == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][25]++;
        }  // G4cout<<" 25 "<<G4endl;}

        //--------25 Lambda K0 Nbar + N K0 Lambda0bar
        Nother = Nproton + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NK0s + NK0l + NLambda
                 + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0
                 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Nneutron == 1) && (NLambdaBar == 1) && (NK0s + NK0l == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][25]++;
          G4cout << " 25 " << G4endl;
          G4int Uzhi;
          G4cin >> Uzhi;
        }  // G4cout<<" 25 "<<G4endl;}

        //--------26 Pi+ Pi- K+ K-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (NKp == 1) && (NKm == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][26]++;
        }  // G4cout<<" 26 "<<G4endl;}

        //--------27 Pi0 Pi- K+ K0s
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npip + NKm + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npi0 == 1) && (Npim == 1) && (NKp == 1) && (NK0s == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][27]++;
        }  // G4cout<<" 27 "<<G4endl;}

        //--------28 Pi+ Pi- K0l K0s
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (NK0l == 1) && (NK0s == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][28]++;
        }  // G4cout<<" 28 "<<G4endl;}

        //--------29 Pi+ Pi- 2 K0s
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (NK0s == 2) && (Nother == 0))
        {
          Xs[Uzhi_run][29]++;
        }  // G4cout<<" 29 "<<G4endl;}

        //--------30 Pi+ PI- Pi0 K+ K-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar +

                 NK0s + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 1) && (NKp == 1) && (NKm == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][30]++;
        }  // G4cout<<" 30 "<<G4endl;}

        //--------31 Pi+ 2Pi- K+ K0s
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 2) && (NKp == 1) && (NK0s == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][31]++;
        }  // G4cout<<" 31 "<<G4endl;}

        //--------31 Pi- 2Pi+ K- K0s  c.c. of 31
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKp + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;
        if ((Npip == 2) && (Npim == 1) && (NKm == 1) && (NK0s == 1) && (Nother == 0))
        {
          Xs[Uzhi_run][31]++;
        }  // G4cout<<" 31 "<<G4endl;}

        //--------32 Pi+ Pi- PI0 Kl Ks
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar +

                 NKm + NKp + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 1) && (NK0l == 1) && (NK0s == 1)
            && (Nother == 0))
        {
          Xs[Uzhi_run][32]++;
        }  // G4cout<<" 32 "<<G4endl;}

        //-------- 33 Pi+ Pi- Pi0 2K0s
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar +

                 NKm + NKp + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p
                 + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta
                 + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 1) && (NK0s == 2) && (Nother == 0))
        {
          Xs[Uzhi_run][33]++;
        }  // G4cout<<" 33 "<<G4endl;}

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if (verbose > 0) G4cout << "Energy/Momentum balance= " << labv << G4endl;

        aChange->Clear();

        if (verbose > 0)
        {
          G4cout << "End event =====================================" << Plab << G4endl;  // Uzhi
          G4int Uzhi_i;  // Uzhi
          G4cin >> Uzhi_i;  // Uzhi
        }
        //

      }  // End of the event loop ------------------------------------

      timer->Stop();
      G4cout << "  " << *timer << G4endl;
      delete timer;

      //++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++
      // sigTot=sigTot; sigEl=sigEl;
      sigTot = sigIn;
      // sigTot=1.;
      if (Ntotal != 0)
      {
        for (G4int ii = 1; ii < 16; ii++)
        {
          TopoUzhi[Uzhi_run][ii] *= sigTot / Ntotal;
        };  // sigIn->sigTot
        for (G4int ii = 4; ii < 34; ii++)
        {
          Xs[Uzhi_run][ii] *= sigTot / Ntotal;
        }
      }

      //-------------------------------------------------------------------
      G4cout << "***********************************************************" << G4endl;

      G4cout << "nevt Ninel " << nevt << " " << Ntotal << G4endl;
      G4cout << "Plab " << Plab / GeV << " SigIn " << sigIn << G4endl;
      //-------------------------------------------------------------------

      if (verbose > 0)
      {
        G4cout << "###### End of run # " << run << "     ######" << G4endl;
      }
      //++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi
      Uzhi_run++;  // Uzhi
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi
      //    G4cerr << "###### End of run # " << run << "     ######" << G4endl;

    } while (end);
  }  // End of job ----------------------------------------------
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //----------------------------- Write distributions------------------

  std::ofstream PbarPtopo("PbarPtopo.dat", std::ios::out);

  // PPtopo<<"################# Topological cross-sections #################"<<G4endl;
  PbarPtopo << G4Version << G4endl;
  PbarPtopo
    << " Plab         S0           S2           S4           S6           S8           S10"
    << "          S12          S14          S16          S18          S20          S22       "
    << "   S24          S26          S28 (mb)" << G4endl;

  for (G4int ii = 0; ii < Uzhi_run; ii++)
  {
    G4double SumXs = 0.;
    for (G4int jj = 0; jj < 16; jj++)
    {
      PbarPtopo << " " << TopoUzhi[ii][jj];
      if (jj != 0)
      {
        SumXs += TopoUzhi[ii][jj];
      };
    };
    PbarPtopo << G4endl;
    if (SumXs != 0.)
    {
      for (G4int jj = 0; jj < 16; jj++)  // Normalization
      {
        if (jj != 0)
        {
          TopoUzhi[ii][jj] /= SumXs;
        };
      };
    };
  };

  std::ofstream PbarPchan1("PbarPchan1.dat", std::ios::out);

  // PbarPchan1 ################### Channels with baryons in final sttates
  // ###################################
  PbarPchan1 << G4Version << G4endl;
  PbarPchan1 << "Plab  Xtot  Xel  Xin  El   NaN   LaL  P_PimNb  N_Pip_Pb  P_Pi0_Pb  P_PipPim_Pb  "
                "P_PipPimPi0_Pb  Lpi0aL  LpppmAL"
             << G4endl;

  for (G4int ii = 0; ii < Uzhi_run; ii++)
  {
    for (G4int jj = 0; jj < 12; jj++)
    {
      PbarPchan1 << " " << Xs[ii][jj];
    }
    PbarPchan1 << " " << Xs[ii][19];
    PbarPchan1 << " " << Xs[ii][20];
    //   G4double Sinv=1.76+1.88*std::sqrt(0.88+Xs[ii][0]*Xs[ii][0]);
    //   G4double SSss=std::sqrt(Sinv);
    //   PbarPchan1<<" "<<Sinv<<" "<<SSss;
    PbarPchan1 << G4endl;
  };
  // PbarPchan1<<G4endl<<G4endl;

  std::ofstream PbarPchan2("PbarPchan2.dat", std::ios::out);
  // PbarPchan2 ################### Channels with only mesons in final sttates
  // ###################################
  PbarPchan2 << G4Version << G4endl;
  PbarPchan2 << " Plab      Pi+Pi-     K+K-    Pi+Pi-Pi0    2Pi+2Pi-   2Pi+2Pi-Pi0 ";
  PbarPchan2 << "  3Pi+3Pi-   3Pi+3Pi-Pi0 " << G4endl;

  for (G4int ii = 0; ii < Uzhi_run; ii++)  //--------------------------------- Uzhi
  {
    PbarPchan2 << " " << Xs[ii][0];
    for (G4int jj = 12; jj < 19; jj++)
    {
      PbarPchan2 << " " << Xs[ii][jj];
    };
    //   G4double Sinv=1.76+1.88*std::sqrt(0.88+Xs[ii][0]*Xs[ii][0]);
    //   G4double SSss=std::sqrt(Sinv);
    //   PbarPchan2<<" "<<Sinv<<" "<<SSss<<" 0";
    PbarPchan2 << G4endl;
  };

  std::ofstream PbarPchan3("PbarPchan3.dat", std::ios::out);
  // PbarPchan3 ################### Channels with baryons in final sttates
  // ###################################
  PbarPchan3 << G4Version << G4endl;
  PbarPchan3 << " Plab      LamLamBar     S+S-Bar  S-S+Bar S0LamBar+cc LamPi0LamBar S+Pi-LamBar+cc "
                "LamK0Nbar+cc "
             << G4endl;
  //             0           6           21        22      23           19          24 25
  for (G4int ii = 0; ii < Uzhi_run; ii++)  //--------------------------------- Uzhi
  {
    PbarPchan3 << " " << Xs[ii][0] << " " << Xs[ii][6];
    for (G4int jj = 21; jj < 24; jj++)
    {
      PbarPchan3 << " " << Xs[ii][jj];
    };

    PbarPchan3 << " " << Xs[ii][19] << " " << Xs[ii][24] << " " << Xs[ii][25];

    PbarPchan3 << G4endl;
  };

  std::ofstream PbarPchan4("PbarPchan4.dat", std::ios::out);
  // PbarPchan4 ################### Channels with mesons in final sttates
  // ###################################
  PbarPchan4 << G4Version << G4endl;
  PbarPchan4 << " Plab   Pi+Pi-K+K-  Pi0Pi-K+K0s  Pi+Pi-K0lK0s  Pi+Pi-2K0s  Pi+PI-Pi0K+K-  "
                "Pi+2Pi-K+K0s  Pi+Pi-PI0KlKs  Pi+Pi-Pi02K0s    "
             << G4endl;
  //                      26           27           28            29          30             31 32
  //                      33
  for (G4int ii = 0; ii < Uzhi_run; ii++)  //--------------------------------- Uzhi
  {
    PbarPchan4 << " " << Xs[ii][0];
    for (G4int jj = 26; jj < 34; jj++)
    {
      PbarPchan4 << " " << Xs[ii][jj];
    };

    PbarPchan4 << G4endl;
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
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"

// Decay all input tracks, put daughters onto end of list

G4DecayKineticTracks::G4DecayKineticTracks(G4KineticTrackVector* tracks)
{
  if (tracks) Decay(tracks);
}

void G4DecayKineticTracks::Decay(G4KineticTrackVector* tracks) const
{
  if (!tracks) return;

  G4KineticTrackVector* daughters = 0;
  for (size_t i = 0; i < tracks->size(); ++i)
  {
    G4KineticTrack* track = (*tracks)[i];
    if (!track) continue;

    // Select decay of current track, put daughters at end of vector
    daughters = track->GetDefinition()->IsShortLived() ? track->Decay() : 0;
    //
    if ((!daughters)
        && ((track->GetDefinition()->GetPDGEncoding() == 221) ||  // Eta
            (track->GetDefinition()->GetPDGEncoding() == 331)))  // Eta_prime
      daughters = track->Decay();
    //
    if (daughters)
    {
      tracks->insert(tracks->end(), daughters->begin(), daughters->end());
      delete track;  // Remove parent track
      delete daughters;
      (*tracks)[i] = NULL;  // Flag parent's slot for removal
    }
  }

  // Find and remove null pointers created by decays above
  for (int j = tracks->size() - 1; j >= 0; --j)
  {
    if (NULL == (*tracks)[j]) tracks->erase(tracks->begin() + j);
  }
}
