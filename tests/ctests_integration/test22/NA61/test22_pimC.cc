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
//               FTF test: Pi- + C --> Pi+/Pi-/K+/K-/Pro/AntiPro at 158 and 350 GeV/c; NA61/SHINE
//               data
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
  // G4double sigTot = 0;
  // G4double sigEl  = 0;
  G4double sigIn = 0;

  //-------------------------- Global histograms  -------------------------
  std::ofstream hAp("PimCna61.dat", std::ios::out);
  // std::ofstream hAy("hAy.dat",std::ios::out);
  //---------------------------------------------------------------------- pi- C 2 P/Pbar
  G4double protP[22] = {1.6,  2.2,  3.1,  3.6,  4.4,  5.5,  7.0,  9.0,  11.0, 14.0,  18.0,
                        22.0, 28.0, 36.0, 46.0, 56.0, 66.0, 76.0, 86.0, 96.0, 106.0, 116.0};
  G4double protDP[22] = {0.3, 0.3, 0.3, 0.4, 0.5, 0.8, 1.0, 1.5, 1.5, 2.0, 2.0,
                         2.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};

  G4double piP[22] = {0.6,  1.4,  2.2,  2.6,  3.0,  3.6,  4.4,  5.5,  7.0,  9.0,   11.0,
                      14.0, 18.0, 22.0, 28.0, 35.0, 45.0, 57.0, 72.0, 92.0, 102.0, 122.0};
  G4double piDP[22] = {0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.4, 0.5, 1.0,  1.0,  1.0,
                       2.0, 2.0, 2.0, 4.0, 5.0, 5.0, 6.0, 6.0, 10.0, 10.0, 10.0};

  G4double kP[22] = {1.0,  1.2,  1.8,  2.6,  3.2,  3.6,  4.3,  5.5,  7.0,  9.0,  11.0,
                     15.0, 18.0, 23.0, 28.0, 36.0, 46.0, 58.0, 72.0, 82.0, 92.0, 102.0};
  G4double kDP[22] = {0.2, 0.2, 0.4, 0.4, 0.2, 0.4, 0.3, 1.0, 1.0, 1.0, 1.5,
                      1.5, 1.5, 2.5, 4.0, 4.0, 6.0, 6.0, 7.0, 5.0, 5.0, 5.0};

  /*
  G4cout<<"Protons"<<G4endl;
  for(G4int i=0; i<22; i++) {G4cout<<protP[i]<<" "<<protDP[i]<<" "<<protP[i]-protDP[i]<<"
  "<<protP[i]+protDP[i]<<G4endl;}

  G4cout<<"Pions"<<G4endl;
  for(G4int i=0; i<22; i++) {G4cout<<piP[i]<<" "<<piDP[i]<<" "<<piP[i]-piDP[i]<<"
  "<<piP[i]+piDP[i]<<G4endl;}

  G4cout<<"Kaons"<<G4endl;
  for(G4int i=0; i<22; i++) {G4cout<<kP[i]<<" "<<kDP[i]<<" "<<kP[i]-kDP[i]<<"
  "<<kP[i]+kDP[i]<<G4endl;}
  */
  //-----------------------------------------------------------------------

  G4double XUzhi[22][7];  // -1 -- +1 x=2 Pz/sqrt(S)
  for (G4int ii = 0; ii < 22; ii++)
  {
    for (G4int j = 0; j < 7; j++)
      XUzhi[ii][j] = 0.;
  }

  // G4double YUzhi[100][7];                 // -5 -- 5 rapidity
  // for(G4int ii=0; ii<100; ii++){for(G4int j=0;j<7;j++) YUzhi[ii][j]=0.;}

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

      G4bool CHIPapplic = true;  // false;   Uzhi 29.01.13
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

        // sigTot=chipsTot;
        // sigEl=chipsEl;
        sigIn = chipsIn;
      }
      else
      {
        // sigTot = cross_sec;
        // sigEl  = cross_secel;
        sigIn = cross_inel;

        G4cout << "Proposed Xs (mb) are used: Tot El In: " << cross_sec << " " << cross_secel << " "
               << cross_inel << G4endl;
      }

      //+++++++++++++++++++++++++++++++++ For each energy +++++++++++++++++++++
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-------------------------------------------------------

      const G4DynamicParticle* sec = 0;
      G4ParticleDefinition* pd;
      G4ThreeVector mom;
      G4LorentzVector labv, fm;
      G4double e(0.);
      G4VParticleChange* aChange = 0;

      // -------- Event loop
      G4cout << "Events start " << nevt << G4endl;
      G4int Ninelast = nevt;  // Uzhi
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
        G4LorentzVector labNN(0.0, 0.0, std::sqrt(e0 * (e0 + 2. * mass)),
                              e0 + mass + 938.);  // amass);
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
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        for (G4int ii = 0; ii < n; ++ii)  // Loop over produced particles
        {
          sec = aChange->GetSecondary(ii)->GetDynamicParticle();
          pd = sec->GetDefinition();
          G4String pname = pd->GetParticleName();

          if (verbose >= 1)
            G4cout << " Part  " << ii << " " << pname << " " << sec->Get4Momentum() / GeV
                   << sec->Get4Momentum().mag() / GeV << G4endl;

          fm = sec->Get4Momentum();
          labv -= fm;  // For checking energy-momentum conservation

          // electron can come only from internal conversion
          // its mass should be added to initial state
          if (pd == electron)
          {
            labv += G4LorentzVector(0.0, 0.0, 0.0, electron_mass_c2);
          }

          mom = fm.vect();

          G4double Pmod = mom.mag() / GeV;
          //+++++++++++++++++ For each particle in the event ++++++++++++++++++++++

          if (pname == "proton")
          {
            for (G4int i = 0; i < 22; i++)
            {
              if (std::abs(Pmod - protP[i]) <= protDP[i]) XUzhi[i][0] += Pmod;
            }
          };

          if (pname == "neutron")
          {};

          if (pname == "pi+")
          {
            for (G4int i = 0; i < 22; i++)
            {
              if (std::abs(Pmod - piP[i]) <= piDP[i]) XUzhi[i][2] += Pmod;
            }
          };

          if (pname == "pi-")
          {
            for (G4int i = 0; i < 22; i++)
            {
              if (std::abs(Pmod - piP[i]) <= piDP[i]) XUzhi[i][3] += Pmod;
            }
          };

          if (pname == "kaon+")
          {
            for (G4int i = 0; i < 22; i++)
            {
              if (std::abs(Pmod - kP[i]) <= kDP[i]) XUzhi[i][4] += Pmod;
            }
          };

          if (pname == "kaon-")
          {
            for (G4int i = 0; i < 22; i++)
            {
              if (std::abs(Pmod - kP[i]) <= kDP[i]) XUzhi[i][5] += Pmod;
            }
          };

          if (pname == "anti_proton")
          {
            for (G4int i = 0; i < 22; i++)
            {
              if (std::abs(Pmod - protP[i]) <= protDP[i]) XUzhi[i][6] += Pmod;
            }
          };

          //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          de += e;

          //	delete sec;
          delete aChange->GetSecondary(ii);

        }  //     end of the loop on particles
        //+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++

        if (n == 0)
        {
          Ninelast--;
          G4cout << "n=0 !!! " << G4endl;
        }
        if (n == 1)
        {
          Ninelast--;
          G4cout << "n=1 !!! " << G4endl;
        }
        if (n == 2)
        {
          Ninelast--;
          G4cout << "n=2 !!! " << G4endl;
        }
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
      G4cout << "***********************************************************" << G4endl;

      G4cout << "nevt Ninel " << nevt << " " << Ninelast << G4endl;
      G4cout << "Plab " << Plab / GeV << " SigIn " << sigIn << G4endl;
      //-------------------------------------------------------------------

      if (verbose > 0)
      {
        G4cout << "###### End of run # " << run << "     ######" << G4endl;
      }
      //++++++++++++++++++++++ After each energy run ++++++++++++++++++++++++++ Uzhi

      // sigTot=sigTot; sigEl=sigEl;

      // ----------------------------- P distributions----------------------
      G4cout << "******** P distr ******* at Plab " << Plab << " Xin " << sigIn << G4endl;

      hAp << G4Version << G4endl;
      hAp << " Ppr  Prot Ppip Pip Ppim Pim Pkp Kp Pkm Km Ppbar Pbar" << G4endl;

      for (G4int ii = 0; ii < 22; ii++)
      {
        XUzhi[ii][0] *= 1. / Ninelast / 2. / protDP[ii];
        XUzhi[ii][6] *= 1. / Ninelast / 2. / protDP[ii];
      }

      for (G4int ii = 0; ii < 22; ii++)
      {
        XUzhi[ii][2] *= 1. / Ninelast / 2. / piDP[ii];
        XUzhi[ii][3] *= 1. / Ninelast / 2. / piDP[ii];
      }

      for (G4int ii = 0; ii < 22; ii++)
      {
        XUzhi[ii][4] *= 1. / Ninelast / 2. / kDP[ii];
        XUzhi[ii][5] *= 1. / Ninelast / 2. / kDP[ii];
      }

      for (G4int ii = 0; ii < 22; ii++)
      {
        hAp << protP[ii] << " " << XUzhi[ii][0] << " " << piP[ii] << " " << XUzhi[ii][2] << " "
            << piP[ii] << " " << XUzhi[ii][3] << " " << kP[ii] << " " << XUzhi[ii][4] << " "
            << kP[ii] << " " << XUzhi[ii][5] << " " << protP[ii] << " " << XUzhi[ii][6] << G4endl;
      }

      //----------------------------- Rapidity distributions------------------// Uzhi ++++
      /*
          G4cout<< "******** Rapidity ******* at Plab "<<Plab<<" Xin " << sigIn<< G4endl;
          hAy<<G4Version<<G4endl;
          hAy<<"  Ycms   P Neut Pip Pim Kp Km Pbar"<< G4endl;

          for(G4int ii=0; ii <100; ii++)
          {for(G4int jj=0;jj < 7; jj++){YUzhi[ii][jj]*= 1./Ninelast/0.1;}} // *sigIn/0.04;};};

          G4double Yu=-5.05;  // For NN cms
          for(G4int ii=0; ii <100; ii++)
          {
           Yu+=0.1;
           hAy<<Yu<<" ";
           for(G4int jj=0;jj < 7; jj++){hAy<<YUzhi[ii][jj]<<" ";}
           hAy<<G4endl;
          }
      */
      //
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi
      //    G4cerr << "###### End of run # " << run << "     ######" << G4endl;

    } while (end);
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
