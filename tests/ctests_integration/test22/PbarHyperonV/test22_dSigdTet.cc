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
#include "G4Version.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>
// #include <iomanip>
#include "G4ChipsComponentXS.hh"  // Uzhi 29.01.13
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "FTFtest1.icc"
#include "UZHI_diffraction.hh"
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <iostream>

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
  G4double sigTot = 0;  // Vova
  // G4double sigEl  = 0;  // Vova
  // G4double sigIn  = 0;  // Vova

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

      G4double CosDi[2][50];
      for (G4int i = 0; i < 50; i++)
      {
        CosDi[0][i] = 0.0;
        CosDi[1][i] = 0.0;
      }

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

        sigTot = chipsTot;  // sigEl=chipsEl; sigIn=chipsIn;    // Vova
      }
      else
      {
        sigTot = cross_sec;  // Vova
        //     sigEl  = cross_secel;   // Vova
        //     sigIn  = cross_inel;    // Vova

        G4cout << "Proposed Xs (mb) are used: Tot El In: " << cross_sec << " " << cross_secel << " "
               << cross_inel << G4endl;
      }

      //+++++++++++++++++++++++++++++++++ For each energy +++++++++++++++++++++
      G4int Ntotal = nevt;
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-------------------------------------------------------

      const G4DynamicParticle* sec = 0;
      G4ParticleDefinition* pd;
      // Vova    G4ThreeVector  mom;
      G4LorentzVector labv, fm;
      G4double e, theta;
      G4VParticleChange* aChange = 0;

      //   G4double E=energy+part->GetPDGMass();                         // Elab Proj
      //   G4double SS=sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E;      // per  Proj+N
      //  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                 //      Proj+N

      // -------- Event loop
      G4cout << "Events start " << nevt << G4endl;

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

        G4LorentzVector Proj4Mom =
          G4LorentzVector(0.0, 0.0, std::sqrt(e0 * (e0 + 2. * mass)), e0 + mass);
        //    G4cout<<Proj4Mom<<G4endl;
        Proj4Mom.boost(-boostNN);
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

        int id;

        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        G4int AntiLam = 0;
        G4int Lam = 0;
        G4int SigmaB = 0;
        G4int Sigma = 0;
        for (G4int i = 0; i < n; ++i)  // Loop over produced particles
        {
          sec = aChange->GetSecondary(i)->GetDynamicParticle();
          pd = sec->GetDefinition();
          if (pd == electron)
          {};  // Vova  to erase warning message
          G4String pname = pd->GetParticleName();

          if (verbose >= 1)
            G4cout << " Part  " << i << " " << pname << " " << sec->Get4Momentum() / GeV
                   << sec->Get4Momentum().mag() / GeV << G4endl;

          id = pd->GetPDGEncoding();

          if (id == 3122) Lam++;
          ;
          if (id == -3122) AntiLam++;
          if (id == -3212) SigmaB++;
          if (id == 3212) Sigma++;
        }  //     end of the loop on particles

        //==============================
        if ((n == 2) && (AntiLam == 1) && (Lam == 1))  // Lambda - AntiLambda
        {
          for (G4int i = 0; i < n; ++i)  // Loop over produced particles
          {
            sec = aChange->GetSecondary(i)->GetDynamicParticle();
            pd = sec->GetDefinition();
            G4String pname = pd->GetParticleName();

            fm = sec->Get4Momentum();
            labv -= fm;  // For checking energy-momentum conservation

            id = pd->GetPDGEncoding();

            if (id == -3122)
            {
              fm.boost(-bst);  // Transformation to CM system
              theta = fm.theta();
              G4double CosTheta = std::cos(theta);
              G4int Ihist = G4int((CosTheta + 1.0) / 0.04);
              if ((Ihist >= 0) && (Ihist <= 49)) CosDi[0][Ihist]++;
            }
          }  //     end of the loop on particles
        }  //     end of if((n == 2) && (AntiLam == 1) && (Lam == 1))

        if ((n == 2)
            && (((Lam == 1) && (SigmaB == 1))
                || ((AntiLam == 1) && (Sigma == 1))))  // Lambda SigmaB+c.c.
        {
          for (G4int i = 0; i < n; ++i)  // Loop over produced particles
          {
            sec = aChange->GetSecondary(i)->GetDynamicParticle();
            pd = sec->GetDefinition();
            G4String pname = pd->GetParticleName();

            fm = sec->Get4Momentum();
            e = fm.e() / GeV;  // - m;

            id = pd->GetPDGEncoding();

            if (id == -3212)
            {
              fm.boost(-bst);  // Transformation to CM system
              theta = fm.theta();
              G4double CosTheta = std::cos(theta);
              G4int Ihist = G4int((CosTheta + 1.0) / 0.04);
              if ((Ihist >= 0) && (Ihist <= 49)) CosDi[1][Ihist]++;
            }
          }  //     end of the loop on particles
        }  //     end of if((n == 2) && (Lam == 1) && (SigmaB == 1))
        de += e;

        for (G4int i = 0; i < n; ++i)  // Loop over produced particles
        {
          sec = aChange->GetSecondary(i)->GetDynamicParticle();
          fm = sec->Get4Momentum();
          e = fm.e() / GeV;  // - m;	//	delete sec;
          delete aChange->GetSecondary(i);
        }  //     end of the loop on particles
        //     }
        //==============================

        if (verbose > 0) G4cout << "Energy/Momentum balance= " << labv << G4endl;

        aChange->Clear();

        if (verbose > 0)
        {
          G4cout << "End event =====================================" << Plab << G4endl;  // Uzhi
          G4int Uzhi_i;  // Uzhi
          G4cin >> Uzhi_i;  // Uzhi
        }

      }  // End of the event loop ------------------------------------
      /*                  // Vova
          h1->Write();
          h2->Write();
      */
      timer->Stop();
      G4cout << "  " << *timer << G4endl;
      delete timer;

      if (verbose > 0)
      {
        G4cout << "###### End of run # " << run << "     ######" << G4endl;
      }
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi

      G4cout << G4Version << G4endl;

      // Vova
      // std::ofstream CosLb1_476("CosDi.dat",std::ios::out);
      std::ofstream CosLb1_476(OutFile, std::ios::out);
      CosLb1_476 << G4Version << G4endl;
      CosLb1_476 << "CosThet  LamBar Sigma" << G4endl;

      G4double X_b = sigTot * 1000.0;  // X_b=15141.; // MicroBarn
      G4double SumLLbar(0.0), SumLbarS(0.0);
      G4double dOmega = 2 * pi * 0.04;

      G4double CosC = -1.02;
      for (G4int i = 0; i < 50; i++)
      {
        CosC += 0.04;
        SumLLbar += CosDi[0][i];
        SumLbarS += CosDi[1][i];
        CosDi[0][i] *= X_b / dOmega / Ntotal;
        CosDi[1][i] *= X_b / dOmega / Ntotal;
        CosLb1_476 << CosC << " " << CosDi[0][i] << " " << CosDi[1][i] << G4endl;
      }
      G4cout << "N LLbar Lbar S " << SumLLbar << " " << SumLbarS << G4endl;
      // Vova

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
  // Vova  f1.Write();

  G4cout << "###### End of test #####" << G4endl;
}
