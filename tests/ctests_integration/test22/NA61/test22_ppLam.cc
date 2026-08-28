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
//               FTF test: P+P interactions; Inclusive
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      File name:     Test30
//      Author:        V.Ivanchenko
//      Creation date: 12 March 2002
// -------------------------------------------------------------------
#include "G4ChipsComponentXS.hh"
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
  G4double SqrtS;

  //-------------------------- Global histograms  -------------------------

  G4double YUzhi[100];
  for (G4int ii = 0; ii < 100; ii++)
  {
    YUzhi[ii] = 0.;
  }

  G4double XfUzhi[100];
  for (G4int ii = 0; ii < 100; ii++)
  {
    XfUzhi[ii] = 0.;
  }

  //------------------------------------------------------------------
  G4double PtY[6][20];
  for (G4int ii = 0; ii < 6; ii++)
  {
    for (G4int jj = 0; jj < 20; jj++)
      PtY[ii][jj] = 0.;
  }

  G4double LlimitY[6] = {-1.75, -1.25, -0.75, -0.25, 0.25, 0.75};
  G4double UlimitY[6] = {-1.25, -0.75, -0.25, 0.25, 0.75, 1.25};
  //--------------------------
  G4double PtX[8][20];
  for (G4int ii = 0; ii < 8; ii++)
  {
    for (G4int jj = 0; jj < 20; jj++)
      PtX[ii][jj] = 0.;
  }

  G4double LlimitX[8] = {-0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3};
  G4double UlimitX[8] = {-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4};
  //--------------------------

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
      G4double E = energy + part->GetPDGMass();
      SqrtS = std::sqrt(sqr(part->GetPDGMass()) + sqr(938.) + 2. * 938. * E);
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-------------------------------------------------------

      const G4DynamicParticle* sec = 0;
      G4ParticleDefinition* pd;
      G4ThreeVector mom;
      G4LorentzVector labv, fm;
      G4double e, px, py, pt;  //, pt2, pz, theta;
      G4VParticleChange* aChange = 0;

      // -------- Event loop
      G4cout << "Events start " << nevt << G4endl;
      G4int Ninelast = nevt;  // Requested number of events
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

        // G4cout<<" Uzhi ------------ N prod. part "<<n<<G4endl;
        //++++++++++++++++ Variables for each event +++++++++++++++++++++++++++++
        if ((verbose > 0) && (n < 2))
        {
          G4cout << "Multiplicity of produced < 2!!!" << G4endl;
        }

        //    G4int Ncharged=0;                          // Uzhi
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        for (G4int i = 0; i < n; ++i)  // Loop over produced particles
        {
          sec = aChange->GetSecondary(i)->GetDynamicParticle();
          pd = sec->GetDefinition();
          G4String pname = pd->GetParticleName();

          // Particles in Lab. sys.
          if (verbose > 0)
            G4cout << " Part  " << i << " " << pname << " " << sec->Get4Momentum() / GeV
                   << sec->Get4Momentum().mag() / GeV << G4endl;

          labv -= fm;  // For checking energy-momentum conservation

          // electron can come only from internal conversion
          // its mass should be added to initial state
          if (pd == electron)
          {
            labv += G4LorentzVector(0.0, 0.0, 0.0, electron_mass_c2);
          }
          if (pname == "lambda")
          {  // +++++++++++++++++++++++++++++++++++++++++++++++
            fm = sec->Get4Momentum();
            fm.boost(-bst);

            //  Particles in CMS sys.

            mom = fm.vect();

            px = mom.x();
            py = mom.y();
            pt = std::sqrt(px * px + py * py) / GeV;
            e = fm.e();
            // Xxxxxxx
            G4double feynmanX = 2 * fm.z() / SqrtS;
            G4double PartWeight = 1.;

            G4int NxUzhi = int((feynmanX + 1.0) / 0.02);
            if (NxUzhi < 0.) NxUzhi = -1;
            if (NxUzhi > 99) NxUzhi = -1;
            // Xxxxxxx
            // Yyyyyyy
            G4double rapidity = fm.rapidity();
            G4int NyUzhi = int((rapidity + 5.) / 0.1);  // For CMS
            if (NyUzhi < 0) NyUzhi = -1;
            if (NyUzhi > 99) NyUzhi = -1;
            // Yyyyyyy

            G4int NPtUzhi = int(pt / 0.1);
            if (NPtUzhi < 0) NPtUzhi = -1;
            if (NPtUzhi > 19) NPtUzhi = -1;

            G4int IinterY = -1;
            for (G4int Iint = 0; Iint < 6; Iint++)
            {
              if ((LlimitY[Iint] <= rapidity) && (rapidity < UlimitY[Iint])) IinterY = Iint;
            }

            G4int IinterX = -1;
            for (G4int Iint = 0; Iint < 8; Iint++)
            {
              if ((LlimitX[Iint] <= feynmanX) && (feynmanX < UlimitX[Iint])) IinterX = Iint;
            }

            if (verbose >= 1)
              G4cout << " Part  " << i << " " << pname << " " << fm / GeV << fm.mag() / GeV << " "
                     << rapidity << " " << NyUzhi << G4endl;
            //

            //+++++++++++++++++ For each particle in the event ++++++++++++++++++++++

            if (NxUzhi >= 0) XfUzhi[NxUzhi] += PartWeight;
            if (NyUzhi >= 0) YUzhi[NyUzhi] += 1.;

            if ((NPtUzhi >= 0) && (IinterY >= 0)) PtY[IinterY][NPtUzhi]++;
            if ((NPtUzhi >= 0) && (IinterX >= 0)) PtX[IinterX][NPtUzhi]++;

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            de += e;

            //	delete sec;
            delete aChange->GetSecondary(i);
          };  // End of if ( pname == "lambda" ) { //
              // +++++++++++++++++++++++++++++++++++++++++++++++
        }  //     end of the loop on particles
        //+++++++++++++++++ Store after each event ++++++++++++++++++++++++++++++

        if (n == 2) Ninelast--;
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

      sigIn = 1.;  // Apr. 2016

      std::ofstream PP158y("pp158yLam.dat", std::ios::out);
      std::ofstream PP158x("pp158xLam.dat", std::ios::out);

      //----------------------------- Rapidity distributions------------------// Uzhi ++++
      G4cout << "******** Rapidity ******* at Plab " << Plab << " Xin " << sigIn << " "
             << SqrtS / 2000. << G4endl;
      PP158y << G4Version << G4endl;
      PP158y << "  Y    Lamb" << G4endl;

      for (G4int ii = 0; ii < 100; ii++)
      {
        YUzhi[ii] *= 10. / Ninelast * sigIn;
      };

      G4double Yu = -0.05 - 5.;  // For CMS                                      // June 2017
      for (G4int ii = 0; ii < 100; ii++)
      {
        Yu += 0.1;
        if (std::abs(Yu) < 3.05)
        {
          PP158y << Yu << " ";
          PP158y << YUzhi[ii] << " ";
          PP158y << G4endl;
        }
      }
      //

      // ----------------------------- Xf distributions----------------------// Uzhi ++++
      G4cout << "******** Xf distr ******* at Plab " << Plab << " Xin " << sigIn << G4endl;
      PP158x << G4Version << G4endl;
      PP158x << "  Xf Lamb " << G4endl;

      for (G4int ii = 0; ii < 100; ii++)
      {
        XfUzhi[ii] *= 1. / Ninelast * sigIn / 0.02;
      }

      G4double Xu = -1.01;
      for (G4int ii = 0; ii < 100; ii++)
      {
        Xu += 0.02;
        PP158x << Xu << " ";
        PP158x << XfUzhi[ii] << " ";
        PP158x << G4endl;
      }

      //----------------------------- Pt -Y distributions----------------------
      std::ofstream PP158pty("pp158ptyLam.dat", std::ios::out);

      G4cout << "******** Pt-Y ******* at Plab " << Plab << " Xin " << sigIn << G4endl;
      PP158pty << G4Version << G4endl;
      PP158pty << "  Pt Ym1p5  Ym1p0 Ym0p5 Y0 Y0p5 Y0p1" << G4endl;

      for (G4int ii = 0; ii < 6; ii++)
      {
        for (G4int jj = 0; jj < 20; jj++)
        {
          PtY[ii][jj] *= 1. / Ninelast * sigIn / 0.1 / 0.5;
        };
      };

      G4double Pt = -0.05;
      for (G4int ii = 0; ii < 20; ii++)
      {
        Pt += 0.1;
        PP158pty << Pt << " ";
        for (G4int jj = 0; jj < 6; jj++)
        {
          PP158pty << PtY[jj][ii] << " ";
        }
        PP158pty << G4endl;
      }

      //----------------------------- Pt -Xf distributions----------------------
      std::ofstream PP158ptx("pp158ptxLam.dat", std::ios::out);

      G4cout << "******** Pt-Xf ******* at Plab " << Plab << " Xin " << sigIn << G4endl;
      PP158ptx << G4Version << G4endl;
      PP158ptx << "  Pt Xm0p35  Xm0p25 Xm0p15 Xm0p05 X0p05 X0p15 X0p25 X0p35" << G4endl;

      for (G4int ii = 0; ii < 8; ii++)
      {
        for (G4int jj = 0; jj < 20; jj++)
        {
          PtX[ii][jj] *= 1. / Ninelast * sigIn / 0.1 / 0.1;
        };
      };

      Pt = -0.05;
      for (G4int ii = 0; ii < 20; ii++)
      {
        Pt += 0.1;
        PP158ptx << Pt << " ";
        for (G4int jj = 0; jj < 8; jj++)
        {
          PP158ptx << PtX[jj][ii] << " ";
        }
        PP158ptx << G4endl;
      }

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
