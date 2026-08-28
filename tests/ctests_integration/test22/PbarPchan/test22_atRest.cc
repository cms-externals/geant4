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
  // G4double sigTot = 0;
  // G4double sigEl  = 0;
  G4double sigIn = 0;

  //-------------------------- Global histograms  -------------------------

  G4double Xs[50];
  for (G4int ii = 0; ii < 50; ii++)
  {
    Xs[ii] = 0.;
  }

  G4double MultDistr[14];
  for (G4int ii = 0; ii < 14; ii++)
  {
    MultDistr[ii] = 0.;
  }

  G4double MomDistr[40];
  for (G4int ii = 0; ii < 40; ii++)
  {
    MomDistr[ii] = 0.;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;

  G4int Ninelast{0};

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

        // sigTot=chipsTot; sigEl=chipsEl;
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
      G4int Ntotal = nevt;
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //-------------------------------------------------------

      const G4DynamicParticle* sec = 0;
      G4ParticleDefinition* pd;
      G4ThreeVector mom;
      G4LorentzVector labv, fm;
      G4double e, px, py, pz, pt, theta;
      G4VParticleChange* aChange = 0;

      //  G4double E=energy+part->GetPDGMass();                                  // Elab Proj
      //  G4double SqrtS=std::sqrt(sqr(part->GetPDGMass())+sqr(938.)+2.*938.*E); // per  Proj+N
      //  G4double Ycms=0.5*std::log((E+Plab)/(E-Plab));                         //      Proj+N

      // -------- Event loop
      G4cout << "Events start " << nevt << G4endl;
      //    G4int
      Ninelast = nevt;  // Uzhi
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

          labv -= fm;  // For checking energy-momentum conservation

          // electron can come only from internal conversion
          // its mass should be added to initial state
          if (pd == electron)
          {
            labv += G4LorentzVector(0.0, 0.0, 0.0, electron_mass_c2);
          }

          px = mom.x();
          py = mom.y();
          pz = mom.z();
          pt = std::sqrt(px * px + py * py);
          e = fm.e() - m;
          theta = mom.theta();

          //        G4double CosTheta=std::cos(theta);

          theta = theta * 180. / pi;

          G4int Imom = 0;
          G4double p = std::sqrt(sqr(pt) + sqr(pz));
          if (pname == "pi-")
          {
            Imom = G4int(p / 25.);
          }
          else if (pname == "pi+")
          {
            Imom = G4int(p / 25.);
          }
          else
          {}

          if ((Imom != 0) && (Imom < 40))
          {
            MomDistr[Imom]++;
          }

          fm.boost(-bst);

          //        G4double costcm = std::cos(fm.theta());
          //+++++++++++++++++ For each particle in the event ++++++++++++++++++++++
          if (pname == "pi-")
          {
            Npim++;
          }
          else if (pname == "pi+")
          {
            Npip++;
          }
          else if (pname == "pi0")
          {
            Npi0++;
          }

          else if (pname == "kaon-")
          {
            NKm++;
          }
          else if (pname == "kaon+")
          {
            NKp++;
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
          }
          else if (pname == "sigma-")
          {
            NSigma_m++;
          }
          else if (pname == "anti_sigma+")
          {
            NSigma_pBar++;
          }
          else if (pname == "anti_sigma-")
          {
            NSigma_mBar++;
          }

          else if (pname == "xi0")
            NXi0++;
          else if (pname == "xi-")
          {
            NXi_m++;
          }
          else if (pname == "anti_xi0")
            NXi0Bar++;
          else if (pname == "anti_xi-")
          {
            NXi_mBar++;
          }

          else if (pname == "omega-")
          {}
          else if (pname == "anti_omega-")
          {}

          else if (pname == "proton")
          {
            Nproton++;
          }
          else if (pname == "anti_proton")
          {
            NprotonBar++;
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
        if ((Npip + Npim + Npi0 + NKm + NKp < 14)
            && (Nproton + Nneutron + NprotonBar + NneutronBar == 0))
        {
          MultDistr[Npip + Npim + Npi0 + NKm + NKp]++;
        }

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
        */
        G4int Nother = Nneutron + NneutronBar + Npim + Npip + Npi0 + NKm + NKp + NK0s + NK0l
                       + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                       + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta
                       + Neta_prime + Ngamma;
        if ((NprotonBar == 1) && (Nproton == 1) && (Nother == 0))
        {
          Xs[30]++;
          Ninelast--;
        }  // G4cout<<" 0 "<<G4endl;}

        //--------1 2 Pi+ Pi-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;
        if ((Npip == 1) && (Npim == 1) && (Nother == 0))
        {
          Xs[1]++;
        }  //  G4cout<<" 1 "<<G4endl;}

        //--------2 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 1) && (Neta == 1) && (Nother == 0))
        {
          Xs[2]++;
        }  //  G4cout<<" 2 "<<G4endl;}

        //---------3 2 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npi0 == 2) && (Nother == 0))
        {
          Xs[3]++;
        }  // G4cout<<" 3 "<<G4endl;}

        //----------4 2 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + Npi0 + NKm + NKp
                 + NK0s + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p
                 + NSigma_pBar + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar
                 + Neta_prime + Ngamma;

        if ((Neta == 2) && (Nother == 0))
        {
          Xs[4]++;
        }  // G4cout<<" 4 "<<G4endl;}

        //-------------5 Pi+ Pi- Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[5]++;
        }  // G4cout<<" 5 "<<G4endl;}

        //-----------6 Pi+ Pi- Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Neta == 1) && (Nother == 0))
        {
          Xs[6]++;
        }  // G4cout<<" 6 "<<G4endl;}

        //-----------7 3 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npi0 == 3) && (Nother == 0))
        {
          Xs[7]++;
        }  // G4cout<<" 7 "<<G4endl;}

        //-----------8 2 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 2) && (Neta == 1) && (Nother == 0))
        {
          Xs[8]++;
        }  // G4cout<<" 8 "<<G4endl;}

        //-----------9 Pi0 2 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 1) && (Neta == 2) && (Nother == 0))
        {
          Xs[9]++;
        }  // G4cout<<" 9 "<<G4endl;}

        //-----------10 Pi+ Pi- 2 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 2) && (Nother == 0))
        {
          Xs[10]++;
        }  // G4cout<<" 10 "<<G4endl;}

        //-----------11 2 Pi+ 2 Pi-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Nother == 0))
        {
          Xs[11]++;
        }  // G4cout<<" 11 "<<G4endl;}

        //-----------12 Pi+ Pi- Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 1) && (Neta == 1) && (Nother == 0))
        {
          Xs[12]++;
        }  // G4cout<<" 12 "<<G4endl;}

        //-----------13 4 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npi0 == 4) && (Nother == 0))
        {
          Xs[13]++;
        }  // G4cout<<" 13 "<<G4endl;}

        //-----------14 3 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 2) && (Neta == 1) && (Nother == 0))  // !!!! *****************
        {
          Xs[14]++;
        }  // G4cout<<" 14 "<<G4endl;} //  G4int Uzhi; G4cin>>Uzhi;

        //----------15  Pi+ Pi- 2 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Neta == 2) && (Nother == 0))
        {
          Xs[15]++;
        }  // G4cout<<" 15 "<<G4endl;} //  G4int Uzhi; G4cin>>Uzhi;

        //----------16 2 Pi0 2 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 2) && (Neta == 2) && (Nother == 0))
        {
          Xs[16]++;
        }  // G4cout<<" 16 "<<G4endl;}

        //-----------17 2 Pi+ 2 Pi- Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Npi0 == 1) && (Nother == 0))
        {
          Xs[17]++;
        }  // G4cout<<" 17 "<<G4endl;}

        //-----------18 Pi+ Pi- 3Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 3) && (Nother == 0))
        {
          Xs[18]++;
        }  // G4cout<<" 18 "<<G4endl;}

        //-----------19 Pi+ Pi- 2 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 2) && (Neta == 1) && (Nother == 0))
        {
          Xs[19]++;
        }  // G4cout<<" 19 "<<G4endl;}

        //-----------20 2 Pi+ 2 Pi- Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Neta == 1) && (Nother == 0))
        {
          Xs[20]++;
        }  // G4cout<<" 20 "<<G4endl;}

        //-----------21 5 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npi0 == 5) && (Nother == 0))
        {
          Xs[21]++;
        }  // G4cout<<" 18 "<<G4endl;}

        //-----------22 4 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 4) && (Neta == 1) && (Nother == 0))
        {
          Xs[22]++;
        }  // G4cout<<" 22 "<<G4endl;}

        //-----------23 2 Pi+ 2 Pi- 2 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Npi0 == 2) && (Nother == 0))
        {
          Xs[23]++;
        }  // G4cout<<" 23 "<<G4endl;}

        //-----------24 3 Pi+ 3 Pi-
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npi0 + NKm + NKp + NK0s + NK0l
                 + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 3) && (Npim == 3) && (Nother == 0))
        {
          Xs[24]++;
        }  // G4cout<<" 24 "<<G4endl;}

        //-----------25 Pi+ Pi- 4 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 4) && (Nother == 0))
        {
          Xs[25]++;
        }  // G4cout<<" 25 "<<G4endl;}

        //-----------26 2 Pi+ 2 Pi- Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 2) && (Npim == 2) && (Npi0 == 1) && (Neta == 1) && (Nother == 0))
        {
          Xs[26]++;
        }  // G4cout<<" 26 "<<G4endl;}

        //-----------27 Pi+ Pi- 3 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + NKm + NKp + NK0s + NK0l + NLambda
                 + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar + NSigma_m
                 + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npip == 1) && (Npim == 1) && (Npi0 == 3) && (Neta == 1) && (Nother == 0))
        {
          Xs[27]++;
        }  // G4cout<<" 27 "<<G4endl;}

        //-----------28 6 Pi0
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta + Neta_prime
                 + Ngamma;

        if ((Npi0 == 6) && (Nother == 0))
        {
          Xs[28]++;
        }  // G4cout<<" 28 "<<G4endl;}

        //-----------29 5 Pi0 Eta
        Nother = Nproton + Nneutron + NprotonBar + NneutronBar + Npim + Npip + NKm + NKp + NK0s
                 + NK0l + NLambda + NLambdaBar + NSigma0 + NSigma0Bar + NSigma_p + NSigma_pBar
                 + NSigma_m + NSigma_mBar + NXi0 + NXi0Bar + NXi_m + NXi_mBar + Neta_prime + Ngamma;

        if ((Npi0 == 5) && (Neta == 1) && (Nother == 0))
        {
          Xs[29]++;
        }  // G4cout<<" 29 "<<G4endl;}

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
      // sigTot=1.;
      G4cout << "***********************************************************" << G4endl;

      G4cout << "nevt Ninel " << nevt << " " << Ntotal << G4endl;
      G4cout << "Plab " << Plab / GeV << " SigIn " << sigIn << G4endl;
      //-------------------------------------------------------------------

      if (verbose > 0)
      {
        G4cout << "###### End of run # " << run << "     ######" << G4endl;
      }

    } while (end);
  }  // End of job ----------------------------------------------
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //----------------------------- Write distributions------------------
  std::ofstream Branch("Branch.dat", std::ios::out);
  std::ofstream PionMult("PionMult.dat", std::ios::out);
  std::ofstream MomDi("MomDistr.dat", std::ios::out);

  //++++++++++++++++++++++++++++++++++++++++++++
  G4cout << "nevt Ninel " << nevt << " " << Ninelast << G4endl;
  Branch << G4Version << G4endl;
  Branch << "Nchan  BraRat nevt Ninel " << nevt << " " << Ninelast << G4endl;

  G4double XchanSum(0.);

  for (G4int ii = 1; ii < 30; ii++)
  {
    Xs[ii] /= Ninelast;
    XchanSum += Xs[ii];
    if (Xs[ii] == 0.)
    {
      Xs[ii] = 1.0e-5;
    }
    Branch << ii - 0.5 << " " << Xs[ii] << G4endl;
    Branch << ii + 0.5 << " " << Xs[ii] << G4endl;
  }

  Xs[30] /= nevt;  // Elastic
  Branch << 30.5 << " " << 1.0e-5 << G4endl;
  Branch << 30.5 << " " << Xs[30] << G4endl;
  Branch << +31.5 << " " << Xs[30] << G4endl;
  Branch << +31.5 << " " << 1.0e-5 << G4endl;

  G4cout << G4endl << "Sum Br " << XchanSum << G4endl;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi
  G4cout << G4endl << "Pion multiplicity distrivution" << G4endl;
  PionMult << G4Version << G4endl;
  PionMult << "Mult Prob " << G4endl;
  G4double AverMult(0.);
  for (G4int ii = 0; ii < 14; ii++)
  {
    MultDistr[ii] /= Ninelast;
    AverMult += ii * MultDistr[ii];
    PionMult << ii << "   " << MultDistr[ii] << G4endl;
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Uzhi
  G4cout << G4endl << "Aver Mult " << AverMult << G4endl;

  G4cout << G4endl << "Momentum distribution" << G4endl;
  MomDi << G4Version << G4endl;
  MomDi << "Pmom   Distr " << G4endl;
  for (G4int ii = 0; ii < 40; ii++)
  {
    MomDistr[ii] /= (Ninelast * 0.025);
    MomDi << 0.025 * (ii + 0.5) << "  " << MomDistr[ii] << G4endl;
  }

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
    // if((!daughters) && ((track->GetDefinition()->GetPDGEncoding() == 221) ||    // Eta
    //                    (track->GetDefinition()->GetPDGEncoding() == 331)   ))  // Eta_prime
    // daughters = track->Decay();
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
