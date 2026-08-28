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

#include "Test19ExecProcessLevel.hh"

#include "G4CascadeInterface.hh"
#include "G4INCLXXInterface.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "FTFPWrapper.hh"
#include "QGSPWrapper.hh"
#include "TstDiscreteProcessReader.hh"
// bertini
#include "G4CascadeParameters.hh"
// preco
#include "G4DeexPrecoParameters.hh"
#include "G4FTFTunings.hh"
#include "G4HadronicDeveloperParameters.hh"
#include "G4NuclearLevelData.hh"

#ifdef G4_USE_FLUKA
// interface to FLUKA.CERN, if specified
#  include "FLUKAParticleTable.hh"
#  include "fluka_interface.hh"
// --> no need here --> #include "FLUKAInelasticScatteringXS.hh"
#  include "FLUKANuclearInelasticModel.hh"
#endif

void Test19ExecProcessLevel::InitProcess(const TstReader* pset)
{
  // G4String name = (dynamic_cast<const TstDiscreteProcessReader*>(pset))->GetProcessName();

  G4String name = (pset)->GetPhysics();
  G4String name_lc = G4StrUtil::to_lower_copy(name);

  ProcessWrapper* pw = 0;

  if (name.find("ftf") != std::string::npos)
  {
    if (name.find("_tune") != std::string::npos)
    {
      // example(s) of setting FTF parameters as tunes
      // (i.e. as groups where selected parameters have been obtained
      //  via collective fits vs thin target data)
      //
      // G4cout << " Is it locked: ? " << G4FTFTunings::Instance()->IsLocked() << std::endl;

      bool isMaster = G4Threading::IsMasterThread();
      G4cout << " Is Master ? " << isMaster << G4endl;

      bool isState_PreInit =
        (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit);
      G4cout << " Is state PreInit ? " << isState_PreInit << G4endl;

      G4ApplicationState currentstate = G4StateManager::GetStateManager()->GetCurrentState();
      bool ok = G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);
      G4cout << " Is PreInit set ? " << ok << G4endl;

      isState_PreInit = (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit);

      bool isLocked = (!isMaster || !isState_PreInit);
      G4cout << " Is locked ?  " << isLocked << G4endl;

      // -->      G4FTFTunings::Instance()->SetTuneApplicabilityState( 1,      1 );
      // --> below is attempt to set "combined tune" (both baryon-2022 and pion-2022)
      // --> int tuneID = 3;
      int tuneID = atoi(&(name.back()));
      G4cout << " tuneID = " << tuneID << G4endl;
      G4FTFTunings::Instance()->SetTuneApplicabilityState(tuneID, 1);

      // NOTE (JVY): For some reason when trying to set tune via UI
      //             the parameter values do not seem to propagate;
      //             maybe it should be at a different stage/state ?
      //
      // -->   G4UImanager* uim = G4UImanager::GetUIpointer();
      // -->   uim->ApplyCommand( "/process/had/models/ftf/selectTuneByIndex 1" );

      ok = G4StateManager::GetStateManager()->SetNewState(currentstate);
    }

    // example(s) of setting FTF parameters one by one

    // G4HadronicDeveloperParameters& HDP = G4HadronicDeveloperParameters::GetInstance();

    // baryon projectile
    //
    // setting a1 and a2 for *both* processes will turn off
    // quark exchange with AND w/o excitation
    // thus only non-diff should remain...
    //
    // setting a1=1 and b1=0 and a2=0 for one or another process
    // whill tuns the process on while while the other will be off
    // this will also turn off non-diff completely because
    // P(ndiff0=1.-P0-P1 (-P2-P3 but bother P2 and P3 are off in 11.1)
    //
    /*
          double br_pr0_a1 = 0.;
          HDP.Set( "FTF_BARYON_PROC0_A1", br_pr0_a1 );
          double br_pr0_b1 = 0.;
          HDP.Set( "FTF_BARYON_PROC0_B1", br_pr0_b1 );
          double br_pr0_a2 = 0.;
          HDP.Set( "FTF_BARYON_PROC0_A2", br_pr0_a2 );
          double br_pr0_a3 = 0.;
          HDP.Set( "FTF_BARYON_PROC0_A3", br_pr0_a3 );
          double br_pr0_ymin = 0.;
          HDP.Set("FTF_BARYON_PROC0_YMIN", br_pr0_ymin );
          double br_pr1_a1 = 0.;
          HDP.Set( "FTF_BARYON_PROC1_A1", br_pr1_a1 );
          double br_pr1_b1 = 0.;
          HDP.Set( "FTF_BARYON_PROC1_B1", br_pr1_b1 );
          double br_pr1_a2 = 0.;
          HDP.Set( "FTF_BARYON_PROC1_A2", br_pr1_a2 );
          double br_pr1_a3 = 0.;
          HDP.Set( "FTF_BARYON_PROC1_A3", br_pr1_a3 );
          double br_pr1_ymin = 0.;
          HDP.Set("FTF_BARYON_PROC1_YMIN", br_pr1_ymin );
          double br_pr4_a1 = 1.;
          HDP.Set( "FTF_BARYON_PROC4_A1", br_pr4_a1 );
          double br_pr4_b1 = 0.;
          HDP.Set( "FTF_BARYON_PROC4_B1", br_pr4_b1 );
          double br_pr4_a2 = 0.;
          HDP.Set( "FTF_BARYON_PROC4_A2", br_pr4_a2 );
          double br_pr4_a3 = 0.;
          HDP.Set( "FTF_BARYON_PROC4_A3", br_pr4_a3 );
          double br_pr4_ymin = 0.;
          HDP.Set("FTF_BARYON_PROC4_YMIN", br_pr4_ymin );
    */
    /*
          // NO big deal for either 5GeV p+Ta or 31GeV p+C
          // if setting these two ("PROJ") to 2.9 which is almost max
          //
          double br_diff_m_pr = 2.9;
    //      HDP.Set( "FTF_BARYON_DIFF_M_PROJ", br_diff_m_pr ); // 1.16-3.0
          double br_ndiff_m_pr = 2.9;
    //      HDP.Set( "FTF_BARYON_NONDIFF_M_PROJ", br_ndiff_m_pr );  // 1.16-3.0

          // Slight change in the right direction, as checked with LA-1 and LA-8,
          // at least for 5GeV p+Ta; but not really for 31GeV p+C
          //
          double br_diff_m_tg = 2.9;
    //      HDP.Set( "FTF_BARYON_DIFF_M_TGT", br_diff_m_tg ); // 1.16-3.0

          // Largely pulls down the "LA" spectra, as least for for 5Gev p+Ta,
          // with the smallest effect being on piminus_LA_8;
          // but not much for 31GeV p+C, except pi+ 300-360 moving a bit up,
          // pi- 0-10 moving somewhat up at 0-10 ("at peak) and pi- 300-360
          // moving slightly up at ~1.5GeV (mid-tail)
          //
          double br_ndiff_m_tg = 2.9;
    //      HDP.Set( "FTF_BARYON_NONDIFF_M_TGT", br_ndiff_m_tg ); // 1.16-3.0

          // At 0.1 there's a mild effect for 5GeV p+Ta but it seems to be
          // slightly re-shaping the pi+ LA_8 spectrum the right way, and
          // the pi- one as well, very mildly
          // Very mild effect for 31GeV p+C but not exactly the right way,
          // at least for pi- at 0-10
          // BTW, at the time of ~10.4 --> D=0.15
          //
          // At 0.95 for the 31GeV p+C it very slightly pulls down (right way)
          // both pi+ & pi- at 0-10 & 300-360
          // For 5GeV p+Ta slight pull down (right) for pi+ LA_1, however
          // for pi+ LA_8 the re-shaping effect of 0.1 seems to be slightly better
          // for pi- slight pull down in both cases which is not ideal in BCK
          //
          double br_avg_pt2 = 0.95; // range : 0.08 - 1.0, D=0.3
    //      HDP.Set( "FTF_BARYON_AVRG_PT2", br_avg_pt2);
    */
    // pion projectile
    //
    // setting a1 and a2 to zero for *both* processes will turn off
    // quark exchange with AND w/o excitation
    // thus only non-diff should remain (unless one turns ON diffraction)
    //
    // setting a1=1 and b1=0 and a2=0 for one or another process
    // whill tuns the process on while while the other will be off
    // this will also turn off non-diff completely because
    // P(ndiff0=1.-P0-P1 (-P2-P3 but both P2 and P3 are off in 11.1)
    //
    /*
          double pi_pr0_a1 = 0.;
          HDP.Set( "FTF_PION_PROC0_A1", pi_pr0_a1 );
          double pi_pr0_b1 = 0.;
          HDP.Set( "FTF_PION_PROC0_B1", pi_pr0_b1 );
          double pi_pr0_a2 = 0.;
          HDP.Set( "FTF_PION_PROC0_A2", pi_pr0_a2 );
          double pi_pr0_a3 = 0.;
          HDP.Set( "FTF_PION_PROC0_A3", pi_pr0_a3 );
          double pi_pr0_ymin = 0.;
          HDP.Set( "FTF_PION_PROC0_YMIN", pi_pr0_ymin );
    */
    /*
          double pi_pr1_a1 = 0.;
          HDP.Set( "FTF_PION_PROC1_A1", pi_pr1_a1 );
          double pi_pr1_b1 = 0.;
          HDP.Set( "FTF_PION_PROC1_B1", pi_pr1_b1 );
          double pi_pr1_a2 = 0.;
          HDP.Set( "FTF_PION_PROC1_A2", pi_pr1_a2 );
          double pi_pr1_b2 = 0.;
          HDP.Set( "FTF_PION_PROC1_B2", pi_pr1_b2 );
          double pi_pr1_a3 = 0.;
          HDP.Set( "FTF_PION_PROC1_A3", pi_pr1_a3 );
          double pi_pr1_ymin = 0.;
          HDP.Set( "FTF_PION_PROC1_YMIN", pi_pr1_ymin );
    */

    /*
          double pi_pr3_a1 = 1.;
          HDP.Set( "FTF_PION_PROC3_A1", pi_pr3_a1 );
          double pi_pr3_b1 = 0.;
          HDP.Set( "FTF_PION_PROC3_B1", pi_pr3_b1 );
          double pi_pr3_a2 = 0.;
          HDP.Set( "FTF_PION_PROC3_A2", pi_pr3_a2 );
          double pi_pr3_a3 = 0.;
          HDP.Set( "FTF_PION_PROC3_A3", pi_pr3_a3 );
          double pi_pr3_ymin = 0.;
          HDP.Set( "FTF_PION_PROC3_YMIN", pi_pr3_ymin );
    */

    /*
          double pi_pr4_a1 = 1.;
          HDP.Set( "FTF_PION_PROC4_A1", pi_pr4_a1 );
          double pi_pr4_b1 = 0.;
          HDP.Set( "FTF_PION_PROC4_B1", pi_pr4_b1 );
          double pi_pr4_a2 = 0.;
          HDP.Set( "FTF_PION_PROC4_A2", pi_pr4_a2 );
          double pi_pr4_a3 = 0.;
          HDP.Set( "FTF_PION_PROC4_A3", pi_pr4_a3 );
          double pi_pr4_ymin = 0.;
          HDP.Set( "FTF_PION_PROC4_YMIN", pi_pr4_ymin );
    */
    /* -->
          // Overall, the "PROJ" ones result in very small variations
          // or NO variations at all (apparently some limits implemented ?)
          // Q : Was the 2.9 upper limit recommended by A.R./V.Uzh. ?
          //     Or wha it an "educated guess" by JY ???
          //
          // At 0.6 slightly pull up, but at 2.9 it's *identical* to default
          // - perhaps maxes out somewhere ?
          //
          double pi_diff_m_pr = 2.9;
    //      HDP.Set( "FTF_PION_DIFF_M_PROJ", pi_diff_m_pr ); // 0.5-3.0, D=1.

          // At 0.6 slightly pulls down at LA_1, at 2.9 identical to default
          //
          double pi_ndiff_m_pr = 2.9;
    //      HDP.Set( "FTF_PION_NONDIFF_M_PROJ", pi_ndiff_m_pr ); // 0.5-3.0, D=1.
    --> */
    /*
          // Very mild fluctuations for 5GeV p+ta
          //
          double pi_diff_m_tg = 2.9;
    //      HDP.Set( "FTF_PION_DIFF_M_TGT", pi_diff_m_tg ); // 1.16-3.0

          // For 5GeV p+Ta largely pills down LA_1 & LA_8 for pi+,
          // also LA_1 for pi-, but only mildly LA_8 for pi-
          //
          double pi_ndiff_m_tg = 2.9;
    //      HDP.Set( "FTF_PION_NONDIFF_M_TGT", pi_ndiff_m_tg ); // 1.16-3.0

          // Effects similar to those for proton beam
          //
          double pi_avg_pt2 = 0.95; // range : 0.08 - 1.0, D=0.3
    //      HDP.Set( "FTF_PION_AVRG_PT2", pi_avg_pt2);

          // NOTE: The following, for both baryon or pion projectile, can NOT change at run-time,
          //       NO API as of now (per A.R. back in 2017)
          //       SetProbLogDistrPrD( 0.3 );
          //       SetProbLogDistr( 0.3 );
          //       From the doc on the matter:
          //       Parameters ... are very important for a correct description of produced
          //       particles multiplicity in non-diffractive interactions.
          //       ...
          //       A variation of the parameters can have an essential influence on multiplicity
          //       of produced particles.
    */
    /*
          // projectile and target diffraction
          //
          // the following probabilities:
          // Probability of quark exchange process without excitation of participants (Fig. 44(b));
       (Proc# 0)
          // Probability of quark exchange process with excitation of participants (Fig. 44(c));
       (Proc# 1)
          // Probability of projectile diffraction dissociation; (Proc# 2)
          // Probability of target diffraction dissociation. (Proc# 3)
          // are described with the formila:
          // P = A1 * exp(-B1*y) + A2 * exp(-B2*y) + A3
          // where y is the rapidity of the projectile in the target rest frame
          //
          // example of setting "proc1" for pion projectile
          // (similar can be done for baryon proj)
          //
          double br_pr1_a1 = 0.;
          HDP.Get( "FTF_PION_PROC1_A1", br_pr1_a1 ); // D=5.77
          std::cout << " br_pr1_a1 = " << br_pr1_a1 << std::endl;
          br_pr1_a1 = 5.836; // 7.05; // 7.43; // 6.55; // 7.87;
          HDP.Set( "FTF_PION_PROC1_A1", br_pr1_a1 );
          //
          double br_pr1_b1 = 0.;
          HDP.Get( "FTF_PION_PROC1_B1", br_pr1_b1 ); // D=0.6
          std::cout << " br_pr1_b1 = " << br_pr1_b1 << std::endl;
          br_pr1_b1 = 0.337; // 0.37; // 0.36; // 0.5; // 0.54;
          HDP.Set( "FTF_PION_PROC1_B1", br_pr1_b1 );
          //
          double br_pr1_a2 = 0.;
          HDP.Get( "FTF_PION_PROC1_A2", br_pr1_a2 ); // D=-5.77
          std::cout << " br_pr1_a2 = " << br_pr1_a2 << std::endl;
          br_pr1_a2 = -7.568; // -8.31; // -7.575; // -0.65; // -4.18;
          HDP.Set( "FTF_PION_PROC1_A2", br_pr1_a2 );
          //
          double br_pr1_b2 = 0.;
          HDP.Get( "FTF_PION_PROC1_B2", br_pr1_b2 ); // D=0.8
          std::cout << " br_pr1_b2 = " << br_pr1_b2 << std::endl;
          br_pr1_b2 = 0.44; // 0.458; // 0.374; // 1.24; // 1.14;
          HDP.Set( "FTF_PION_PROC1_B2", br_pr1_b2 );
          //
    */
    /*
          // examples of switching ON/OFF proc2 & proc3 for baryon and/or pion projectile
          //
          bool br_disso_proj = true;
          // HDP.Set( "FTF_BARYON_DIFF_DISSO_PROJ", br_disso_proj ); // D=false
          bool br_disso_tgt  = true;
          HDP.Set( "FTF_BARYON_DIFF_DISSO_TGT",  br_disso_tgt  ); // D=false
          //
          bool pi_disso_proj = true;
          // HDP.Set( "FTF_PION_DIFF_DISSO_PROJ", pi_disso_proj ); // D=false
          bool pi_disso_tgt  = true;
          HDP.Set( "FTF_PION_DIFF_DISSO_TGT",  pi_disso_tgt  ); // D=false

          // example of setting selected nuclear target destruction parameters
          // for baryon projectile
          // (similar can be done for pion)
          //
    */
    /*
          double br_nd_p1_tgt = 0.00173;
          HDP.Set( "FTF_BARYON_NUCDESTR_P1_TGT", br_nd_p1_tgt );
          bool br_nd_p1_adep_tgt = true;
          HDP.Set( "FTF_BARYON_NUCDESTR_P1_ADEP_TGT", br_nd_p1_adep_tgt );
          double br_exci = 26.1;
          HDP.Set( "FTF_BARYON_EXCI_E_PER_WNDNUCLN", br_exci );
    */
    // fProcWrapper = new FTFPWrapper();
    pw = new FTFPWrapper();
  }
  else if (name.find("qgs") != std::string::npos)
  {
    // fProcWrapper = new QGSPWrapper();
    pw = new QGSPWrapper();
  }
  else if (name.find("bertini") != std::string::npos)
  {
    //
    /*
          G4CascadeParameters::Instance();
          G4UImanager* uim = 0;
          G4ApplicationState currentstate = G4StateManager::GetStateManager()->GetCurrentState();
          bool ok = G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);
          if ( !ok )
          {
             G4cout << " G4StateManager PROBLEM: can NOT change state to G4State_Idle !" << G4endl;
          }
          else
          {
             uim = G4UImanager::GetUIpointer();
             // uim->ApplyCommand( "/process/had/cascade/nuclearRadiusScale 1.5" );
             // std::cout << " Updated: Bertini RadiusScale = " <<
       G4CascadeParameters::radiusScale() << std::endl; std::cout << " Bertini UsePreCo = " <<
       G4CascadeParameters::usePreCompound() << std::endl; uim->ApplyCommand(
       "/process/had/cascade/usePreCompound 1" ); std::cout << " Update: Bertini UsePreCo = " <<
       G4CascadeParameters::usePreCompound() << std::endl;

             // G4DeexPrecoParameters* precoparams =
       G4NuclearLevelData::GetInstance()->GetParameters();
       // precoparams->SetLevelDensity( 0.01/CLHEP::MeV );
       // std::cout << " Update: PreCo LevelDens = " << precoparams->GetLevelDensity() << std::endl;
       // precoparams->SetPrecoModelType(1);
       // precoparams->SetDeexModelType(1);
       // precoparams->SetUseAngularGen(false);
       // precoparams->SetDeexChannelsType(fGEM); // it's an enum, with fEvaporation=0, etc.
          }
          ok = G4StateManager::GetStateManager()->SetNewState( currentstate );
    */
    //
    pw = new ProcessWrapper("BertiniProcessWrapper");
    G4CascadeInterface* bert = new G4CascadeInterface();
    bert->SetMaxEnergy(15. * GeV);
    pw->RegisterMe(bert);
  }
  else if (name.find("incl++") != std::string::npos || name.find("inclxx") != std::string::npos)
  {
    //      std::cout << " HyperNuclei enabled : "
    //                << G4HadronicParameters::Instance()->EnableHyperNuclei() << std::endl;
    //      G4HadronicParameters::Instance()->SetEnableHyperNuclei(false);
    pw = new ProcessWrapper("INCLXXProcessWrapper");
    G4INCLXXInterface* inclxx = new G4INCLXXInterface();
    //      std::cout << " HyperNuclei enabled : "
    //                << G4HadronicParameters::Instance()->EnableHyperNuclei() << std::endl;
    inclxx->SetMinEnergy(1. * MeV);  // this is what stands in the INCL-based PL for protons
                                     // it's 0. for pions, but let's make it uniform
    inclxx->SetMaxEnergy(20. * GeV);  // this is what stands in the INCL-based PL for p,n,pi
                                      // overlap is typically set at 15-20GeV with a string model
    // (higher energy model)
    //
    // in principle INCL++ is good up to 3AGeV for p,n,pi, A(<18)
    pw->RegisterMe(inclxx);
  }
#ifdef G4_USE_FLUKA
  else if (name_lc.find("fluka") != std::string::npos)
  {
    const G4bool activateCoalescence = true;
    const G4bool activateHeavyFragmentsEvaporation = true;
    fluka_interface::initialize(activateCoalescence, activateHeavyFragmentsEvaporation);
    fluka_particle_table::initialize();
    pw = new ProcessWrapper("FLUKA.CERN.ProcessQrapper");
    FLUKANuclearInelasticModel* flukaHI = new FLUKANuclearInelasticModel();
    flukaHI->SetMaxEnergy(1. * TeV);
    pw->RegisterMe(flukaHI);
  }
#endif

  // if (!fProcWrapper)
  if (!pw)
  {
    G4cout << " generator " << name << " is unavailable" << G4endl;
    exit(1);
  }

  // std::cout << " process = " << fProcWrapper->GetProcessName() << std::endl;

  if (name.find("lund-str-fragm") != std::string::npos)
  {
    // fProcWrapper->UseG4LundStringFragm(true);
    pw->UseG4LundStringFragm(true);
  }
  else if (name.find("qgsm-str-fragm") != std::string::npos)
  {
    pw->UseG4LundStringFragm(false);
  }

  // fProcWrapper->Compose();
  pw->Compose();

  fProcWrapper = pw;

  return;
}
