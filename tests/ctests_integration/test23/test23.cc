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

#include "G4Timer.hh"

#include "Tst23Histo.hh"

// random engine/seed settings
#include "G4GeometryManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "FTFP_BERT.hh"
#include "FTF_BIC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_FTFP_BERT.hh"
#include "Tst23ActionInit.hh"
#include "Tst23SteppingAction.hh"
#include "TstPhysListReader.hh"
#include "TstPrimaryGeneratorAction.hh"
#include "TstTarget.hh"

// experimental lists !!!
#include "NuBeam.hh"
// for cross-checks
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT_HP.hh"

//
// available staring 4.10.1.b01
//
#include "G4ParticleTable.hh"

#include "Shielding.hh"
#include "ShieldingLEND.hh"
//
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4HadronicProcessStore.hh"

//
// G4 HAD model(s) parameters variation
//
#include "G4StateManager.hh"
#include "G4UImanager.hh"
// bertini
#include "G4CascadeParameters.hh"
// preco
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

#ifdef G4_USE_FLUKA
// FLUKA Interface if required
#  include "G4_HP_CernFLUKAHadronInelastic_PhysicsList.hh"

#  include "FLUKAParticleTable.hh"
#endif

using namespace std;

int main(int argc, char** argv)
{
  TstPhysListReader* theConfigReader = new TstPhysListReader();

  // Seed the Random engine
  // NOTE: Use the default random engine;
  //       do NOT replace default engine with anything else
  //
  // --> CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(1234);  // just something... to be on the safe side

  // Control on input
  if (argc < 2)
  {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  // std::string input = argv[1];
  // G4STring input = argv[1];
  theConfigReader->OpenAppConfig(argv[1]);
  theConfigReader->Help();

  // Note-1: Run manager once and for all;
  //       I can NOT create/delete run manager because
  //       once it's deleated, it'll also cause deletion
  //       of several other singleton - some in the way that
  //       will prevent their re-creation in the same job;
  //       however, run manager can be somewhat re-configured
  //
  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  // Note-2: This is related to an attempt to run a loop
  //       over physics lists, changing them as we go.
  //       It turned out to be "possible with restrictions",
  //       but is very much a headache.
  //       The reason I keep this piece of code, commented,
  //       is just a written record for myself.
  //       The reason why I can NOT create/delete a list
  //       is that when a list is deleted, it EMPTIES
  //       the particle table - and the table can NOT
  //       be re-populated because instances of particles
  //       are already there in the memory, so the system
  //       will check on that and wont move any further
  //       ... thus we'll have particles everywhere in
  //       the memory but not in the table... oops !
  //
  // std::vector<G4VModularPhysicsList*> listCollection;

  // Note-3: All other elements can be re-created in the loop,
  //         although a run with just one phys.list is best,
  //         thus there's no real need to recreate anything...
  //
  TstTarget* geom = 0;
  TstPrimaryGeneratorAction* beam = 0;
  TstHisto* histo = 0;
  Tst23SteppingAction* stepping = 0;
  //
  Tst23ActionInit* action = 0;

  do
  {
    theConfigReader->ProcessConfig();
    if (theConfigReader->IsDone())
      break;  // double protection because
              // in the previous cycle it stopped at #run

    // -->      G4cout << " Check Random Seed: " << theConfigReader->GetRndmSeed() << G4endl;

    CLHEP::HepRandom::setTheSeed(theConfigReader->GetRndmSeed());  // reset the seed

    // -->      G4cout << " Re-Check Random Seed: " << CLHEP::HepRandom::getTheSeed() << G4endl;

    if (!G4StateManager::GetStateManager()->SetNewState(G4State_PreInit))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    if (geom) delete geom;
    geom = new TstTarget();
    G4ThreeVector targetSize = theConfigReader->GetTargetSize();
    geom->SetDimentions(targetSize.x(), targetSize.y(), targetSize.z());
    geom->SetShape(theConfigReader->GetTargetShape());
    geom->ResetMaterial(theConfigReader->GetTargetMaterial());

    runManager->SetUserInitialization(geom);
    runManager->GeometryHasBeenModified();

    // Again, this commented potion below is related to an attempt
    // to run a loop over physics lists.
    // See a note earlier in the code.
    /*
          if ( theConfigReader->GetPhysics() == "ftfp_bert" )
          {
       listCollection.push_back( new FTFP_BERT() );
          }
          else if ( theConfigReader->GetPhysics() == "qgsp_ftfp_bert" )
          {
       listCollection.push_back( new QGSP_FTFP_BERT() );
          }
          else if ( theConfigReader->GetPhysics() == "qgsp_bert" )
          {
       listCollection.push_back( new QGSP_BERT() );
          }
          else if (  theConfigReader->GetPhysics() == "ftf_bic" )
          {
             listCollection.push_back( new FTF_BIC() );
          }
          else if (  theConfigReader->GetPhysics() ==  "NuBeam" )
          {
             listCollection.push_back( new TentativeNuMI() ) ;
          }

          int listCounter = listCollection.size() - 1;
          runManager->SetUserInitialization( listCollection[listCounter] );
          // runManager->PhysicsHasBeenModified();
    */

    /* ---> FOR INTERNAL USE ONLY !!!
          //
          // G4 HAD model(s) parameters variation
          //
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
       G4CascadeParameters::usePreCompound() << std::endl; G4DeexPrecoParameters* precoparams =
       G4NuclearLevelData::GetInstance()->GetParameters();
       precoparams->SetLevelDensity( 1.00/CLHEP::MeV );

          }
          ok = G4StateManager::GetStateManager()->SetNewState( currentstate );
    */
    G4VModularPhysicsList* plist = 0;
    if (theConfigReader->GetPhysics() == "ftfp_bert")
    {
      plist = new FTFP_BERT();
    }
    if (theConfigReader->GetPhysics() == "ftfp_bert_hp")
    {
      plist = new FTFP_BERT_HP();
    }
    else if (theConfigReader->GetPhysics() == "qgsp_ftfp_bert")
    {
      plist = new QGSP_FTFP_BERT();
    }
    else if (theConfigReader->GetPhysics() == "qgsp_bert")
    {
      plist = new QGSP_BERT();
    }
    else if (theConfigReader->GetPhysics() == "qgsp_bert_hp")
    {
      plist = new QGSP_BERT_HP();
    }
    else if (theConfigReader->GetPhysics() == "ftf_bic")
    {
      plist = new FTF_BIC();
    }
    else if (theConfigReader->GetPhysics() == "NuBeam")
    {
      plist = new NuBeam();
    }
    else if (theConfigReader->GetPhysics() == "Shielding")
    {
      plist = new Shielding();
    }
    else if (theConfigReader->GetPhysics() == "ShieldingLEND")
    {
      plist = new ShieldingLEND();
    }
    // --> leave as example --> #if (USE_G4REF==0 && USE_G4PUBLIC==0) || USE_G4REF>4100006 ||
    // USE_G4PUBLIC>41000
    else if (theConfigReader->GetPhysics() == "ShieldingM")
    {
      plist = new Shielding(1, "HP", "M");
    }
#ifdef G4_USE_FLUKA
    else if (theConfigReader->GetPhysics().find("G4_HP_CFLUKAHI") != std::string::npos)
    {
      plist = new G4_HP_CernFLUKAHadronInelastic_PhysicsList();
    }
#endif

    runManager->SetUserInitialization(plist);

#ifdef G4_USE_FLUKA
    if (theConfigReader->GetPhysics().find("G4_HP_CFLUKAHI") != std::string::npos)
    {
      fluka_particle_table::initialize();
    }
#endif

    if (beam) delete beam;
    beam = new TstPrimaryGeneratorAction();
    beam->InitBeam(theConfigReader);
    // up to 4.9.6.p02
    runManager->SetUserAction(beam);

    if (histo) delete histo;
    histo = new Tst23Histo(theConfigReader);

    if (stepping) delete stepping;
    stepping = new Tst23SteppingAction(histo);

    // up to 4.9.6.p02
    // runManager->SetUserAction( stepping );

    // from 4.10.b01 onwards
    if (action) delete action;
    action = new Tst23ActionInit();
    action->SetAct(beam);
    action->SetAct(stepping);

    runManager->SetUserInitialization(action);

    if (!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    runManager->InitializeGeometry();
    runManager->InitializePhysics();
    runManager->Initialize();

    // G4double xsec = beam->GetXSecOnTarget() / millibarn;
    // G4cout << " xsec = " << xsec << G4endl;

    runManager->ConfirmBeamOnCondition();
    runManager->RunInitialization();  // this is part of BeamOn
                                      // and needs be done (at least) to set GeomClosed status

    stepping->SetTargetPtr(geom->GetTarget());

    G4Timer timer;
    timer.Start();

    for (G4int iter = 0; iter < theConfigReader->GetNEvents(); ++iter)
    {
      /*
         if ( (iter%10) == 0 )
         {
            G4cout << " event # " << iter << G4endl;
         }
      */
      runManager->ProcessOneEvent(iter);
      // G4Event* event = runManager->GenerateEvent( iter );
    }

    timer.Stop();
    G4cout << " CPU = " << timer.GetUserElapsed() << G4endl;
    G4cout << " Real Time = " << timer.GetRealElapsed() << G4endl;

    // G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* partDef =
      G4ParticleTable::GetParticleTable()->FindParticle(theConfigReader->GetBeamParticle());
    G4double partMass = partDef->GetPDGMass();
    G4double partMom = theConfigReader->GetBeamMomentum();
    G4double partEnergy =
      std::sqrt(partMom * partMom + partMass * partMass);  // this is total energy
    G4DynamicParticle dParticle(partDef, theConfigReader->GetDirection(), partEnergy - partMass);
    // dParticle.DumpInfo();
    G4Material* mat = geom->GetCurrentMaterial();
    const G4Element* elm = mat->GetElement(0);

    G4HadronicProcessStore* const store = G4HadronicProcessStore::Instance();
    // store->PrintInfo(partDef);
    // store->Dump(1);
    // NOTE: Technically speaking, this method below does extract xsec *without* material (nullptr)
    //       but it produces a warning of not being able to calculate something unless mat is given
    double xsec_from_store =
      store->GetInelasticCrossSectionPerAtom(partDef, dParticle.GetKineticEnergy(), elm, mat);

    std::cout << " xsec_from_beam = " << beam->GetXSecOnTarget() / millibarn << std::endl;
    std::cout << " xsec_from_store = " << xsec_from_store / millibarn << std::endl;

    histo->Write(theConfigReader->GetNEvents(), (beam->GetXSecOnTarget() / millibarn));

    runManager->RunTermination();

  } while (!theConfigReader->IsDone());

  G4cout << "###### End of test #####" << G4endl;

  // See notes earlier in the code, on phys.lists business
  /*
     int NLists = listCollection.size();

     // remove all physics lists but the last one
     // as it HAS to be removed by run manager dtor
     //
     for ( int i=0; i<NLists-1; i++ )
     {
        if ( listCollection[i] )
        {
           delete listCollection[i];
     listCollection[i] = 0;
        }
     }
  */
  delete runManager;

  delete theConfigReader;

  return 0;
}
