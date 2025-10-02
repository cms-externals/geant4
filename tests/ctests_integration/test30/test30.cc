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

#include <fstream>
#include <iomanip>

#include "globals.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"

#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ElementVector.hh"
#include "Test30Material.hh"
#include "Test30Physics.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleChange.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4HadronicProcess.hh"
#include "G4ChargeExchangeProcess.hh"
#include "G4VCrossSectionDataSet.hh"

#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionElastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4GammaNuclearXS.hh"

#include "G4ElasticHadrNucleusHE.hh"
#include "G4ChargeExchange.hh"
#include "G4HadronElastic.hh"
#include "G4DiffuseElastic.hh"
#include "G4ChipsElasticModel.hh"
#include "G4LowEHadronElastic.hh"
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4HadronElasticProcess.hh"
#include "G4WentzelVIModel.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DynamicParticle.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Alpha.hh"
#include "G4He3.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4GenericIon.hh"
#include "G4DecayPhysics.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4EmParameters.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

#include "Histo.hh"
#include "G4Timer.hh"

#include "G4ForceCondition.hh"
#include "G4TouchableHistory.hh"

#include "G4NucleiProperties.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4TransportationManager.hh"
#include "G4ExceptionHandler.hh"

int main(int argc, char** argv)
{
  G4cout << "========================================================" << G4endl;
  G4cout << "======        Cascade Test (test30) Start       ========" << G4endl;
  G4cout << "========================================================" << G4endl;

  G4VParticleChange* aChange = nullptr;
  G4StateManager* g4State=G4StateManager::GetStateManager();
  g4State->SetNewState(G4State_PreInit);
  g4State->SetExceptionHandler(new G4ExceptionHandler());

  //-----------------------------------------------------------------------------
  // ------- Initialisation 
  Histo     histo;
  G4String  namePart = "proton";
  G4bool    ionParticle = false;
  G4int     ionZ(0), ionA(0);
  G4String  nameMat  = "G4_Al";
  G4String  nameGen  = "binary";
  G4bool    usehis   = false;
  G4bool    isInitH  = false;
  G4bool    inclusive= true;
  G4bool    saverand = false;
  G4bool    applyMSC = false;
  G4int     verbose  = 0;
  G4int     verboseD = 0;
  G4double  energy   = 100.*MeV;
  G4double  sigmae   = 0.0;
  G4double  elim     = 30.*MeV;
  G4double  energyg  = 40.*MeV;
  G4double  balance  = 5.0*MeV;
  G4double  dangl    = 5.0;
  G4int     nevt     = 1000;
  G4int     nbins    = 100;
  G4int     nbinsa   = 40;
  G4int     nbinse   = 80;
  G4int     nbinsd   = 20;
  G4int     nbinspi  = 20;
  G4int     nangl    = 0;
  G4int     nanglpr  = 0;
  G4int     nanglpi  = 0;
  G4int     modu     = 20000;
  G4int     targetA  = 0;
  G4String hFile     = "test30";
  G4double aStep     = 0.01*micrometer;
  G4double bStep     = 0.01*micrometer;
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double  emax     = 160.*MeV;
  G4double  emaxpi   = 200.*MeV;
  G4double ebinlog   = 2.0*MeV;
  G4double eBound    = 70.*MeV;
  G4double tetmax    = 180.;
  G4double costmax   = -1.0;
  G4double xxl       = 100.;
  G4double kBound    = 0.2;

  G4double rmsProton = 0.0;
  G4double rmsNeutron= 0.0;
  G4double rmsPion   = 0.0;

  G4Material* material = nullptr;
  G4NuclearLevelData* fLevelData = G4NuclearLevelData::GetInstance(); 
  G4DeexPrecoParameters* fParam = fLevelData->GetParameters();
  
  G4bool xsbgg = true;
  G4bool xschips = false;

  G4double ang[20] = {0.0};
  G4double bng1[20] = {0.0};
  G4double bng2[20] = {0.0};
  G4double cng[20] = {0.0};
  G4double angpr[20] = {0.0};
  G4double bng1pr[20] = {0.0};
  G4double bng2pr[20] = {0.0};
  G4double cngpr[20] = {0.0};
  G4double angpi[20] = {0.0};
  G4double bngpi1[20] = {0.0};
  G4double bngpi2[20] = {0.0};
  G4double cngpi[20] = {0.0};
  float bestZ[250] = {
    0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0,   //0
    4.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 9.0, 0.0,   //10
    10.0, 10.0, 11.0, 0.0, 11.0, 0.0, 13.0, 0.0, 0.0, 0.0,   //20
    0.0, 0.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //30
    0.0, 0.0, 0.0, 0.0, 21.0, 0.0, 21.0, 21.0, 23.0, 0.0,   //40
    0.0, 24.0, 25.0, 0.0, 25.0, 0.0, 27.0, 27.0, 27.0, 26.0,   //50
    27.0, 0.0, 0.0, 0.0, 0.0, 30.0, 31.0, 31.0, 32.0, 32.0,   //60
    33.0, 33.0, 33.0, 34.0, 33.0, 34.0, 35.0, 35.0, 0.0, 36.0,   //70
    36.0, 37.0, 36.0, 37.0, 37.0, 38.0, 39.0, 39.0, 41.0, 40.0,   //80
    41.0, 39.0, 41.0, 38.0, 39.0, 40.0, 40.0, 39.0, 40.0, 0.0,   //90
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //100
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //110
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //120
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //130
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //140
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //150
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //160
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //170
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //180
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //190
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //200
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //210
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //220
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //230
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };        //240

  // Track
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4double      aTime     = 0. ;
  G4ThreeVector aDirection= G4ThreeVector(0.0,0.0,1.0);
  G4ThreeVector bDirection= G4ThreeVector(0.0,0.0,1.0);
  G4double nx(0.0), ny(0.0), nz(0.0);

  G4cout.setf( std::ios::scientific, std::ios::floatfield );

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

  // -------------------------------------------------------------------
  //--------- Materials definition ---------
  Test30Material*  mate = new Test30Material();
  G4NistManager::Instance()->SetVerbose(0);

  //--------- Particles definition ---------
  Test30Physics*   phys = new Test30Physics();

  const G4ParticleDefinition* gamma = G4Gamma::Gamma();
  const G4ParticleDefinition* electron = G4Electron::Electron();
  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4AntiProton::AntiProton();
  G4AntiNeutron::AntiNeutron();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pi0 = G4PionZero::PionZero();
  const G4ParticleDefinition* deu = G4Deuteron::DeuteronDefinition();
  const G4ParticleDefinition* tri = G4Triton::TritonDefinition();
  const G4ParticleDefinition* he3 = G4He3::He3Definition();
  const G4ParticleDefinition* alp = G4Alpha::AlphaDefinition();
  const G4double massalp = alp->GetPDGMass();
  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));

  G4DecayPhysics decays;
  decays.ConstructParticle();  

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ions = partTable->GetIonTable();

  partTable->SetReadiness();
  ions->CreateAllIon();
  ions->CreateAllIsomer();

  //--------- Geometry definition

  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);
  G4TransportationManager::GetTransportationManager()
    ->SetWorldForTracking(pFrame);

  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;

  for(G4int run=0; run<100; run++) {

    // ---- Read input file
    do {
      if( fin->eof() ) {
        end = false;
        break;
      }
      (*fin) >> line;
      if(verbose > 0) G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
        ionParticle= false;
      } else if(line == "#ion") {
        ionParticle= true;
	namePart="GenericIon";
        (*fin) >> ionZ >> ionA;
      } else if(line == "#energy(MeV)") {
        (*fin) >> energy;
        energy *= MeV;
        emax    = energy;
      } else if(line == "#sigmae(MeV)") {
        (*fin) >> sigmae;
        sigmae *= MeV;
      } else if(line == "#emax(MeV)") {
        (*fin) >> emax;
        emax *= MeV;
      } else if(line == "#emaxg(MeV)") {
        (*fin) >> energyg;
        energyg *= MeV;
      } else if(line == "#emaxpi(MeV)") {
        (*fin) >> emaxpi;
        emaxpi *= MeV;
      } else if(line == "#elim(MeV)") {
        (*fin) >> elim;
        elim *= MeV;
      } else if(line == "#ebalance(MeV)") {
        (*fin) >> balance;
        balance *= MeV;
      } else if(line == "#ebinlog(MeV)") {
        (*fin) >> ebinlog;
	if (ebinlog < 1.1) ebinlog = 1.1;
        ebinlog *= MeV;
      } else if(line == "#events") {
        (*fin) >> nevt;
	char* st = std::getenv("statistic");
	if(st) {
	  G4int st0 = atol(st);
	  G4cout<<"Number of events according to macro "<< nevt << G4endl;
	  if (st0 > 0){  
	    nevt = st0;
	    G4cout<<"Number of events is forced to "<< nevt << G4endl;
	  }
	}
      } else if(line == "#exclusive") {
        inclusive = false;
      } else if(line == "#inclusive") {
        inclusive = true;
      } else if(line == "#nbins") {
        (*fin) >> nbins;
      } else if(line == "#nbinse") {
        (*fin) >> nbinse;
      } else if(line == "#nbinsa") {
        (*fin) >> nbinsa;
      } else if(line == "#nbinsd") {
        (*fin) >> nbinsd;
      } else if(line == "#nbinspi") {
        (*fin) >> nbinspi;
      } else if(line == "#nangle") {
        (*fin) >> nangl;
      } else if(line == "#nanglepr") {
        (*fin) >> nanglpr;
      } else if(line == "#nanglepi") {
        (*fin) >> nanglpi;
      } else if(line == "#anglespi") {
        for(int k=0; k<nanglpi; k++) {(*fin) >> angpi[k];}
      } else if(line == "#dangle") {
        (*fin) >> dangl;
      } else if(line == "#angles") {
        for(int k=0; k<nangl; k++) {(*fin) >> ang[k];}
      } else if(line == "#anglespr") {
        for(int k=0; k<nanglpr; k++) {(*fin) >> angpr[k];}
      } else if(line == "#thetamax") {
        (*fin) >> tetmax;
        if(tetmax <= 0.0 || tetmax >= 180.) tetmax = 180.;
        costmax = std::cos(tetmax*degree); 
      } else if(line == "#range(mm)") {
        (*fin) >> range;
        range *= mm;
      } else if(line == "#step(mm)") {
        (*fin) >> bStep;
        bStep *= mm;
      } else if(line == "#step(g/cm2)") {
        (*fin) >> aStep;
        aStep *= g/cm2;
        applyMSC = true;
      } else if(line == "#print") {
        (*fin) >> modu;
      } else if(line == "#eBound(MeV)") {
        (*fin) >> eBound;
        eBound *= MeV;
      } else if(line == "#kBound") {
        (*fin) >> kBound;
      } else if(line == "#material") {
        (*fin) >> nameMat;
      } else if(line == "#targetA") {
        (*fin) >> targetA;
      } else if(line == "#generator") {
	nameGen = "";
        (*fin) >> nameGen;
	if (nameGen == "" || nameGen == "chips" || 
	    nameGen == "lepar" || nameGen == "rpg") {
	  G4cout << "Generator " << nameGen << " is not allowed" << G4endl; 
          nameGen = "";
	  continue;
	}
        usehis = true;
	hFile = nameGen;
	char* cc = std::getenv(nameGen);
        if(nullptr == cc) {
	  G4cout << "Generator <" << nameGen << "> is not included in the "
		 << " list SetModels.csh, so is ignored!" 
		 << G4endl; 
	  G4cout << "Use #run command to overcome this limitation " << G4endl;
	  continue;
	}
	G4String ss(cc);
	if(ss=="1") { break; }
      } else if(line == "#run") { 
	if("" != nameGen) { break; }
      } else if(line == "#verbose") {
        (*fin) >> verbose;
	G4NistManager::Instance()->SetVerbose(verbose);
        if(2 < verbose) fParam->SetVerbose(verbose - 2);
      } else if(line == "#verboseD") {
        (*fin) >> verboseD;
	G4NuclearLevelData::GetInstance()->GetParameters()->SetVerbose(verboseD);
      } else if(line == "#position(mm)") {
        (*fin) >> nx >> ny >> nz;
        aPosition = G4ThreeVector(nx*mm, ny*mm, nz*mm);
      } else if(line == "#direction") {
        (*fin) >> nx >> ny >> nz;
        if(nx*nx + ny*ny + nz*nz > 0.0) {
          aDirection = G4ThreeVector(nx, ny, nz);
          aDirection = aDirection.unit();
	}
      } else if(line == "#time(ns)") {
        (*fin) >> aTime;
        aTime *= ns;
      } else if(line == "#saverand") {
        saverand = true;
      } else if(line == "#random") {
        G4String sss("");
        (*fin) >> sss;
        CLHEP::HepRandom::restoreEngineStatus(sss);
        if(verbose>0) G4cout << "Random Engine restored from file <"
                             << sss << ">" << G4endl;
        CLHEP::HepRandom::showEngineStatus();
      } else if(line == "#xs_ghad") {
	xsbgg = false;
	xschips = false;
      } else if(line == "#xs_chips") {
	xsbgg = false;
	xschips = true;
      } else if(line == "#xs_bgg") {
	xsbgg = true;
	xschips = false;
      } else if(line == "#exit") {
        end = false;
        break;
      } else if(line == "#rmsEp") {
        (*fin) >> rmsProton;
      } else if(line == "#rmsEn") {
        (*fin) >> rmsNeutron;
      } else if(line == "#rmsEpi") {
        (*fin) >> rmsPion;
      } else if(line == "#HETCEmission") {
        fParam->SetUseHETC(true);
      } else if(line == "#PrecoAngularON") {
	fParam->SetUseAngularGen(true);
      } else if(line == "#PrecoAngularOFF") {
	fParam->SetUseAngularGen(false);
      } else if(line == "#GEMEvaporation") {
        G4cout<<"### GEM evaporation is set"<<G4endl;
        fParam->SetDeexChannelsType(fGEM);
      } else if(line == "#DefGEMEvaporation") {
        G4cout<<"### Combined Default+GEM evaporation is set"<<G4endl;
        fParam->SetDeexChannelsType(fCombined);
      } else if(line == "#Evaporation") {
        G4cout<<"### Default evaporation is set"<<G4endl;
        fParam->SetDeexChannelsType(fEvaporation);
      } else if(line == "#XSoption") {
        G4int OPTxs;
        (*fin)>>OPTxs;
	if (OPTxs< 0 || OPTxs >4  ){
	  G4cout << "### Warning: BAD CROSS SECTION OPTION for PreCompound model " 
		 << OPTxs << " ignored" << G4endl;
	} else {

          fParam->SetPrecoModelType(OPTxs);
	  fParam->SetDeexModelType(OPTxs);
	  G4cout<<" Option for inverse cross sections : OPTxs="<<OPTxs<<G4endl;
	}

      } else if(line == "#UseNeverGoBack") {
        G4cout<<" Never Go Back hypothesis has been assumed at preequilibrium"
	      <<G4endl;
	fParam->SetNeverGoBack(true);

      } else if(line == "#UseSoftCutOff") {
        G4cout<<" Soft Cut Off  hypothesis has been assumed at preequilibrium"
	      <<G4endl;
	fParam->SetUseSoftCutoff(true);

      } else if(line == "#UseCEMTransitions") {
        G4cout<<" Transition probabilities at preequilibrium based on CEM model"
	      <<G4endl;
	fParam->SetUseCEM(true);

      } else if(line == "#MFenergy(MeV)") {
        G4double tmin;
        (*fin) >> tmin;
        G4cout<<" Min energy for multi-fragmentation is set to " << tmin 
	      << " MeV" << G4endl;
	fParam->SetMinExPerNucleounForMF(tmin*MeV);
      }      
    } while(end);

    if(!end) { break; }

    // -------------------------------------------------------------------
    // -------- Start run processing
    G4cout << "###### Start new run # " << run << "   for "
	   << nevt << " events  #####" << G4endl;

    // -------- Target 
    material = mate->GetMaterial(nameMat);
    if(nullptr == material) {
      G4cout << "Material <" << nameMat << "> is not found out"
	     << G4endl;
      exit(1);
    }
    G4ProductionCuts* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(1.0, 0);
    cuts->SetProductionCut(1.0, 1);
    cuts->SetProductionCut(1.0, 2);
    cuts->SetProductionCut(0.0, 3);

    G4MaterialCutsCouple* couple = new G4MaterialCutsCouple(material,cuts);
    couple->SetIndex(0);
    G4ProductionCutsTable* pkt = G4ProductionCutsTable::GetProductionCutsTable();
    std::vector<G4double>* vc = 
      const_cast<std::vector<G4double>*>(pkt->GetEnergyCutsVector(3));
    vc->push_back(0.0); 
    const G4Element* elm = material->GetElement(0); 
    theStep = bStep;

    // -------- Projectile
    G4ParticleDefinition* part = nullptr;
    if (!ionParticle) { part = partTable->FindParticle(namePart); }
    else { part = partTable->GetIonTable()->GetIon(ionZ,ionA,0.0); }
    if (nullptr == part) {
      G4cout << " Sorry, No definition for particle" << namePart 
	     << " found" << G4endl;
      exit(1);  
    }
    G4DynamicParticle dParticle(part,aDirection,energy);

    // ------- Select model
    G4HadronicProcess* extraproc = nullptr;
    G4VProcess* proc = nullptr;
    G4String namegen1 = nameGen.substr(0, 4);
    G4cout << "<" << namegen1 << ">" << G4endl; 
    if(namegen1 == "Elas") { 
      extraproc = new G4HadronElasticProcess(); 
    } else if (nameGen == "chargeex") { 
      extraproc = new G4ChargeExchangeProcess(); 
    } else { 
      proc = phys->GetProcess(nameGen, part, material); 
    }

    if(!proc && !extraproc) {
      if("" != nameGen) {
	G4cout << "For particle: " << part->GetParticleName()
	       << " generator <" << nameGen << " is unavailable"<< G4endl;
      }
      continue;
    }

    // ------- Define target A
    G4int Z = elm->GetZasInt();
    G4int A = G4lrint(elm->GetN());
    if(targetA > 0) { A = targetA; }
    if(nullptr != proc) { phys->SetA(A); }

    // ------- Binning 

    G4int maxn = A + 1;

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    if(verbose > 0) {
      G4cout << "The material:  " << material->GetName() 
	     << " Z= " << Z << " A= " << A
	     << "  Amax= " << maxn << G4endl;
      G4cout << "The step:      " << theStep/mm << " mm" << G4endl;
      G4cout << "The position:  " << aPosition/mm << " mm" << G4endl;
      G4cout << "The direction: " << aDirection << G4endl;
      G4cout << "The time:      " << aTime/ns << " ns" << G4endl;
    }
    // linear histograms
    G4double mass = part->GetPDGMass();
    G4double pmax = std::sqrt(energy*(energy + 2.0*mass));
    G4double bine = emax/(G4double)nbinse;
    G4double bind = emax/(G4double)nbinsd;

    G4cout << "energy = " << energy/MeV << " MeV   " 
	   << " RMS(MeV)= " << sigmae/MeV << G4endl;
    if(verbose > 0) {
      G4cout << "emax   = " << emax/MeV << " MeV" << G4endl;
      G4cout << "pmax   = " << pmax/MeV << " MeV" << G4endl;
    }
    // logarithmic histograms
    G4double logmax = 2.0;
    G4double logmin = -2.0;
    G4int nbinlog   = 80;
    if(emax > 100.*MeV) {
      logmax  = 3.0;
      logmin  = -1.0;
      nbinlog = 80;
    }
    if(emax > 1000.*MeV) {
      logmax  = 4.0;
      logmin  = -1.0;
      nbinlog = 100;
    }
    G4double binlog = (logmax - logmin)/G4double(nbinlog);
    G4cout << "### Log10 scale from " << logmin << " to " << logmax 
	   << " in " << nbinlog << " bins" << G4endl; 

    // ------- Histograms
    if(usehis && !isInitH) {

      isInitH = true;
    
      histo.Add1D("1","Number of Secondaries",100,-0.5,99.5);
      histo.Add1D("2","Type of secondary",10,-0.5,9.5);
      histo.Add1D("3","Phi(degrees) of Secondaries",90,-180.0,180.0);

      histo.Add1D("4","ds/dE for protons at theta = 0",nbinsd,0.,emax);
      histo.Add1D("5","ds/dE for protons at theta = 1",nbinsd,0.,emax);
      histo.Add1D("6","ds/dE for protons at theta = 2",nbinsd,0.,emax);
      histo.Add1D("7","ds/dE for protons at theta = 3",nbinsd,0.,emax);
      histo.Add1D("8","ds/dE for protons at theta = 4",nbinsd,0.,emax);
      histo.Add1D("9","ds/dE for protons at theta = 5",nbinsd,0.,emax);
      histo.Add1D("10","ds/dE for protons at theta = 6",nbinsd,0.,emax);
      histo.Add1D("11","ds/dE for protons at theta = 7",nbinsd,0.,emax);
      histo.Add1D("12","ds/dE for protons at theta = 8",nbinsd,0.,emax);
      histo.Add1D("13","ds/dE for protons at theta = 9",nbinsd,0.,emax);
      histo.Add1D("14","ds/dE for protons at theta = 10",nbinsd,0.,emax);

      G4int i;
      // proton double differential histograms are active by request
      for(i=nanglpr; i<11; ++i) {histo.Activate(3+i, false);}

      histo.Add1D("15","E(MeV) for gamma",200,0.,energyg);
      histo.Add1D("16","delta E (MeV) ",100,-balance,balance);
      histo.Add1D("17","delta Pz (MeV/c)",100,-balance,balance);
      histo.Add1D("18","delta Pt (MeV/c)",100,-balance,balance);

      histo.Add1D("19","E(MeV) for pi0",nbinse,0.,emax);
      histo.Add1D("20","E(MeV) for pi+",nbinse,0.,emax);
      histo.Add1D("21","E(MeV) for pi-",nbinse,0.,emax);
      histo.Add1D("22","E(MeV) protons",nbinse,0.,emax);
      histo.Add1D("23","E(MeV) neutrons",nbinse,0.,emax);

      histo.Add1D("24","Phi(degrees) of neutrons",90,-180.0,180.0);

      histo.Add1D("25","cos(theta) protons",nbinsa,costmax,1.);
      histo.Add1D("26","cos(theta) neutrons",nbinsa,costmax,1.);

      histo.Add1D("27","Baryon number (mbn)",maxn,-0.5,(G4double)maxn + 0.5);

      histo.Add1D("28","ds/dE for neutrons at theta = 0",nbinsd,0.,emax);
      histo.Add1D("29","ds/dE for neutrons at theta = 1",nbinsd,0.,emax);
      histo.Add1D("30","ds/dE for neutrons at theta = 2",nbinsd,0.,emax);
      histo.Add1D("31","ds/dE for neutrons at theta = 3",nbinsd,0.,emax);
      histo.Add1D("32","ds/dE for neutrons at theta = 4",nbinsd,0.,emax);
      histo.Add1D("33","ds/dE for neutrons at theta = 5",nbinsd,0.,emax);
      histo.Add1D("34","ds/dE for neutrons at theta = 6",nbinsd,0.,emax);
      histo.Add1D("35","ds/dE for neutrons at theta = 7",nbinsd,0.,emax);
      histo.Add1D("36","ds/dE for neutrons at theta = 8",nbinsd,0.,emax);
      histo.Add1D("37","ds/dE for neutrons at theta = 9",nbinsd,0.,emax);
      histo.Add1D("38","ds/dE for neutrons at theta = 10",nbinsd,0.,emax);
      histo.Add1D("39","ds/dE for neutrons at theta = 11",nbinsd,0.,emax);
      histo.Add1D("40","ds/dE for neutrons at theta = 12",nbinsd,0.,emax);

      // neutron double differential histograms are active by request
      for(i=nangl; i<13; ++i) {histo.Activate(27+i, false);}
 
      histo.Add1D("41","ds/dE for pi- at theta = 0",nbinspi,0.,emaxpi);
      histo.Add1D("42","ds/dE for pi- at theta = 1",nbinspi,0.,emaxpi);
      histo.Add1D("43","ds/dE for pi- at theta = 2",nbinspi,0.,emaxpi);
      histo.Add1D("44","ds/dE for pi- at theta = 3",nbinspi,0.,emaxpi);
      histo.Add1D("45","ds/dE for pi- at theta = 4",nbinspi,0.,emaxpi);
      histo.Add1D("46","ds/dE for pi+ at theta = 0",nbinspi,0.,emaxpi);
      histo.Add1D("47","ds/dE for pi+ at theta = 1",nbinspi,0.,emaxpi);
      histo.Add1D("48","ds/dE for pi+ at theta = 2",nbinspi,0.,emaxpi);
      histo.Add1D("49","ds/dE for pi+ at theta = 3",nbinspi,0.,emaxpi);
      histo.Add1D("50","ds/dE for pi+ at theta = 4",nbinspi,0.,emaxpi);

      // pion double differential histograms are active by request
      for(i=nanglpi; i<5; ++i) {
	histo.Activate(40+i, false);
	histo.Activate(45+i, false);
      }
      histo.Add1D("51","E(MeV) neutrons",nbinlog,logmin,logmax);
      histo.Add1D("52","ds/dE for neutrons at theta = 0",nbinlog,logmin,logmax);
      histo.Add1D("53","ds/dE for neutrons at theta = 1",nbinlog,logmin,logmax);
      histo.Add1D("54","ds/dE for neutrons at theta = 2",nbinlog,logmin,logmax);
      histo.Add1D("55","ds/dE for neutrons at theta = 3",nbinlog,logmin,logmax);
      histo.Add1D("56","ds/dE for neutrons at theta = 4",nbinlog,logmin,logmax);
      histo.Add1D("57","ds/dE for neutrons at theta = 5",nbinlog,logmin,logmax);

      // neutron double differential histograms are active by request
      for(i=nangl; i<6; ++i) {histo.Activate(51+i, false);}

      // elastic
      histo.Add1D("58","Ekin (MeV) for primary particle",120,0.,energy*1.2/MeV);
      histo.Add1D("59","cos(Theta) for recoil particle in Lab.Sys.",nbinsa,costmax,1.);
      histo.Add1D("60","cos(Theta) for primary particle in Lab.Sys.",nbinsa,costmax,1.);
      histo.Add1D("61","cos(Theta) for recoil particle in CM.Sys.",nbinsa,costmax,1.);
      histo.Add1D("62","cos(Theta) for primary particle in CM.Sys.",nbinsa,costmax,1.);
      histo.Add1D("63","Ekin (MeV) for recoil particle",120,0.,energy*1.2/MeV);
      histo.Add1D("64","Theta (degree) for primary particle in Lab.Sys.",nbinsa,0.0,tetmax);
      G4double x2 = std::log10(tetmax);
      G4double x1 = x2 - std::log10((double)nbinsa);
      xxl = x2 - x1;
      histo.Add1D("65","log10(theta (degree)) for primary particle in Lab.Sys.",nbinsa,x1,x2);
      histo.Add1D("66","Theta (degree) for primary particle in CM.Sys.",nbinsa,0.0,tetmax);
      for(i=57; i<67; ++i) {histo.Activate(i, extraproc);}

      // alpha spectra
      histo.Add1D("67","E(MeV) alpha",nbinse,0.,emax);
      // energy of the relative motion of two alpha particles
      histo.Add1D("68","Erel(MeV) of alpha pair",nbinse,0.,emax);

      if(extraproc) {

	// desactivate not needed hist for elastic
	histo.Activate(50, false);
	for(i=0; i<13; ++i) {histo.Activate(14+i, false);}
      }
    }

    if(usehis) {
      histo.SetFileName(hFile);
      histo.Book();
      G4cout << "Histograms are booked output file <" << hFile << ".root> "
	     << G4endl;
    }

    // -------- Normalisation

    G4VCrossSectionDataSet* cs = nullptr;
    G4WentzelVIModel* emcscat  = nullptr;
    G4ParticleChangeForMSC* gParticleChange = nullptr;

    G4DataVector emcuts;
    emcuts.resize(1, 0.0);
    G4double cross_sec = 0.0;

    if(extraproc && namegen1 == "Elas") {

      if(nameGen == "ElasticHP") {
	cs = new G4ParticleHPElasticData(); 
      } else if(part == proton) {
	if(xschips)      { cs = new G4ChipsProtonElasticXS(); }
	else             { cs = new G4BGGNucleonElasticXS(part); }
      } else if(part == neutron) {
	if(xschips)      { cs = new G4ChipsNeutronElasticXS(); }
        else             { cs = new G4NeutronElasticXS(); } 
      } else if(part == pip) {
	if(xschips)      { cs = new G4ChipsPionPlusElasticXS(); }
	else             { cs = new G4BGGPionElasticXS(part); }
      } else if(part == pin) {
	if(xschips)      { cs = new G4ChipsPionMinusElasticXS(); }
	else             { cs = new G4BGGPionElasticXS(part); }
      } else {
	cs = new G4CrossSectionElastic(new G4ComponentGGHadronNucleusXsc());
      }
      extraproc->AddDataSet(cs);
      G4ProcessManager* man = new G4ProcessManager(part);
      if(nullptr == man) {  
	G4cout << " Sorry, No manager available for particle" 
	       << namePart << G4endl;
	exit(1);  
      }

      G4HadronicInteraction* els = nullptr;
      if(nameGen == "Elastic")           { els = new G4HadronElastic(); }
      else if(nameGen == "ElasticLHEP")  { els = new G4HadronElastic(); }
      else if(nameGen == "ElasticHE")    { els = new G4ElasticHadrNucleusHE(); }
      else if(nameGen == "ElasticDIFF")  { els = new G4DiffuseElastic(); }
      else if(nameGen == "ElasticCHIPS") { els = new G4ChipsElasticModel(); }
      else if(nameGen == "ElasticHP")    { els = new G4ParticleHPElastic(); }
      else if(nameGen == "ElasticLE")    { els = new G4LowEHadronElastic(); }
      extraproc->RegisterMe(els);
      man->AddDiscreteProcess(extraproc); 

      if(applyMSC) {
	theStep = aStep/material->GetDensity();
	emcscat = new G4WentzelVIModel(false);
	G4EmParameters::Instance()->SetMscThetaLimit(pi);
	G4EmParameters::Instance()->SetVerbose(0);
	G4EmParameters::Instance()->SetMuHadLateralDisplacement(false);
	gParticleChange = new G4ParticleChangeForMSC();
	emcscat->SetParticleChange(gParticleChange, nullptr);
	emcscat->SetCurrentCouple(couple);
	emcscat->Initialise(part, emcuts);
      }
    } else if( part == proton ) {
      if(xsbgg)        { cs = new G4BGGNucleonInelasticXS(part); }
      else             { cs = new G4ParticleInelasticXS(part); }
    } else if( part == neutron ) {
      cs = new G4NeutronInelasticXS();
    } else if( part == pin || part == pip ) {
      cs = new G4BGGPionInelasticXS(part);
    } else if( part == gamma ) { 
      cs = new G4GammaNuclearXS();  
    } else if( part == deu || part == tri ||
               part == he3 || part == alp) {
      cs = new G4ParticleInelasticXS(part);
    } else if( ionParticle ) {
      cs = new G4CrossSectionInelastic(new G4ComponentGGNuclNuclXsc());
    } else {
      cs = new G4CrossSectionInelastic(new G4ComponentGGHadronNucleusXsc());
    }

    // -------- Track

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition);
    G4TouchableHandle fpTouchable(new G4TouchableHistory());
    gTrack->SetTouchableHandle(fpTouchable);
    std::vector<G4DynamicParticle*> vdp; 

    // -------- Step
    G4Step* step;
    step = new G4Step();
    step->SetTrack(gTrack);
    gTrack->SetStep(step);
    
    G4StepPoint *aPoint, *bPoint;
    aPoint = new G4StepPoint();
    aPoint->SetPosition(aPosition);
    aPoint->SetMaterial(material);
    aPoint->SetMaterialCutsCouple(couple);
    G4double safety = 10000.*cm;
    aPoint->SetSafety(safety);
    step->SetPreStepPoint(aPoint);

    bPoint = aPoint;
    G4ThreeVector bPosition = aDirection*theStep;
    bPosition += aPosition;
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);
    step->SetStepLength(theStep);
    G4cout << "Step Length(mm)= " << theStep << " applyMSC: " << applyMSC << G4endl; 

    if(extraproc) {
      extraproc->PreparePhysicsTable(*part);
      extraproc->BuildPhysicsTable(*part);
      cross_sec = extraproc->GetElementCrossSection(&dParticle,elm,material);

    } else {
      cs->BuildPhysicsTable(*part);
      cross_sec = cs->GetCrossSection(&dParticle, elm);
    }

    G4double factor  = 
      cross_sec*MeV*1000.0*(G4double)nbinse/(energy*barn*(G4double)nevt);
    G4double factora = 
      cross_sec*1000.0*(G4double)nbinsa/(twopi*(1.0 - costmax)*barn*(G4double)nevt);
    G4double factoraa = 
      cross_sec*1000.0*(G4double)nbinsa/(twopi*tetmax*degree*barn*(G4double)nevt);
    G4double factoral = 
      cross_sec*1000.0*(G4double)nbinsa/(twopi*xxl*std::log(10.)*barn*(G4double)nevt);
    G4double factorb= cross_sec*1000.0/(barn*(G4double)nevt);
    G4cout << "### factor  = " << factor
           << "### factora = " << factora
           << "### factorb = " << factorb
           << "    cross(mb)= " << cross_sec*1000./barn << G4endl;

    // -------- Limit number of angles

    if(nangl   > 13) nangl   = 13;
    if(nanglpr > 11) nanglpr = 11;
    if(nanglpi >  5) nanglpi = 5;

    if(nangl > 0) {
      for(G4int k=0; k<nangl; k++) {

        if(nangl == 1) {
          bng1[0] = std::max(0.0,ang[0] - dangl);
          bng2[0] = std::min(180., ang[0] + dangl);
        } else if(k == 0) {
          bng1[0] = std::max(0.0,ang[0] - dangl);
          bng2[0] = std::min(0.5*(ang[0] + ang[1]), ang[0] + dangl);
        } else if(k < nangl-1) {
          bng1[k] = std::max(bng2[k-1], ang[k]-dangl);
          bng2[k] = std::min(0.5*(ang[k] + ang[k+1]), ang[k] + dangl);
        } else {
          bng1[k] = std::max(bng2[k-1], ang[k]-dangl);
          bng2[k] = std::min(180., ang[k] + dangl);
        }

        cng[k] = cross_sec*MeV*1000.0*(G4double)nbinsd/
         (twopi*(std::cos(degree*bng1[k]) - std::cos(degree*bng2[k]))*
                barn*emax*(G4double)nevt);
      }
    }

    if(nanglpr > 0) {
      for(G4int k=0; k<nanglpr; k++) {

        if(nanglpr == 1) {
          bng1pr[0] = std::max(0.0,angpr[0] - dangl);
          bng2pr[0] = std::min(180., angpr[0] + dangl);
        } else if(k == 0) {
          bng1pr[0] = std::max(0.0,angpr[0] - dangl);
          bng2pr[0] = std::min(0.5*(angpr[0] + angpr[1]), angpr[0] + dangl);
        } else if(k < nanglpr-1) {
          bng1pr[k] = std::max(bng2pr[k-1], angpr[k]-dangl);
          bng2pr[k] = std::min(0.5*(angpr[k] + angpr[k+1]), angpr[k] + dangl);
        } else {
          bng1pr[k] = std::max(bng2pr[k-1], angpr[k]-dangl);
          bng2pr[k] = std::min(180., angpr[k] + dangl);
        }

        cngpr[k] = cross_sec*MeV*1000.0*(G4double)nbinsd/
         (twopi*(std::cos(degree*bng1pr[k]) - std::cos(degree*bng2pr[k]))*
                barn*emax*(G4double)nevt);
      }
    }

    if(nanglpi > 0) {
      for(G4int k=0; k<nanglpi; k++) {

        if(nangl == 1) {
          bngpi1[0] = std::max(0.0,angpi[0] - dangl);
          bngpi2[0] = std::min(180., angpi[0] + dangl);
        } else if(k == 0) {
          bngpi1[0] = std::max(0.0,angpi[0] - dangl);
          bngpi2[0] = std::min(0.5*(angpi[0] + angpi[1]), angpi[0] + dangl);
        } else if(k < nanglpi-1) {
          bngpi1[k] = std::max(bngpi2[k-1], angpi[k]-dangl);
          bngpi2[k] = std::min(0.5*(angpi[k] + angpi[k+1]), angpi[k] + dangl);
        } else {
          bngpi1[k] = std::max(bngpi2[k-1], angpi[k]-dangl);
          bngpi2[k] = std::min(180., angpi[k] + dangl);
        }

        cngpi[k] = cross_sec*MeV*1000.0*(G4double)nbinspi/
         (twopi*(std::cos(degree*bngpi1[k]) - std::cos(degree*bngpi2[k]))*
                 barn*emaxpi*(G4double)nevt);
      }
    }

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      { G4cout << "G4StateManager PROBLEM! " << G4endl; }
    G4Timer* timer = new G4Timer();
    timer->Start();
    const G4DynamicParticle* sec = nullptr;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4ThreeVector  momprim;
    std::vector<G4ThreeVector> momalpha;
    G4LorentzVector labv, fm;
    G4double e, p, px, py, pz, pt, theta;
    G4double Epair;
    G4ThreeVector momtmp;
    G4int warn = 0;

    // -------- Event loop
    momalpha.reserve(10);

    for (G4int iter=0; iter<nevt; ++iter) {
    
      momalpha.clear();
      
      if(verbose>=1 || iter == modu*(iter/modu)) { 
        G4cout << "###=== " << iter << "-th event start " << G4endl;
      }
      if(saverand) { CLHEP::HepRandom::saveEngineStatus("random.txt"); }

      G4double e0 = energy;
      if(sigmae > 0.0)  {
	do {
	  e0 = G4RandGauss::shoot(energy,sigmae);
	} while (e0 < 0.0);
      }
      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      G4double amass = 0.0;
      
      if(extraproc) {
        if(emcscat) {
	  // EM scattering
          gParticleChange->Initialize(*gTrack);
          emcscat->StartTracking(gTrack);
          G4double zzz = theStep;
	  G4double xxx = emcscat->ComputeTruePathLengthLimit(*gTrack, zzz);
          G4double yyy = emcscat->ComputeTrueStepLength(zzz);
	  if(verbose > 1) {
	    G4cout << "MSC: xxx= " << xxx << "  yyy= " << yyy << " zzz= " << zzz 
		   << " step= " << theStep << G4endl;
	  }
          emcscat->SampleScattering(aDirection, 0.0);
	  bDirection = *(gParticleChange->GetProposedMomentumDirection());
	}
	// elastic and change exchenage processes
	// for extra processes isotope is selected by the process
	aChange = extraproc->PostStepDoIt(*gTrack,*step); 
	amass = G4NucleiProperties::GetNuclearMass(A, Z);

      } else { 
	// test30 special process, isotope may be set in macro (Z and A)
	// or is sampled in GetNuclearMass method
        amass = phys->GetNucleusMass();
	aChange = proc->PostStepDoIt(*gTrack,*step); 
      }
      /*
      G4cout << "Mnuc= " << amass/GeV << " GeV"
	     << " From prop: " 
	     << G4NucleiProperties::GetNuclearMass(A, Z)/GeV << " GeV" 
	     << G4endl;
      */
      labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)), 
			     e0 + mass + amass);
      G4ThreeVector bst = labv.boostVector();

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de); 
      labv -= dee;
      //G4cout << " deltaE= " << de/MeV << G4endl;

      G4int n = aChange->GetNumberOfSecondaries();

      G4int nbar = 0;

      for(G4int j=0; j<n; ++j) {

        sec = aChange->GetSecondary(j)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        if(pd->GetPDGMass() > 100.*MeV) { ++nbar; }
      }

      for(G4int i=0; i<n; ++i) {

        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        fm  = sec->Get4Momentum();
        mom = sec->GetMomentum();

	// for exclusive reaction 2 particles in final state
        if(!inclusive && nbar != 2) { break; }

        G4double mass1 = pd->GetPDGMass();
	p = mom.mag();
        labv -= fm;

	// electron can come only from internal conversion
	// its mass should be added to initial state
        if(pd == electron) { 
	  //G4cout << "+e- " << G4endl;
	  labv += G4LorentzVector(0.0,0.0,0.0,electron_mass_c2); 
	}

        px = mom.x();
        py = mom.y();
        pz = mom.z();
        pt = std::sqrt(px*px +py*py);
        e  = fm.e() - mass1;

        theta = mom.theta();
        G4double cost  = std::cos(theta);
        G4double thetad = theta/degree;

        fm.boost(-bst);
        G4double tetcm  = fm.theta();
	//        G4double tetcmd = tetcm/degree;
        G4double costcm = std::cos(tetcm);

	if(usehis) {
          if(extraproc) {
	    if(i==0)  {
	      // fill recoil parameters
	      histo.Fill(62,e/MeV,1.0);
	      histo.Fill(58,cost,factora);
	      histo.Fill(60,costcm,factora);
	    }
	  }

          histo.Fill(2,mom.phi()/degree,1.0);
          if(pd == neutron) histo.Fill(23,mom.phi()/degree,1.0);
	}

	if( (e > e0 + MeV || e == 0.0 || pt == 0.0) && warn < 100 ) {
          warn++;
          G4cout << "Warning! evt# " << iter 
	         << "  " << i << "-th sec  "
		 << pd->GetParticleName() << "   Ekin(MeV)= "
                 << e/MeV
                 << " Pt(MeV/c)= " << pt/MeV
		 << " Ebeam(MeV)= " << e0/MeV
		 << G4endl;
	}
	de += e;
        if(verbose>1) {
          G4cout /*<< "Warning! evt# " << iter*/ 
                 << "  " << i << "-th sec  "
		 << pd->GetParticleName() << "  Ekin(MeV)= "
                 << e/MeV
		 << "  p(MeV)= " << mom/MeV
		 << "  m(MeV)= " << mass1/MeV
		 << "  Etot(MeV)= " << (e+mass1)/MeV
		 << "  pt(MeV)= " << pt/MeV
                 << "  std::sin(tet)= " << pt/p
                 << "  phi(deg)= " << mom.phi()/degree
                 << G4endl;
        }

	if(usehis) {

          if(pd) {
            G4double N = pd->GetBaryonNumber();
            G4double Z1= pd->GetPDGCharge()/eplus;
            G4double Z0= bestZ[(G4int)N];
            if(std::abs(Z0 - Z1) < 0.1 || Z0 == 0.0) 
	      histo.Fill(26, N, factorb);
	  }

          if(pd == proton) {

            if(rmsProton > 0.0) e += e*rmsProton*G4RandGauss::shoot(0.0,1.0);
            histo.Fill(1,1.0, 1.0);
	    histo.Fill(21,e/MeV, factor);
	    histo.Fill(24,cost, factora);
            for(G4int kk=0; kk<nanglpr; kk++) {
              if(bng1pr[kk] <= thetad && thetad <= bng2pr[kk]) {
                histo.Fill(3+kk,e/MeV, cngpr[kk]);
                break;
	      }
	    }

          } else if(pd == pin) {

            if(rmsPion > 0.0) e += e*rmsPion*G4RandGauss::shoot(0.0,1.0);
	    histo.Fill(1,4.0, 1.0);
            histo.Fill(20,e/MeV, 1.0);
            for(G4int kk=0; kk<nanglpi; kk++) {
              if(bngpi1[kk] <= thetad && thetad <= bngpi2[kk]) {
                histo.Fill(40+kk,e/MeV, cngpi[kk]);
                break;
	      }
	    }

          } else if(pd == pip) {

            if(rmsPion > 0.0) e += e*rmsPion*G4RandGauss::shoot(0.0,1.0);
	    histo.Fill(1,3.0, 1.0);
            histo.Fill(19,e/MeV, 1.0);
            for(G4int kk=0; kk<nanglpi; kk++) {
              if(bngpi1[kk] <= thetad && thetad <= bngpi2[kk]) {
                histo.Fill(45+kk,e/MeV, cngpi[kk]);
                break;
	      }
	    }

	  } else if(pd == pi0) {

	    histo.Fill(1,5.0, 1.0);
	    histo.Fill(18,e/MeV, 1.0);

	  } else if(pd == neutron) {

            if(rmsNeutron > 0.0) e += e*rmsNeutron*G4RandGauss::shoot(0.0,1.0);
	    histo.Fill(1,2.0, 1.0);
	    histo.Fill(22,e/MeV, factor);
            G4double ee = std::log10(e/MeV);
	    G4double e2 = ee;
            G4bool islog= false;
            if(ee >= logmin && ee <= logmax) { 
	      islog = true;
	      G4int nbb = G4int(((ee - logmin)/binlog));
              G4double e1 = logmin + binlog*nbb;
	      e2 = std::pow(10.,e1 + binlog) - std::pow(10.,e1);
	      histo.Fill(50, ee, factor*bine/e2);
	    } 
	    if(e >= elim) histo.Fill(25, cost, factora);
            for(G4int kk=0; kk<nangl; kk++) {
              if(bng1[kk] <= thetad && thetad <= bng2[kk]) {
                histo.Fill(27+kk,e/MeV, cng[kk]);
                if(islog && kk < 6) histo.Fill(51+kk,ee,cng[kk]*bind/e2);
                break;
	      }
	    }

	  } else if(pd == gamma) {
	    histo.Fill(14,e/MeV, 1.0);
	  } else if(pd == deu) {
	    histo.Fill(1,6.0, 1.0);
	  } else if(pd == tri) {
	    histo.Fill(1,7.0, 1.0);
	  } else if(pd == alp) {
	    histo.Fill(1,8.0, 1.0);
	    histo.Fill(66,e/MeV, 1.0);
	    momalpha.push_back(mom);
	  } else {
	    histo.Fill(1,9.0, 1.0);
	  }
	}
	
	// calculation of energy of the relative motion of two alpha particles
	G4int alpSizeM = (G4int) momalpha.size();
    	if (alpSizeM > 1) {
      	  for (G4int itr = 0; itr < alpSizeM-1; itr++) {
            for(G4int jtr = itr+1; jtr < alpSizeM; jtr++) {
              momtmp=momalpha[itr]-momalpha[jtr];
              Epair = (momtmp.mag2())/(4*massalp);
              histo.Fill(67, Epair/MeV, factor);
            }
       	  }
    	}	
    
	//	delete sec;       	 
        delete aChange->GetSecondary(i);
      }
      if(aChange->GetTrackStatus() == fAlive && extraproc) {
        G4double ekin;
        G4ThreeVector dir;
	G4ParticleChange* bChange = static_cast<G4ParticleChange*>(aChange);
	ekin = bChange->GetEnergy();
	dir = *(bChange->GetMomentumDirection());
        if(applyMSC) { dir.rotateUz(bDirection); }
	
        G4double mom1 = std::sqrt(ekin*(ekin + 2*mass));

        theta = dir.theta();
        G4double cost  = dir.z();
        G4double thetad = theta/degree;
        G4LorentzVector fm1(dir.x()*mom1,dir.y()*mom1,dir.z()*mom1,ekin + mass);

        labv -= fm1;

        fm1.boost(-bst);
        G4double tetcm  = fm1.theta();
        G4double tetcmd = tetcm/degree;
        G4double costcm = std::cos(tetcm);
       
	histo.Fill(57,ekin/MeV,1.0);
	histo.Fill(59,cost,factora);
	histo.Fill(61,costcm,factora);
	histo.Fill(63,thetad,factoraa/std::sin(theta));
	histo.Fill(64,std::log10(thetad),factoral*theta/std::sin(theta));
	histo.Fill(65,tetcmd,factoraa/std::sin(tetcm));
        if(verbose>1) {
          G4cout /*<< "Warning! evt# " << iter*/ 
                 << "primary  "
		 << namePart << "  Ekin(MeV)= "
                 << ekin/MeV
		 << "  p(MeV)= " << mom1/MeV
		 << "  m(MeV)= " << mass/MeV
		 << "  Etot(MeV)= " << (ekin+mass)/MeV
                 << "  std::sin(tet)= " << std::sin(theta)
                 << "  phi(deg)= " << dir.phi()/degree
                 << G4endl;
        }
      }

      px = labv.px();
      py = labv.py();
      pz = labv.pz();
      p  = std::sqrt(px*px +py*py + pz*pz);
      pt = std::sqrt(px*px +py*py);
      if(verbose > 1 || 
	 (verbose > 0 && (std::abs(labv.e()) > 0.1*MeV || p > 100*MeV))) {
        G4cout << iter << "-event Energy/Momentum balance= " << labv 
	       << " p=" << p << " pt= " << pt << G4endl;
      }

      if(usehis) {
        histo.Fill(0,(G4double)n,1.0);

        G4double ex = labv.e()/MeV;
        if(ex >= balance) { ex = balance - 0.0001; }
        else if(ex <= -balance) { ex = -balance + 0.0001; }
	histo.Fill(15,ex, 1.0);

        ex = pz/MeV;
        if(ex >= balance) { ex = balance - 0.0001; }
        else if(ex <= -balance) { ex = -balance + 0.0001; }
	histo.Fill(16,ex, 1.0);

        ex = pt/MeV;
        if(ex >= balance) { ex = balance - 0.0001; }
	histo.Fill(17,ex, 1.0);
      }
      aChange->Clear();
    }

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // -------- Committing the transaction with the tree

    if(usehis) {
      if(verbose > 0) { G4cout << "###### Save histograms" << G4endl; }
      histo.Save();
    }
    if(verbose > 0) {
      G4cout << "###### End of run # " << run << "     ######" << G4endl;
    }
    emcscat = nullptr;
    extraproc = nullptr;
  }

  delete pFrame;
  delete lFrame;
  delete sFrame;
  delete mate;
  delete fin;
  delete phys;
  partTable->DeleteAllParticles();

  G4cout << "###### End of test #####" << G4endl;
}
