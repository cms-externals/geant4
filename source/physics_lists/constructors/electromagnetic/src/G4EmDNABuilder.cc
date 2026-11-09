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
// Geant4 class G4EmDNABuilder
//
// Author V.Ivanchenko 22.05.2020
//

#include "G4EmDNABuilder.hh"

#include "G4SystemOfUnits.hh"

// particles
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4Deuteron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

// utilities
#include "G4EmBuilder.hh"
#include "G4EmParameters.hh"
#include "G4LowEnergyEmProcessSubType.hh"
#include "G4PhysListUtil.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"

// standard processes
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4NuclearStopping.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

// standard models
#include "G4BetheBlochModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BraggModel.hh"
#include "G4DummyModel.hh"
#include "G4Generator2BS.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LindhardSorensenIonModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4LowEWentzelVIModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4UrbanMscModel.hh"

// DNA models
#include "G4DNABornExcitationModel.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNABornIonisationModel1.hh"
#include "G4DNACPA100ElasticModel.hh"
#include "G4DNACPA100ExcitationModel.hh"
#include "G4DNACPA100IonisationModel.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"
#include "G4DNADingfelderChargeIncreaseModel.hh"
#include "G4DNADoubleIonisationModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"
#include "G4DNAGeneralIonIonisationModel.hh"
#include "G4DNAIonElasticModel.hh"
#include "G4DNAMeltonAttachmentModel.hh"
#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4DNAOneStepThermalizationModel.hh"
#include "G4DNAQuadrupleIonisationModel.hh"
#include "G4DNARPWBAExcitationModel.hh"
#include "G4DNARPWBAIonisationModel.hh"
#include "G4DNARuddIonisationDynamicModel.hh"
#include "G4DNARuddIonisationExtendedModel.hh"
#include "G4DNARuddIonisationModel.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNATripleIonisationModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"

// Energy limits of models:
//   1. lowEnergyRPWBA is the low energy limit of proton RPWA models
//      G4DNARPWBAExcitationModel and G4DNARPWBAIonisationModel.
//   2. lowEnergyBorn is a recommended low energy limit of proton
//      Born models G4DNABornExcitationModel and G4DNABornIonisationModel.
//      It is highEnergy limit for proton and hydrogen ion excitation model
//      G4DNAMillerGreenExcitationModel, for alpha ions this model is used
//      up to maximum ion energy.
//   3. lowEnergyBetheBloch is the low energy limit of standard
//      G4BetheBlochModel for protons, this model is used for all heavy
//      charged particles. For energies below G4BraggModel is used for
//      all hadrons. These standard models and transition between them
//      are needed for consistency of computation of range, which is needed
//      to compare standard/DNA models and for multiple scattering models.
//      The real transition energy from DNA to standard ionisation models
//      is above 100 MeV.
//   4. highEnergyElastic - is upper limit for DNA models of elastic
//      scattering. Above this energy multiple scattering models for electrons
//      protons, and ions are applied.
//   5. highEnergyChargeExchange - is the upper energy of charge increase
//      and decrease models for protons, neutral atoms, and light ions.
//   6. highEnergyCPA100 is the high energy limit of all CPA100 electron
//      models.
//   7. highEnergyEmfietzoglou is the upper energy limit of electron models
//      G4DNAEmfietzoglouExcitationModel and G4DNAEmfietzoglouIonisationModel.

namespace
{
constexpr G4double lowEnergyRPWBA = 100 * CLHEP::MeV;
constexpr G4double lowEnergyBorn = 0.5 * CLHEP::MeV;
constexpr G4double lowEnergyBetheBloch = 2 * CLHEP::MeV;
// SI: increased high energylimit for DNA Option4 only, temporary till G4 12.0
constexpr G4double highEnergyOpt4 = 10 * CLHEP::MeV;
//
constexpr G4double highEnergyElastic = 1 * CLHEP::MeV;
constexpr G4double highEnergyChargeExchange = 100 * CLHEP::MeV;
constexpr G4double highEnergyCPA100 = 250 * CLHEP::keV;
};  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructDNAParticles()
{
  // standard particles
  G4EmBuilder::ConstructMinimalEmSet();

  // DNA ions
  G4DNAGenericIonsManager* genericIonsManager = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructStandardEmPhysics(const G4double emin_elec,
                                                const G4double emin_proton, const G4double emin_ion,
                                                const G4EmDNAMscModelType mscType, const G4bool)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();
  const G4double emax = param->MaxKinEnergy();
  G4EmBuilder::PrepareEMPhysics();

  // gamma
  G4ParticleDefinition* gamma = G4Gamma::Gamma();

  // photoelectric effect - Livermore model
  auto thePEEffect = new G4PhotoElectricEffect();
  thePEEffect->SetEmModel(new G4LivermorePhotoElectricModel());
  ph->RegisterProcess(thePEEffect, gamma);

  // Compton scattering - Klein-Nishina
  auto theComptonScattering = new G4ComptonScattering();
  theComptonScattering->SetEmModel(new G4KleinNishinaModel());
  auto cModel = new G4LowEPComptonModel();
  cModel->SetHighEnergyLimit(20 * CLHEP::MeV);
  theComptonScattering->AddEmModel(0, cModel);
  ph->RegisterProcess(theComptonScattering, gamma);

  // gamma conversion - 5D model
  auto theGammaConversion = new G4GammaConversion();
  ph->RegisterProcess(theGammaConversion, gamma);

  // Rayleigh scattering - Livermore model
  auto theRayleigh = new G4RayleighScattering();
  ph->RegisterProcess(theRayleigh, gamma);

  // electron
  if (emin_elec < emax)
  {
    G4ParticleDefinition* elec = G4Electron::Electron();
    auto msc_el = new G4eMultipleScattering();
    G4VMscModel* msc_model_el;
    if (mscType == dnaWVI)
    {
      msc_model_el = new G4LowEWentzelVIModel();
    }
    else if (mscType == dnaGS)
    {
      msc_model_el = new G4GoudsmitSaundersonMscModel();
    }
    else
    {
      msc_model_el = new G4UrbanMscModel();
    }
    msc_model_el->SetActivationLowEnergyLimit(highEnergyElastic);

    msc_el->SetEmModel(msc_model_el);
    ph->RegisterProcess(msc_el, elec);

    auto ioni = new G4eIonisation();
    auto mb_el = new G4MollerBhabhaModel();
    mb_el->SetActivationLowEnergyLimit(emin_elec);
    ioni->SetEmModel(mb_el);
    ph->RegisterProcess(ioni, elec);

    auto brem = new G4eBremsstrahlung();
    auto sb_el = new G4SeltzerBergerModel();
    sb_el->SetActivationLowEnergyLimit(emin_elec);
    sb_el->SetHighEnergyLimit(emax);
    sb_el->SetAngularDistribution(new G4Generator2BS());
    brem->SetEmModel(sb_el);
    ph->RegisterProcess(brem, elec);
  }

  // positron
  G4ParticleDefinition* posi = G4Positron::Positron();
  auto msc_pos = new G4eMultipleScattering();
  G4VMscModel* msc_model_pos;
  if (mscType == dnaWVI)
  {
    msc_model_pos = new G4LowEWentzelVIModel();
  }
  else if (mscType == dnaGS)
  {
    msc_model_pos = new G4GoudsmitSaundersonMscModel();
  }
  else
  {
    msc_model_pos = new G4UrbanMscModel();
  }
  msc_pos->SetEmModel(msc_model_pos);
  ph->RegisterProcess(msc_pos, posi);
  ph->RegisterProcess(new G4eIonisation(), posi);

  auto brem = new G4eBremsstrahlung();
  auto sb = new G4SeltzerBergerModel();
  sb->SetHighEnergyLimit(emax);
  sb->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(sb);
  ph->RegisterProcess(brem, posi);
  ph->RegisterProcess(new G4eplusAnnihilation(), posi);

  // proton
  if (emin_proton < emax)
  {
    G4ParticleDefinition* part = G4Proton::Proton();
    StandardHadronPhysics(part, highEnergyElastic, emin_proton, emax, mscType, false);
  }

  // GenericIon
  if (emin_ion < emax)
  {
    G4ParticleDefinition* ion = G4GenericIon::GenericIon();
    StandardHadronPhysics(ion, 0.0, emin_ion, emax, dnaUrban, true);
    // deuteron
    ion = G4Deuteron::Deuteron();
    StandardHadronPhysics(ion, highEnergyElastic, emin_ion, emax, dnaUrban, false);
    // alpha
    G4ParticleDefinition* part = G4Alpha::Alpha();
    StandardHadronPhysics(part, highEnergyElastic, emin_ion, emax, dnaUrban, true);

    // alpha+
    G4DNAGenericIonsManager* genericIonsManager = G4DNAGenericIonsManager::Instance();
    part = genericIonsManager->GetIon("alpha+");
    StandardHadronPhysics(part, highEnergyElastic, emin_ion, emax, dnaUrban, false);
  }
  // list of main standard charged particles, which have no DNA models
  const std::vector<G4int> chargedParticles = {13,   -13,   211,        -211,      321,
                                               -321, -2212, 1000010030, 1000020030};
  auto msc = new G4hMultipleScattering();
  msc->SetEmModel(new G4WentzelVIModel());
  G4EmBuilder::ConstructBasicEmPhysics(msc, chargedParticles);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::StandardHadronPhysics(G4ParticleDefinition* part,
                                           const G4double lowELimitForMSC,
                                           const G4double lowELimitForIoni,
                                           const G4double maxEnergy,
                                           const G4EmDNAMscModelType mscType, const G4bool isIon)
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4hMultipleScattering* msc = new G4hMultipleScattering();
  G4VMscModel* msc_model = nullptr;
  if (mscType == dnaWVI)
  {
    msc_model = new G4LowEWentzelVIModel();
  }
  else
  {
    msc_model = new G4UrbanMscModel();
  }
  msc_model->SetActivationLowEnergyLimit(lowELimitForMSC);
  msc_model->SetLowEnergyLimit(lowELimitForMSC);
  msc_model->SetHighEnergyLimit(maxEnergy);
  msc->SetEmModel(msc_model);
  ph->RegisterProcess(msc, part);

  G4VEnergyLossProcess* ioni = nullptr;
  G4VEmModel* mod1 = nullptr;
  G4VEmModel* mod2 = nullptr;
  G4int pdg = part->GetPDGEncoding();
  G4double emin = lowELimitForIoni * part->GetPDGMass() / CLHEP::proton_mass_c2;

  // Low energy models are needed only to compute correctly range
  // tables for standard models. DNA models are used for sampling
  // of final state at low energies.
  G4double mRatio = 1.0;
  if (isIon)
  {
    ioni = new G4ionIonisation();
    if (pdg == 1000020040)
    {
      mod1 = new G4BraggIonModel();
      mod2 = new G4BetheBlochModel();
      mRatio = part->GetPDGMass() / CLHEP::proton_mass_c2;
    }
    else
    {
      mod1 = new G4LindhardSorensenIonModel();
    }
  }
  else
  {
    ioni = new G4hIonisation();
    if (pdg != 2212)
    {
      ioni->SetBaseParticle(G4Proton::Proton());
      mRatio = part->GetPDGMass() / CLHEP::proton_mass_c2;
    }
    mod1 = new G4BraggModel();
    mod2 = new G4BetheBlochModel();
  }
  ioni->SetEmModel(mod1);
  if (nullptr == mod2)
  {
    mod1->SetActivationLowEnergyLimit(emin);
  }
  else
  {
    G4double eth = lowEnergyBetheBloch * mRatio;
    mod1->SetHighEnergyLimit(eth);
    mod1->SetActivationLowEnergyLimit(emin);
    mod2->SetActivationLowEnergyLimit(emin);
    mod2->SetLowEnergyLimit(eth);
    ioni->SetEmModel(mod2);
  }
  ph->RegisterProcess(ioni, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructDNAElectronPhysics(const G4double emaxDNA, const G4int opt,
                                                 const G4bool fast, const G4bool stationary,
                                                 const G4Region* reg)
{
  G4ParticleDefinition* part = G4Electron::Electron();
  G4EmParameters* param = G4EmParameters::Instance();

  // old limit of the Emfietzoglou models
  // G4double emaxE = (4 == opt) ? highEnergyEmfietzoglou : 0.0;
  // limit of the elastic and solvation models
  G4double emaxT = param->MinDNAElectronEnergy();

  // *** Solvation ***
  G4DNAElectronSolvation* pSolvation = FindOrBuildElectronSolvation();
  auto therm = G4DNASolvationModelFactory::GetMacroDefinedModel();
  therm->SetHighEnergyLimit(emaxT);
  pSolvation->AddEmModel(-1, therm, reg);

  // *** Elastic scattering ***
  auto pElasticProcess = FindOrBuildElastic(part, "e-_G4DNAElastic");
  G4VEmModel* elast;
  G4VEmModel* elast2 = nullptr;
  if (4 == opt)
  {
    // SI: choice of default elastic model for option4
    elast = new G4DNAUeharaScreenedRutherfordElasticModel();
    // elast = new G4DNAELSEPAElasticModel();
    //
  }
  else if (6 <= opt)
  {
    auto mod = new G4DNACPA100ElasticModel();
    mod->SelectStationary(stationary);
    elast = mod;
    elast2 = new G4DNAChampionElasticModel();
  }
  else
  {
    elast = new G4DNAChampionElasticModel();
  }
  elast->SetLowEnergyLimit(emaxT);
  // SI temporary (1/3)
  // elast->SetHighEnergyLimit(highEnergyElastic);
  elast->SetHighEnergyLimit((opt == 4) ? highEnergyOpt4 : highEnergyElastic);
  //

  pElasticProcess->AddEmModel(-2, elast, reg);

  if (nullptr != elast2)
  {
    elast->SetHighEnergyLimit(highEnergyCPA100);
    elast2->SetLowEnergyLimit(highEnergyCPA100);
    elast2->SetHighEnergyLimit(highEnergyElastic);
    pElasticProcess->AddEmModel(-3, elast2, reg);
  }

  // BEGIN SI
  // it is possible to have 1 or 2 DNA models
  // Emfietzoglou models up to 10 MeV for opt4

  // *** Excitation ***
  auto theDNAExc = FindOrBuildExcitation(part, "e-_G4DNAExcitation");
  G4VEmModel* modB;
  G4VEmModel* modB2 = nullptr;
  if (4 == opt)
  {
    auto modE = new G4DNAEmfietzoglouExcitationModel();
    modE->SelectStationary(stationary);
    modB = modE;
  }
  else if (6 == opt)
  {
    auto mod = new G4DNACPA100ExcitationModel();
    mod->SelectStationary(stationary);
    modB = mod;
    auto mod1 = new G4DNABornExcitationModel();
    mod1->SelectStationary(stationary);
    modB2 = mod1;
  }
  else
  {
    auto mod = new G4DNABornExcitationModel();
    mod->SelectStationary(stationary);
    modB = mod;
  }
  theDNAExc->AddEmModel(-2, modB, reg);
  if (nullptr != modB2)
  {
    modB->SetHighEnergyLimit(highEnergyCPA100);
    modB2->SetLowEnergyLimit(highEnergyCPA100);
    modB2->SetHighEnergyLimit(emaxDNA);
    theDNAExc->AddEmModel(-3, modB2, reg);
  }
  else
  {
    modB->SetHighEnergyLimit(emaxDNA);
  }

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, "e-_G4DNAIonisation");
  G4VEmModel* modI;
  G4VEmModel* modI2 = nullptr;
  if (4 == opt)
  {
    auto modE = new G4DNAEmfietzoglouIonisationModel();
    modE->SelectFasterComputation(fast);
    modE->SelectStationary(stationary);
    modI = modE;
  }
  else if (6 == opt)
  {
    auto mod = new G4DNACPA100IonisationModel();
    mod->SelectStationary(stationary);
    mod->SelectFasterComputation(fast);
    modI = mod;
    modI2 = new G4DNABornIonisationModel1();
  }
  else
  {
    modI = new G4DNABornIonisationModel1();
  }
  theDNAIoni->AddEmModel(-2, modI, reg);
  if (nullptr != modI2)
  {
    modI->SetHighEnergyLimit(highEnergyCPA100);
    modI2->SetLowEnergyLimit(highEnergyCPA100);
    modI2->SetHighEnergyLimit(emaxDNA);
    theDNAIoni->AddEmModel(-3, modI2, reg);
  }
  else
  {
    modI->SetHighEnergyLimit(emaxDNA);
  }
  // END SI

  if (4 != opt && 6 != opt)
  {
    // *** Vibrational excitation ***
    auto theDNAVibExc = FindOrBuildVibExcitation(part, "e-_G4DNAVibExcitation");
    auto modS = new G4DNASancheExcitationModel();
    theDNAVibExc->AddEmModel(-1, modS, reg);
    modS->SelectStationary(stationary);

    // *** Attachment ***
    auto theDNAAttach = FindOrBuildAttachment(part, "e-_G4DNAAttachment");
    auto modM = new G4DNAMeltonAttachmentModel();
    theDNAAttach->AddEmModel(-1, modM, reg);
    modM->SelectStationary(stationary);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructDNAProtonPhysics(const G4double emaxProtonDNA, const G4int opt,
                                               const G4bool fast, const G4bool stationary,
                                               const G4Region* reg)
{
  G4EmParameters* param = G4EmParameters::Instance();
  G4ParticleDefinition* part = G4Proton::Proton();
  G4double e2DNA = (8 == opt) ? lowEnergyRPWBA : lowEnergyBorn;

  // *** Elastic scattering ***
  auto pElasticProcess = FindOrBuildElastic(part, "proton_G4DNAElastic");
  auto modE = new G4DNAIonElasticModel();
  modE->SetHighEnergyLimit(highEnergyElastic);
  modE->SelectStationary(stationary);
  pElasticProcess->AddEmModel(-1, modE, reg);

  // *** Excitation ***
  auto theDNAExc = FindOrBuildExcitation(part, "proton_G4DNAExcitation");
  auto modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SetHighEnergyLimit(lowEnergyBorn);
  modMGE->SelectStationary(stationary);
  theDNAExc->AddEmModel(-1, modMGE, reg);

  auto modB = new G4DNABornExcitationModel();
  modB->SelectStationary(stationary);
  modB->SetLowEnergyLimit(lowEnergyBorn);
  modB->SetHighEnergyLimit(lowEnergyRPWBA);
  theDNAExc->AddEmModel(-2, modB, reg);

  if (lowEnergyRPWBA < emaxProtonDNA)
  {
    auto modC = new G4DNARPWBAExcitationModel();
    modC->SelectStationary(stationary);
    modC->SetLowEnergyLimit(lowEnergyRPWBA);
    modC->SetHighEnergyLimit(emaxProtonDNA);
    theDNAExc->AddEmModel(-3, modC, reg);
  }

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, "proton_G4DNAIonisation");
  G4VEmModel* modRI;
  if (8 == opt)
  {
    modRI = new G4DNARuddIonisationDynamicModel();
  }
  else
  {
    modRI = new G4DNARuddIonisationExtendedModel();
  }
  modRI->SetHighEnergyLimit(e2DNA);
  theDNAIoni->AddEmModel(-1, modRI, reg);

  if (e2DNA < lowEnergyRPWBA)
  {
    G4VEmModel* modI = new G4DNABornIonisationModel1();
    modI->SetLowEnergyLimit(e2DNA);
    modI->SetHighEnergyLimit(lowEnergyRPWBA);
    theDNAIoni->AddEmModel(-2, modI, reg);
  }
  if (lowEnergyRPWBA < emaxProtonDNA)
  {
    auto modJ = new G4DNARPWBAIonisationModel();
    modJ->SelectFasterComputation(fast);
    modJ->SelectStationary(stationary);
    modJ->SetLowEnergyLimit(lowEnergyRPWBA);
    modJ->SetHighEnergyLimit(emaxProtonDNA);
    theDNAIoni->AddEmModel(-3, modJ, reg);
  }

  if (param->DNAMultipleIonisation())
  {
    // *** Double Ionisation ***
    FindOrBuildDoubleIonisation(part, "proton_G4DNADoubleIonisation", reg);

    // *** Triple Ionisation ***
    FindOrBuildTripleIonisation(part, "proton_G4DNATripleIonisation", reg);

    // *** Quadruple Ionisation ***
    FindOrBuildQuadrupleIonisation(part, "proton_G4DNAQuadrupleIonisation", reg);
  }

  // *** Charge decrease ***
  auto theDNAChargeDecreaseProcess = FindOrBuildChargeDecrease(part, "proton_G4DNAChargeDecrease");
  auto modDCD = new G4DNADingfelderChargeDecreaseModel();
  modDCD->SelectStationary(stationary);
  modDCD->SetHighEnergyLimit(highEnergyChargeExchange);
  theDNAChargeDecreaseProcess->AddEmModel(-1, modDCD, reg);

  // *** Tracking cut ***
  G4double cut = param->MinDNAProtonEnergy();
  FindOrBuildCapture(cut, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructDNAIonPhysics(const G4double emaxIonDNA, const G4int,
                                            const G4Region* reg)
{
  G4ParticleDefinition* part = G4GenericIon::GenericIon();
  G4EmParameters* param = G4EmParameters::Instance();

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, "GenericIon_G4DNAIonisation");
  G4VEmModel* mod = new G4DNAGeneralIonIonisationModel();
  mod->SetHighEnergyLimit(emaxIonDNA);
  theDNAIoni->AddEmModel(-1, mod, reg);

  if (param->DNAMultipleIonisation())
  {
    // *** Double Ionisation ***
    FindOrBuildDoubleIonisation(part, "GenericIon_G4DNADoubleIonisation", reg);

    // *** Triple Ionisation ***
    FindOrBuildTripleIonisation(part, "GenericIon_G4DNATripleIonisation", reg);

    // *** Quadruple Ionisation ***
    FindOrBuildQuadrupleIonisation(part, "GenericIon_G4DNAQuadrupleIonisation", reg);
  }

  // *** NIEL ***
  FindOrBuildNuclearStopping(part, CLHEP::MeV);

  // *** Tracking cut ***
  G4double cut = param->MinDNAIonEnergy();
  FindOrBuildCapture(cut, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::ConstructDNALightIonPhysics(G4ParticleDefinition* part, const G4int charge,
                                                 const G4int opt, const G4double emaxIonDNA,
                                                 const G4bool, const G4bool stationary,
                                                 const G4Region* reg)
{
  G4EmParameters* param = G4EmParameters::Instance();
  const G4String& name = part->GetParticleName();
  G4double emax = emaxIonDNA * part->GetPDGMass() / CLHEP::proton_mass_c2;
  G4double emax1 = ("hydrogen" == name) ? lowEnergyBorn : emax;

  // *** Elastic ***
  auto theDNAElastic = FindOrBuildElastic(part, name + "_G4DNAElastic");
  auto modEI = new G4DNAIonElasticModel();
  modEI->SelectStationary(stationary);
  modEI->SetHighEnergyLimit(highEnergyElastic);
  theDNAElastic->AddEmModel(-1, modEI, reg);

  // *** Excitation ***
  auto theDNAExc = FindOrBuildExcitation(part, name + "_G4DNAExcitation");
  auto modMGE = new G4DNAMillerGreenExcitationModel();
  modMGE->SelectStationary(stationary);
  modMGE->SetHighEnergyLimit(emax1);
  theDNAExc->AddEmModel(-1, modMGE, reg);

  // *** Ionisation ***
  auto theDNAIoni = FindOrBuildIonisation(part, name + "_G4DNAIonisation");
  G4VEmModel* modRI;
  if (8 == opt && name != "deuteron")
  {
    modRI = new G4DNARuddIonisationDynamicModel();
  }
  else
  {
    modRI = new G4DNARuddIonisationExtendedModel();
  }
  modRI->SetHighEnergyLimit(emax);
  theDNAIoni->AddEmModel(-2, modRI, reg);

  if (G4EmParameters::Instance()->DNAMultipleIonisation())
  {
    // *** Double Ionisation ***
    FindOrBuildDoubleIonisation(part, name + "_G4DNADoubleIonisation", reg);

    // *** Triple Ionisation ***
    FindOrBuildTripleIonisation(part, name + "_G4DNATripleIonisation", reg);

    // *** Quadruple Ionisation ***
    FindOrBuildQuadrupleIonisation(part, name + "_G4DNAQuadrupleIonisation", reg);
  }

  // *** Charge increase ***
  if (2 > charge && name != "deuteron")
  {
    auto theDNAChargeIncrease = FindOrBuildChargeIncrease(part, name + "_G4DNAChargeIncrease");
    auto modDCI = new G4DNADingfelderChargeIncreaseModel();
    modDCI->SelectStationary(stationary);
    modDCI->SetHighEnergyLimit(highEnergyChargeExchange);
    theDNAChargeIncrease->AddEmModel(-1, modDCI, reg);
  }

  // *** Charge decrease ***
  if (0 < charge && name != "deuteron")
  {
    auto theDNAChargeDecrease = FindOrBuildChargeDecrease(part, name + "_G4DNAChargeDecrease");
    auto modDCD = new G4DNADingfelderChargeDecreaseModel();
    modDCD->SelectStationary(stationary);
    modDCD->SetHighEnergyLimit(highEnergyChargeExchange);
    theDNAChargeDecrease->AddEmModel(-1, modDCD, reg);
  }

  // *** Tracking cut ***
  G4double cut = param->MinDNAProtonEnergy();
  FindOrBuildCapture(cut, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAElectronSolvation* G4EmDNABuilder::FindOrBuildElectronSolvation()
{
  auto elec = G4Electron::Electron();
  auto* p = G4PhysListUtil::FindProcess(elec, fLowEnergyElectronSolvation);
  G4DNAElectronSolvation* ptr = dynamic_cast<G4DNAElectronSolvation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, elec);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAElastic* G4EmDNABuilder::FindOrBuildElastic(G4ParticleDefinition* part, const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyElastic);
  G4DNAElastic* ptr = dynamic_cast<G4DNAElastic*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAElastic(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAExcitation* G4EmDNABuilder::FindOrBuildExcitation(G4ParticleDefinition* part,
                                                       const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyExcitation);
  G4DNAExcitation* ptr = dynamic_cast<G4DNAExcitation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAExcitation(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAVibExcitation* G4EmDNABuilder::FindOrBuildVibExcitation(G4ParticleDefinition* part,
                                                             const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyVibrationalExcitation);
  G4DNAVibExcitation* ptr = dynamic_cast<G4DNAVibExcitation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAVibExcitation(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAIonisation* G4EmDNABuilder::FindOrBuildIonisation(G4ParticleDefinition* part,
                                                       const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyIonisation);
  G4DNAIonisation* ptr = dynamic_cast<G4DNAIonisation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAIonisation(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAAttachment* G4EmDNABuilder::FindOrBuildAttachment(G4ParticleDefinition* part,
                                                       const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyAttachment);
  G4DNAAttachment* ptr = dynamic_cast<G4DNAAttachment*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAAttachment(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAChargeDecrease* G4EmDNABuilder::FindOrBuildChargeDecrease(G4ParticleDefinition* part,
                                                               const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyChargeDecrease);
  G4DNAChargeDecrease* ptr = dynamic_cast<G4DNAChargeDecrease*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAChargeDecrease(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNAChargeIncrease* G4EmDNABuilder::FindOrBuildChargeIncrease(G4ParticleDefinition* part,
                                                               const G4String& name)
{
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyChargeIncrease);
  G4DNAChargeIncrease* ptr = dynamic_cast<G4DNAChargeIncrease*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAChargeIncrease(name);
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LowECapture* G4EmDNABuilder::FindOrBuildCapture(const G4double elim, G4ParticleDefinition* part)
{
  auto p = G4PhysListUtil::FindProcess(part, -1);
  G4LowECapture* ptr = dynamic_cast<G4LowECapture*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4LowECapture(elim);
    auto mng = part->GetProcessManager();
    mng->AddDiscreteProcess(ptr);
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::FindOrBuildNuclearStopping(G4ParticleDefinition* part, const G4double elim)
{
  auto p = G4PhysListUtil::FindProcess(part, fNuclearStopping);
  auto ptr = dynamic_cast<G4NuclearStopping*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4NuclearStopping();
  }
  ptr->SetMaxKinEnergy(elim);
  auto ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(ptr, part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::FindOrBuildDoubleIonisation(G4ParticleDefinition* part, const G4String& name,
                                                 const G4Region* reg)
{
  auto kind = part->GetParticleName();
  if (kind != "proton" && kind != "alpha" && kind != "GenericIon")
  {
    return;
  }
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyDoubleIonisation);
  G4DNADoubleIonisation* ptr = dynamic_cast<G4DNADoubleIonisation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNADoubleIonisation(name);
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  ptr->AddEmModel(-1, new G4DNADoubleIonisationModel(), reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::FindOrBuildTripleIonisation(G4ParticleDefinition* part, const G4String& name,
                                                 const G4Region* reg)
{
  auto kind = part->GetParticleName();
  if (kind != "proton" && kind != "alpha" && kind != "GenericIon")
  {
    return;
  }
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyTripleIonisation);
  G4DNATripleIonisation* ptr = dynamic_cast<G4DNATripleIonisation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNATripleIonisation(name);
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  ptr->AddEmModel(-1, new G4DNATripleIonisationModel(), reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNABuilder::FindOrBuildQuadrupleIonisation(G4ParticleDefinition* part,
                                                    const G4String& name, const G4Region* reg)
{
  auto kind = part->GetParticleName();
  if (kind != "proton" && kind != "alpha" && kind != "GenericIon")
  {
    return;
  }
  auto p = G4PhysListUtil::FindProcess(part, fLowEnergyQuadrupleIonisation);
  G4DNAQuadrupleIonisation* ptr = dynamic_cast<G4DNAQuadrupleIonisation*>(p);
  if (nullptr == ptr)
  {
    ptr = new G4DNAQuadrupleIonisation(name);
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(ptr, part);
    ptr->SetEmModel(new G4DummyModel());
  }
  ptr->AddEmModel(-1, new G4DNAQuadrupleIonisationModel(), reg);
}
