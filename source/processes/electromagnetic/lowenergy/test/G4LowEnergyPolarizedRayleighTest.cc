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
// --------------------------------------------------------------
//
// File name:     G4LowEnergyPolarizedRayleighTest.cc
//
// Author:        Capra Riccardo
//
// Creation date: May 2005
//
// History:
// -----------
// 03 May 2005  R. Capra         1st implementation
//
//----------------------------------------------------------------

//! \file    G4LowEnergyPolarizedRayleighTest.cc
//! \brief   Tests G4LowEnergyPolarizedRayleigh process
//! \author  Capra Riccardo
//! \date    May 2005
//! \par     History:
//! <TABLE>
//!  <TR><TD> 03 May 2005 </TD><TD> R. Capra	</TD><TD> 1<SUP>st</SUP> implementation </TD></TR>
//! </TABLE>
//! \sa      G4LowEnergyPolarizedRayleigh.hh         

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <memory>
#include <cstdlib>

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4VDiscreteProcess.hh"
#include "G4VLowEnergyDiscretePhotonProcess.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"

#include "G4LowEnergyPolarizedRayleigh.hh"
#include "G4LowEnergyRayleigh.hh"

#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ForceCondition.hh"

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4RunManager.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"

#include "G4UnitsTable.hh"

#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"

//! \brief Options structure
struct Options
{
 //! \brief Mean free path test
 bool meanFreePathTest;
 //! \brief Post step do it test
 bool postStepDoItTest;
 //! \brief Post step do it test
 bool randomEnergy;
 //! \brief Output file name
 const char *outputFileName;
 //! \brief Material name
 const char *material;
 //! \brief Process name
 const char *process;
 //! \brief Minimum energy
 G4double minEnergy;
 //! \brief Maximum energy
 G4double maxEnergy;
 //! \brief Number of energy step
 G4int nEnergySteps;
 //! \brief Number of interactions
 G4int nIterations;
};

//! \brief Default output file name
struct Options defaultOptions = { false, false, false, "G4LowEnergyPolarizedRayleighTest.hbook", "Iron", "polarLowEnRayleigh", 100*eV, 100*GeV, 300, 1 };

//! \brief Creates some materials
void CreateMaterials(void)
{
 G4Element * H       = new G4Element ("Hydrogen", "H",   1.,   1.01*g/mole);
 G4Element * O       = new G4Element ("Oxygen",   "O",   8.,  16.00*g/mole);
 G4Element * C       = new G4Element ("Carbon",   "C",   6.,  12.00*g/mole);
 G4Element * Cs      = new G4Element ("Cesium",   "Cs", 55., 132.905*g/mole);
 G4Element * I       = new G4Element ("Iodine",   "I",  53., 126.9044*g/mole);

 G4Material * Si     = new G4Material("Silicon",   14., 28.055*g/mole,  2.33*g/cm3);
 G4Material * Fe     = new G4Material("Iron",      26.,  55.85*g/mole,  7.87*g/cm3);
 G4Material * Cu     = new G4Material("Copper",    29.,  63.55*g/mole,  8.96*g/cm3);
 G4Material * W      = new G4Material("Tungsten",  74., 183.85*g/mole, 19.30*g/cm3);
 G4Material * Pb     = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
 G4Material * U      = new G4Material("Uranium",   92., 238.03*g/mole, 18.95*g/cm3);
 G4Material * maO    = new G4Material("Oxygen",     8.,  16.00*g/mole,   1.1*g/cm3);
 G4Material * water  = new G4Material ("Water",  1.*g/cm3,     2);
 water->AddElement(H, 2);
 water->AddElement(O, 1);

 G4Material* ethane = new G4Material ("Ethane", 0.4241*g/cm3, 2);
 ethane->AddElement(H, 6);
 ethane->AddElement(C, 2);

 G4Material* csI    = new G4Material ("CsI",    4.53*g/cm3,   2);
 csI->AddElement(Cs, 1);
 csI->AddElement(I, 1);
 
 // This is needed to suppress some warnings. These lines can be deleted;
 Si->GetTemperature();
 Fe->GetTemperature();
 Cu->GetTemperature();
 W->GetTemperature();
 Pb->GetTemperature();
 U->GetTemperature();
 maO->GetTemperature();
 water->GetTemperature();
 ethane->GetTemperature();
 csI->GetTemperature();
}

//! \brief Process the options arguments
//! \param argc Number of arguments
//! \param argv Pointer to the arguments
//! \param options Structure to fill-in
void processOptions(int argc, char ** argv, struct Options * options)
{
 options->meanFreePathTest = defaultOptions.meanFreePathTest;
 options->postStepDoItTest = defaultOptions.meanFreePathTest;
 options->randomEnergy     = defaultOptions.randomEnergy;
 options->outputFileName   = defaultOptions.outputFileName;
 options->material         = defaultOptions.material;
 options->process          = defaultOptions.process;
 options->minEnergy        = defaultOptions.minEnergy;
 options->maxEnergy        = defaultOptions.maxEnergy;
 options->nEnergySteps     = defaultOptions.nEnergySteps;
 options->nIterations      = defaultOptions.nIterations;
 
 int i(1);
 
 while (i<argc)
 {
  if (argv[i][0]=='-' && argv[i][2]==0)
  {
   switch(argv[i][1])
   {
    case 'h':
    case '?':
     G4cout << argv[0] << " [-h|-?] [-a] [-b] [-r] [-o <file name>] [-m <material name>] [-p <process name>] [-e <min energy in MeV>] [-E <max energy in MeV>] [-s <energy steps>] [-n <iterations>] " << G4endl
            << G4endl
            << "-h|-?     Shows this help" << G4endl
            << "-a        Enables mean free path test" << G4endl
            << "-b        Enables post step do it test" << G4endl
            << "-r        Energy is choosen at random within the range" << G4endl
            << "-o <arg>  Set the output file name (default: \"" << defaultOptions.outputFileName << "\")" << G4endl
            << "-m <arg>  Set the material (default: \"" << defaultOptions.material << "\")" << G4endl
            << "-p <arg>  Set the process (default: \"" << defaultOptions.process << "\")" << G4endl
            << "-e <arg>  Set the low energy range in MeV (default: " << defaultOptions.minEnergy/MeV << " MeV)" << G4endl
            << "-E <arg>  Set the high energy range in MeV (default: " << defaultOptions.maxEnergy/MeV << " MeV)" << G4endl
            << "-s <arg>  Set the energy range step (default: " << defaultOptions.nEnergySteps << ")"<< G4endl
            << "-n <arg>  Set the number of iterations for the post step do it (default: " << defaultOptions.nIterations << ")" << G4endl;
     exit(0);
     break;
     
    case 'a':
     options->meanFreePathTest=true;
     break;
     
    case 'b':
     options->postStepDoItTest=true;
     break;
     
    case 'r':
     options->randomEnergy=true;
     break;
     
    case 'o':
     i++;
     if (i<argc)
     {
      options->outputFileName = argv[i];
      break;
     }

    case 'm':
     i++;
     if (i<argc)
     {
      options->material       = argv[i];
      break;
     }
     
    case 'p':
     i++;
     if (i<argc)
     {
      options->process        = argv[i];
      break;
     }

    case 'e':
     i++;
     if (i<argc)
     {
      options->minEnergy      = std::atof(argv[i])*MeV;
      if (options->minEnergy <= 0.)
      {
       G4cout << argv[0] << ": Energy must be > 0." << G4endl;
       exit(-1);
      }

      break;
     }

    case 'E':
     i++;
     if (i<argc)
     {
      options->maxEnergy      = std::atof(argv[i])*MeV;
      if (options->maxEnergy <= 0.)
      {
       G4cout << argv[0] << ": Energy must be > 0." << G4endl;
       exit(-1);
      }

      break;
     }

    case 's':
     i++;
     if (i<argc)
     {
      options->nEnergySteps   = atoi(argv[i]);
      if (options->nEnergySteps <= 1)
      {
       G4cout << argv[0] << ": Expected at least two steps." << G4endl;
       exit(-1);
      }
      
      break;
     }
     
     G4cout << argv[0] << ": Expected one more parameter in " << argv[i] << " option. Use -h option for help." << G4endl;
     exit(-1);
     
    case 'n':
     i++;
     if (i<argc)
     {
      options->nIterations    = atoi(argv[i]);
      if (options->nIterations <= 0)
      {
       G4cout << argv[0] << ": Expected at least one iteration." << G4endl;
       exit(-1);
      }
      
      break;
     }
     
     G4cout << argv[0] << ": Expected one more parameter in " << argv[i] << " option. Use -h option for help." << G4endl;
     exit(-1);
     
    default:
     G4cout << argv[0] << ": Unknown " << argv[i] << " option. Use -h option for help." << G4endl;
     exit(-1);
   }
  }
  else
  {
   G4cout << argv[0] << ": Bad arguments. Use -h option for help." << G4endl;
   exit(-1);
  }
  
  i++;
 }
 
 if (options->minEnergy >= options->maxEnergy)
 {
  G4cout << argv[0] << ": Mininum energy is higher than maximum energy" << G4endl;
  exit(-1);
 }

 G4cout << "Mean free path test:      ";
 if (options->meanFreePathTest)
  G4cout << "On";
 else
  G4cout << "Off";
 G4cout << G4endl << "Post step do it test:     ";
 if (options->postStepDoItTest)
  G4cout << "On";
 else
  G4cout << "Off";
 G4cout << G4endl << "Random energy generation: ";
 if (options->randomEnergy)
  G4cout << "On";
 else
  G4cout << "Off";
 G4cout << G4endl << "Output file:              " << options->outputFileName << G4endl;
 G4cout << "Material:                 " << options->material << G4endl;
 G4cout << "Process:                  " << options->process << G4endl;
 G4cout << "Min energy:               " << options->minEnergy/MeV << " MeV" << G4endl;
 G4cout << "Max energy:               " << options->maxEnergy/MeV << " MeV" << G4endl;
 G4cout << "N energy steps:           " << options->nEnergySteps << G4endl;
 G4cout << "N iterations:             " << options->nIterations << G4endl;
}

//! \brief Return the selected material
//! \param options Options for the material choice
//! \return The material
G4Material * GetSelectedMaterial(const struct Options & options)
{
 const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();

 G4int i(G4Material::GetNumberOfMaterials());
 
 while (i>0)
 {
  i--;
  
  if ((*theMaterialTable)[i]->GetName()==options.material)
   return (*theMaterialTable)[i];
 }
 
 i=G4Material::GetNumberOfMaterials();
 
 G4cout << "Available materials are: " << G4endl;
 while (i>0)
 {
  i--;
  G4cout << (*theMaterialTable)[i]->GetName();

  if (i>0)
   G4cout << ", ";
 }
 
 G4cout << G4endl;
 
 exit(-2);
 return 0;
}

//! \brief Creates the geometry
//! \param options Options for the material choice
//! \return The world volume
G4PVPlacement * CreateGeometry(const struct Options & options)
{
 G4Box* theFrame = new G4Box ("Frame", 1*mm, 1*mm, 1*mm);
  
 G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame, GetSelectedMaterial(options), "LFrame", 0, 0, 0);
  
 G4PVPlacement * placement = new G4PVPlacement(0, G4ThreeVector(), "PFrame", logicalFrame, 0, false, 0);
   
 G4cout << "[OK] Geometry built" << G4endl;
 return placement;
}

//! \brief Get process from options
//! \param options Options for the process choice
//! \return The choosen process
G4VLowEnergyTestableDiscreteProcess * GetSelectedProcess(const struct Options & options)
{
 static G4VLowEnergyTestableDiscreteProcess ** processes=0;
 if (!processes)
 {
  processes=new G4VLowEnergyTestableDiscreteProcess *[3];
  processes[0]=new G4LowEnergyPolarizedRayleigh;
  processes[1]=reinterpret_cast<G4VLowEnergyTestableDiscreteProcess *>(new G4LowEnergyRayleigh);
  processes[2]=0;
 }
 
 unsigned long i(0);
 while (processes[i])
 {
  if (processes[i]->GetProcessName()==options.process)
   return processes[i];
   
  i++;
 }
 
 G4cout << "Available processes are: " << G4endl;
 i=0;
 while (processes[i])
 {
  G4cout << processes[i]->GetProcessName();
  i++;
  
  if (processes[i])
   G4cout << ", ";
 }
 
 G4cout << G4endl;
 
 exit(-2);
 return 0;
}

//! \brief Setup processes
//! \param options Options for the process choice
void SetPhysics(const struct Options & options)
{
 G4ParticleDefinition * gamma(G4Gamma::GammaDefinition());
 G4ParticleDefinition * electron(G4Electron::ElectronDefinition());
 G4ParticleDefinition * positron(G4Positron::PositronDefinition());
  
 G4ProductionCutsTable * cutsTable(G4ProductionCutsTable::GetProductionCutsTable());
 G4ProductionCuts * cuts(cutsTable->GetDefaultProductionCuts());
 G4double cutG(1*micrometer);
 G4double cutE(1*micrometer);
 cuts->SetProductionCut(cutG, gamma);
 cuts->SetProductionCut(cutE, electron);
 cuts->SetProductionCut(cutE, positron);
 cutsTable->UpdateCoupleTable();
 G4cout << "[OK] Cuts are defined " << G4endl;

 G4VProcess * gammaProcess=GetSelectedProcess(options);
 if (! (gammaProcess->IsApplicable(*gamma)))
 {
  G4cout<< "Process " << gammaProcess->GetProcessName() << " is not applicable to photons" << G4endl;
  exit(0);
  return;
 }
 
 G4cout<< "[OK] Process " << gammaProcess->GetProcessName() << " is applicable to photons" << G4endl;
  
 G4ProcessManager * gProcessManager(new G4ProcessManager(gamma));
 gamma->SetProcessManager(gProcessManager);
 gProcessManager->AddDiscreteProcess(gammaProcess);

/* G4ProcessManager * eProcessManager(new G4ProcessManager(electron));
 G4VProcess * theEMinusMultipleScattering(new G4MultipleScattering());
 G4VProcess * theEMinusIonisation(new G4eIonisation());
 G4VProcess * theEMinusBremsstrahlung(new G4eBremsstrahlung());
 electron->SetProcessManager(eProcessManager);
 eProcessManager->AddProcess(theEMinusMultipleScattering);
 eProcessManager->AddProcess(theEMinusIonisation);
 eProcessManager->AddProcess(theEMinusBremsstrahlung);
 eProcessManager->SetProcessOrdering(theEMinusMultipleScattering, idxAlongStep, 1);
 eProcessManager->SetProcessOrdering(theEMinusIonisation,         idxAlongStep, 2);
 eProcessManager->SetProcessOrdering(theEMinusMultipleScattering, idxPostStep,  1);
 eProcessManager->SetProcessOrdering(theEMinusIonisation,         idxPostStep,  2);
 eProcessManager->SetProcessOrdering(theEMinusBremsstrahlung,     idxPostStep,  3);
  
 G4ProcessManager * pProcessManager(new G4ProcessManager(positron));
 G4VProcess * theEPlusMultipleScattering(new G4MultipleScattering());
 G4VProcess * theEPlusIonisation(new G4eIonisation());
 G4VProcess * theEPlusBremsstrahlung(new G4eBremsstrahlung());
 G4VProcess * theEPlusAnnihilation(new G4eplusAnnihilation());
 positron->SetProcessManager(pProcessManager);
 pProcessManager->AddProcess(theEPlusMultipleScattering);
 pProcessManager->AddProcess(theEPlusIonisation);
 pProcessManager->AddProcess(theEPlusBremsstrahlung);
 pProcessManager->AddProcess(theEPlusAnnihilation);
 pProcessManager->SetProcessOrderingToFirst(theEPlusAnnihilation, idxAtRest);
 pProcessManager->SetProcessOrdering(theEPlusMultipleScattering,  idxAlongStep, 1);
 pProcessManager->SetProcessOrdering(theEPlusIonisation,          idxAlongStep, 2);
 pProcessManager->SetProcessOrdering(theEPlusMultipleScattering,  idxPostStep,  1);
 pProcessManager->SetProcessOrdering(theEPlusIonisation,          idxPostStep,  2);
 pProcessManager->SetProcessOrdering(theEPlusBremsstrahlung,      idxPostStep,  3);
 pProcessManager->SetProcessOrdering(theEPlusAnnihilation,        idxPostStep,  4);*/
 G4cout << "[OK] Processes are defined " << G4endl;
 

 G4cout << "[OK] Building physics tables" << G4endl;
 gammaProcess->BuildPhysicsTable(* gamma);

/* theEMinusMultipleScattering->BuildPhysicsTable(* electron);
 theEMinusIonisation->BuildPhysicsTable(* electron);        
 theEMinusBremsstrahlung->BuildPhysicsTable(* electron);
 theEPlusMultipleScattering->BuildPhysicsTable(* positron);
 theEPlusIonisation->BuildPhysicsTable(* positron);
 theEPlusBremsstrahlung->BuildPhysicsTable(* positron);     
 theEPlusAnnihilation->BuildPhysicsTable(* positron);*/
 G4cout << "[OK] Physics tables built" << G4endl;
}

//! \brief Generates the step
//! \param options Options related to the track generation
//! \return The generated track
G4Step * GenerateStep(const struct Options & options)
{
 G4ThreeVector momentumDirection;
 
 momentumDirection.setRThetaPhi(1., std::acos(2.*G4UniformRand()-1.), twopi * G4UniformRand());

 G4ThreeVector vecA(momentumDirection.orthogonal());
 G4ThreeVector vecB(vecA.cross(momentumDirection));
 G4double beta(twopi * G4UniformRand());

 G4ThreeVector polarizationDirection(vecA * std::cos(beta)+ vecB * std::sin(beta));
 
 G4double lnEnergyMin=G4Log(options.minEnergy);
 G4double lnEnergyMax=G4Log(options.maxEnergy); 
 G4DynamicParticle * dynamicPhoton(new G4DynamicParticle(G4Gamma::Gamma(), momentumDirection, G4Exp(lnEnergyMin+(lnEnergyMax-lnEnergyMin)*G4UniformRand())));
 dynamicPhoton->SetPolarization(polarizationDirection.getX(), polarizationDirection.getY(), polarizationDirection.getZ());
  
 G4Track * aTrack(new G4Track(dynamicPhoton, 0., G4ThreeVector(0., 0., 0.)));

 G4Step* aStep(new G4Step());  
 aStep->SetTrack(aTrack);
 aTrack->SetStep(aStep);

 G4Material * material(GetSelectedMaterial(options));
 G4ProductionCutsTable * cutsTable(G4ProductionCutsTable::GetProductionCutsTable());
 const G4MaterialCutsCouple * theCouple(cutsTable->GetMaterialCutsCouple(material, cutsTable->GetDefaultProductionCuts()));

 G4StepPoint * aPoint(new G4StepPoint());
 aPoint->SetPosition(G4ThreeVector(0., 0., 0.));
 aPoint->SetMaterial(material);
 aPoint->SetMaterialCutsCouple(theCouple);
 aPoint->SetSafety(10000.*cm);

 aStep->SetPreStepPoint(aPoint);  
 
 return aStep;
}

void ProgressBar(G4int remainingIterations)
{
 static time_t startingTime;
 static time_t nextDumpTime;
 static G4int startingIteration(0);
 time_t now;
 
 if (remainingIterations==0)
 {
  startingIteration=0;
 }
 else if (startingIteration==0)
 {
  startingTime=time(0);
  nextDumpTime=startingTime+3;
  startingIteration=remainingIterations;
 }
 else
 {
  now=time(0);
  if (now>nextDumpTime)
  {
   nextDumpTime=now+10;
   G4double time;
   G4double perc;
   
   time=std::floor(static_cast<G4double>(now-startingTime)/static_cast<G4double>(startingIteration-remainingIterations)*static_cast<G4double>(remainingIterations)+0.5);
   perc=std::floor(static_cast<G4double>(remainingIterations)/static_cast<G4double>(startingIteration)*200.+.5)/2.;
   
   G4cout << "  " << perc << " % Remaining time: " << time << " s        \r";
   G4cout.flush();
  }
 }
}

//! \brief Test the mean free path table
//! \param tupleFactory The tuple factory
//! \param options Options related to the mean free path test
void MeanFreePathTest(AIDA::ITupleFactory * tupleFactory, const struct Options & options)
{
 AIDA::ITuple* iTuple = tupleFactory->create("1", "Mean Free Path Ntuple", "double k, log_k, mfp, log_mfp, cpu_time");
 
 G4double energy(options.minEnergy);
 G4double stpEnergy(std::pow(options.maxEnergy/energy, 1./static_cast<G4double>(options.nEnergySteps-1)));
 G4int step(options.nEnergySteps);
 
 G4ForceCondition condition;
 G4VLowEnergyTestableDiscreteProcess * process(GetSelectedProcess(options));
 
 G4double mfp;
 clock_t time;
 
 ProgressBar(0);
 while (step>0)
 {
  G4Step * aStep(GenerateStep(options));
  G4Track * aTrack(aStep->GetTrack());

  if (!options.randomEnergy)
  {
   aTrack->SetKineticEnergy(energy);
   energy*=stpEnergy;
  }
  ProgressBar(step);
  step--;

  time=clock();
  mfp=process->DumpMeanFreePath(*aTrack, 1.*mm, &condition)/cm;
  time=clock()-time;

  iTuple->fill(iTuple->findColumn("k"), aTrack->GetKineticEnergy()/eV);
  iTuple->fill(iTuple->findColumn("log_k"), G4Log(aTrack->GetKineticEnergy()/eV)/g4pow->logZ(10));
  iTuple->fill(iTuple->findColumn("mfp"), mfp);
  iTuple->fill(iTuple->findColumn("log_mfp"), G4Log(mfp)/g4pow->logZ(10));
  iTuple->fill(iTuple->findColumn("cpu_time"), static_cast<G4double>(time)/static_cast<G4double>(CLOCKS_PER_SEC));
  iTuple->addRow();
 
  delete aTrack;
  delete aStep;
 }
}

//! \brief Test the post step do it
//! \param tupleFactory The tuple factory
//! \param options Options related to the post step do it test
void PostStepDoItTest(AIDA::ITupleFactory * tupleFactory, const struct Options & options)
{
 AIDA::ITuple* iTuple = tupleFactory->create("2", "Post Step Do It Test", "double iteration, step, in_k, log_in_k, in_theta, in_phi, in_pol_theta, in_pol_phi, e_deposit, log_e_deposit, trk_status, out_k, log_out_k, out_theta, out_phi, out_pol_theta, out_pol_phi, cpu_time");
 
 G4double energy(options.minEnergy);
 G4double stpEnergy(std::pow(options.maxEnergy/energy, 1./static_cast<G4double>(options.nEnergySteps-1)));
 G4int step(options.nEnergySteps);
 
 G4VDiscreteProcess * process(GetSelectedProcess(options));
 clock_t time;
 
 ProgressBar(0);
 while (step>0)
 {
  G4Step * aStep(GenerateStep(options));
  G4Track * aTrack(aStep->GetTrack());
  const G4DynamicParticle * aParticle(aTrack->GetDynamicParticle());
  G4ThreeVector vector;
   
  if (!options.randomEnergy)
  {
   aTrack->SetKineticEnergy(energy);
   energy*=stpEnergy;
  }
  ProgressBar(step);
  step--;

  G4int iteration(options.nIterations);
  
  while (iteration>0)
  {
   iteration--;
   
   aStep->SetStepLength(1*micrometer);
      
   iTuple->fill(iTuple->findColumn("iteration"), iteration);
   iTuple->fill(iTuple->findColumn("step"), aStep->GetStepLength()/cm);

   iTuple->fill(iTuple->findColumn("in_k"), aParticle->GetKineticEnergy()/eV);
   iTuple->fill(iTuple->findColumn("log_in_k"), G4Log(aParticle->GetKineticEnergy()/eV)/g4pow->logZ(10));
   vector=aParticle->GetMomentumDirection();
   iTuple->fill(iTuple->findColumn("in_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("in_phi"), vector.phi());
   vector=aParticle->GetPolarization();
   iTuple->fill(iTuple->findColumn("in_pol_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("in_pol_phi"), vector.phi());
   
   time=clock();
   G4ParticleChange * particleChange(dynamic_cast<G4ParticleChange *>(process->PostStepDoIt(*aTrack, *aStep)));
   time=clock()-time;
   
   aTrack->SetKineticEnergy(particleChange->GetEnergy());
   aTrack->SetMomentumDirection(*particleChange->GetMomentumDirection());
   aTrack->SetPolarization(*particleChange->GetPolarization());
   
   iTuple->fill(iTuple->findColumn("e_deposit"), particleChange->GetLocalEnergyDeposit()/eV);
   iTuple->fill(iTuple->findColumn("log_e_deposit"), G4Log(particleChange->GetLocalEnergyDeposit()/eV)/g4pow->logZ(10));
   iTuple->fill(iTuple->findColumn("trk_status"), particleChange->GetTrackStatus());
   
   iTuple->fill(iTuple->findColumn("out_k"), aParticle->GetKineticEnergy()/eV);
   iTuple->fill(iTuple->findColumn("log_out_k"), G4Log(aParticle->GetKineticEnergy()/eV)/g4pow->logZ(10));
   vector=aParticle->GetMomentumDirection();
   iTuple->fill(iTuple->findColumn("out_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("out_phi"), vector.phi());
   vector=aParticle->GetPolarization();
   iTuple->fill(iTuple->findColumn("out_pol_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("out_pol_phi"), vector.phi());
   iTuple->fill(iTuple->findColumn("cpu_time"), static_cast<G4double>(time)/static_cast<G4double>(CLOCKS_PER_SEC));
   
   iTuple->addRow();
   
   particleChange->Clear();
  }

  delete aTrack;
  delete aStep;
 }
}

//! \brief Main function
//! \param argc Number of arguments
//! \param argv Pointer to the arguments
//! \return The exit value 
int main(int argc, char ** argv)
{
 G4Pow* g4pow;
 g4pow = G4Pow::GetInstance();

 struct Options options;
 processOptions(argc, argv, &options);
 
 CreateMaterials();
 
 GetSelectedProcess(options);
 GetSelectedMaterial(options);
 
 G4RunManager* rm = new G4RunManager();
 rm->GeometryHasBeenModified();
 rm->DefineWorldVolume(CreateGeometry(options));
 G4cout << "[OK] World is defined " << G4endl;
 
 SetPhysics(options);
 
 if (!(options.meanFreePathTest || options.postStepDoItTest))
 {
  G4cout << "[OK] Program completed" << G4endl;
  return 0;
 }
 
 // HBOOK initialization 
 AIDA::IAnalysisFactory * analysisFactory(AIDA_createAnalysisFactory());
 AIDA::ITreeFactory * treeFactory(analysisFactory->createTreeFactory());
 AIDA::ITree * tree(treeFactory->create(options.outputFileName, "hbook", false, true));
 G4cout << "[OK] Tree store: " << tree->storeName() << G4endl;
 
 AIDA::ITupleFactory * tupleFactory(analysisFactory->createTupleFactory(*tree));
 
 // Mean free path test
 if (options.meanFreePathTest)
 {
  G4cout << "[OK] Mean free path test started" << G4endl;
  MeanFreePathTest(tupleFactory, options);
  G4cout << "[OK] Mean free path test completed" << G4endl;
 }
 
 // Post step do it test
 if (options.postStepDoItTest)
 {
  G4cout << "[OK] Post step do it test started" << G4endl;
  PostStepDoItTest(tupleFactory, options);
  G4cout << "[OK] Post step do it test completed" << G4endl;
 }
 
 G4cout << "[OK] Storing analysis data" << G4endl;
 tree->commit();
 tree->close();
 
 G4cout << "[OK] Deleting analysis data" << G4endl;
 delete tupleFactory;
 delete tree;
 delete treeFactory;
 delete analysisFactory;

 G4cout << "[OK] Program completed" << G4endl;
 return 0;
}
