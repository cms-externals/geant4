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
// Implementation of the G4FieldSetup class
//
// Author: Ivana Hrivnacova (IJClab, Orsay), 2024.
// --------------------------------------------------------------------

#include "G4FieldSetup.hh"

#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4CachedMagneticField.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4ConstRK4.hh"
#include "G4DoLoMcPriRK34.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4ElectroMagneticField.hh"
#include "G4EqEMFieldWithEDM.hh"
#include "G4EqEMFieldWithSpin.hh"
#include "G4EqMagElectricField.hh"
#include "G4ExactHelixStepper.hh"
#include "G4Exception.hh"
#include "G4ExplicitEuler.hh"
#include "G4FSALIntegrationDriver.hh"
#include "G4FieldManager.hh"
#include "G4FieldSetupMessenger.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixMixedStepper.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4ImplicitEuler.hh"
#include "G4IntegrationDriver.hh"
#include "G4LogicalVolume.hh"
#include "G4MagErrorStepper.hh"
#include "G4MagHelicalStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4Mag_EqRhs.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagneticField.hh"
#include "G4NystromRK4.hh"
#include "G4RK547FEq1.hh"
#include "G4RK547FEq2.hh"
#include "G4RK547FEq3.hh"
#include "G4RKG3_Stepper.hh"
#include "G4SimpleHeum.hh"
#include "G4SimpleRunge.hh"
#include "G4SystemOfUnits.hh"
#include "G4TsitourasRK45.hh"
#include "G4VIntegrationDriver.hh"

namespace
{
template<class Stepper, class Equation>
G4VIntegrationDriver* CreateIntegrationDriver(Equation* equation, G4double minStep)
{
  auto* stepper = new Stepper(equation);
  return new G4IntegrationDriver<Stepper>(minStep, stepper, stepper->GetNumberOfVariables());
}

template<class Stepper>
G4VIntegrationDriver* CreateFSALIntegrationDriver(G4EquationOfMotion* equation, G4double minStep)
{
  auto* stepper = new Stepper(equation);
  return new G4FSALIntegrationDriver<Stepper>(minStep, stepper, stepper->GetNumberOfVariables());
}
}  // namespace

//_____________________________________________________________________________
G4FieldSetup::G4FieldSetup(const G4FieldParameters& parameters, G4Field* field, G4LogicalVolume* lv)
  : fParameters(parameters), fG4Field(field), fLogicalVolume(lv)
{
  // Standard constructor

  fMessenger = new G4FieldSetupMessenger(this);

  // Get or create field manager
  if (fLogicalVolume == nullptr)
  {
    // global field
    fFieldManager = G4FieldManager::GetGlobalFieldManager();
  }
  else
  {
    // local field
    fFieldManager = new G4FieldManager();
    G4bool overwriteDaughtersField = true;
    // TO DO: this parameter should be made optional for users
    fLogicalVolume->SetFieldManager(fFieldManager, overwriteDaughtersField);
  }
}

//_____________________________________________________________________________
G4FieldSetup::~G4FieldSetup()
{
  // Destructor

  delete fG4Field;
  
  // G4ChordFinder owns and deletes its integration driver.
  // Delete fDriver directly only if no chord finder took ownership.
  if (fChordFinder != nullptr)
  {
    delete fChordFinder;
  }
  else
  {
    delete fDriver;
  }

  // Integration drivers do not own their steppers.
  delete fStepper;
}

//
// private methods
//

//_____________________________________________________________________________
G4Field* G4FieldSetup::CreateCachedField(const G4FieldParameters& parameters, G4Field* field)
{
  // Create cached magnetic field if const distance is set > 0.
  // and field is of G4MagneticField.
  // Return the input field otherwise.

  if (parameters.GetConstDistance() > 0.)
  {
    auto magField = dynamic_cast<G4MagneticField*>(field);
    if (magField == nullptr)
    {
      G4Exception("G4FieldSetup::CreateCachedField:", "GeomFieldParameters0001", JustWarning,
                  "Incompatible field type.");
      return field;
    }
    return new G4CachedMagneticField(magField, parameters.GetConstDistance());
  }

  return field;
}

//_____________________________________________________________________________
G4EquationOfMotion* G4FieldSetup::CreateEquation(G4EquationType equation)
{
  // Set the equation of motion of a particle in a field

  // magnetic fields
  G4MagneticField* magField = nullptr;
  if (equation == kEqMagnetic || equation == kEqMagneticWithSpin)
  {
    magField = dynamic_cast<G4MagneticField*>(fG4Field);
    if (magField == nullptr)
    {
      G4Exception("G4FieldSetup::CreateEquation:", "GeomFieldParameters0001", FatalErrorInArgument,
                  "Incompatible field and equation.\n"
                  "The field type must be set explicitly for other than magnetic field.");
      return nullptr;
    }
  }

  // electromagnetic fields
  G4ElectroMagneticField* elMagField = nullptr;
  if (equation >= kEqElectroMagnetic && equation <= kEqEMfieldWithEDM)
  {
    elMagField = dynamic_cast<G4ElectroMagneticField*>(fG4Field);
    if (elMagField == nullptr)
    {
      G4Exception("G4FieldSetup::CreateEquation:", "GeomFieldParameters0001", FatalErrorInArgument,
                  "Incompatible field and equation.\n"
                  "The field type must be set explicitly for other than magnetic field.");
      return nullptr;
    }
  }

  // electromagnetic fields
  switch (equation)
  {
    case kEqMagnetic:
      return new G4Mag_UsualEqRhs(magField);
      break;

    case kEqMagneticWithSpin:
      return new G4Mag_SpinEqRhs(magField);
      break;

    case kEqElectroMagnetic:
      return new G4EqMagElectricField(elMagField);
      break;

    case kEqEMfieldWithSpin:
      return new G4EqEMFieldWithSpin(elMagField);
      break;

    case kEqEMfieldWithEDM:
      return new G4EqEMFieldWithEDM(elMagField);
      break;

    // other fields
    case kEqGravity:
    case kEqMonopole:
    case kEqReplate:
      G4Exception(
        "G4FieldSetup::CreateEquation:", "GeomFieldParameters0001", FatalErrorInArgument,
        "Limitation: Equation not supported in G4FieldBuilder.\n"
        "Only magnetic and electromagnetic field can be constructed with G4FieldBuilder.");
      return nullptr;
      break;
    case kUserEquation:
      // nothing to be done
      return nullptr;
      break;
  }

  G4Exception("G4FieldSetup::CreateEquation:", "GeomFieldParameters0001", FatalErrorInArgument,
              "Unknown equation type.");
  return nullptr;
}

//_____________________________________________________________________________
G4VIntegrationDriver* G4FieldSetup::CreateStepperAndDriver(G4EquationOfMotion* equation,
                                                           G4StepperType stepper, G4double minStep)
{
  // Create the integrator and preserve its concrete type in the driver.

  // Check steppers which require equation of motion of G4Mag_EqRhs type
  auto eqRhs = dynamic_cast<G4Mag_EqRhs*>(equation);
  if ((eqRhs == nullptr) && (stepper > kTsitourasRK45) && (stepper < kUserStepper))
  {
    G4Exception("G4FieldSetup::CreateStepperAndDriver:", "GeomFieldParameters0001",
                FatalErrorInArgument,
                "The stepper type requires equation of motion of G4Mag_EqRhs type.");
    return nullptr;
  }

  switch (stepper)
  {
    case kBogackiShampine23:
      return CreateIntegrationDriver<G4BogackiShampine23>(equation, minStep);

    case kBogackiShampine45:
      return CreateIntegrationDriver<G4BogackiShampine45>(equation, minStep);

    case kCashKarpRKF45:
      return CreateIntegrationDriver<G4CashKarpRKF45>(equation, minStep);

    case kClassicalRK4:
      return CreateIntegrationDriver<G4ClassicalRK4>(equation, minStep);

    case kDoLoMcPriRK34:
      return CreateIntegrationDriver<G4DoLoMcPriRK34>(equation, minStep);

    case kDormandPrince745:
      return CreateIntegrationDriver<G4DormandPrince745>(equation, minStep);

    case kDormandPrinceRK56:
      return CreateIntegrationDriver<G4DormandPrinceRK56>(equation, minStep);

    case kDormandPrinceRK78:
      return CreateIntegrationDriver<G4DormandPrinceRK78>(equation, minStep);

    case kExplicitEuler:
      return CreateIntegrationDriver<G4ExplicitEuler>(equation, minStep);

    case kImplicitEuler:
      return CreateIntegrationDriver<G4ImplicitEuler>(equation, minStep);

    case kSimpleHeum:
      return CreateIntegrationDriver<G4SimpleHeum>(equation, minStep);

    case kSimpleRunge:
      return CreateIntegrationDriver<G4SimpleRunge>(equation, minStep);

    case kTsitourasRK45:
      return CreateIntegrationDriver<G4TsitourasRK45>(equation, minStep);

    case kConstRK4:
      return CreateIntegrationDriver<G4ConstRK4>(eqRhs, minStep);

    case kExactHelixStepper:
      return CreateIntegrationDriver<G4ExactHelixStepper>(eqRhs, minStep);

    case kHelixExplicitEuler:
      return CreateIntegrationDriver<G4HelixExplicitEuler>(eqRhs, minStep);

    case kHelixHeum:
      return CreateIntegrationDriver<G4HelixHeum>(eqRhs, minStep);

    case kHelixImplicitEuler:
      return CreateIntegrationDriver<G4HelixImplicitEuler>(eqRhs, minStep);

    case kHelixMixedStepper:
      return CreateIntegrationDriver<G4HelixMixedStepper>(eqRhs, minStep);

    case kHelixSimpleRunge:
      return CreateIntegrationDriver<G4HelixSimpleRunge>(eqRhs, minStep);

    case kNystromRK4:
      return CreateIntegrationDriver<G4NystromRK4>(eqRhs, minStep);

    case kRKG3Stepper:
      return CreateIntegrationDriver<G4RKG3_Stepper>(eqRhs, minStep);

    case kRK547FEq1:
      return CreateFSALIntegrationDriver<G4RK547FEq1>(equation, minStep);

    case kRK547FEq2:
      return CreateFSALIntegrationDriver<G4RK547FEq2>(equation, minStep);

    case kRK547FEq3:
      return CreateFSALIntegrationDriver<G4RK547FEq3>(equation, minStep);

    case kUserStepper:
      // nothing to be done
      return nullptr;

    // templated steppers
    case kTCashKarpRKF45:
    case kTDormandPrince45:
    case kTMagErrorStepper:
    case kQSStepper:
      G4Exception("G4FieldSetup::CreateStepperAndDriver:", "GeomFieldParameters0001",
                  FatalErrorInArgument,
                  "Limitation: Templated steppers not supported in G4FieldBuilder");
      return nullptr;

    default:
      G4Exception("G4FieldSetup::CreateStepperAndDriver:", "GeomFieldParameters0001",
                  FatalErrorInArgument, "Incorrect stepper type.");
      return nullptr;
  }
}

//_____________________________________________________________________________
void G4FieldSetup::CreateCachedField()
{
  // Create cached field (if ConstDistance is set)
  fG4Field = CreateCachedField(fParameters, fG4Field);
}

//_____________________________________________________________________________
void G4FieldSetup::CreateStepper()
{
  // Create equation of motion (or get the user one if defined)
  if (fParameters.GetEquationType() == kUserEquation)
  {
    fEquation = fParameters.GetUserEquationOfMotion();
  }
  else
  {
    delete fEquation;
    fEquation = nullptr;
    fEquation = CreateEquation(fParameters.GetEquationType());
  }
  fEquation->SetFieldObj(fG4Field);

  // Create stepper  (or get the user one if defined)
  if (fParameters.GetStepperType() == kUserStepper)
  {
    // User stepper
    fDriver = nullptr;
    fStepper = fParameters.GetUserStepper();
  }
  else
  {
    // Built-in stepper and type-specific driver
    delete fStepper;
    fStepper = nullptr;
    fDriver =
      CreateStepperAndDriver(fEquation, fParameters.GetStepperType(), fParameters.GetMinimumStep());
    if (fDriver != nullptr)
    {
      fStepper = fDriver->GetStepper();
    }
  }
}

//_____________________________________________________________________________
void G4FieldSetup::CreateChordFinder()
{
  // Chord finder
  if (fDriver != nullptr)
  {
    fChordFinder = new G4ChordFinder(fDriver);
  }
  else if (fParameters.GetFieldType() == kMagnetic)
  {
    // Polymorphic fallback for a user-provided stepper
    fChordFinder = new G4ChordFinder(static_cast<G4MagneticField*>(fG4Field),
                                     fParameters.GetMinimumStep(), fStepper);
  }
  else if (fParameters.GetFieldType() == kElectroMagnetic)
  {
    // Legacy polymorphic fallback for a user-provided stepper
    auto intDriver =
      new G4MagInt_Driver(fParameters.GetMinimumStep(), fStepper, fStepper->GetNumberOfVariables());
    if (intDriver != nullptr)
    {
      // Chord finder
      fChordFinder = new G4ChordFinder(intDriver);
    }
  }
  if (fChordFinder != nullptr)
  {
    fChordFinder->SetDeltaChord(fParameters.GetDeltaChord());
  }
}

//_____________________________________________________________________________
void G4FieldSetup::UpdateFieldManager()
{
  fFieldManager->SetChordFinder(fChordFinder);
  fFieldManager->SetDetectorField(fG4Field);

  fFieldManager->SetMinimumEpsilonStep(fParameters.GetMinimumEpsilonStep());
  fFieldManager->SetMaximumEpsilonStep(fParameters.GetMaximumEpsilonStep());
  fFieldManager->SetDeltaOneStep(fParameters.GetDeltaOneStep());
  fFieldManager->SetDeltaIntersection(fParameters.GetDeltaIntersection());
}

//
// public methods
//

//_____________________________________________________________________________
void G4FieldSetup::Clear()
{
  // First clean up previous state.
  if (fChordFinder != nullptr)
  {
    delete fChordFinder;
  }
  else
  {
    delete fDriver;
  }
  fChordFinder = nullptr;
  fDriver = nullptr;

  if (fG4Field == nullptr)
  {
    delete fEquation;
    delete fStepper;
    fEquation = nullptr;
    fStepper = nullptr;
    fFieldManager->SetChordFinder(fChordFinder);
    fFieldManager->SetDetectorField(fG4Field);
  }
}

//_____________________________________________________________________________
void G4FieldSetup::Update()
{
  // Update field with new field parameters

  Clear();
  if (fG4Field == nullptr)
  {
    // No further update needed
    return;
  }

  CreateCachedField();
  CreateStepper();
  CreateChordFinder();
  UpdateFieldManager();
}

//_____________________________________________________________________________
void G4FieldSetup::PrintInfo(G4int verboseLevel, const G4String& about)
{
  if (verboseLevel == 0)
  {
    return;
  }

  auto fieldType = G4FieldParameters::FieldTypeName(fParameters.GetFieldType());
  auto isCachedMagneticField = (fParameters.GetConstDistance() > 0.);
  if (fLogicalVolume == nullptr)
  {
    fieldType = "Global";
  }
  else
  {
    fieldType = "Local (in ";
    fieldType.append(fLogicalVolume->GetName());
    fieldType.append(")");
  }
  if (isCachedMagneticField)
  {
    fieldType.append(" cached");
  }

  G4cout << fieldType << " field " << about << " with stepper ";
  G4cout << G4FieldParameters::StepperTypeName(fParameters.GetStepperType()) << G4endl;

  if (verboseLevel > 1)
  {
    fParameters.PrintParameters();
  }
}
