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
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldTrack.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4FieldUtils.hh"

#include "G4ClassicalRK4.hh"
#include "G4BogackiShampine23.hh"
#include "G4BogackiShampine45.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ConstRK4.hh"
#include "G4DoLoMcPriRK34.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4ExactHelixStepper.hh"
#include "G4FSALBogackiShampine45.hh"
#include "G4FSALDormandPrince745.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixMixedStepper.hh"
#include "G4MagHelicalStepper.hh"
#include "G4NystromRK4.hh"
#include "G4RK547FEq1.hh"
#include "G4RK547FEq2.hh"
#include "G4RK547FEq3.hh"
#include "G4RKG3_Stepper.hh"
#include "G4TsitourasRK45.hh"


#include <memory>
#include <cstdio>
#include <assert.h>
#include <stdexcept>

const size_t TOTAL_COMPONENTS = G4FieldTrack::ncompSVEC;
const size_t NUMBER_OF_INTEGRATION_STEPS = 10;
const G4double STEP_LENGTH = 1. * CLHEP::mm;
const G4double LAB_TIME = 42;

using State = double[TOTAL_COMPONENTS];

std::shared_ptr<G4Mag_EqRhs> EQUATION = nullptr;


class TestField : public G4QuadrupoleMagField {
public:

    TestField(G4double pGradient):
        G4QuadrupoleMagField(pGradient)
    {}

    virtual void  GetFieldValue(const  double Point[4],
                                double *fieldArr) const
    {
        G4double time = Point[3];
        if (time != LAB_TIME) {
            std::ostringstream what;
            what << "time = " << time;
            throw(std::invalid_argument(what.str()));
        }
        G4QuadrupoleMagField::GetFieldValue(Point, fieldArr);
    }
};



template <typename Stepper>
std::shared_ptr<G4MagIntegratorStepper> make()
{
    assert(EQUATION);
    return std::make_shared<Stepper>(EQUATION.get());
}

void test(G4MagIntegratorStepper* stepper, const State& state)
{
    State yIn, dydxIn, yOut, yError;
    field_utils::copy(yIn, state);

    for (size_t i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        EQUATION->RightHandSide(yIn, dydxIn);
        stepper->Stepper(yIn, dydxIn, STEP_LENGTH, yOut, yError);
        field_utils::copy(yOut, yIn);
    }
}

int main()
{
    G4DynamicParticle dynParticle(
            G4Proton::Definition(),
            G4ThreeVector(1, 0, 2).unit(),
            0.01 * CLHEP::GeV);

    auto track =
        std::make_shared<G4FieldTrack>(
            G4ThreeVector(0., 1., 0.), // start position
            LAB_TIME,                  // LaboratoryTimeOfFlight
            dynParticle.GetMomentumDirection(),
            dynParticle.GetKineticEnergy(),
            dynParticle.GetMass(),
            dynParticle.GetCharge(),
            dynParticle.GetPolarization());

    auto field = std::make_shared<TestField>(1 * CLHEP::tesla);
    
    EQUATION = std::make_shared<G4Mag_UsualEqRhs>(field.get());
    EQUATION->SetChargeMomentumMass(
            {
                dynParticle.GetCharge(),
                dynParticle.GetSpin(),
                dynParticle.GetMagneticMoment()
            },
            dynParticle.GetMomentum().mag(),
            dynParticle.GetMass());

    State state;
    track->DumpToArray(state);
    std::vector<std::pair<std::shared_ptr<G4MagIntegratorStepper>, std::string>> steppers = {
        {make<G4ClassicalRK4>(), "G4ClassicalRK4"},
        {make<G4BogackiShampine23>(), "G4BogackiShampine23"},
        {make<G4BogackiShampine45>(), "G4BogackiShampine45"},
        {make<G4CashKarpRKF45>(), "G4CashKarpRKF45"},
        //{make<G4ConstRK4>(), "G4ConstRK4"},
        {make<G4DoLoMcPriRK34>(), "G4DoLoMcPriRK34"},
        {make<G4DormandPrince745>(), "G4DormandPrince745"},
        {make<G4DormandPrinceRK56>(), "G4DormandPrinceRK56"},
        {make<G4DormandPrinceRK78>(), "G4DormandPrinceRK78"},
        {make<G4ExactHelixStepper>(), "G4ExactHelixStepper"},
        //make<G4FSALBogackiShampine45>(),
        //make<G4FSALDormandPrince745>(),
        {make<G4HelixExplicitEuler>(), "G4HelixExplicitEuler"},
        {make<G4HelixMixedStepper>(), "G4HelixMixedStepper"},
        //{make1<G4MagHelicalStepper>(), "G4MagHelicalStepper"},
        {make<G4NystromRK4>(), "G4NystromRK4"},
        {make<G4RK547FEq1>(), "G4RK547FEq1"},
        {make<G4RK547FEq2>(), "G4RK547FEq2"},
        {make<G4RK547FEq3>(), "G4RK547FEq3"},
        //make<G4RKG3_Stepper>(),
        {make<G4TsitourasRK45>(), "G4TsitourasRK45"}
    };
    for (auto& stepper : steppers) {
        G4cout << "testing " << stepper.second << " ... ";
        try {
            test(stepper.first.get(), state);
            G4cout << " SUCCESS" << G4endl;
        } catch (std::invalid_argument e) {
            G4cout << " FAIL: " << e.what() << G4endl;
        }
    }

    return 0;
}
