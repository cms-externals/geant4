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
#include "G4ModifiedMidpoint.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldTrack.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4FieldUtils.hh"

#include "boost/numeric/odeint/stepper/modified_midpoint.hpp"

#include <memory>
#include <cstdio>

const size_t INTEGRATED_COMPONENTS = 6;
const size_t TOTAL_COMPONENTS = G4FieldTrack::ncompSVEC;
const size_t NUMBER_OF_INTEGRATION_STEPS = 1000;
const double DIFF_THRESHOLD = 1e-7;
const std::string VERBOSE_MODE = "verbose";

using State = double[TOTAL_COMPONENTS];

using BoostModifiedMidpoint = typename boost::numeric::odeint::modified_midpoint<State>;
using BoostEquation = std::function<void(const State&, State&, double)>;

using GeantEquation = std::shared_ptr<G4EquationOfMotion>;

using namespace field_utils;

double diff(const State& state1, const State& state2)
{
    double result = 0;
    for (size_t i = 0; i < INTEGRATED_COMPONENTS; ++i) {
        result = std::max(result, std::abs(state1[i] - state2[i]));
    }
    return result;
}

void report(size_t iter, const State& geantY, const State& boostY)
{
    auto geantPosition = makeVector(geantY, Value3D::Position);
    auto boostPosition = makeVector(boostY, Value3D::Position);
    auto diffPosition =
        (geantPosition - boostPosition).mag() / boostPosition.mag();
    
    auto geantMomentum = makeVector(geantY, Value3D::Momentum);
    auto boostMomentum = makeVector(boostY, Value3D::Momentum);
    auto diffMomentum =
        (geantMomentum - boostMomentum).mag() / boostMomentum.mag();

    if (diffPosition > DIFF_THRESHOLD || diffMomentum > DIFF_THRESHOLD) {
        fprintf(stderr, "STEP: %zd    ", iter);
        fprintf(stderr, "pos diff: %.3e    ", diffPosition);
        fprintf(stderr, "mom diff: %.3e    \n", diffMomentum);
    }
}

void test(G4ModifiedMidpoint& geantMethod, const GeantEquation& geantEquation, const State& geantState,
          BoostModifiedMidpoint& boostMethod, const BoostEquation& boostEquation, const State& boostState)
{
    State geantY, geantDydx;
    copy(geantY, geantState);

    State boostY, boostDydx;
    copy(boostY, boostState);

    double curveLength = 0;
    double stepLength = 1. * CLHEP::mm;

    for (size_t i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {

        geantEquation->RightHandSide(geantY, geantDydx);
        boostEquation(boostY, boostDydx, curveLength);

        geantMethod.DoStep(geantY, geantDydx, geantY, stepLength);
        boostMethod.do_step_impl(boostEquation, boostY, boostDydx, curveLength,
                                 boostY, stepLength);

        report(i, geantY, boostY);
        
        curveLength += stepLength;
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
            45,                         // LaboratoryTimeOfFlight
            dynParticle.GetMomentumDirection(),
            dynParticle.GetKineticEnergy(),
            dynParticle.GetMass(),
            dynParticle.GetCharge(),
            dynParticle.GetPolarization());

    auto field = std::make_shared<G4QuadrupoleMagField>(1 * CLHEP::tesla);
    auto equation = std::make_shared<G4Mag_UsualEqRhs>(field.get());

    equation->SetChargeMomentumMass(
            {
                dynParticle.GetCharge(),
                dynParticle.GetSpin(),
                dynParticle.GetMagneticMoment()
            },
            dynParticle.GetMomentum().mag(),
            dynParticle.GetMass());

    auto system = [equation](const State& yIn, State& dydxOut, double /*t*/)
    {
        equation->RightHandSide(yIn, dydxOut);
    };

    G4ModifiedMidpoint g4modifiedMidpoint(equation.get(), INTEGRATED_COMPONENTS);
    State geantState;
    track->DumpToArray(geantState);

    BoostModifiedMidpoint boostModifiedMidpoint;
    State boostState;
    track->DumpToArray(boostState);

    test(g4modifiedMidpoint, equation, geantState,
         boostModifiedMidpoint, system, boostState);

    return 0;
}
