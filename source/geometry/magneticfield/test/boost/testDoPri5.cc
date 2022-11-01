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
#include "G4DormandPrince745.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldTrack.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4FieldUtils.hh"
#include "G4LineSection.hh"

#include "boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp"

#include <memory>
#include <cstdio>

const size_t INTEGRATED_COMPONENTS = 6;
const size_t TOTAL_COMPONENTS = G4FieldTrack::ncompSVEC;
const size_t NUMBER_OF_INTEGRATION_STEPS = 1000;
const size_t NUMBER_OF_INTERPOLATION_STEPS = 10;
const double DIFF_THRESHOLD = 1e-7;
const std::string VERBOSE_MODE = "verbose";

using State = double[TOTAL_COMPONENTS];

using BoostDoPri5 = typename boost::numeric::odeint::runge_kutta_dopri5<State>;
using BoostEquation = std::function<void(const State&, State&, double)>;

using GeantDoPri5 = G4DormandPrince745;
using GeantEquation = std::shared_ptr<G4EquationOfMotion>;

using namespace field_utils;


enum class Verbosity {
    Silent,
    Verbose
} verbosity;

double diff(const State& state1, const State& state2)
{
    double result = 0;
    for (size_t i = 0; i < INTEGRATED_COMPONENTS; ++i) {
        result = std::max(result, std::abs(state1[i] - state2[i]));
    }
    return result;
}

void report(size_t iter, size_t subiter,
            const State& geantY, const State& geantError,
            const State& boostY, const State& boostError)
{
    auto geantPosition = makeVector(geantY, Value3D::Position);
    auto boostPosition = makeVector(boostY, Value3D::Position);
    auto diffPosition =
        (geantPosition - boostPosition).mag() / boostPosition.mag();
    
    auto geantMomentum = makeVector(geantY, Value3D::Momentum);
    auto boostMomentum = makeVector(boostY, Value3D::Momentum);
    auto diffMomentum =
        (geantMomentum - boostMomentum).mag() / boostMomentum.mag();
    
    auto geantPositionError = makeVector(geantError, Value3D::Position);
    auto boostPositionError = makeVector(boostError, Value3D::Position);
    auto diffPositionError =
        (geantPositionError - boostPositionError).mag() / boostPositionError.mag();
    
    auto geantMomentumError = makeVector(geantError, Value3D::Momentum);
    auto boostMomentumError = makeVector(boostError, Value3D::Momentum);
    auto diffMomentumError =
        (geantMomentumError - boostMomentumError).mag() / boostMomentumError.mag();

    fprintf(stderr, "STEP: %zd    ", iter);
    fprintf(stderr, "SUBSTEP: %zd    ", subiter);
    fprintf(stderr, "pos diff: %.3e    ", diffPosition);
    fprintf(stderr, "mom diff: %.3e    ", diffMomentum);
    fprintf(stderr, "pos error diff: %.3e    ", diffPositionError);
    fprintf(stderr, "mom error diff: %.3e\n", diffMomentumError);
}

void test(GeantDoPri5& geantDoPri5, const GeantEquation& geantEquation, const State& geantState,
          BoostDoPri5& boostDoPri5, const BoostEquation& boostEquation, const State& boostState)
{
    State geantY0, geantY1, geantDydx, geantError;
    copy(geantY0, geantState);

    State boostY0, boostDydx0, boostError, boostY1, boostDydx1;
    copy(boostY0, boostState);

    double curveLength = 0;
    double stepLength = 1. * CLHEP::mm;

    // boost version is FSAL!
    boostEquation(boostY0, boostDydx0, curveLength);

    for (size_t i = 0; i < NUMBER_OF_INTEGRATION_STEPS; ++i) {
        geantEquation->RightHandSide(geantY0, geantDydx);

        geantDoPri5.Stepper(geantY0, geantDydx, stepLength, geantY1, geantError);
        boostDoPri5.do_step_impl(boostEquation, boostY0, boostDydx0, curveLength,
                                 boostY1, boostDydx1, stepLength, boostError);
        
        G4cout << "DistChord: " <<geantDoPri5.DistChord() << G4endl;

        State geantY2;
        geantDoPri5.SetupInterpolation_low();
        geantDoPri5.Interpolate_low(geantY2, 0.5);
        G4cout << "InterpLow: " <<G4LineSection::Distline(
            makeVector(geantY2, Value3D::Position), 
            makeVector(geantY0, Value3D::Position), 
            makeVector(geantY1, Value3D::Position)) << G4endl;

        geantDoPri5.SetupInterpolation_high();
        geantDoPri5.Interpolate_high(geantY2, 0.5);
        G4cout << "InterpHigh: " <<G4LineSection::Distline(
            makeVector(geantY2, Value3D::Position),
            makeVector(geantY0, Value3D::Position),  
            makeVector(geantY1, Value3D::Position)) << G4endl;

        G4cout << G4endl << G4endl;


        for (size_t j = 1; j <= NUMBER_OF_INTERPOLATION_STEPS; ++j) {
            const double tau = j * 1. / NUMBER_OF_INTERPOLATION_STEPS;
            State geantY;
            geantDoPri5.Interpolate_low(geantY, tau);
            
            State boostY;
            boostDoPri5.calc_state(curveLength + tau * stepLength, boostY,
                                   boostY0, boostDydx0, curveLength, boostY1,
                                   boostDydx1, curveLength + stepLength);
            
            if (diff(geantY, boostY) > DIFF_THRESHOLD || verbosity == Verbosity::Verbose) {
                report(i, j, geantY, geantError, boostY, boostError);
            }
        }
        
        copy(geantY0, geantY1);
        copy(boostY0, boostY1);
        copy(boostDydx0, boostDydx1);
        
        curveLength += stepLength;
    }
}

int main(int argc, const char* argv[])
{
    verbosity = Verbosity::Silent;
    if (argc > 1 && argv[1] == VERBOSE_MODE) {
        verbosity = Verbosity::Verbose;
    }
    
    G4DynamicParticle dynParticle(
            G4Proton::Definition(),
            G4ThreeVector(1, 0, 2).unit(),
            0.01 * CLHEP::GeV);

    auto track =
        std::make_shared<G4FieldTrack>(
            G4ThreeVector(0., 1., 0.), // start position
            0,                         // LaboratoryTimeOfFlight
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

    GeantDoPri5 geantDoPri5(equation.get(), INTEGRATED_COMPONENTS);
    State geantState;
    track->DumpToArray(geantState);

    BoostDoPri5 boostDoPri5;
    State boostState;
    track->DumpToArray(boostState);

    test(geantDoPri5, equation, geantState,
         boostDoPri5, system, boostState);

    return 0;
}
