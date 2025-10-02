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
// created by jrmadsen on Sun Oct 21 21:31:53 2018
//
//
//

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include <limits>
#include <thread>

#include "G4ConvergenceTester.hh"
#include "G4StatAnalysis.hh"
#include "Randomize.hh"

#include "testing.hh"

//============================================================================//

static intmax_t             nrand = 1000000;
static intmax_t             nzero = 1000;
static G4ConvergenceTester  conv_obj("test76");
static G4StatAnalysis       stat_obj;
static const G4double       epsilon = std::numeric_limits<G4float>::epsilon();

//============================================================================//

void check_mean();
void check_relative_error();
void check_standard_deviation();
void check_variance();
void check_efficiency();
void check_r2int();
void check_r2eff();

//============================================================================//

int main(int argc, char** argv)
{
    G4StatAnalysis::ResetCpuClock();

    if(argc > 1)
        nrand = atol(argv[1]);
    if(argc > 2)
        nzero = atol(argv[2]);

    std::cout << "\nExecuting " << argv[0] << ". # random values = "
              << nrand << ", # zero values = " << nzero << "..." << std::endl;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values
    // from the mean
    std::normal_distribution<> d{5,2};

    for(intmax_t i = 0; i < nrand - nzero; ++i)
    {
        G4double rand = d(gen);
        // G4ConvergenceTester does not accept negative values
        if(rand < 0.0)
            rand -= 2.0*rand;
        conv_obj += rand;
        stat_obj += rand;
    }

    // so that we have efficiency less than 1.0
    for(intmax_t i = 0; i < nzero; ++i)
    {
        conv_obj += 0.0;
        stat_obj += 0.0;
    }

    // ensure statistics are computed
    conv_obj.ComputeStatistics();

    // for reference
    std::cout << "\nG4StatAnalysis: " << stat_obj << std::endl;
    conv_obj.ShowResult(std::cout);

    int num_test = 0;
    int num_fail = 0;

    // run the tests
    RUN_TEST(check_mean, num_test, num_fail);
    RUN_TEST(check_relative_error, num_test, num_fail);
    RUN_TEST(check_standard_deviation, num_test, num_fail);
    RUN_TEST(check_variance, num_test, num_fail);
    RUN_TEST(check_efficiency, num_test, num_fail);
    RUN_TEST(check_r2int, num_test, num_fail);
    RUN_TEST(check_r2eff, num_test, num_fail);

    // print how many tests passed or failed
    TEST_SUMMARY(argv[0], num_test, num_fail);

    return (num_fail > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}

//============================================================================//

void check_mean()
{
    ASSERT_CLOSE(conv_obj.GetMean(), stat_obj.GetMean(), epsilon);
}

//----------------------------------------------------------------------------//

void check_relative_error()
{
    ASSERT_CLOSE(conv_obj.GetR(), stat_obj.GetRelativeError(), epsilon);
}

//----------------------------------------------------------------------------//

void check_standard_deviation()
{
    ASSERT_CLOSE(conv_obj.GetStandardDeviation(), stat_obj.GetStdDev(), epsilon);
}

//----------------------------------------------------------------------------//

void check_variance()
{
    ASSERT_CLOSE(conv_obj.GetVariance(), stat_obj.GetVariance(), epsilon);
}

//----------------------------------------------------------------------------//

void check_efficiency()
{
    ASSERT_CLOSE(conv_obj.GetEfficiency(), stat_obj.GetEfficiency(), epsilon);
}

//----------------------------------------------------------------------------//

void check_r2int()
{
    ASSERT_CLOSE(conv_obj.GetR2int(), stat_obj.GetR2Int(), epsilon);
}

//----------------------------------------------------------------------------//

void check_r2eff()
{
    ASSERT_CLOSE(conv_obj.GetR2eff(), stat_obj.GetR2Eff(), epsilon);
}

//============================================================================//
