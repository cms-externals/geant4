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
//-------------------------------------------------------------------//
// Unit test for G4FPEDetection
//-------------------------------------------------------------------//
#include "G4FPEDetection.hh"

#include <gtest/gtest.h>

#include <cmath>

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
class G4FPEDetectionTest : public ::testing::Test
{
  protected:

    // Use volatile variables to prevent compile-time constant folding:
    // the compiler must emit runtime FP instructions.
    volatile double zero = 0.0;
    volatile double max = DBL_MAX;
    volatile double min = DBL_MIN;

    void SetUp() override
    {
// Skip if we don't have FPE handling
#ifndef G4FPE_SIGNAL
      GTEST_SKIP() << "No FPE handling support in this build, skipping";
#endif
      InvalidOperationDetection();
    }
};

TEST_F(G4FPEDetectionTest, InvalidOperation)
{
  EXPECT_DEATH(
    {
      double y = zero / zero;
      std::cout << y << "\n";
    },
    ".*Floating point invalid operation.*");  // IOC trap
}

TEST_F(G4FPEDetectionTest, DivideByZero)
{
  EXPECT_DEATH(
    {
      double y = max / zero;
      std::cout << y << "\n";
    },
    ".*Floating point divide by zero.*");  // DZC trap
}

TEST_F(G4FPEDetectionTest, Overflow)
{
  // Enable explicitly
#ifdef G4FPE_SIGNAL
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_OVERFLOW);
  EXPECT_DEATH(
    {
      double y = max * 2.0;
      std::cout << y << "\n";
    },
    ".*Floating point overflow.*");  // OFC trap

  // Disable and try again
  fedisableexcept(FE_OVERFLOW);
#endif
  EXPECT_TRUE(std::isinf(max * 2.0));
}

TEST_F(G4FPEDetectionTest, Underflow)
{
  // Enable explicitly
#ifdef G4FPE_SIGNAL
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_UNDERFLOW);
  EXPECT_DEATH(
    {
      double y = min / 2.0;
      std::cout << y << "\n";
    },
    ".*Floating point underflow.*");  // UFC trap

  // Disable and try again: result is a subnormal (finite, nonzero, < DBL_MIN)
  fedisableexcept(FE_UNDERFLOW);
#endif
  volatile double y = min / 2.0;
  EXPECT_FALSE(std::isnormal(y));
  EXPECT_GT(y, 0.0);  // subnormal, not zero or negative
}
