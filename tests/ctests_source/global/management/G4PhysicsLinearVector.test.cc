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

#include "G4ExceptionHelper.hh"
#include "G4Filesystem.hh"

#include <G4PhysicsLinearVector.hh>
#include <G4PhysicsVector.hh>
#include <gtest/gtest.h>

#include <sstream>

// This defines a class derived from G4PhysicsLinearVector
// to enable testing of protected properties such as the derivatives

struct G4PhysicsTestLinearVector : public G4PhysicsLinearVector
{
    G4PhysicsTestLinearVector(G4double Emin, G4double Emax, std::size_t Nbin, G4bool spline = false)
      : G4PhysicsLinearVector(Emin, Emax, Nbin, spline)
    {
      ;
    }
    auto& getBinRef() { return binVector; }
    auto& getDataRef() { return dataVector; }
    auto& getDerivRef() { return secDerivative; }
};

// Fixture for tests which need to write to a temporary file
// This ensures a suitable directory is created and removed after
// all tests have run
// If a test does not delete its temp files, they will get destroyed
// at the end
class TemporaryDirectory : public testing::Test
{
  public:

    TemporaryDirectory() { ; }
    static auto getPath()
    {
      auto filepath = G4fs::temp_directory_path();
      const std::string subdir = "geant4-tests";
      filepath /= subdir;
      return filepath;
    }

    // Handles creating and removing the temporary directory
    // Use the getPath method so that we don't have an empty path
    // even if SetUp hasn't been called
    static void SetUpTestSuite() { G4fs::create_directory(getPath()); }
    static void TearDownTestSuite() { G4fs::remove_all(getPath()); }
};

TEST(G4PhysicsLinearVector, ConstructionTest)
{
  G4PhysicsLinearVector vec;
  EXPECT_FALSE(vec.GetSpline());
  EXPECT_EQ(vec.GetVectorLength(), 0);
  EXPECT_EQ(vec.GetType(), G4PhysicsVectorType::T_G4PhysicsLinearVector);
}

/*
NOTE: use 'bins' to refer to the energy bin centres - this implies n+1 nodes or edges.
The values of energy and data are to be specified at the nodes

*/
TEST(G4PhysicsLinearVector, ConstructionWithEnergyRange)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  EXPECT_FALSE(vec.GetSpline());
  EXPECT_EQ(vec.GetVectorLength(), 100 + 1);
  EXPECT_FLOAT_EQ(vec.GetLowEdgeEnergy(0), 1.0);
  EXPECT_FLOAT_EQ(vec.GetLowEdgeEnergy(100), 1001.0);
}

// Force check of the energy bins
//  IMPORTANT: this is skirting along checking the implementation details instead of the
//  contracts, BUT closes some important gaps in how things work
TEST(G4PhysicsLinearVector, ConstructedBins)
{
  G4PhysicsTestLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  auto bins = vec.getBinRef();
  EXPECT_EQ(bins.size(), 101);
  for (std::size_t i = 0; i < 101; i++)
  {
    EXPECT_EQ(bins[i], vec.GetLowEdgeEnergy(i));
  }
}

TEST(G4PhysicsLinearVector, ConstructionUndersized)
{
  // Attempting to create at size 0
  //  Create exception handler - this also sets it as in use
  G4ExceptionHelper::TestExceptionHandler X;
  EXPECT_THROW(G4PhysicsLinearVector vec(1.0, 1001.0, 0),
               G4ExceptionHelper::specifiedException<FatalException>);
  EXPECT_THROW(G4PhysicsLinearVector vec(1.0, 1001.0, 0), std::runtime_error);
}
TEST(G4PhysicsLinearVector, ConstructionWithBadEnergies)
{
  // Attempting to create with negative energy or max < min
  G4ExceptionHelper::TestExceptionHandler X;
  // EXPECT_THROW( G4PhysicsLinearVector vec(-1.0, 1001.0, 10),
  // G4ExceptionHelper::specifiedException<FatalException> );
  EXPECT_THROW(G4PhysicsLinearVector vec(10.0, 1.0, 10),
               G4ExceptionHelper::specifiedException<FatalException>);
}
TEST(G4PhysicsLinearVector, CheckingEnergyRange)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  EXPECT_FLOAT_EQ(vec.GetMinEnergy(), 1.0);
  EXPECT_FLOAT_EQ(vec.GetMaxEnergy(), 1001.0);
}

TEST(G4PhysicsLinearVector, GettingBinIndex)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  std::size_t bin = 0;
  //"No" initial guess
  // Value is in range -> found
  EXPECT_TRUE(vec.CheckIndex(20.0, bin));
  EXPECT_EQ(bin, 1);
  EXPECT_TRUE(vec.CheckIndex(501.0, bin));
  EXPECT_EQ(bin, 50);
  EXPECT_TRUE(vec.CheckIndex(999.0, bin));
  EXPECT_EQ(bin, 99);
}

/*//Direct tests of dependent functions - GetBin and Interpolation - can't as these are private...
TEST(G4PhysicsLinearVector, DirectBinIndex){
    G4PhysicsLinearVector vec{1.0, 1001.0, 100}; // 1 to 1001 with 100 bins
    EXPECT_EQ(vec.GetBin(40.0), 2);

}*/

TEST(G4PhysicsLinearVector, GettingBinIndexWithGuess)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  std::size_t bin = 12;
  // Value is in range -> found
  EXPECT_TRUE(vec.CheckIndex(200.0, bin));
  EXPECT_EQ(bin, 19);
  bin = 27;
  EXPECT_TRUE(vec.CheckIndex(200.0, bin));
  EXPECT_EQ(bin, 19);
}
TEST(G4PhysicsLinearVector, ManuallySettingValues)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }
  EXPECT_FLOAT_EQ(vec.Value(1.0), 0.0);  // Fetching value at _Energy_
  EXPECT_FLOAT_EQ(vec.Value(101.0), 22.0);
}

TEST(G4PhysicsLinearVector, GettingValueAndIndex)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }
  std::size_t bin = 0;
  EXPECT_FLOAT_EQ(vec.Value(1.0, bin), 0.0);  // Fetching value at _Energy_
  EXPECT_EQ(bin, 0);
  EXPECT_FLOAT_EQ(vec.Value(101.0, bin), 22.0);
  EXPECT_EQ(bin, 10);
  EXPECT_FLOAT_EQ(vec.Value(206.0, bin), 45.1);  //(Bin 20: [201.0, 211.0] -> [44.0, 46.2])
  EXPECT_EQ(bin, 20);
}

TEST(G4PhysicsLinearVector, MinAndMaxValue)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {  // Last Value is at index 100
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }
  EXPECT_FLOAT_EQ(vec.GetMinValue(), 0.0);
  EXPECT_FLOAT_EQ(vec.GetMaxValue(), 220.0);
}

TEST(G4PhysicsLinearVector, CopyAssign)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }
  auto vec2 = vec;
  EXPECT_EQ(vec.GetVectorLength(), vec2.GetVectorLength());
  EXPECT_FLOAT_EQ(vec.Energy(10), vec2.Energy(10));
  EXPECT_FLOAT_EQ(vec.Energy(100), vec2.Energy(100));
  EXPECT_NE(vec2.Energy(99), vec2.Energy(100));
}
TEST(G4PhysicsLinearVector, CopyConstruct)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }
  auto vec2{vec};
  EXPECT_EQ(vec.GetVectorLength(), vec2.GetVectorLength());
  EXPECT_FLOAT_EQ(vec.Energy(10), vec2.Energy(10));
  EXPECT_FLOAT_EQ(vec.Energy(100), vec2.Energy(100));
  EXPECT_NE(vec2.Energy(99), vec2.Energy(100));
}

TEST(G4PhysicsLinearVector, Accessors)
{
  // () and []
  // Need idiom for throw/no-throw
}

// TODO - test Retrieve from know file

TEST_F(TemporaryDirectory, G4PhysicsTestLinearVectorStoreRetrieve)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }

  // Simple subdir of temp called geant4-tests/
  // TODO - make this more general?
  // TODO - make this a fixture as may be needed in several places
  auto filename = TemporaryDirectory::getPath();
  filename /= "LinearStoreTest.dat";
  {
    std::ofstream strm;
    strm.open(filename, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    vec.Store(strm);
  }  // strm closes now
  G4PhysicsLinearVector vec_in;
  vec_in.SetDataLength(101);
  {
    std::ifstream strm;
    strm.open(filename, std::ios_base::in | std::ios_base::binary);
    vec_in.Retrieve(strm);
  }
  EXPECT_EQ(vec.GetVectorLength(), vec_in.GetVectorLength());
  EXPECT_EQ(vec_in.GetVectorLength(), 101);
  // Check the restored energy bins
  for (std::size_t i = 0; i < 101; i++)
  {
    EXPECT_FLOAT_EQ(vec.GetLowEdgeEnergy(i), vec_in.GetLowEdgeEnergy(i));
  }
  // Now check the restored values
  //  TODO - exact eq, or float eq?
  for (std::size_t i = 0; i < 101; i++)
  {
    EXPECT_FLOAT_EQ(vec.Value(vec.GetLowEdgeEnergy(i)), vec_in.Value(vec_in.GetLowEdgeEnergy(i)));
  }
  // Finally check a few interpolated values to ensure any flags are correct
  for (double x : {24.0, 72.1, 189.5, 592.4, 766.6})
  {
    EXPECT_FLOAT_EQ(vec.Value(x), vec_in.Value(x));
  }
}

TEST_F(TemporaryDirectory, G4PhysicsTestLinearVectorStoreRetrieveAscii)
{
  G4PhysicsLinearVector vec{1.0, 1001.0, 100};  // 1 to 1001 with 100 bins
  for (std::size_t i = 0; i < 101; i++)
  {
    vec.PutValue(i, (G4float)i * 2.2);  // Assigning monotonic values
  }

  // Simple subdir of temp called geant4-tests/
  // TODO - make this more general?
  // TODO - make this a fixture as may be needed in several places
  auto filename = TemporaryDirectory::getPath();
  filename /= "LinearStoreTest.dat";
  {
    std::ofstream strm;
    strm.open(filename, std::ios_base::out | std::ios_base::trunc);
    vec.Store(strm, true);
  }  // strm closes now
  G4PhysicsLinearVector vec_in;
  vec_in.SetDataLength(101);
  {
    std::ifstream strm;
    strm.open(filename, std::ios_base::in);
    vec_in.Retrieve(strm, true);
  }
  EXPECT_EQ(vec.GetVectorLength(), vec_in.GetVectorLength());
  EXPECT_EQ(vec_in.GetVectorLength(), 101);
  // Check the restored energy bins
  for (std::size_t i = 0; i < 101; i++)
  {
    EXPECT_FLOAT_EQ(vec.GetLowEdgeEnergy(i), vec_in.GetLowEdgeEnergy(i));
  }
  // Now check the restored values
  //  TODO - exact eq, or float eq?
  for (std::size_t i = 0; i < 101; i++)
  {
    EXPECT_FLOAT_EQ(vec.Value(vec.GetLowEdgeEnergy(i)), vec_in.Value(vec_in.GetLowEdgeEnergy(i)));
  }
  // Finally check a few interpolated values to ensure any flags are correct
  for (double x : {24.0, 72.1, 189.5, 592.4, 766.6})
  {
    EXPECT_FLOAT_EQ(vec.Value(x), vec_in.Value(x));
  }
}

// Checking derivs - only used in case of spline
TEST(G4PhysicsLinearVector, DerivsNoSpline)
{
  G4PhysicsTestLinearVector vec(1.0, 1001.0, 100, false);  // 1 to 1001 with 100 bins

  auto derivs = vec.getDerivRef();
  EXPECT_EQ(derivs.size(), 0);
}
// Included for completeness, however for LinearVector all Second
//  derivatives are in fact 0
TEST(G4PhysicsLinearVector, Derivs)
{
  G4PhysicsTestLinearVector vec(1.0, 1001.0, 100, true);  // 1 to 1001 with 100 bins
  vec.FillSecondDerivatives();
  auto derivs = vec.getDerivRef();
  EXPECT_EQ(derivs.size(), 101);
  // Linear, thus all derivs should be equal, and second-derivs are 0
  EXPECT_FLOAT_EQ(derivs[0], 0.0);
  for (std::size_t i = 1; i < 101; i++)
  {
    EXPECT_FLOAT_EQ(derivs[i], derivs[i - 1]);
  }
}
