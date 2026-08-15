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

#include <G4PhysicsVector.hh>
#include <gtest/gtest.h>

/** Helper structure to access protected fields
 *
 * This uses all implementations from G4PhysicsVector but exfiltrates the protected vectors by
 * reference for direct access inside tests. Instantiate a G4PhysicsTestVector instead of a plain
 * G4PhysicsVector only for tests which need this.
 */
struct G4PhysicsTestVector : public G4PhysicsVector
{
    G4PhysicsTestVector() : G4PhysicsVector() { ; }
    auto& getBinRef() { return binVector; }
    auto& getDataRef() { return dataVector; }
    auto& getDerivRef() { return secDerivative; }
};

TEST(G4PhysicsVector, ConstructionTest)
{
  G4PhysicsVector vec;
  EXPECT_FALSE(vec.GetSpline());
  EXPECT_EQ(vec.GetType(), G4PhysicsVectorType::T_G4PhysicsFreeVector);
}

TEST(G4PhysicsVector, Resize)
{
  G4PhysicsVector vec;
  const std::size_t len = 10;
  vec.SetDataLength(len);
  EXPECT_EQ(vec.GetVectorLength(), len);
}
TEST(G4PhysicsVector, ResizeTwice)
{
  G4PhysicsVector vec;
  const std::size_t len = 10;
  vec.SetDataLength(len);
  EXPECT_EQ(vec.GetVectorLength(), len);
  vec.SetDataLength(173);  // EXPECT no-op
  EXPECT_EQ(vec.GetVectorLength(), len);
}
TEST(G4PhysicsVector, ResizeCheckDataLength)
{
  G4PhysicsTestVector vec;
  auto& bins = vec.getBinRef();
  auto& data = vec.getDataRef();
  const std::size_t len = 10;
  vec.SetDataLength(len);
  EXPECT_EQ(vec.GetVectorLength(), len);
  EXPECT_EQ(bins.size(), len);
  EXPECT_EQ(data.size(), len);
  vec.SetDataLength(173);  // EXPECT no-op
  // Check all the items:
  EXPECT_EQ(vec.GetVectorLength(), len);
  EXPECT_EQ(bins.size(), len);
  EXPECT_EQ(data.size(), len);
}

/*TEST(G4PhysicsVector, SetAndRetrieve){
    G4PhysicsVector vec;
    const std::size_t len=10;
    vec.SetDataLength(len);
    vec.

}*/