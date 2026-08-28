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
//---------------------------------------------------------------------------//
// Unit tests for G4MuonDecayChannel
//---------------------------------------------------------------------------//

#include "G4MuonDecayChannel.hh"

#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4DecayProducts.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4Positron.hh"

// TODO: Factor the fixed random engine into a global helper library
#include "../../event/G4FixedRandomEngine.hh"
#include "CLHEP/Random/Random.h"
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

//---------------------------------------------------------------------------//
// TEST FIXTURES
//---------------------------------------------------------------------------//
namespace
{
// Tolerances somewhat arbitrary, but larger than that actually used
// for endpoint singularity removal in G4MuonDecay::DecayIt.
// Note that this corresponds to an energy tolerance 10^{-3}eV, so
// extremely tight already.
constexpr G4double kEnergyTolerance = 1e-9;
constexpr G4double kMomentumTolerance = 1e-9;
constexpr G4double kMassShellTolerance = 1e-9;

/** Fixture for testing G4MuonDecayChannel
 *
 * Ensures all needed particle definitions are available.
 * Provides helpers to set up forced RNG sequences for testing kinematic endpoints.
 */
class G4MuonDecayChannelTest : public ::testing::Test
{
  protected:

    void SetUp() override
    {
      previousEngine = CLHEP::HepRandom::getTheEngine();

      // Ensure all particles needed by the channel are defined in the table.
      G4MuonMinus::Definition();
      G4MuonPlus::Definition();
      G4Electron::Definition();
      G4Positron::Definition();
      G4NeutrinoE::Definition();
      G4AntiNeutrinoE::Definition();
      G4NeutrinoMu::Definition();
      G4AntiNeutrinoMu::Definition();
    }

    void TearDown() override
    {
      CLHEP::HepRandom::setTheEngine(previousEngine);
      delete fixedEngine;
      fixedEngine = nullptr;
    }

    /** Set sequence of random numbers to values in input vector
     *
     * For G4MuonDecayChannel, we generally need six random numbers:
     * 1. Sample the electron-neutrino energy
     * 2. Accept/Reject electron-neutrino energy
     * 3. Sample electron energy
     * 4. G4RandomDirection: u value
     * 5. G4RandomDirection: v value
     * 5. Sample azimuth angle for electron-neutrino plane
     *
     * \warning Does not check that values in sequence are in [0,1]
     */
    void SetSequence(const std::vector<double>& sequence)
    {
      delete fixedEngine;
      fixedEngine = new G4FixedRandomEngine(sequence);
      CLHEP::HepRandom::setTheEngine(fixedEngine);
    }

    CLHEP::HepRandomEngine* previousEngine = nullptr;
    G4FixedRandomEngine* fixedEngine = nullptr;
};

/** Return minimum kinematically allowed energy for electron/positron in muon decay
 *
 * Accounts for electron mass and electron neutrino energy
 *
 * \param electron_neutrino_energy energy of the corresponding electron neutrino
 */
G4double MinimumElectronEnergy(G4double electron_neutrino_energy)
{
  const G4double muon_mass = G4MuonMinus::Definition()->GetPDGMass();
  const G4double electron_mass = G4Electron::Definition()->GetPDGMass();

  return 0.5 * muon_mass - electron_neutrino_energy
         + electron_mass * electron_mass / (2.0 * (muon_mass - 2.0 * electron_neutrino_energy));
}

/** Return maximum kinematically allowed energy for electron/positron in muon decay
 *
 * Accounts for electron mass
 */
G4double MaximumElectronEnergy()
{
  const G4double muon_mass = G4MuonMinus::Definition()->GetPDGMass();
  const G4double electron_mass = G4Electron::Definition()->GetPDGMass();

  return (muon_mass * muon_mass + electron_mass * electron_mass) / (2.0 * muon_mass);
}

/** Assert that decay products satisfy kinematic constraints for Muon decay
 *
 * \param products
 * \param parentMass
 */
void CheckDecayKinematicsAtRest(const G4DecayProducts* products, G4double parentMass)
{
  ASSERT_NE(products, nullptr);
  ASSERT_EQ(products->entries(), 3);

  // Total E/p conservation
  G4double totalEnergy = 0.0;
  G4ThreeVector totalMomentum(0.0, 0.0, 0.0);

  // Store electron/electron neutrino energy for Dalitz check
  G4double electronEnergy = 0.0;
  G4double electronNeutrinoEnergy = 0.0;

  for (G4int i = 0; i < products->entries(); ++i)
  {
    const G4DynamicParticle* daughter = (*products)[i];
    ASSERT_NE(daughter, nullptr);

    const G4double e = daughter->GetTotalEnergy();
    const G4double p = daughter->GetTotalMomentum();
    const G4double m = daughter->GetMass();

    ASSERT_TRUE(std::isfinite(e));
    ASSERT_TRUE(std::isfinite(p));
    ASSERT_TRUE(std::isfinite(m));

    EXPECT_GE(e, m);

    const G4double massShellResidual = e * e - p * p - m * m;
    EXPECT_NEAR(massShellResidual, 0.0, kMassShellTolerance);

    totalEnergy += e;
    totalMomentum += daughter->GetMomentum();

    if (daughter->GetDefinition() == G4Electron::Definition()) electronEnergy = e;
    if (daughter->GetDefinition() == G4AntiNeutrinoE::Definition()) electronNeutrinoEnergy = e;
  }

  EXPECT_NEAR(totalEnergy, parentMass, kEnergyTolerance);
  EXPECT_NEAR(totalMomentum.mag(), 0.0, kMomentumTolerance);

  // Check Dalitz constraint, accounting for FP rounding.
  const G4double eMin = MinimumElectronEnergy(electronNeutrinoEnergy);
  const G4double eMax = MaximumElectronEnergy();
  EXPECT_LE(eMin - electronEnergy, kEnergyTolerance);
  EXPECT_LE(electronEnergy - eMax, kEnergyTolerance);
}
}  // namespace

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(G4MuonDecayChannelTest, GlobalFourMomentumConservation)
{
  // Run across N arbitrary iterations to scan a good fraction of the interior phase space
  G4MuonDecayChannel channel("mu-", 1.0);

  for (int i = 0; i < 1000000; ++i)
  {
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, G4MuonMinus::Definition()->GetPDGMass());
    delete products;
  }
}

TEST_F(G4MuonDecayChannelTest, MaximumElectronEnergy)
{
  const G4double muon_mass = G4MuonMinus::Definition()->GetPDGMass();
  const G4double electron_mass = G4Electron::Definition()->GetPDGMass();
  const G4double expected_energy =
    (muon_mass * muon_mass + electron_mass * electron_mass) / (2.0 * muon_mass);

  G4MuonDecayChannel channel("mu-", 1.0);

  // There are two endpoints to check - one for each neutrino
  auto check_electron_maximum = [&]() {
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);
    const G4DynamicParticle* electron = (*products)[0];
    EXPECT_NEAR(electron->GetTotalEnergy(), expected_energy, kEnergyTolerance);
    delete products;
  };

  // 1. Electron anti-neutrino takes its min energy
  SetSequence({0.0, 0.0, 0.0, 0.1, 0.2, 0.3});
  check_electron_maximum();

  // 2. Electron antineutrino takes its max energy
  SetSequence({1.0, 0.0, 0.0, 0.4, 0.3, 0.2});
  check_electron_maximum();
}

TEST_F(G4MuonDecayChannelTest, MinimumElectronEnergy)
{
  const G4double muon_mass = G4MuonMinus::Definition()->GetPDGMass();
  const G4double electron_mass = G4Electron::Definition()->GetPDGMass();
  const G4double expected_energy = electron_mass;
  const G4double nu_bias = 1. / (1 + electron_mass / muon_mass);

  G4MuonDecayChannel channel("mu-", 1.0);
  SetSequence({nu_bias, 0.0, 0.0, 0.4, 0.3, 0.2});
  G4DecayProducts* products = channel.DecayIt(0.0);

  // General on-mass-shell
  CheckDecayKinematicsAtRest(products, muon_mass);

  // Electron must be at rest
  const G4DynamicParticle* electron = (*products)[0];
  EXPECT_NEAR(electron->GetTotalEnergy(), expected_energy, kEnergyTolerance);
  EXPECT_NEAR(electron->GetTotalMomentum(), 0.0, kMomentumTolerance);

  // Neutrinos must be back to back with equal momenta
  const G4DynamicParticle* neutrino_e = (*products)[1];
  const G4DynamicParticle* neutrino_mu = (*products)[2];
  EXPECT_NEAR((neutrino_e->GetMomentum() + neutrino_mu->GetMomentum()).mag(), 0.0,
              kMomentumTolerance);

  delete products;
}

TEST_F(G4MuonDecayChannelTest, EnergyMomentumConservationInEndpoints)
{
  const G4double muon_mass = G4MuonMinus::Definition()->GetPDGMass();
  const G4double electron_mass = G4Electron::Definition()->GetPDGMass();
  const G4double nu_bias = 1. / (1 + electron_mass / muon_mass);

  G4MuonDecayChannel channel("mu-", 1.0);

  // There are three cases where the energy of one daughter can go to zero
  // A threshold is applied here to avoid divisions by zero, giving a small
  // window in which E/p conservation must be carefully preserved

  // 1. Electron momentum < threshold (practically never reachable)
  {
    SetSequence({nu_bias, 0.0, 0.0, 0.4, 0.3, 0.2});
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);
    delete products;
  }

  // 2. Electron Neutrino momentum < threshold
  {
    SetSequence({1e-11, 0.0, 1 - 1e-11, 0.4, 0.3, 0.2});
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);
    delete products;
  }

  // 3. Muon Neutrino momentum < threshold
  {
    SetSequence({1.0 - 1e-11, 0.0, 1 - 1e-11, 0.4, 0.3, 0.2});
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);
    delete products;
  }

  // 4. General case : tests tolerance on cos/sintheta
  {
    // Pathological/regression case identified by Codex:
    SetSequence({0.9999994999883047, 0.0, 0.9999761104444773, 0.4, 0.3, 0.2});
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);
    delete products;
  }

  // 5. Interior point: small electron-neutrino momentum (just above threshold)
  {
    // Keep U_accept = 0.0 so the first y proposal is always accepted
    SetSequence({3e-11, 0.0, 0.5, 0.4, 0.3, 0.2});
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);

    const G4DynamicParticle* neutrino_e = (*products)[1];
    ASSERT_NE(neutrino_e, nullptr);
    const G4double p_nue = neutrino_e->GetTotalMomentum();

    const G4DynamicParticle* electron = (*products)[0];
    const G4DynamicParticle* neutrino_mu = (*products)[2];
    ASSERT_NE(electron, nullptr);
    ASSERT_NE(neutrino_mu, nullptr);
    const G4double p_e = electron->GetTotalMomentum();
    const G4double p_numu = neutrino_mu->GetTotalMomentum();

    // Above singular branch threshold, but still in a narrow near-threshold interior window
    EXPECT_GT(p_nue, kMomentumTolerance);
    // Endpoint windows are very tight in this corner. Keep a small upper bound to
    // ensure this stays close to the singularity without collapsing to the endpoint.
    EXPECT_LT(p_nue, 1e-8);

    // Ensure this remains in the general branch region (none of the daughters in singular windows)
    EXPECT_GT(p_e, kMomentumTolerance);
    EXPECT_GT(p_numu, kMomentumTolerance);

    delete products;
  }

  // 6. Interior point: small electron momentum (just above threshold)
  {
    // U_x = 0 puts x at xmin; small offset from nu_bias lifts Pe above zero
    SetSequence({nu_bias + 4.9e-11, 0.0, 0.0, 0.4, 0.3, 0.2});
    G4DecayProducts* products = channel.DecayIt(0.0);
    CheckDecayKinematicsAtRest(products, muon_mass);

    const G4DynamicParticle* electron = (*products)[0];
    ASSERT_NE(electron, nullptr);
    const G4double p_e = electron->GetTotalMomentum();

    const G4DynamicParticle* neutrino_e = (*products)[1];
    const G4DynamicParticle* neutrino_mu = (*products)[2];
    ASSERT_NE(neutrino_e, nullptr);
    ASSERT_NE(neutrino_mu, nullptr);
    const G4double p_nue = neutrino_e->GetTotalMomentum();
    const G4double p_numu = neutrino_mu->GetTotalMomentum();

    // Above singular branch threshold, but still near-threshold
    EXPECT_GT(p_e, kMomentumTolerance);
    // Near this endpoint the sampled momentum jumps discretely between exact endpoint
    // (p_e = 0) and ~1.0544e-8 MeV for this deterministic stream.
    // Keep the upper bound tight but above that quantized interior value.
    EXPECT_LT(p_e, 1.1e-8);

    // Ensure this remains in the general branch region (none of the daughters in singular windows)
    EXPECT_GT(p_nue, kMomentumTolerance);
    EXPECT_GT(p_numu, kMomentumTolerance);

    delete products;
  }

  // General on-mass-shell
}
