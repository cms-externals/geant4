#-----------------------------------------------------------------------------#
# - Setup
#-----------------------------------------------------------------------------#
geant4_setup_unit_tests(LINK_LIBRARIES G4globman G4hepnumerics G4heprandom)

#-----------------------------------------------------------------------------#
# - Tests
#-----------------------------------------------------------------------------#
# - Management
geant4_add_unit_test(management/G4Allocator.test.cc)
geant4_add_unit_test(management/G4AllocatorDealloc.test.cc)
geant4_add_unit_test(management/G4AllocatorInSTL.test.cc)
geant4_add_unit_test(management/G4AutoDelete.test.cc)
geant4_add_unit_test(management/G4Cache.test.cc)
geant4_add_unit_test(management/G4OrderedTable.test.cc)
geant4_add_unit_test(management/G4PhysicsTable.test.cc)
geant4_add_unit_test(management/G4ReferenceCountedHandle.test.cc)
geant4_add_unit_test(management/G4ThreadLocalSingleton.test.cc)
geant4_add_unit_test(management/G4UnitDefinition.test.cc)
geant4_add_unit_test(management/G4UnitsTable.test.cc)

# - Numerics
geant4_add_unit_test(HEPNumerics/G4AnalyticalPolSolver.test.cc)
geant4_add_unit_test(HEPNumerics/G4ChebyshevApproximation.test.cc)
geant4_add_unit_test(HEPNumerics/G4DataInterpolation.test.cc)

# Fails with too many iterations...
#geant4_add_unit_test(HEPNumerics/G4GaussianQuadrature.test.cc)
# Doesn't compile
#geant4_add_unit_test(HEPNumerics/G4Integrator.test.cc)

geant4_add_unit_test(HEPNumerics/G4JTPolynomialSolver.test.cc)
geant4_add_unit_test(HEPNumerics/G4PolynomialSolver.test.cc)
geant4_add_unit_test(HEPNumerics/G4SimpleIntegration.test.cc)

# - Random
geant4_add_unit_test(HEPRandom/G4RandomDirection.test.cc)
# Is interactive...
#geant4_add_unit_test(HEPRandom/G4RandomDistribution.test.cc)
# Is interactive...
#geant4_add_unit_test(HEPRandom/G4RandomGlobal.test.cc)

geant4_add_unit_test(HEPRandom/G4RandomMarsaglia.test.cc)

# This is quite long running, so may be more benchmark than unit test
#geant4_add_unit_test(HEPRandom/G4RandomMT.test.cc)
# Same as above
#geant4_add_unit_test(HEPRandom/G4RandomPerformance.test.cc)
# So is this, but though runtime is more reasonable
#geant4_add_unit_test(HEPRandom/G4RandomPool.test.cc)

# Temporary disable this test on Windows when a system CLHEP is used.
# See Merge Request !4864 for details, this MR requiring to be upstreamed
# to CLHEP and Geant4 to minimum require that version before we can remove
# this protection.
if(NOT (WIN32 AND GEANT4_USE_SYSTEM_CLHEP))
  geant4_add_unit_test(HEPRandom/G4ShootBit.test.cc)
endif()

# These tests only run if multithreading is enabled
if(GEANT4_BUILD_MULTITHREADED)
  geant4_add_unit_test(HEPRandom/testRngMT/DualRand.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/HepJamesRandom.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/MixMaxRng.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/MTwistEngine.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/RanecuEngine.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/Ranlux64Engine.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/RanluxEngine.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/RanluxppEngine.test.cc)
  geant4_add_unit_test(HEPRandom/testRngMT/RanshiEngine.test.cc)
endif()
