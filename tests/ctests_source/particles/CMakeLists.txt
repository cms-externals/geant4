#-----------------------------------------------------------------------------#
# - Setup
#-----------------------------------------------------------------------------#
geant4_setup_unit_tests(LINK_LIBRARIES G4partman G4leptons)

#-----------------------------------------------------------------------------#
# - Tests
#-----------------------------------------------------------------------------#
# - Muon decay
geant4_add_unit_test(management/G4MuonDecayChannel.test.cc)

#-----------------------------------------------------------------------------#
# - Micro Benchmarks
#-----------------------------------------------------------------------------#
if(GEANT4_USE_benchmark)
  # We don't add these as a test just yet.
  add_executable(G4MuonDecayChannel_bench management/G4MuonDecayChannel.bench.cc)
  geant4_test_link_libraries(G4MuonDecayChannel_bench
    PRIVATE
      benchmark::benchmark
      G4partman
      G4leptons)
endif()

#-----------------------------------------------------------------------------#
# - Validation Utility
#-----------------------------------------------------------------------------#
find_package(ROOT QUIET COMPONENTS Core Hist Graf Gpad RIO)
if(ROOT_FOUND)
  # Visual validation utility: not registered as a test.
  add_executable(G4MuonDecayChannel_validation management/G4MuonDecayChannel.validation.cc)
  geant4_test_link_libraries(G4MuonDecayChannel_validation
    PRIVATE
      G4partman
      G4leptons
      ROOT::Core
      ROOT::Hist
      ROOT::Graf
      ROOT::Gpad
      ROOT::RIO)
endif()
