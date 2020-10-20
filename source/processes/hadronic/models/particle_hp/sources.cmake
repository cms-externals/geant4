#-----------------------------------------------------------------------
# sources.cmake
# Module : G4had_par_hp
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_par_hp
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#-----------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4had_par_hp
  HEADERS
    G4InterpolationIterator.hh
    G4InterpolationManager.hh
    G4InterpolationScheme.hh
    G4NRESP71M03.hh
    G4ParticleHPList.hh
    G4ParticleHPIsoData.hh
    G4ParticleHPLevel.hh
    G4ParticleHP2AInelasticFS.hh
    G4ParticleHPNames.hh
    G4ParticleHP2N2AInelasticFS.hh
    G4ParticleHPPartial.hh
    G4ParticleHP2NAInelasticFS.hh
    G4ParticleHPProduct.hh
    G4ParticleHP2NDInelasticFS.hh
    G4ParticleHPVector.hh
    G4ParticleHP2NInelasticFS.hh
    G4VParticleHPEDis.hh
    G4ParticleHP2NPInelasticFS.hh
    G4ParticleHP2PInelasticFS.hh
    G4ParticleHP3AInelasticFS.hh
    G4ParticleHP3NAInelasticFS.hh
    G4ParticleHP3NInelasticFS.hh
    G4ParticleHP3NPInelasticFS.hh
    G4ParticleHP4NInelasticFS.hh
    G4ParticleHPAInelasticFS.hh
    G4ParticleHPAngular.hh
    G4ParticleHPAngularP.hh
    G4ParticleHPArbitaryTab.hh
    G4ParticleHPCapture.hh
    G4ParticleHPCaptureData.hh
    G4ParticleHPCaptureFS.hh
    G4ParticleHPChannel.hh
    G4ParticleHPChannelList.hh
    G4ParticleHPContAngularPar.hh
    G4ParticleHPContEnergyAngular.hh
    G4ParticleHPD2AInelasticFS.hh
    G4ParticleHPDAInelasticFS.hh
    G4ParticleHPDInelasticFS.hh
    G4ParticleHPData.hh
    G4ParticleHPDataPoint.hh
    G4ParticleHPDataUsed.hh
    G4ParticleHPDeExGammas.hh
    G4ParticleHPDiscreteTwoBody.hh
    G4ParticleHPElastic.hh
    G4ParticleHPElasticData.hh
    G4ParticleHPElasticFS.hh
    G4ParticleHPElementData.hh
    G4ParticleHPEnAngCorrelation.hh
    G4ParticleHPEnergyDistribution.hh
    G4ParticleHPEvapSpectrum.hh
    G4ParticleHPFCFissionFS.hh
    G4ParticleHPFFFissionFS.hh
    G4ParticleHPFSFissionFS.hh
    G4ParticleHPFastLegendre.hh
    G4ParticleHPField.hh
    G4ParticleHPFieldPoint.hh
    G4ParticleHPFinalState.hh
    G4ParticleHPFission.hh
    G4ParticleHPFissionBaseFS.hh
    G4ParticleHPFissionData.hh
    G4ParticleHPFissionERelease.hh
    G4ParticleHPFissionFS.hh
    G4ParticleHPFissionSpectrum.hh
    G4ParticleHPGamma.hh
    G4ParticleHPHash.hh
    G4ParticleHPHe3InelasticFS.hh
    G4ParticleHPInelastic.hh
    G4ParticleHPInelasticBaseFS.hh
    G4ParticleHPInelasticCompFS.hh
    G4ParticleHPInelasticData.hh
    G4ParticleHPInterpolator.hh
    G4ParticleHPIsotropic.hh
    G4ParticleHPJENDLHEData.hh
    G4ParticleHPJENDLHEElasticData.hh
    G4ParticleHPJENDLHEInelasticData.hh
    G4ParticleHPKallbachMannSyst.hh
    G4ParticleHPLCFissionFS.hh
    G4ParticleHPLabAngularEnergy.hh
    G4ParticleHPLegendreStore.hh
    G4ParticleHPLegendreTable.hh
    G4ParticleHPMadlandNixSpectrum.hh
    G4ParticleHPN2AInelasticFS.hh
    G4ParticleHPN2PInelasticFS.hh
    G4ParticleHPN3AInelasticFS.hh
    G4ParticleHPNAInelasticFS.hh
    G4ParticleHPNBodyPhaseSpace.hh
    G4ParticleHPND2AInelasticFS.hh
    G4ParticleHPNDInelasticFS.hh
    G4ParticleHPNHe3InelasticFS.hh
    G4ParticleHPNInelasticFS.hh
    G4ParticleHPNPAInelasticFS.hh
    G4ParticleHPNPInelasticFS.hh
    G4ParticleHPNT2AInelasticFS.hh
    G4ParticleHPNTInelasticFS.hh
    G4ParticleHPNXInelasticFS.hh
    G4ParticleHPParticleYield.hh
    G4ParticleHPPAInelasticFS.hh
    G4ParticleHPPDInelasticFS.hh
    G4ParticleHPPInelasticFS.hh
    G4ParticleHPPTInelasticFS.hh
    G4ParticleHPPhotonDist.hh
    G4ParticleHPPolynomExpansion.hh
    G4ParticleHPSCFissionFS.hh
    G4ParticleHPSimpleEvapSpectrum.hh
    G4ParticleHPT2AInelasticFS.hh
    G4ParticleHPTCFissionFS.hh
    G4ParticleHPTInelasticFS.hh
    G4ParticleHPThermalBoost.hh
    G4ParticleHPThermalScattering.hh
    G4ParticleHPThermalScatteringData.hh
    G4ParticleHPThermalScatteringNames.hh
    G4ParticleHPWattSpectrum.hh
    G4VParticleHPEnergyAngular.hh
    G4ParticleHPBGGNucleonInelasticXS.hh
    G4ParticleHPManager.hh
    G4ParticleHPThreadLocalManager.hh
    G4ParticleHPReactionWhiteBoard.hh
    G4ParticleHPMessenger.hh
### Fission Fragment Generator - start
	G4ArrayOps.hh
	G4ENDFTapeRead.hh
	G4ENDFYieldDataContainer.hh
	G4FFGDebuggingMacros.hh
	G4FFGDefaultValues.hh
	G4FFGEnumerations.hh
	G4FFGVerboseMacros.hh
	G4FissionFragmentGenerator.hh
	G4FissionProductYieldDist.hh
	G4FPYBiasedLightFragmentDist.hh
	G4FPYNormalFragmentDist.hh
	G4FPYNubarValues.hh
	G4FPYSamplingOps.hh
	G4FPYTreeStructures.hh
	G4ShiftedGaussian.hh
	G4TableTemplate.hh
	G4WendtFissionFragmentGenerator.hh
	G4WattFissionSpectrumValues.hh
### FissionFragment Generator - end
### Headers of NeutronHP for backward compatibility - start
    G4NeutronHPList.hh
    G4NeutronHPIsoData.hh
    G4NeutronHPLevel.hh
    G4NeutronHP2AInelasticFS.hh
    G4NeutronHPNames.hh
    G4NeutronHP2N2AInelasticFS.hh
    G4NeutronHPPartial.hh
    G4NeutronHP2NAInelasticFS.hh
    G4NeutronHPProduct.hh
    G4NeutronHP2NDInelasticFS.hh
    G4NeutronHPVector.hh
    G4NeutronHP2NInelasticFS.hh
    G4VNeutronHPEDis.hh
    G4NeutronHP2NPInelasticFS.hh
    G4NeutronHP2PInelasticFS.hh
    G4NeutronHP3AInelasticFS.hh
    G4NeutronHP3NAInelasticFS.hh
    G4NeutronHP3NInelasticFS.hh
    G4NeutronHP3NPInelasticFS.hh
    G4NeutronHP4NInelasticFS.hh
    G4NeutronHPAInelasticFS.hh
    G4NeutronHPAngular.hh
    G4NeutronHPAngularP.hh
    G4NeutronHPArbitaryTab.hh
    G4NeutronHPCapture.hh
    G4NeutronHPCaptureData.hh
    G4NeutronHPCaptureFS.hh
    G4NeutronHPChannel.hh
    G4NeutronHPChannelList.hh
    G4NeutronHPContAngularPar.hh
    G4NeutronHPContEnergyAngular.hh
    G4NeutronHPD2AInelasticFS.hh
    G4NeutronHPDAInelasticFS.hh
    G4NeutronHPDInelasticFS.hh
    G4NeutronHPData.hh
    G4NeutronHPDataPoint.hh
    G4NeutronHPDataUsed.hh
    G4NeutronHPDeExGammas.hh
    G4NeutronHPDiscreteTwoBody.hh
    G4NeutronHPElastic.hh
    G4NeutronHPElasticData.hh
    G4NeutronHPElasticFS.hh
    G4NeutronHPElementData.hh
    G4NeutronHPEnAngCorrelation.hh
    G4NeutronHPEnergyDistribution.hh
    G4NeutronHPEvapSpectrum.hh
    G4NeutronHPFCFissionFS.hh
    G4NeutronHPFFFissionFS.hh
    G4NeutronHPFSFissionFS.hh
    G4NeutronHPFastLegendre.hh
    G4NeutronHPField.hh
    G4NeutronHPFieldPoint.hh
    G4NeutronHPFinalState.hh
    G4NeutronHPFission.hh
    G4NeutronHPFissionBaseFS.hh
    G4NeutronHPFissionData.hh
    G4NeutronHPFissionERelease.hh
    G4NeutronHPFissionFS.hh
    G4NeutronHPFissionSpectrum.hh
    G4NeutronHPGamma.hh
    G4NeutronHPHash.hh
    G4NeutronHPHe3InelasticFS.hh
    G4NeutronHPInelastic.hh
    G4NeutronHPInelasticBaseFS.hh
    G4NeutronHPInelasticCompFS.hh
    G4NeutronHPInelasticData.hh
    G4NeutronHPInterpolator.hh
    G4NeutronHPIsotropic.hh
    G4NeutronHPJENDLHEData.hh
    G4NeutronHPJENDLHEElasticData.hh
    G4NeutronHPJENDLHEInelasticData.hh
    G4NeutronHPKallbachMannSyst.hh
    G4NeutronHPLCFissionFS.hh
    G4NeutronHPLabAngularEnergy.hh
    G4NeutronHPLegendreStore.hh
    G4NeutronHPLegendreTable.hh
    G4NeutronHPMadlandNixSpectrum.hh
    G4NeutronHPN2AInelasticFS.hh
    G4NeutronHPN2PInelasticFS.hh
    G4NeutronHPN3AInelasticFS.hh
    G4NeutronHPNAInelasticFS.hh
    G4NeutronHPNBodyPhaseSpace.hh
    G4NeutronHPND2AInelasticFS.hh
    G4NeutronHPNDInelasticFS.hh
    G4NeutronHPNHe3InelasticFS.hh
    G4NeutronHPNInelasticFS.hh
    G4NeutronHPNPAInelasticFS.hh
    G4NeutronHPNPInelasticFS.hh
    G4NeutronHPNT2AInelasticFS.hh
    G4NeutronHPNTInelasticFS.hh
    G4NeutronHPNXInelasticFS.hh
    G4NeutronHPNeutronYield.hh
    G4NeutronHPPAInelasticFS.hh
    G4NeutronHPPDInelasticFS.hh
    G4NeutronHPPInelasticFS.hh
    G4NeutronHPPTInelasticFS.hh
    G4NeutronHPPhotonDist.hh
    G4NeutronHPPolynomExpansion.hh
    G4NeutronHPSCFissionFS.hh
    G4NeutronHPSimpleEvapSpectrum.hh
    G4NeutronHPT2AInelasticFS.hh
    G4NeutronHPTCFissionFS.hh
    G4NeutronHPTInelasticFS.hh
    G4NeutronHPThermalBoost.hh
    G4NeutronHPThermalScattering.hh
    G4NeutronHPThermalScatteringData.hh
    G4NeutronHPThermalScatteringNames.hh
    G4NeutronHPWattSpectrum.hh
    G4VNeutronHPEnergyAngular.hh
    G4NeutronHPBGGNucleonInelasticXS.hh
    G4NeutronHPManager.hh
    G4NeutronHPThreadLocalManager.hh
    G4NeutronHPReactionWhiteBoard.hh
    G4NeutronHPMessenger.hh
### Headers of NeutronHP for backward compatibility - end
  SOURCES
    G4InterpolationManager.cc
    G4NRESP71M03.cc
    G4ParticleHPIsoData.cc
    G4ParticleHPLevel.cc
    G4ParticleHP2AInelasticFS.cc
    G4ParticleHPList.cc
    G4ParticleHP2N2AInelasticFS.cc
    G4ParticleHPPartial.cc
    G4ParticleHP2NAInelasticFS.cc
    G4ParticleHPNames.cc
    G4ParticleHP2NDInelasticFS.cc
    G4ParticleHP2NInelasticFS.cc
    G4ParticleHPProduct.cc
    G4ParticleHP2NPInelasticFS.cc
    G4ParticleHPVector.cc
    G4ParticleHP2PInelasticFS.cc
    G4ParticleHP3AInelasticFS.cc
    G4ParticleHP3NAInelasticFS.cc
    G4ParticleHP3NInelasticFS.cc
    G4ParticleHP3NPInelasticFS.cc
    G4ParticleHP4NInelasticFS.cc
    G4ParticleHPAInelasticFS.cc
    G4ParticleHPAngular.cc
    G4ParticleHPArbitaryTab.cc
    G4ParticleHPCapture.cc
    G4ParticleHPCaptureData.cc
    G4ParticleHPCaptureFS.cc
    G4ParticleHPChannel.cc
    G4ParticleHPChannelList.cc
    G4ParticleHPContAngularPar.cc
    G4ParticleHPContEnergyAngular.cc
    G4ParticleHPD2AInelasticFS.cc
    G4ParticleHPDAInelasticFS.cc
    G4ParticleHPDInelasticFS.cc
    G4ParticleHPData.cc
    G4ParticleHPDeExGammas.cc
    G4ParticleHPDiscreteTwoBody.cc
    G4ParticleHPElastic.cc
    G4ParticleHPElasticData.cc
    G4ParticleHPElasticFS.cc
    G4ParticleHPElementData.cc
    G4ParticleHPEnAngCorrelation.cc
    G4ParticleHPFCFissionFS.cc
    G4ParticleHPFFFissionFS.cc
    G4ParticleHPFSFissionFS.cc
    G4ParticleHPFastLegendre.cc
    G4ParticleHPFastLegendre_14.cc
    G4ParticleHPFastLegendre_18.cc
    G4ParticleHPFastLegendre_21.cc
    G4ParticleHPFastLegendre_24.cc
    G4ParticleHPFastLegendre_26.cc
    G4ParticleHPFastLegendre_28.cc
    G4ParticleHPFastLegendre_30.cc
    G4ParticleHPField.cc
    G4ParticleHPFieldPoint.cc
    G4ParticleHPFinalState.cc
    G4ParticleHPFission.cc
    G4ParticleHPFissionBaseFS.cc
    G4ParticleHPFissionData.cc
    G4ParticleHPFissionFS.cc
    G4ParticleHPGamma.cc
    G4ParticleHPHe3InelasticFS.cc
    G4ParticleHPInelastic.cc
    G4ParticleHPInelasticBaseFS.cc
    G4ParticleHPInelasticCompFS.cc
    G4ParticleHPInelasticData.cc
    G4ParticleHPInterpolator.cc
    G4ParticleHPIsotropic.cc
    G4ParticleHPJENDLHEData.cc
    G4ParticleHPJENDLHEElasticData.cc
    G4ParticleHPJENDLHEInelasticData.cc
    G4ParticleHPKallbachMannSyst.cc
    G4ParticleHPLCFissionFS.cc
    G4ParticleHPLabAngularEnergy.cc
    G4ParticleHPLegendreStore.cc
    G4ParticleHPMadlandNixSpectrum.cc
    G4ParticleHPN2AInelasticFS.cc
    G4ParticleHPN2PInelasticFS.cc
    G4ParticleHPN3AInelasticFS.cc
    G4ParticleHPNAInelasticFS.cc
    G4ParticleHPNBodyPhaseSpace.cc
    G4ParticleHPND2AInelasticFS.cc
    G4ParticleHPNDInelasticFS.cc
    G4ParticleHPNHe3InelasticFS.cc
    G4ParticleHPNInelasticFS.cc
    G4ParticleHPNPAInelasticFS.cc
    G4ParticleHPNPInelasticFS.cc
    G4ParticleHPNT2AInelasticFS.cc
    G4ParticleHPNTInelasticFS.cc
    G4ParticleHPNXInelasticFS.cc
    G4ParticleHPPAInelasticFS.cc
    G4ParticleHPPDInelasticFS.cc
    G4ParticleHPPInelasticFS.cc
    G4ParticleHPPTInelasticFS.cc
    G4ParticleHPPhotonDist.cc
    G4ParticleHPSCFissionFS.cc
    G4ParticleHPT2AInelasticFS.cc
    G4ParticleHPTCFissionFS.cc
    G4ParticleHPTInelasticFS.cc
    G4ParticleHPThermalScattering.cc
    G4ParticleHPThermalScatteringData.cc
    G4ParticleHPThermalScatteringNames.cc
    G4ParticleHPWattSpectrum.cc
    G4ParticleHPBGGNucleonInelasticXS.cc
    G4ParticleHPManager.cc
    G4ParticleHPThreadLocalManager.cc
    G4ParticleHPReactionWhiteBoard.cc
    G4ParticleHPMessenger.cc
### Fission Fragment Generator - start
	G4ENDFTapeRead.cc
	G4ENDFYieldDataContainer.cc
	G4FFGDebuggingMacros.cc
	G4FFGVerboseMacros.cc
	G4FissionFragmentGenerator.cc
	G4FissionProductYieldDist.cc
	G4FPYBiasedLightFragmentDist.cc
	G4FPYNormalFragmentDist.cc
	G4FPYSamplingOps.cc
	G4ShiftedGaussian.cc
	G4WendtFissionFragmentGenerator.cc
### Fission Fragment Generator - end
  GRANULAR_DEPENDENCIES
    G4baryons
    G4bosons
    G4geometrymng
    G4globman
    G4had_mod_man
    G4had_mod_util
    G4hadronic_deex_management
    G4hadronic_deex_photon_evaporation
    G4hadronic_deex_util
    G4hadronic_mgt
    G4hadronic_proc
    G4hadronic_util
    G4hadronic_xsect
    G4ions
    G4leptons
    G4materials
    G4mesons
    G4partman
    G4procman
    G4track
    G4volumes
  GLOBAL_DEPENDENCIES
    G4geometry
    G4global
    G4materials
    G4particles
    G4track
  LINK_LIBRARIES
    ${ZLIB_LIBRARIES}
  SOURCES_EXCLUDE_FORMAT
    G4ParticleHPFastLegendre.cc
    G4ParticleHPFastLegendre_14.cc
    G4ParticleHPFastLegendre_18.cc
    G4ParticleHPFastLegendre_21.cc
    G4ParticleHPFastLegendre_24.cc
    G4ParticleHPFastLegendre_26.cc
    G4ParticleHPFastLegendre_28.cc
    G4ParticleHPFastLegendre_30.cc
  )

# List any source specific properties here
if(GEANT4_BUILD_PHP_AS_HP)
  add_definitions(-DPHP_AS_HP)
endif()



