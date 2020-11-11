#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_xsect
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_xsect
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
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4hadronic_xsect
    HEADERS
	G4BGGNucleonElasticXS.hh
	G4BGGNucleonInelasticXS.hh
	G4BGGPionElasticXS.hh
	G4BGGPionInelasticXS.hh
	G4BarashenkovData.hh
	G4ChipsAntiBaryonElasticXS.hh
	G4ChipsAntiBaryonInelasticXS.hh
	G4ChipsComponentXS.hh
	G4ChipsHyperonElasticXS.hh
	G4ChipsHyperonInelasticXS.hh
	G4ChipsKaonMinusElasticXS.hh
	G4ChipsKaonMinusInelasticXS.hh
	G4ChipsKaonPlusElasticXS.hh
	G4ChipsKaonPlusInelasticXS.hh
	G4ChipsKaonZeroElasticXS.hh
	G4ChipsKaonZeroInelasticXS.hh
	G4ChipsNeutronElasticXS.hh
	G4ChipsNeutronInelasticXS.hh
	G4ChipsPionMinusElasticXS.hh
	G4ChipsPionMinusInelasticXS.hh
	G4ChipsPionPlusElasticXS.hh
	G4ChipsPionPlusInelasticXS.hh
	G4ChipsProtonElasticXS.hh
	G4ChipsProtonInelasticXS.hh
	G4ComponentAntiNuclNuclearXS.hh
	G4ComponentBarNucleonNucleusXsc.hh
	G4ComponentGGHadronNucleusXsc.hh
	G4ComponentGGNuclNuclXsc.hh
	G4ComponentSAIDTotalXS.hh
	G4CrossSectionDataSetRegistry.hh
	G4CrossSectionDataStore.hh
	G4CrossSectionElastic.hh
	G4CrossSectionFactory.hh
	G4CrossSectionInelastic.hh
	G4CrossSectionPairGG.hh
	G4ElectroNuclearCrossSection.hh
        G4ElNeutrinoNucleusTotXsc.hh
        G4DiffElasticRatio.hh
	G4EMDissociationCrossSection.hh
	G4EMDissociationSpectrum.hh
	G4GammaNuclearXS.hh
	G4GeneralSpaceNNCrossSection.hh
	G4HadronCaptureDataSet.hh
	G4HadronCrossSections.hh
	G4HadronElasticDataSet.hh
	G4HadronFissionDataSet.hh
	G4HadronInelasticDataSet.hh
	G4HadronNucleonXsc.hh
        G4HadronXSDataTable.hh
	G4IonProtonCrossSection.hh
	G4IonsKoxCrossSection.hh
	G4IonsShenCrossSection.hh
	G4IonsSihverCrossSection.hh
	G4IsotopeList.hh
	G4KokoulinMuonNuclearXS.hh
        G4NeutrinoElectronCcXsc.hh
        G4NeutrinoElectronNcXsc.hh
        G4NeutrinoElectronTotXsc.hh
	G4NeutronCaptureXS.hh
	G4NeutronElasticXS.hh
        G4NeutronElectronElXsc.hh
	G4NeutronInelasticCrossSection.hh
	G4NeutronInelasticXS.hh
	G4NucleonNuclearCrossSection.hh
	G4ParticleInelasticXS.hh
	G4PhotoNuclearCrossSection.hh
	G4PiData.hh
	G4PiNuclearCrossSection.hh
	G4ProtonInelasticCrossSection.hh
	G4TripathiCrossSection.hh
	G4TripathiLightCrossSection.hh
	G4UPiNuclearCrossSection.hh
	G4VComponentCrossSection.hh
	G4VCrossSectionDataSet.hh
	G4VCrossSectionRatio.hh
        G4ZeroXS.hh
	G4CrossSectionFactoryRegistry.hh
	G4FastPathHadronicCrossSection.hh
	G4MuNeutrinoNucleusTotXsc.hh
    SOURCES
	G4BGGNucleonElasticXS.cc
	G4BGGNucleonInelasticXS.cc
	G4BGGPionElasticXS.cc
	G4BGGPionInelasticXS.cc
	G4ChipsAntiBaryonElasticXS.cc
	G4ChipsAntiBaryonInelasticXS.cc
	G4ChipsComponentXS.cc
	G4ChipsHyperonElasticXS.cc
	G4ChipsHyperonInelasticXS.cc
	G4ChipsKaonMinusElasticXS.cc
	G4ChipsKaonMinusInelasticXS.cc
	G4ChipsKaonPlusElasticXS.cc
	G4ChipsKaonPlusInelasticXS.cc
	G4ChipsKaonZeroElasticXS.cc
	G4ChipsKaonZeroInelasticXS.cc
	G4ChipsNeutronElasticXS.cc
	G4ChipsNeutronInelasticXS.cc
	G4ChipsPionMinusElasticXS.cc
	G4ChipsPionMinusInelasticXS.cc
	G4ChipsPionPlusElasticXS.cc
	G4ChipsPionPlusInelasticXS.cc
	G4ChipsProtonElasticXS.cc
	G4ChipsProtonInelasticXS.cc
	G4ComponentAntiNuclNuclearXS.cc
	G4ComponentBarNucleonNucleusXsc.cc
	G4ComponentGGHadronNucleusXsc.cc
	G4ComponentGGNuclNuclXsc.cc
	G4ComponentSAIDTotalXS.cc
	G4CrossSectionDataSetRegistry.cc
	G4CrossSectionDataStore.cc
	G4CrossSectionElastic.cc
	G4CrossSectionInelastic.cc
	G4CrossSectionPairGG.cc
        G4DiffElasticRatio.cc
	G4ElectroNuclearCrossSection.cc
        G4ElNeutrinoNucleusTotXsc.cc
	G4EMDissociationCrossSection.cc
	G4EMDissociationSpectrum.cc
	G4GammaNuclearXS.cc
	G4GeneralSpaceNNCrossSection.cc
	G4HadronCaptureDataSet.cc
	G4HadronCrossSections.cc
	G4HadronElasticDataSet.cc
	G4HadronFissionDataSet.cc
	G4HadronInelasticDataSet.cc
	G4HadronNucleonXsc.cc
        G4HadronXSDataTable.cc
	G4IonProtonCrossSection.cc
	G4IonsKoxCrossSection.cc
	G4IonsShenCrossSection.cc
	G4IonsSihverCrossSection.cc
	G4KokoulinMuonNuclearXS.cc
        G4NeutrinoElectronCcXsc.cc
        G4NeutrinoElectronNcXsc.cc
        G4NeutrinoElectronTotXsc.cc
	G4NeutronCaptureXS.cc
	G4NeutronElasticXS.cc
        G4NeutronElectronElXsc.cc
	G4NeutronInelasticCrossSection.cc
	G4NeutronInelasticXS.cc
	G4NucleonNuclearCrossSection.cc
	G4ParticleInelasticXS.cc
	G4PhotoNuclearCrossSection.cc
	G4PiData.cc
	G4PiNuclearCrossSection.cc
	G4ProtonInelasticCrossSection.cc
	G4TripathiCrossSection.cc
	G4TripathiLightCrossSection.cc
	G4UPiNuclearCrossSection.cc
	G4VComponentCrossSection.cc
	G4VCrossSectionDataSet.cc
	G4VCrossSectionRatio.cc
        G4ZeroXS.cc
	G4CrossSectionFactoryRegistry.cc
	G4FastPathHadronicCrossSection.cc
	G4MuNeutrinoNucleusTotXsc.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_util
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
)

# List any source specific properties here

