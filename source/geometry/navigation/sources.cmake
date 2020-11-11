#------------------------------------------------------------------------------
# Module : G4navigation
# Package: Geant4.src.G4geometry.G4navigation
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4navigation
  HEADERS
    G4AuxiliaryNavServices.hh
    G4AuxiliaryNavServices.icc
    G4BrentLocator.hh
    G4DrawVoxels.hh
    G4ErrorPropagationNavigator.hh
    G4GeomTestVolume.hh
    G4GeometryMessenger.hh
    G4GlobalMagFieldMessenger.hh
    G4LocatorChangeRecord.hh
    G4LocatorChangeLogger.hh	    
    G4MultiLevelLocator.hh
    G4MultiNavigator.hh
    G4NavigationLogger.hh
    G4Navigator.hh
    G4Navigator.icc
    G4NormalNavigation.hh
    G4NormalNavigation.icc
    G4ParameterisedNavigation.hh
    G4ParameterisedNavigation.icc
    G4PartialPhantomParameterisation.hh
    G4PathFinder.hh
    G4PhantomParameterisation.hh
    G4PhantomParameterisation.icc
    G4PropagatorInField.hh
    G4PropagatorInField.icc
    G4RegularNavigation.hh
    G4RegularNavigationHelper.hh
    G4ReplicaNavigation.hh
    G4ReplicaNavigation.icc
    G4SafetyHelper.hh
    G4SimpleLocator.hh
    G4TransportationManager.hh
    G4TransportationManager.icc
    G4VExternalNavigation.hh
    G4VIntersectionLocator.hh
    G4VIntersectionLocator.icc
    G4VoxelNavigation.hh
    G4VoxelNavigation.icc
    G4VoxelSafety.hh
  SOURCES
    G4AuxiliaryNavServices.cc
    G4BrentLocator.cc
    G4DrawVoxels.cc
    G4ErrorPropagationNavigator.cc
    G4GeomTestVolume.cc
    G4GeometryMessenger.cc
    G4GlobalMagFieldMessenger.cc
    G4LocatorChangeRecord.cc
    G4LocatorChangeLogger.cc    
    G4MultiLevelLocator.cc
    G4MultiNavigator.cc
    G4NavigationLogger.cc
    G4Navigator.cc
    G4NormalNavigation.cc
    G4ParameterisedNavigation.cc
    G4PartialPhantomParameterisation.cc
    G4PathFinder.cc
    G4PhantomParameterisation.cc
    G4PropagatorInField.cc
    G4RegularNavigation.cc
    G4RegularNavigationHelper.cc
    G4ReplicaNavigation.cc
    G4SafetyHelper.cc
    G4SimpleLocator.cc
    G4TransportationManager.cc
    G4VExternalNavigation.cc
    G4VIntersectionLocator.cc
    G4VoxelNavigation.cc
    G4VoxelSafety.cc
  GRANULAR_DEPENDENCIES
    G4geometrymng
    G4globman
    G4graphics_reps
    G4intercoms
    G4magneticfield
    G4materials
    G4volumes
  GLOBAL_DEPENDENCIES
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
)

# List any source specific properties here
