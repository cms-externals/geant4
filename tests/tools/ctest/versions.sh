# functions returning version numbes
get_CMake_version()
{
 # Needed to support latest CDash instance
 local CMAKE_version=3.30.6
 echo $CMAKE_version
}

get_CLHEP_version()
{
 local CLHEP_version=2.4.7.1
 echo $CLHEP_version
}

get_EXTERNALS()
{
 local EXTERNALS=LCG_109a_geant4ext20260605
 echo $EXTERNALS
}

get_LCG_version()
{
 local LCG_version=LCG_109a
 echo $LCG_version
}
