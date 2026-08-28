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

#include <vector>
#include <list>
#include <memory>

#include <G4Types.hh>

#ifndef G4GIDI_hh_included
#define G4GIDI_hh_included 1

#define channelID std::string

class G4GIDI_Product {

    public:
        int A, Z, m;
        double kineticEnergy, px, py, pz;
        double birthTimeSec;
};

class G4GIDI_target {

    private:
        struct Impl;
        std::unique_ptr<Impl> m_impl;

        G4GIDI_target( std::unique_ptr<Impl> a_impl );

        friend class G4GIDI;

    public:
        ~G4GIDI_target( );

        std::string const *getName( ) const ;
        std::string const *getFilename( ) const ;
        std::string const *getEvaluation( ) const ;
        int getZ( ) const ;
        int getA( ) const ;
        int getM( ) const ;
        double getMass( ) const ;

//        int getTemperatures( double *a_temperatures ) const ;
//        int readTemperature( int index );

//        std::string getEqualProbableBinSampleMethod( );
//        int setEqualProbableBinSampleMethod( std::string const &a_method );

        std::vector<int> const &elasticIndices( ) ;
        std::vector<int> const &captureIndices( ) ;
        std::vector<int> const &fissionIndices( ) ;
        std::vector<int> const &othersIndices( ) ;

        int getNumberOfChannels( ) const ;
        int getNumberOfProductionChannels( ) const ;
        channelID getChannelsID( int channelIndex ) const ;
        std::vector<channelID> *getChannelIDs( ) const ;
        std::vector<channelID> *getProductionChannelIDs( ) const ;

//        std::vector<double> *getEnergyGridAtTIndex( int index );

        double getTotalCrossSectionAtE( double a_energy, double a_temperature ) const ;
        double getElasticCrossSectionAtE( double a_energy, double a_temperature ) const ;
        double getCaptureCrossSectionAtE( double a_energy, double a_temperature ) const ;
        double getFissionCrossSectionAtE( double a_energy, double a_temperature ) const ;
        double getOthersCrossSectionAtE( double a_energy, double a_temperature ) const ;
        double sumChannelCrossSectionAtE( std::vector<int> const &a_indices, double a_energy, double a_temperature ) const ;
        double sumChannelCrossSectionAtE( int a_nIndices, int const *a_indices, double a_energy, double a_temperature ) const ;
        int sampleChannelCrossSectionAtE( std::vector<int> const &a_indices, double a_energy, double a_temperature, 
                        double (*a_rng)( void * ), void *a_rngState ) const ;
        int sampleChannelCrossSectionAtE( int a_nIndices, int const *a_indices, double a_energy, double a_temperature, 
                        double (*a_rng)( void * ), void *a_rngState ) const ;

        double getElasticFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const ;
        std::vector<G4GIDI_Product> *getCaptureFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const ;
        std::vector<G4GIDI_Product> *getFissionFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const ;
        std::vector<G4GIDI_Product> *getOthersFinalState( double a_energy, double a_temperature, double (*a_rng)( void * ), void *a_rngState ) const ;
        std::vector<G4GIDI_Product> *getFinalState( std::vector<int> const &a_indices, double a_energy, double a_temperature, 
                        double (*a_rng)( void * ), void *a_rngState ) const ;
        std::vector<G4GIDI_Product> *getFinalState( int a_nIndices, int const *a_indices, double a_energy, double a_temperature, 
                        double (*a_rng)( void * ), void *a_rngState ) const ;

//        double getReactionsThreshold( int a_index ) const ;
//        void getReactionsDomain( int a_index, double *a_EMin, double *a_EMax ) const ;
};

class G4GIDI {

    private:
        struct Impl;

        std::unique_ptr<Impl> m_impl;

    public:
        G4GIDI( G4int a_ip, std::string const &a_dataDirectory );
        G4GIDI( G4int a_ip, std::list<std::string> const &a_dataDirectory );
        ~G4GIDI( );

        G4int projectileIP( ) const ;

        G4int numberOfDataDirectories( ) const ;
        G4int addDataDirectory( std::string const &a_dataDirectory );
        G4int removeDataDirectory( std::string const &a_dataDirectory );
        std::string const getDataDirectoryAtIndex( G4int a_index ) const ;
        std::vector<std::string> *getDataDirectories( ) const ;

        bool isThisDataAvailable( std::string const &a_lib_name, G4int a_Z, G4int a_A, G4int a_M = 0 ) const ;
        bool isThisDataAvailable( std::string const &a_lib_name, std::string const &a_targetName ) const ;

        std::string dataFilename( std::string const &lib_name, G4int a_Z, G4int a_A, G4int a_M = 0 ) const ;
        std::string dataFilename( std::string const &lib_name, std::string const &a_targetName ) const ;

        std::vector<std::string> *getNamesOfAvailableLibraries( G4int a_Z, G4int a_A, G4int a_M = 0 ) const ;
        std::vector<std::string> *getNamesOfAvailableLibraries( std::string const &a_targetName ) const ;

        std::vector<std::string> *getNamesOfAvailableTargets( ) const ;

        G4GIDI_target *readTarget( std::string const &lib_name, G4int a_Z, G4int a_A, G4int a_M = 0, bool a_bind = true );
        G4GIDI_target *readTarget( std::string const &lib_name, std::string const &a_targetName, bool a_bind = true );

        G4GIDI_target *getAlreadyReadTarget( G4int a_Z, G4int a_A, G4int a_M = 0 );
        G4GIDI_target *getAlreadyReadTarget( std::string const &a_targetName );

        G4int freeTarget( G4int a_Z, G4int a_A, G4int a_M = 0 );
        G4int freeTarget( std::string const &a_targetSymbol );
        G4int freeTarget( G4GIDI_target *a_target );

        std::vector<std::string> *getListOfReadTargetsNames( );
};

std::string G4GIDI_version( );
int G4GIDI_versionMajor( );
int G4GIDI_versionMinor( );
int G4GIDI_versionPatchLevel( );
std::string G4GIDI_GitHash( );
void G4GIDI_initialize( std::string const &a_dataPath );
std::string G4GIDI_Misc_Z_toSymbol( int a_Z );
std::string G4GIDI_Misc_Z_A_m_ToName( int a_Z, int a_A, int a_M );

#endif      // End of G4GIDI_hh_included
