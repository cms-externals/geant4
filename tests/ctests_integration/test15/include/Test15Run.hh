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
/// \file biasing/Test15/include/Test15Run.hh
/// \brief Definition of the Test15Run class
//
//
//
//---------------------------------------------------------------------
// (Purpose)
//    Example implementation for multi-functional-detector and
//   primitive scorer.
//    This Test15Run class has collections which accumulate
//   a event information into a run information.
//
//---------------------------------------------------------------------

#ifndef Test15Run_h
#  define Test15Run_h 1

#  include "G4Event.hh"
#  include "G4Run.hh"
#  include "G4THitsMap.hh"
#  include "globals.hh"

#  include <vector>
//
class Test15Run : public G4Run
{
  public:

    // constructor and destructor.
    //  vector of multifunctionaldetector name has to given to constructor.
    // Test15Run(const std::vector<G4String> mfdName);
    Test15Run();
    virtual ~Test15Run();

  public:

    // virtual method from G4Run.
    // The method is overriden in this class for scoring.
    virtual void RecordEvent(const G4Event*);

    // Access methods for scoring information.
    // - Number of HitsMap for this RUN.
    //   This is equal to number of collections.
    G4int GetNumberOfHitsMap() const { return (G4int)fRunMap.size(); }
    // - Get HitsMap of this RUN.
    //   by sequential number, by multifucntional name and collection name,
    //   and by collection name with full path.
    G4THitsMap<G4double>* GetHitsMap(G4int idx) { return fRunMap[idx]; }
    G4THitsMap<G4double>* GetHitsMap(const G4String& detName, const G4String& colName);
    G4THitsMap<G4double>* GetHitsMap(const G4String& fullName);
    // - Dump All HitsMap of this RUN.
    //   This method calls G4THisMap::PrintAll() for individual HitsMap.
    void DumpAllScorer();

    virtual void Merge(const G4Run*);

    void FillPerEvent(G4double);
    void AddExitingFlux(G4double);
    void AddExitingGrichineFlux();
    void AddExitingCheckFlux();

    void AddFlux(const G4String&);

    void analyseNeutronFlux(G4double energy, G4double time, G4double startEnergy, G4int TrackID,
                            G4int ParentID, G4double zMomentum, G4double startTime, G4double radius,
                            G4double zPos, G4double parentEnergy, G4String parentParticle,
                            G4double cos_angle, G4int number_generations, G4String Particle,
                            G4bool reduced_tally);

    void analyseNeutronShellFluence(G4double energy, G4double time, G4double startEnergy,
                                    G4int TrackID, G4int ParentID, G4double zMomentum,
                                    G4double startTime, G4double radius, G4double zPos,
                                    G4double parentEnergy, G4String parentParticle,
                                    G4double steplength, G4bool enter_sph, G4bool enter_cyl,
                                    G4bool exit_sph, G4bool exit_cyl, G4String Volume,
                                    G4bool enter_sph_front, G4bool exit_sph_front,
                                    G4int preParentReplica, G4int postParentReplica,
                                    G4int preReplica, G4int postReplica);

    void analyseNeutronRadialFluence(G4double, G4double, G4double, G4int);

    G4int GetExitingFlux() const { return exiting_flux; }
    G4int GetExitingGrichineFlux() const { return exitinggrichine_flux; }
    G4int GetExitingCheckFlux() const { return exiting_check_flux; }

    G4double GetExitingEnergy() const { return exiting_energy; }
    G4double GetIntegralFlux_46cm() const { return integral_flux_46cm; }
    G4double GetIntegralEFlux_46cm() const { return integral_Eflux_46cm; }

  private:

    std::vector<G4String> fCollName;
    std::vector<G4int> fCollID;
    std::vector<G4THitsMap<G4double>*> fRunMap;

    //--------------- Analysis Variables - should be RunAction or Run?

  public:

    G4int total_flux;
    G4double flux_energy[33];
    G4double fluence_spectrum[1000];
    G4int n_max;
    G4double radii[10];
    G4double radii_energies[10];
    G4double flux_radius[10][10];
    G4double fine_energy[65];
    G4double flux_data[32];
    G4double eflux_data[32];
    G4double flux_stat_error[32];
    G4double flux_syst_error[32];
    G4double flux[32];
    G4double cos_flux[32];
    G4double fluence[32];
    G4double fluence_step[32];
    G4double fluence_front_step[32];
    G4double fluence_cyl[32];
    G4double fluence_step_cyl[32];
    G4double fluence_step_shell[32];
    G4double eflux[32];
    G4double fine_eflux[64];
    G4double energy_integral[4];
    G4double enflux[4];
    G4int neutflux[4];

    //  G4double low_energy[79];
    G4double low_energy[101];
    G4double low_flux_data[100];
    //  G4double low_flux[78];
    G4double low_flux[100];
    G4double cos_low_flux[100];
    G4double low_fluence[100];
    G4double low_fluence_step[100];
    G4double low_fluence_front_step[100];
    G4double low_fluence_cyl[100];
    G4double low_fluence_step_cyl[100];
    G4double low_fluence_step_shell[100];
    G4double low_stat[100];
    G4double low_syst[100];

    G4double lithium_radial_energy_lower[10];
    G4double lithium_radial_energy_upper[10];
    G4double lithium_radial_mean[10];
    G4double lithium_radial_true_mean[10];
    G4double radial_fluence_step[26][10];

    G4double lithium_energy[101];
    G4double lithium_flux_data[100];
    G4double cos_lithium_flux[100];
    G4double lithium_flux[100];
    G4double lithium_fluence[100];
    G4double lithium_fluence_step[100];
    G4double lithium_fluence_front_step[100];
    G4double lithium_fluence_cyl[100];
    G4double lithium_fluence_step_cyl[100];
    G4double lithium_fluence_step_shell[100];
    G4double lithium_Zflux[100];
    G4double lithium_flux_5cm[100];
    G4double lithium_stat[100];
    G4double lithium_syst[100];

    G4double eflux_integral;
    G4double lithium_integral_data;
    G4double lithium_Eintegral_data;
    G4double integral_flux_5cm;
    G4double integral_flux_10cm;
    G4double integral_flux_46cm;
    G4double integral_Eflux_46cm;
    G4double integral_Eflux_46cm_restricted;
    G4double integral_Zflux_46cm;
    G4double integral_flux_70cm;
    G4double integral_flux_100cm;
    G4double integral_flux_120cm;

    // G4double energy; G4double name;

    G4double time;

    G4double neutron_energy, neutron_time;

    G4int gamma_flux, neutron_flux, neutron_check, electron_flux, piminus_flux, piplus_flux,
      pizero_flux, positron_flux, proton_flux, muon_flux, other_flux, exiting_flux,
      exitinggrichine_flux, exiting_check_flux, neutron_fluence, neutron_fluence_46cm,
      neutron_fluence_cyl;

    G4double exiting_energy;

    G4int integral_scintillation;
    G4double integral_scintillation_E;
    G4int integral_lithium;
    G4double integral_lithium_E;
    G4int integral_helium;
    G4double integral_helium_E;

    G4int duplicate_neutron;
    G4int oldTrackID;
    G4int duplicate_neutron2;
    G4int oldTrackID2;

    G4double fractional_bin_width;
};

inline void Test15Run::FillPerEvent(G4double) {}

inline void Test15Run::AddExitingFlux(G4double ex_energy)
{
  exiting_flux++;
  exiting_energy += ex_energy;
}

inline void Test15Run::AddExitingGrichineFlux()
{
  exitinggrichine_flux++;
}

inline void Test15Run::AddExitingCheckFlux()
{
  exiting_check_flux++;
}

inline void Test15Run::AddFlux(const G4String& particleName)
{
  if (particleName == "gamma") gamma_flux++;
  if (particleName == "neutron") neutron_flux++;
  if (particleName == "e-") electron_flux++;
  if (particleName == "pi-") piminus_flux++;
  if (particleName == "pi+") piplus_flux++;
  if (particleName == "pi0") pizero_flux++;
  if (particleName == "e+") positron_flux++;
  if (particleName == "proton") proton_flux++;
  if (particleName == "mu-") muon_flux++;
  if (particleName == "mu+") muon_flux++;
  if (particleName == "other") other_flux++;
  if (particleName == "neutron_check") neutron_check++;
  if (particleName == "neutron_fluence") neutron_fluence++;
}

//

#endif
