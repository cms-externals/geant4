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
// G4VSIntegration
//
// Created 03.03.2025 V.Ivanchenko
//
// Class description:
//
// Numerical algorithm for integration of a 1-D probability density function
// and sampling of a value. Parameters of the algorithm should be defined by
// the consumer class via the InitialiseIntegrator(..) method. The algorithm
// is effective for the case of functions with a peak and long falling tail.
// Tunning of parameters is needed to gurantee efficiency of the algorithm.
// The default set of parameters is optimized for the pre-compound model.
//
// The method subdivide energy csale on 3 energy intervals: the 1st expect
// to have maximum (Pmax) of the probability density function inside; the 2nd
// interval starts from energy, at which the probability density function
// become less than the maximum value miltiplied by a factor; the 3d energy
// interval starts from energy, at which the probability density function is less
// than the maximal value multiplied by the factor squared. 
//
// Class methods:
//   ProbabilityDensityFunction(G4double e) - pure virtual mandatory user method 
//
//   InitialiseIntegrator(G4double acc, G4double fact1,
//                        G4double fact2, G4double de,
//                        G4double dmin, G4double dmax) - optinal method,
//       which allows optimisation of integration and sampling.
//       Parameters:
//         acc      - (1.e-8 < acc < 0.2) required accuracy of numerical integration
//         fact1    - (0.01 < fact1 < 1) defines the end energy E1 of the first
//                    energy interval: Pmax*fact1 > P(E1)
//         fact2    - (1 < fact2 < 5) provides a margine for max cross section in
//                    each energy interval
//         deltaE   - (deltaE > 0) the initial step in energy for integration
//         dmin     - (0 < dmin < deltaE) - minimim value of energy step for
//                    integration procedure
//         dmax     - (deltaE < dmax) - maximum value of energy step for
//                    integration procedure
//      Integration step is dynamic to improve method performance.
//      The integrator may be initialized many times in run time.
//
//   ComputeIntegral(const G4double emin, const G4double emax) is the MANDATORY
//      method, which should be performed before each new sampling.
//
//   SampleValue() return random value according to probability density
//      functions.
//
//   The usage of this utility is an advantage if in run time parameters of the
//   probability density function are permanetly different.
//
//   The method CANNOT BE APPLIED on an arbitrary probability density function.
//   If this function has narrow peaks, then the preformance of the method is not
//   guranteed. For each case tuning of parameters of the method is required.
//
//   Computations in this integrator do not assume any physical units, consumer
//   code is responsible to provide all values coherently.
//
// --------------------------------------------------------------------

#ifndef G4VSINTEGRATION_HH
#define G4VSINTEGRATION_HH

#include "globals.hh"

class G4VSIntegration
{
  public:

    G4VSIntegration() = default;

    virtual ~G4VSIntegration() = default;

    virtual G4double ProbabilityDensityFunction(G4double) = 0;

    virtual const G4String& ModelName() const;

    void InitialiseIntegrator(G4double accuracy, G4double fact1, G4double fact2,
                              G4double de, G4double dmin, G4double dmax);

    G4double ComputeIntegral(const G4double emin, const G4double emax);

    G4double SampleValue();

    G4VSIntegration(const G4VSIntegration&) = delete;
    G4VSIntegration& operator=(const G4VSIntegration&) = delete;
    G4bool operator==(const G4VSIntegration& right) const = delete;
    G4bool operator!=(const G4VSIntegration& right) const = delete;

    void SetVerbose(G4int verb) { fVerbose = verb; }

  private:

    G4double fAcc{0.001};  // accuracy of integration
    G4double fMinDelta{0.1};  // minimal step integration
    G4double fMaxDelta{2.0};  // maximal step integration
    G4double fDelta{1.0};  // the default step
    G4double fFactor1{0.25};
    G4double fFactor2{1.05};

    // parameters describing function
    G4double fEmin{0.0};
    G4double fEmax{0.0};
    G4double fE1{0.0};
    G4double fP1{0.0};
    G4double fE2{0.0};
    G4double fP2{0.0};
    G4double fE3{0.0};
    G4double fP3{0.0};
    G4double fPmax{0.0};

    G4int fVerbose{0};
    G4int fWarnLimit{4};
    G4int fnWarn{0};

    G4String dummy{""};
};

#endif
