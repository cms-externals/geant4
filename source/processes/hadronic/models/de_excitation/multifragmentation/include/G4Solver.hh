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
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modification: 13.08.2025 V.Ivanchenko rewrite

#ifndef G4Solver_h
#define G4Solver_h 1

#include "globals.hh"

template <class T_Function> class G4Solver 
{
public:
	
  G4Solver(T_Function* ff, const G4int iterations, const G4double tol)
    : maxIter(iterations), tolerance(tol), tF(ff) {};

  // copy constructor	
  G4Solver(const G4Solver& right) = delete;

  // destructor
  ~G4Solver() = default;
	
  // operators
  G4Solver& operator=(const G4Solver& right) = delete;
  G4bool operator==(const G4Solver& right) const = delete;
  G4bool operator!=(const G4Solver& right) const = delete;

  void SetMaxIterations(const G4int iterations)
  { maxIter = iterations;}

  void SetTolerance(const G4double epsilon)
  { tolerance = epsilon; }

  void SetIntervalLimits(const G4double Limit1, const G4double Limit2)
  {
    aa = std::min(Limit1, Limit2);
    bb = std::max(Limit1, Limit2);
  }

  // Calculates the root of equation Function(x)=0
  // by the Regula-Falsi method
  G4bool FindRoot(G4double& x)
  {
    G4double a = aa;
    G4double b = bb;

    // define accuracy
    G4double epsX = tolerance*(b - a);

    // Check the interval before start
    if (std::abs(a-b) <= epsX) { return true; }

    G4double fa = tF->Function(a);
    if (0.0 == fa) { x = a; return true; }
    G4double fb = tF->Function(b);
    if (0.0 == fb) { x = b; return true; }

    // root should be inside interval
    if (fa*fb > 0.0) { return false; }

    G4double epsY = tolerance*std::min(std::abs(fa), std::abs(fb));

    // Finding the root
    for (G4int i = 0; i < maxIter; ++i)
    {
      cc = (a*fb - b*fa)/(fb - fa);
      G4double fc = tF->Function(cc);
      if (std::abs(fc) < epsY) { x = cc; return true; }

      G4double delta = std::min((cc - a), (b - cc));
      
      if (delta < epsX) { x = cc; return true; }
      else if (fa*fc < 0.0) { b = cc; fb = fc;}
      else { a = cc; fa = fc; }
    }
    x = cc; 
    return false;
  }

private:

  // Maximum number of iterations
  G4int maxIter;
  G4double tolerance;

  // interval limits [a,b] which should bracket the root
  G4double aa{0.0};
  G4double bb{0.0};
  G4double cc{0.0};

  T_Function* tF;
};

#endif
