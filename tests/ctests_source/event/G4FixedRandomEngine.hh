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

// G4FixedRandomEngine.hh
// A CLHEP::HepRandomEngine that returns a fixed sequence of numbers for
// deterministic testing when sequence/coverage needs to be explicitly set.

#ifndef G4FIXED_RANDOM_ENGINE_HH
#define G4FIXED_RANDOM_ENGINE_HH

#include "CLHEP/Random/RandomEngine.h"

#include <stdexcept>
#include <vector>

class G4FixedRandomEngine : public CLHEP::HepRandomEngine
{
  public:

    G4FixedRandomEngine(const std::vector<double>& sequence) : sequence_(sequence), index_(0) {}

    double flat() override
    {
      if (index_ >= sequence_.size())
      {
        throw std::out_of_range("G4FixedRandomEngine: Ran out of random numbers");
      }
      return sequence_[index_++];
    }

    void flatArray(int size, double* vect) override
    {
      for (int i = 0; i < size; ++i)
      {
        vect[i] = flat();
      }
    }

    void setSeed(long, int) override {}
    void setSeeds(const long*, int) override {}
    void saveStatus(const char*) const override {}
    void restoreStatus(const char*) override {}
    void showStatus() const override {}
    std::string name() const override { return "G4FixedHepRandomEngine"; }

  private:

    std::vector<double> sequence_;
    size_t index_;
};

#endif  // G4FIXED_RANDOM_ENGINE_HH
