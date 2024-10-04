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

// Class template for accumulable vectors handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 12/07/2024

#ifndef G4AccVector_h
#define G4AccVector_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

#include <vector>

// using vector_std::size_t = std::vector::std::size_t;

template <class T, class Allocator = std::allocator<T>>
class G4AccVector : public G4VAccumulable,
                    public std::vector<T, Allocator>
{
  public:
    // ctors to be supported
    // (https://en.cppreference.com/w/cpp/container/vector/vector)
    // 1) Default constructor. Constructs an empty container with a default-constructed allocator.
    // vector() noexcept(noexcept(Allocator()));
    // 3) Constructs the container with count copies of elements with value value.
    // vector( std::size_t count,
    //         const T& value,
    //         const Allocator& alloc = Allocator() );
    // 4) Constructs the container with count default-inserted instances of T. No copies are made.
    // vector( std::size_t count,
    //         const Allocator& alloc = Allocator() );
    // 6) Copy constructor. Constructs the container with the copy of the contents of other.
    // vector( const vector& other );
    // 8) Move constructor. Constructs the container with the contents of other using move semantics.
    // vector( vector&& other );
    // 10) Constructs the container with the contents of the initializer list init.
    // vector( std::initializer_list<T> init,
    //    const Allocator& alloc = Allocator() );

    // Default constructor (1)
    // Constructs an empty container with a default-constructed allocator.
    G4AccVector(const G4String& name = "",
                G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (2)
    // Constructs an empty container with the given allocator alloc.
    G4AccVector(const Allocator& alloc,
                G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (2) with name
    // Constructs an empty container with the given allocator alloc.
    G4AccVector(const G4String& name,
                const Allocator& alloc,
                G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (3)
    // Constructs the container with count copies of elements with value value
    // with a default-constructed allocator.
    // G4AccVector(std::size_t count, const T& value,
    //                     const Allocator& alloc = Allocator());
    G4AccVector(std::size_t count, const T& value,
                G4MergeMode mergeMode = G4MergeMode::kAddition,
                const Allocator& allocator = Allocator());

    // Constructor (3) with name
    // Constructs the container with count copies of elements with value value
    // with a default-constructed allocator.
    // G4AccVector(std::size_t count, const T& value,
    //                     const Allocator& alloc = Allocator());
    G4AccVector(const G4String& name,
                std::size_t count, const T& value,
                G4MergeMode mergeMode = G4MergeMode::kAddition,
                const Allocator& allocator = Allocator());
    // Constructor (4)
    // Constructs the container with count default-inserted instances of T.
    // No copies are made.
    // G4AccVector(std::size_t count,
    //                     const Allocator& alloc = Allocator() );
    G4AccVector(std::size_t count,
                G4MergeMode mergeMode = G4MergeMode::kAddition,
                const Allocator& allocator = Allocator());

    // Constructor (4) with name
    // Constructs the container with count default-inserted instances of T.
    // No copies are made.
    // G4AccVector(std::size_t count,
    //                     const Allocator& alloc = Allocator() );
    G4AccVector(const G4String& name,
                std::size_t count,
                G4MergeMode mergeMode = G4MergeMode::kAddition,
                const Allocator& allocator = Allocator());

    // Constructor (10)
    // Constructs the container with the contents of the initializer list init.
    // G4AccVector(std::initializer_list<T> init,
    //                     const Allocator& alloc = Allocator() );
    G4AccVector(std::initializer_list<T> init,
                G4MergeMode mergeMode = G4MergeMode::kAddition,
                const Allocator& allocator = Allocator());

    // Constructor (10) with name
    // Constructs the container with the contents of the initializer list init.
    // G4AccVector(std::initializer_list<T> init,
    //                     const Allocator& alloc = Allocator() );
    G4AccVector(const G4String& name,
                std::initializer_list<T> init,
                G4MergeMode mergeMode = G4MergeMode::kAddition,
                const Allocator& allocator = Allocator());

    // Copy constructor
    G4AccVector(const G4AccVector& rhs) = default;
    G4AccVector(const G4AccVector& rhs, const Allocator& allocator);
    // Move constructor
    G4AccVector(G4AccVector&& rhs) = default;
    G4AccVector(G4AccVector&& rhs, const Allocator& allocator);

    // Destructor
    ~G4AccVector() override = default;

    // Methods
    void Merge(const G4VAccumulable& other) final;
    void Reset() final;
    void Print(G4PrintOptions options = G4PrintOptions()) const final;
    void SetMergeMode(G4MergeMode value) final;

    // Get methods
    G4AccType GetType() const final { return G4AccType::kVector; }

  private:
    // Data members
    T fInitValue = 0;
    G4MergeFunction<T> fMergeFunction;
 };

// inline functions

#include "G4AccVector.icc"

#endif
