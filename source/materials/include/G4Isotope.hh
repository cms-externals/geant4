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

// class description
//
// An isotope is a chemical isotope defined by its name,
//                                                 Z: atomic number
//                                                 N: number of nucleons
//                                                 A: mass of a mole (optional)
//                                                 m: isomer state (optional)
// If A is not defined it is taken from Geant4 database
//
// The class contains as a private static member the table of defined
// isotopes (an ordered vector of isotopes).
//
// Isotopes can be assembled into elements via the G4Element class.
//
// 20.08.11: Add flag fm for isomer level (mma)
// 15.11.05: GetIsotope(isotopeName, G4bool warning=false)
// 31.03.05: A becomes optional. Taken from Nist data base by default (mma)
// 26.02.02: fIndexInTable renewed
// 14.09.01: fCountUse: nb of elements which use this isotope
// 13.09.01: stl migration. Suppression of the data member fIndexInTable
// 30.03.01: suppression of the warning message in GetIsotope
// 04.08.98: new method GetIsotope(isotopeName) (mma)
// 17.01.97: aesthetic rearrangement (mma)

#ifndef G4ISOTOPE_HH
#define G4ISOTOPE_HH

#include "G4ios.hh"
#include "globals.hh"

#include <vector>

class G4Isotope;
using G4IsotopeTable = std::vector<G4Isotope*>;

/** A chemical isotope, defined by its name, Z: atomic number and N: number of nucleons
 * 
 * Optionally, one can define A, the mass of a mole, otherwise the Geant4 database value is used; and m, the isomer level.
 * 
 * Once created, an Isotope is added to the internal isotope table.
 * 
 * Isotopes can be assembled into elements via the G4Element class.
 * 
 * IMPORTANT: as with G4Element, user code should not delete isotopes once created, as this will invalidate their G4IsotopeTable entry. Isotopes are deleted at the end of a Geant4 run.
 */
class G4Isotope
{
  public:  // with description

    /**
     * Create an isotope
     * @param name Isotope name
     * @param z Atomic number
     * @param n Number of nucleons
     * @param a Mass of mole (optional, uses database if not provided)
     * @param mlevel Isomer level (optional)
     */
    G4Isotope(const G4String& name,  // its name
              G4int z,  // atomic number
              G4int n,  // number of nucleons
              G4double a = 0.,  // mass of mole
              G4int mlevel = 0);  // isomer level

    /** Destructor */
    ~G4Isotope();

    /** Deleted copy constructor */
    G4Isotope(const G4Isotope&) = delete;
    /** Deleted assignment operator */
    G4Isotope& operator=(const G4Isotope&) = delete;

    /**
     * Get the isotope name
     * @return Reference to the isotope name
     */
    const G4String& GetName() const { return fName; }

    /**
     * Get the atomic number (Z)
     * @return Atomic number
     */
    G4int GetZ() const { return fZ; }

    /**
     * Get the number of nucleons (N)
     * @return Number of nucleons
     */
    G4int GetN() const { return fN; }

    /**
     * Get the atomic mass of mole in Geant4 units with electron shell
     * @return Atomic mass
     */
    G4double GetA() const { return fA; }

    /**
     * Get the isomer level
     * @return Isomer level
     */
    G4int Getm() const { return fm; }

    /**
     * Get an isotope by name from the static table
     * @param name Isotope name
     * @param warning If true, print warning if isotope not found
     * @return Pointer to the isotope, or nullptr if not found
     */
    static G4Isotope* GetIsotope(const G4String& name, G4bool warning = false);

    /**
     * Get the static table of all defined isotopes
     * @return Pointer to the static isotope table
     */
    static const G4IsotopeTable* GetIsotopeTable();

    /**
     * Get the total number of defined isotopes
     * @return Number of isotopes in the static table
     */
    static std::size_t GetNumberOfIsotopes();

    /**
     * Get the index of this isotope in the static table
     * @return Index in the isotope table
     */
    std::size_t GetIndex() const { return fIndexInTable; }

    /**
     * Stream insertion operator for pointer to isotope
     * @param os Output stream
     * @param isotope Pointer to isotope to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4Isotope* isotope);

    /**
     * Stream insertion operator for isotope reference
     * @param os Output stream
     * @param isotope Reference to isotope to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4Isotope& isotope);

    /**
     * Stream insertion operator for isotope table
     * @param os Output stream
     * @param table Reference to isotope table to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4IsotopeTable& table);

    /**
     * Equality operator
     * @param other Isotope to compare with
     * @return True if isotopes are equal
     */
    G4bool operator==(const G4Isotope& other) const;
    /**
     * Inequality operator
     * @param other Isotope to compare with
     * @return True if isotopes are not equal
     */
    G4bool operator!=(const G4Isotope& other) const;

    /**
     * Set the isotope name
     * @param name New name for the isotope
     */
    void SetName(const G4String& name) { fName = name; }

  private:

    G4String fName;  // name of the Isotope
    G4int fZ;  // atomic number
    G4int fN;  // number of nucleons
    G4double fA;  // atomic mass of a mole
    G4int fm;  // isomer level

    /**
     * Get mutable reference to the static isotope table
     * @return Reference to the static isotope table
     */
    static G4IsotopeTable& GetIsotopeTableRef();

    std::size_t fIndexInTable;  // index in the Isotope table
};

#endif
