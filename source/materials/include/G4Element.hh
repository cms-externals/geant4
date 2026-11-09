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

//---------------------------------------------------------------------------
//
// ClassName:   G4Element
//
// Description: Contains element properties
//
// Class description:
//
// An element is a chemical element either directly defined in terms of
// its characteristics: its name, symbol,
//                      Z (effective atomic number)
//                      N (effective number of nucleons)
//                      A (effective mass of a mole)
// or in terms of a collection of constituent isotopes with specified
// relative abundance (i.e. fraction of nb of atoms per volume).
//
// Quantities, with physical meaning or not, which are constant in a given
// element are computed and stored here as Derived data members.
//
// The class contains as a private static member the table of defined
// elements (an ordered vector of elements).
//
// Elements can be assembled singly or in mixtures into materials used
// in volume definitions via the G4Material class.
//
// It is strongly recommended do not delete G4Element instance in the
// user code. All G4Elements will be automatically deleted at the end
// of Geant4 session

// 09-07-96, new data members added by L.Urban
// 17-01-97, aesthetic rearrangement, M.Maire
// 20-01-97, Tsai formula for the rad length, M.Maire
// 21-01-97, remove mixture flag, M.Maire
// 24-01-97, new data member: fTaul
//           new method: ComputeIonisationPara, M.Maire
// 20-03-97, corrected initialization of pointers, M.Maire
// 27-06-97, new function GetIsotope(int), M.Maire
// 24-02-98, fWeightVector becomes fRelativeAbundanceVector
// 27-04-98, atomic shell stuff, V. Grichine
// 09-07-98, Ionisation parameters removed from the class, M.Maire
// 04-08-98, new method GetElement(elementName), M.Maire
// 16-11-98, Subshell -> Shell, mma
// 30-03-01, suppression of the warning message in GetElement
// 17-07-01, migration to STL, M. Verderi
// 13-09-01, stl migration. Suppression of the data member fIndexInTable
// 14-09-01, fCountUse: nb of materials which use this element
// 26-02-02, fIndexInTable renewed
// 01-04-05, new data member fIndexZ to count the number of elements with same Z
// 17-10-06: Add Get/Set fNaturalAbundance (V.Ivanchenko)
// 17.09.09, add fNbOfShellElectrons and methods (V. Grichine)

#ifndef G4ELEMENT_HH
#define G4ELEMENT_HH

#include "G4ElementTable.hh"
#include "G4ElementVector.hh"
#include "G4IonisParamElm.hh"
#include "G4Isotope.hh"
#include "G4IsotopeVector.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <vector>

/** A chemical element
 * 
 * An element is a chemical element either directly defined in terms of its characteristics: its name, symbol, Z (effective atomic number), N (effective number of nucleons), A (effective mass of a mole) or in terms of a collection of constituent isotopes with specified relative abundance (i.e. fraction of nb of atoms per volume).
 *
 * Quantities, with physical meaning or not, which are constant in a given element are computed and stored here as Derived data members.
 * 
 * Once an element is created it is added to an internal table of defined elements.
 * 
 * Elements can be assembled singly or in mixtures into materials used in volume definitions via the G4Material class.
 * 
 * IMPORTANT: It is strongly recommended NOT to delete G4Element instances in the user code, as this will invalidate their G4ElementTable entry. All G4Elements will be automatically deleted at the end of Geant4 session.
 */
class G4Element
{
  public:  // with description

    /**
     * Constructor to build an element directly without reference to isotopes
     * @param name Element name
     * @param symbol Element symbol
     * @param Zeff Effective atomic number
     * @param Aeff Effective mass of a mole
     * @pre Zeff >= 1.0
     * @post A warning is generated if Z is not close to an integer
     */
    G4Element(const G4String& name,  // its name
              const G4String& symbol,  // its symbol
              G4double Zeff,  // atomic number
              G4double Aeff);  // mass of mole

    /**
     * Constructor to build an element from constituent isotopes
     * @param name Element name
     * @param symbol Element symbol
     * @param nbIsotopes Number of isotopes to be added via AddIsotope
     * @pre nIsotopes is strictly > 0
     * @post The element is ready for isotopes to be added
     */
    G4Element(const G4String& name,  // its name
              const G4String& symbol,  // its symbol
              G4int nbIsotopes);  // nb of isotopes

    /** Destructor */
    virtual ~G4Element();

    /** Deleted copy constructor */
    G4Element(G4Element&) = delete;
    /** Deleted assignment operator */
    const G4Element& operator=(const G4Element&) = delete;

    /**
     * Add an isotope to the element
     * @param isotope Pointer to the isotope to add
     * @param RelativeAbundance Fraction of atoms per volume for this isotope
     * @pre The Element was constructed ready to add isotopes
     * @pre The given isotope is not null
     * @pre The given isotope has the same Z as the element
     * @check No check is done for negative abundances
     */
    void AddIsotope(G4Isotope* isotope,  // isotope
                    G4double RelativeAbundance);  // fraction of nb of
                                                  // atomes per volume

    /**
     * Get the name of the element
     * @return Reference to the element name
     */
    inline const G4String& GetName() const { return fName; }
    /**
     * Get the symbol of the element
     * @return Reference to the element symbol
     */
    inline const G4String& GetSymbol() const { return fSymbol; }

    /**
     * Get the effective atomic number (Z)
     * @return Effective atomic number
     */
    inline G4double GetZ() const { return fZeff; }
    /**
     * Get the atomic number as integer
     * @return Atomic number (integer)
     */
    inline G4int GetZasInt() const { return fZ; }

    /**
     * Get the effective number of nucleons (N)
     * @return Effective number of nucleons in atomic units
     */
    inline G4double GetN() const { return fNeff; }
    /**
     * Get the atomic mass in atomic mass units
     * @return Atomic mass (same as GetN()) in AMU
     */
    inline G4double GetAtomicMassAmu() const { return fNeff; }

    /**
     * Get the effective mass of a mole in Geant4 units
     * @return Effective mass of a mole for atoms with atomic shell
     */
    inline G4double GetA() const { return fAeff; }

    inline G4bool GetNaturalAbundanceFlag() const;

    inline void SetNaturalAbundanceFlag(G4bool);

    /**
     * Get the number of atomic shells in this element
     * @return Number of atomic shells
     */
    inline G4int GetNbOfAtomicShells() const { return fNbOfAtomicShells; }

    /**
     * Get the binding energy of an atomic shell
     * @param index Shell index (ground shell is index 0)
     * @return Binding energy of the specified shell
     */
    G4double GetAtomicShell(G4int index) const;

    /**
     * Get the number of electrons in an atomic shell
     * @param index Shell index (ground shell is index 0)
     * @return Number of electrons in the specified shell
     */
    G4int GetNbOfShellElectrons(G4int index) const;

    /**
     * Get the number of isotopes constituting this element
     * @return Number of isotopes
     */
    inline std::size_t GetNumberOfIsotopes() const { return fNumberOfIsotopes; }

    /**
     * Get the vector of isotopes constituting this element
     * @return Pointer to the isotope vector
     */
    inline G4IsotopeVector* GetIsotopeVector() const { return theIsotopeVector; }

    /**
     * Get the vector of relative abundances for each isotope
     * @return Pointer to the relative abundance vector
     */
    inline G4double* GetRelativeAbundanceVector() const { return fRelativeAbundanceVector; }

    /**
     * Get a specific isotope from the element
     * @param iso Isotope index
     * @return Pointer to the isotope at the specified index
     */
    inline const G4Isotope* GetIsotope(G4int iso) const { return (*theIsotopeVector)[iso]; }

    /**
     * Get the static table of all defined elements
     * @return Pointer to the static element table
     */
    static const G4ElementTable* GetElementTable();

    /**
     * Get the total number of defined elements
     * @return Number of elements in the static table
     */
    static std::size_t GetNumberOfElements();

    /**
     * Get the index of this element in the static table
     * @return Index in the element table
     */
    inline std::size_t GetIndex() const { return fIndexInTable; }

    /**
     * Get a pointer to an element by name from the static table
     * @param name Name of the element to retrieve
     * @param warning If true, print warning if element not found
     * @return Pointer to the element, or nullptr if not found
     */
    static G4Element* GetElement(const G4String& name, G4bool warning = true);

    /**
     * Get the Coulomb correction factor
     * @return Coulomb correction factor
     */
    inline G4double GetfCoulomb() const { return fCoulomb; }

    /**
     * Get the radiation length factor from Tsai formula
     * @return Radiation length parameter
     */
    inline G4double GetfRadTsai() const { return fRadTsai; }

    /**
     * Get the ionisation parameters for this element
     * @return Pointer to ionisation parameters
     */
    inline G4IonisParamElm* GetIonisation() const { return fIonisation; }

    /**
     * Stream insertion operator for pointer to element
     * @param os Output stream
     * @param element Pointer to element to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4Element* element);
    /**
     * Stream insertion operator for element reference
     * @param os Output stream
     * @param element Reference to element to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4Element& element);
    /**
     * Stream insertion operator for element table
     * @param os Output stream
     * @param table Reference to element table to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4ElementTable& table);
    /**
     * Stream insertion operator for element vector
     * @param os Output stream
     * @param vector Reference to element vector to print
     * @return Reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const G4ElementVector& vector);

    /**
     * Set the name of the element
     * @param name New name for the element
     */
    inline void SetName(const G4String& name) { fName = name; }

    /** Deleted equality operator */
    G4bool operator==(const G4Element&) const = delete;
    /** Deleted inequality operator */
    G4bool operator!=(const G4Element&) const = delete;

  private:

    void InitializePointers();
    void ComputeDerivedQuantities();
    void ComputeCoulombFactor();
    void ComputeLradTsaiFactor();
    
    /** Populate isotopes from NIST data
     *
     * Adds the naturally occuring set of isotopes from Nist data, and their relative abundances
     * @pre Atomic number fZeff must have been set
     * @post fSymbol is set iff it was not already
     * @post Isotopic abundances are set for the element
    */
    void AddNaturalIsotopes();

    // Mutable access to the element table.
    static G4ElementTable& GetElementTableRef();

    // Basic data members (which define an Element)

    G4String fName;  // name
    G4String fSymbol;  // symbol
    G4double fZeff;  // Effective atomic number
    G4double fNeff;  // Effective number of nucleons
    G4double fAeff;  // Effective mass of a mole
    G4int fZ;

    G4int fNbOfAtomicShells;  // number  of atomic shells
    G4double* fAtomicShells;  // Pointer to atomic shell binding energies
    G4int* fNbOfShellElectrons;  // Pointer to the number of subshell electrons

    G4int fNumberOfIsotopes;  // Number of isotopes added to the element
    G4IsotopeVector* theIsotopeVector;  // vector of constituent isotopes of the element
    G4double* fRelativeAbundanceVector;  // Fraction nb of atomes per volume
                                         // for each constituent

    // Set up the static Table of Elements
    std::size_t fIndexInTable;
    G4bool fNaturalAbundance;

    // Derived data members (computed from the basic data members)

    G4double fCoulomb;  // Coulomb correction factor
    G4double fRadTsai;  // Tsai formula for the radiation length
    G4IonisParamElm* fIonisation;  // Pointer to ionisation parameters
};

/**
 * Get the natural abundance flag (inline definition)
 * @return True if element uses natural isotope abundances
 */
inline G4bool G4Element::GetNaturalAbundanceFlag() const
{
  return fNaturalAbundance;
}

/** Set the natural abundance flag
 * @param val Flag to enable/disable natural isotope abundances
 */
inline void G4Element::SetNaturalAbundanceFlag(G4bool val)
{
  fNaturalAbundance = val;
}

#endif
