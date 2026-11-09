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

#ifndef G4NISTELEMENTBUILDER_HH
#define G4NISTELEMENTBUILDER_HH

//---------------------------------------------------------------------------
//
// ClassName:   G4NistElementBuilder
//
// Description: Utility class to hold and manipulate G4Elements defined from
//              Nist data base
//
// Author:      V.Ivanchenko 21.11.2004
//
// Modifications:
// 27.02.06 V.Ivanchenko Return m=0 if Z&N combination is out of NIST
// 27.02.06 V.Ivanchenko add GetAtomicMassAmu
// 17.10.06 V.Ivanchenko add GetAtomicMass and GetNistElementNames methods
// 02.05.07 V.Ivanchenko add GetNistFirstIsotopeN and GetNumberOfNistIsotopes
// 06.08.08 V.Ivanchenko add binding energy parameterisation and use isotope
//                       mass in G4 units
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//

#include "G4Element.hh"
#include "globals.hh"

#include <CLHEP/Units/PhysicalConstants.h>

#include <vector>

const G4int maxNumElements = 108;
const G4int maxAbundance = 3500;

/** Create and hold G4Elements defined from NIST data
 *
 *  Creates elements as requested using the tabulated data. These can then be fetched and used.
 */
class G4NistElementBuilder
{
  public:

    /**
     * Constructor for G4NistElementBuilder
     * @param vb Verbosity level
     */
    explicit G4NistElementBuilder(G4int vb);
    /** Destructor */
    ~G4NistElementBuilder() = default;

    /** Get the NIST element with Atomic Number Z
     *
     * @param Z Atomic Number
     * @return Pointer to the element, if it exists, else nullptr
    */
    inline G4Element* FindElement(G4int Z) const;
    /**
     * Find or build an element by atomic number
     * @param Z Atomic number
     * @param buildIsotopes If true, build with isotope composition
     * @return Pointer to the element
     */
    G4Element* FindOrBuildElement(G4int Z, G4bool buildIsotopes = true);

    /**
     * Find or build an element by symbol
     * @param symb Element symbol
     * @param buildIsotopes If true, build with isotope composition
     * @return Pointer to the element
     */
    G4Element* FindOrBuildElement(const G4String& symb, G4bool buildIsotopes = true);
    /**
     * Print element information
     * @param Z Atomic number
     */
    void PrintElement(G4int Z) const;

    /**
     * Get the vector of Geant4 predefined element names
     * @return Reference to vector of element symbols
     */
    const std::vector<G4String>& GetElementNames() const;

    /**
     * Get the atomic number from element symbol
     * @param symb Element symbol
     * @return Atomic number
     */
    G4int GetZ(const G4String& symb) const;

    /**
     * Get the atomic mass in atomic mass units by element symbol
     * @param symb Element symbol
     * @return Atomic mass in AMU
     */
    G4double GetAtomicMassAmu(const G4String& symb) const;

    /**
     * Get the atomic mass in atomic mass units by atomic number
     * Mean mass of an atom with electron shell for natural isotope composition
     * @param Z Atomic number
     * @return Atomic mass in AMU
     */
    inline G4double GetAtomicMassAmu(G4int Z) const;

    /**
     * Get the mass of an isotope without electron shell in Geant4 energy units
     * @param Z Atomic number
     * @param N Number of nucleons
     * @return Isotope mass
     */
    inline G4double GetIsotopeMass(G4int Z, G4int N) const;

    /**
     * Get the atomic mass in Geant4 energy units with electron shell
     * @param Z Atomic number
     * @param N Number of nucleons
     * @return Atomic mass
     */
    inline G4double GetAtomicMass(G4int Z, G4int N) const;

    /**
     * Get the total ionisation energy of an atom
     * @param Z Atomic number
     * @return Total electron binding energy
     */
    inline G4double GetTotalElectronBindingEnergy(G4int Z) const;

    /**
     * Get the natural isotope abundance
     * @param Z Atomic number
     * @param N Number of nucleons
     * @return Isotope abundance fraction
     */
    inline G4double GetIsotopeAbundance(G4int Z, G4int N) const;

    /**
     * Get the neutron number (N) for the first natural isotope
     * @param Z Atomic number
     * @return Number of neutrons
     */
    inline G4int GetNistFirstIsotopeN(G4int Z) const;

    /**
     * Get the number of natural isotopes for an element
     * @param Z Atomic number
     * @return Number of isotopes
     */
    inline G4int GetNumberOfNistIsotopes(G4int Z) const;

    /**
     * Get the maximum number of elements in the NIST database
     * @return Maximum atomic number supported
     */
    inline G4int GetMaxNumElements() const;

    /**
     * Set the verbosity level
     * @param val Verbosity level
     */
    inline void SetVerbose(G4int);

  private:

    void Initialise();

    // Add element parameters to internal G4 database:
    // Z - atomic number, N - number of nucleons, A - atomic mass (amu),
    // sigmaA - accuracy of mass in last digits, W - natural abundances (percent)
    void AddElement(const G4String& symbol, G4int Z, G4int NumberOfIsotopes, const G4int& N,
                    const G4double& A, const G4double& sigmaA, const G4double& W);

    // Build a G4Element from the G4 dataBase
    G4Element* BuildElement(G4int Z);

  private:

    G4double atomicMass[maxNumElements];  // amu
    G4double bindingEnergy[maxNumElements];
    G4int nIsotopes[maxNumElements];
    G4int nFirstIsotope[maxNumElements];
    G4int idxIsotopes[maxNumElements];

    G4int elmIndex[maxNumElements];

    G4double massIsotopes[maxAbundance];  // G4 units
    G4double sigMass[maxAbundance];  // G4 units
    G4double relAbundance[maxAbundance];

    G4int index;
    G4int verbose;

    std::vector<G4String> elmSymbol;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetAtomicMassAmu(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? atomicMass[Z] : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetIsotopeMass(G4int Z, G4int N) const
{
  G4double mass = 0.0;
  if (Z > 0 && Z < maxNumElements)
  {
    G4int i = N - nFirstIsotope[Z];
    if (i >= 0 && i < nIsotopes[Z])
    {
      mass = massIsotopes[i + idxIsotopes[Z]];
    }
  }
  return mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetAtomicMass(G4int Z, G4int N) const
{
  G4double mass = 0.0;
  if (Z > 0 && Z < maxNumElements)
  {
    G4int i = N - nFirstIsotope[Z];
    if (i >= 0 && i < nIsotopes[Z])
    {
      mass = massIsotopes[i + idxIsotopes[Z]] + Z * CLHEP::electron_mass_c2 - bindingEnergy[Z];
    }
  }
  return mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetTotalElectronBindingEnergy(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? bindingEnergy[Z] : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistElementBuilder::GetIsotopeAbundance(G4int Z, G4int N) const
{
  G4double x = 0.0;
  if (Z > 0 && Z < maxNumElements)
  {
    G4int i = N - nFirstIsotope[Z];
    if (i >= 0 && i < nIsotopes[Z])
    {
      x = relAbundance[i + idxIsotopes[Z]];
    }
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistElementBuilder::GetNistFirstIsotopeN(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? nFirstIsotope[Z] : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistElementBuilder::GetNumberOfNistIsotopes(G4int Z) const
{
  return (Z > 0 && Z < maxNumElements) ? nIsotopes[Z] : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Get the list of element symbols defined in NIST database
 * @return Reference to vector of element symbols
 */
inline const std::vector<G4String>& G4NistElementBuilder::GetElementNames() const
{
  return elmSymbol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Get the maximum number of elements supported in NIST database
 * @return Maximum element count (maximum Z + 1)
 */
inline G4int G4NistElementBuilder::GetMaxNumElements() const
{
  return maxNumElements - 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * Set the verbosity level
 * @param val Verbosity level
 */
inline void G4NistElementBuilder::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/** Get the NIST element with Atomic Number Z
*
* @param Z Atomic Number
* @return Pointer to the element, if it exists, else nullptr
*/
inline G4Element* G4NistElementBuilder::FindElement(G4int Z) const
{
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  return (Z > 0 && Z < maxNumElements && elmIndex[Z] >= 0) ? (*theElementTable)[elmIndex[Z]]
                                                           : nullptr;
}

#endif
