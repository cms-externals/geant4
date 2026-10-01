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

//---------------------------------------------------------------
//
// G4AtomicBond
//
// Class Description:

#ifndef G4ATOMICBOND_HH
#define G4ATOMICBOND_HH

#include "G4Element.hh"
#include "globals.hh"

/**
 * Represents an atomic bond between two atoms with specified type and properties.
 */
class G4AtomicBond
{
  public:

    /**
     * Enumeration of bond types.
     */
    enum theBondType
    {
      Ionic = 0,       /**< Ionic bond type. */
      Covalent = 1,    /**< Covalent bond type. */
      Metallic = 2,    /**< Metallic bond type. */
      NA = -1          /**< Not applicable/undefined bond type. */
    };

    /**
     * Constructor for G4AtomicBond.
     * Initializes a bond between two atoms of specified types with given indices.
     * @param aType The type of bond (Ionic, Covalent, or Metallic).
     * @param firstAtomKind Pointer to the first atom element.
     * @param firstAtomNumber Index of the first atom.
     * @param secondAtomKind Pointer to the second atom element.
     * @param secondAtomNumber Index of the second atom.
     */
    G4AtomicBond(theBondType aType, G4Element* firstAtomKind, G4int firstAtomNumber,
                 G4Element* secondAtomKind, G4int secondAtomNumber)
      : theFirstAtomKind(firstAtomKind),
        theFirstAtomNumber(firstAtomNumber),
        theSecondAtomKind(secondAtomKind),
        theSecondAtomNumber(secondAtomNumber),
        theType(aType)
    {}

    /**
     * Destructor.
     */
    virtual ~G4AtomicBond() = default;

    /**
     * Gets the element type of the first atom.
     * @return Pointer to the first atom element.
     */
    inline G4Element* GetFirstAtomKind() const { return theFirstAtomKind; };
    /**
     * Sets the element type of the first atom.
     * @param aElement Pointer to the element to set.
     */
    void SetFirstAtomKind(G4Element* aElement) { theFirstAtomKind = aElement; };

    /**
     * Gets the index number of the first atom.
     * @return The index of the first atom.
     */
    inline G4int GetFirstAtomNumber() const { return theFirstAtomNumber; };
    /**
     * Sets the index number of the first atom.
     * @param aInt The index to set.
     */
    void SetFirstAtomNumber(G4int aInt) { theFirstAtomNumber = aInt; };

    /**
     * Gets the element type of the second atom.
     * @return Pointer to the second atom element.
     */
    inline G4Element* GetSecondAtomKind() const { return theSecondAtomKind; };
    /**
     * Sets the element type of the second atom.
     * @param aElement Pointer to the element to set.
     */
    void SetSecondAtomKind(G4Element* aElement) { theSecondAtomKind = aElement; };

    /**
     * Gets the index number of the second atom.
     * @return The index of the second atom.
     */
    inline G4int GetSecondAtomNumber() const { return theSecondAtomNumber; };
    /**
     * Sets the index number of the second atom.
     * @param aInt The index to set.
     */
    void SetSecondAtomNumber(G4int aInt) { theSecondAtomNumber = aInt; };

    /**
     * Gets the type of the bond.
     * @return The bond type.
     */
    inline theBondType GetType() const { return theType; };
    /**
     * Sets the type of the bond.
     * @param aType The bond type to set.
     */
    void SetType(theBondType aType) { theType = aType; };

    /**
     * Gets the aromaticity status of the bond.
     * @return True if the bond is aromatic, false otherwise.
     */
    inline G4bool GetAromaticity() const { return theAromaticity; };
    /**
     * Sets the aromaticity status of the bond.
     * @param aBool True to mark as aromatic, false otherwise.
     */
    void SetAromaticity(G4bool aBool) { theAromaticity = aBool; };

  private:

    // The indexes of the two atoms connected by the bond are stored
    // theFirstAtomKind is the index in the vector of G4Elements in the array
    // and theFirstAtomNumber is the index in the vector of G4ThreeVector in
    // the G4CrystalLattice. For the second atom theSecondAtomKind and
    // theSecondAtomNumber are used.
    G4Element* theFirstAtomKind;
    G4int theFirstAtomNumber;

    G4Element* theSecondAtomKind;
    G4int theSecondAtomNumber;

    // theType stores which kind of bond is the G4AtomicBond object
    // there are three types: Ionic, Covalent and Metallic
    // and they are enumerate in the theBondType enum.
    theBondType theType;

    // theAromaticity stores if the bond is aromatic.
    G4bool theAromaticity{false};
};

#endif
