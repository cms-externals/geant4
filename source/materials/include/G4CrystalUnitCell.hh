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
// G4CrystalUnitCell
//
// Class Description:
//

#ifndef G4CRYSTALUNITCELL_HH
#define G4CRYSTALUNITCELL_HH

#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4CrystalBravaisLattices.h"
#include "G4CrystalLatticeSystems.h"

#include <vector>

/**
 * Represents a crystalline unit cell with lattice parameters.
 * Stores cell dimensions, angles, basis vectors, and provides methods for crystallographic calculations.
 */
class G4CrystalUnitCell
{
  public:

    /**
     * Constructor for a crystalline unit cell.
     * Initializes the cell with specified lattice parameters and space group.
     * @param sizeA Length of cell edge A in angstroms.
     * @param sizeB Length of cell edge B in angstroms.
     * @param sizeC Length of cell edge C in angstroms.
     * @param alpha Angle alpha in degrees.
     * @param beta Angle beta in degrees.
     * @param gamma Angle gamma in degrees.
     * @param spacegroup The space group number.
     * @check Is this actually Angstroms?
     */
    G4CrystalUnitCell(G4double sizeA, G4double sizeB, G4double sizeC, G4double alpha, G4double beta,
                      G4double gamma, G4int spacegroup);

    /**
     * Destructor.
     */
    virtual ~G4CrystalUnitCell() = default;

    /**
     * Gets the space group number.
     * @return The space group number.
     */
    inline G4int GetSpaceGroup() const { return theSpaceGroup; };
    /**
     * Sets the space group number.
     * @param aInt The space group number to set.
     */
    inline void SetSpaceGroup(G4int aInt) { theSpaceGroup = aInt; };

    /**
     * Gets the crystal lattice system type.
     * @return The lattice system type.
     */
    theLatticeSystemType GetLatticeSystem() { return GetLatticeSystem(theSpaceGroup); }
    /**
     * Gets the Bravais lattice type.
     * @return The Bravais lattice type.
     */
    theBravaisLatticeType GetBravaisLattice() { return GetBravaisLattice(theSpaceGroup); }

    /**
     * Gets the basis vector at the specified index.
     * @param idx The index of the basis vector (0, 1, or 2).
     * @return The basis vector.
     */
    const G4ThreeVector& GetBasis(G4int idx) const;
    /**
     * Gets the unit basis vector at the specified index.
     * @param idx The index of the unit basis vector (0, 1, or 2).
     * @return The unit basis vector.
     */
    const G4ThreeVector& GetUnitBasis(G4int idx) const;
    /**
     * Gets the cell dimensions.
     * @return The cell size as a G4ThreeVector (a, b, c).
     */
    inline G4ThreeVector GetSize() const { return theSize; }
    /**
     * Gets the cell angles.
     * @return The cell angles as a G4ThreeVector (alpha, beta, gamma).
     */
    inline G4ThreeVector GetAngle() const { return theAngle; }

    // return theUnitBase[2] vector
    /**
     * Gets the unit basis trigonal vector.
     * Returns the third unit basis vector (theUnitBase[2]).
     * @return The third unit basis vector.
     */
    G4ThreeVector GetUnitBasisTrigonal();

    /**
     * Gets the reciprocal basis vector at the specified index.
     * @param idx The index of the reciprocal basis vector (0, 1, or 2).
     * @return The reciprocal basis vector.
     */
    const G4ThreeVector& GetRecBasis(G4int idx) const;
    /**
     * Gets the reciprocal unit basis vector at the specified index.
     * @param idx The index of the reciprocal unit basis vector (0, 1, or 2).
     * @return The reciprocal unit basis vector.
     */
    const G4ThreeVector& GetRecUnitBasis(G4int idx) const;
    /**
     * Gets the reciprocal cell dimensions.
     * @return The reciprocal cell size.
     */
    inline G4ThreeVector GetRecSize() const { return theRecSize; }
    /**
     * Gets the reciprocal cell angles.
     * @return The reciprocal cell angles.
     */
    inline G4ThreeVector GetRecAngle() const { return theRecAngle; }

    /**
     * Fills atomic positions in unit cell coordinates from a basis vector.
     * Populates the output vector with atom positions based on the unit basis.
     * @param pos The input position vector.
     * @param vecout The output vector to fill with computed positions.
     * @return True if successful, false otherwise.
     */
    G4bool FillAtomicUnitPos(G4ThreeVector& pos, std::vector<G4ThreeVector>& vecout);
    /**
     * Fills atomic positions in cell coordinates from a basis vector.
     * Populates the output vector with atom positions based on the direct basis.
     * @param pos The input position vector.
     * @param vecout The output vector to fill with computed positions.
     * @return True if successful, false otherwise.
     */
    G4bool FillAtomicPos(G4ThreeVector& pos, std::vector<G4ThreeVector>& vecout);

    /**
     * Fills the reduced elasticity tensor for this crystal.
     * Populates a 6x6 elasticity tensor based on the crystal system.
     * @param Cij The 6x6 array to fill with elasticity coefficients.
     * @return True if successful, false otherwise.
     */
    G4bool FillElReduced(G4double Cij[6][6]);

    /**
     * Computes and stores the unit cell volume.
     * @return The computed cell volume.
     */
    G4double ComputeCellVolume();

    /**
     * Gets the stored unit cell volume.
     * @return The cell volume.
     */
    inline G4double GetVolume() const { return theVolume; }
    /**
     * Gets the stored reciprocal cell volume.
     * @return The reciprocal cell volume.
     */
    inline G4double GetRecVolume() const { return theRecVolume; }

    /**
     * Computes the squared direct interplanar spacing.
     * @param h The Miller index h.
     * @param k The Miller index k.
     * @param l The Miller index l.
     * @return The squared interplanar spacing.
     */
    G4double GetIntSp2(G4int h, G4int k,
                       G4int l);
    /**
     * Computes the squared reciprocal interplanar spacing.
     * @param h The Miller index h.
     * @param k The Miller index k.
     * @param l The Miller index l.
     * @return The squared reciprocal interplanar spacing.
     */
    G4double GetRecIntSp2(G4int h, G4int k,
                          G4int l);
    /**
     * Computes the cosine of the angle between two planes.
     * @param h1 The first plane Miller index h.
     * @param k1 The first plane Miller index k.
     * @param l1 The first plane Miller index l.
     * @param h2 The second plane Miller index h.
     * @param k2 The second plane Miller index k.
     * @param l2 The second plane Miller index l.
     * @return The cosine of the angle between the two planes.
     */
    G4double GetIntCosAng(G4int h1, G4int k1, G4int l1, G4int h2, G4int k2,
                          G4int l2);

  protected:

    G4ThreeVector nullVec;  // Size and angles of the crystalline unit cell

    G4ThreeVector theSize;  // cell sizes
    G4ThreeVector theAngle;  // cell angles
    G4ThreeVector theUnitBasis[3];  // Basis unit vectors in direct orientation
    G4ThreeVector theBasis[3];  // Basis vectors in direct orientation
    //
    // Reciprocal size and angles of the crystalline unit cell
    //
    G4ThreeVector theRecSize;  // reciprocal cell sizes
    G4ThreeVector theRecAngle;  // reciprocal cell angles
    G4ThreeVector theRecUnitBasis[3];  // Basis unit vectors in reciprocal orientation
    G4ThreeVector theRecBasis[3];  // Basis vectors in reciprocal orientation

  private:

    theLatticeSystemType GetLatticeSystem(G4int aGroup);
    theBravaisLatticeType GetBravaisLattice(G4int aGroup);
    G4bool FillAmorphous(G4double Cij[6][6]) const;
    G4bool FillCubic(G4double Cij[6][6]) const;
    G4bool FillTetragonal(G4double Cij[6][6]) const;
    G4bool FillOrthorhombic(G4double Cij[6][6]) const;
    G4bool FillRhombohedral(G4double Cij[6][6]) const;
    G4bool FillMonoclinic(G4double Cij[6][6]) const;
    G4bool FillTriclinic(G4double Cij[6][6]) const;
    G4bool FillHexagonal(G4double Cij[6][6]) const;

    G4bool ReflectElReduced(G4double Cij[6][6]) const;

  private:

    G4int theSpaceGroup;  //

    G4double cosa, cosb, cosg;
    G4double sina, sinb, sing;
    G4double cosar, cosbr, cosgr;

    G4double theVolume;  // the cell volume
    G4double theRecVolume;  // the cell volume
};

#endif
