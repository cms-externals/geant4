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
// ClassName:   G4CrystalExtension
//
// Description: Contains crystal properties
//
// Class description:
//
// Extension of G4Material for the management of a crystal
// structure. It has to be attached to a G4ExtendedMaterial
// in order to instantiate a G4LogicalCrystalVolume.
//
// 21-04-16, created by E.Bagli

#ifndef G4CRYSTALEXTENSION_HH
#define G4CRYSTALEXTENSION_HH

#include "G4AtomicBond.hh"
#include "G4CrystalAtomBase.hh"
#include "G4CrystalUnitCell.hh"
#include "G4NistManager.hh"
#include "G4VMaterialExtension.hh"

#include <vector>

/**
 * Extension of G4Material for crystal structure management.
 * Stores crystal lattice properties, elasticity tensors, atomic positions, and bonds.
 */
class G4CrystalExtension : public G4VMaterialExtension
{
  public:

    /**
     * Full 4-dimensional elasticity tensor.
     */
    using Elasticity = G4double[3][3][3][3];
    /**
     * Reduced 2-dimensional elasticity tensor.
     */
    using ReducedElasticity = G4double[6][6];

  public:  // with description

    /**
     * Constructor for crystal extension.
     * Attaches crystal properties to a G4Material.
     * @param aMaterial Pointer to the material to extend.
     * @param name Name of the crystal extension (default: "crystal").
     */
    G4CrystalExtension(G4Material*, const G4String& name = "crystal");

    /**
     * Destructor.
     */
    ~G4CrystalExtension() override = default;

    /**
     * Print the crystal extension properties.
     */
    void Print() const override { ; };

    /**
     * Gets the material associated with this crystal extension.
     * @return Pointer to the G4Material.
     */
    G4Material* GetMaterial() { return fMaterial; };
    /**
     * Sets the material for this crystal extension.
     * @param aMat Pointer to the material to set.
     */
    void SetMaterial(G4Material* aMat) { fMaterial = aMat; };

    /**
     * Sets the unit cell for this crystal.
     * @param aUC Pointer to the G4CrystalUnitCell.
     */
    inline void SetUnitCell(G4CrystalUnitCell* aUC) { theUnitCell = aUC; }
    /**
     * Gets the unit cell for this crystal.
     * @return Pointer to the G4CrystalUnitCell.
     */
    inline G4CrystalUnitCell* GetUnitCell() const { return theUnitCell; }

    /**
     * Gets the full elasticity tensor.
     * @return Reference to the 4-dimensional elasticity tensor.
     */
    const Elasticity& GetElasticity() const { return fElasticity; }
    /**
     * Gets the reduced elasticity tensor.
     * @return Reference to the 2-dimensional reduced elasticity tensor.
     */
    const ReducedElasticity& GetElReduced() const { return fElReduced; }

    /**
     * Gets a specific elasticity tensor component Cijkl.
     * @param i First index.
     * @param j Second index.
     * @param k Third index.
     * @param l Fourth index.
     * @return The elasticity value at position [i][j][k][l].
     */
    G4double GetCijkl(G4int i, G4int j, G4int k, G4int l) const { return fElasticity[i][j][k][l]; }

    // Reduced elasticity tensor: C11-C66 interface for clarity
    /**
     * Sets the reduced elasticity tensor.
     * Converts the reduced 6x6 tensor to the full 4D elasticity tensor.
     * @param mat The 6x6 reduced elasticity tensor to set.
     */
    void SetElReduced(const ReducedElasticity& mat);

    /**
     * Sets a reduced elasticity tensor component Cpq.
     * @param p First reduced index (1-6).
     * @param q Second reduced index (1-6).
     * @param value The elasticity value to set.
     */
    void SetCpq(G4int p, G4int q, G4double value);
    /**
     * Gets a reduced elasticity tensor component Cpq.
     * @param p First reduced index (1-6).
     * @param q Second reduced index (1-6).
     * @return The elasticity value at position [p-1][q-1].
     */
    G4double GetCpq(G4int p, G4int q) const { return fElReduced[p - 1][q - 1]; }

    /**
     * Gets the atomic base for a given element.
     * @param anElement Pointer to the element.
     * @return Pointer to the G4CrystalAtomBase for this element
     * @post If the AtomBase is defined, it is returned. If it is not, then a default Base is registered, a warning is raised and the new base is returned
     */
    G4CrystalAtomBase* GetAtomBase(const G4Element* anElement);
    /**
     * Adds an atomic base for a given element.
     * @param anElement Pointer to the element.
     * @param aBase Pointer to the G4CrystalAtomBase to add.
     */
    void AddAtomBase(const G4Element* anElement, G4CrystalAtomBase* aBase)
    {
      theCrystalAtomBaseMap.insert(
        std::pair<const G4Element*, G4CrystalAtomBase*>(anElement, aBase));
    }

    /**
     * Gets the atomic base for a given element index.
     * @param anElIdx The element index in the material.
     * @return Pointer to the G4CrystalAtomBase for this element.
     * @post If the AtomBase is defined, it is returned. If it is not, then a default Base is registered, a warning is raised and the new base is returned
     */
    G4CrystalAtomBase* GetAtomBase(G4int anElIdx)
    {
      return GetAtomBase(fMaterial->GetElement(anElIdx));
    }

    /**
     * Adds an atomic base for a given element index.
     * @param anElIdx The element index in the material.
     * @param aLattice Pointer to the G4CrystalAtomBase to add.
     */
    void AddAtomBase(G4int anElIdx, G4CrystalAtomBase* aLattice)
    {
      AddAtomBase(fMaterial->GetElement(anElIdx), aLattice);
    }

    /**
     * Gets the positions of atoms in the unit cell for a specific element.
     * @param anEl Pointer to the element.
     * @param vecout Output vector to fill with atom positions.
     * @return True if successful, false otherwise.
     */
    G4bool GetAtomPos(const G4Element* anEl, std::vector<G4ThreeVector>& vecout);
    /**
     * Gets the positions of all atoms in the unit cell for all elements.
     * @param vecout Output vector to fill with atom positions.
     * @return True if successful, false otherwise.
     */
    G4bool GetAtomPos(std::vector<G4ThreeVector>& vecout);

    /**
     * Gets the positions of atoms in the unit cell for a given element index.
     * @param anElIdx The element index in the material.
     * @param vecout Output vector to fill with atom positions.
     * @return True if successful.
     */
    G4bool GetAtomPos(G4int anElIdx, std::vector<G4ThreeVector>& vecout)
    {
      GetAtomPos(fMaterial->GetElement(anElIdx), vecout);
      return true;
    }

    /**
     * Computes the structure factor for a given scattering vector and Miller indices.
     * Uses equation 46 from Chapter 2 of Introduction to Solid State Physics by C. Kittel.
     * @param kScatteringVector The magnitude of the scattering vector.
     * @param h Miller index h.
     * @param k Miller index k.
     * @param l Miller index l.
     * @return The complex structure factor.
     */
    G4complex ComputeStructureFactor(G4double kScatteringVector, G4int h, G4int k, G4int l);
    /**
     * Computes the geometrical structure factor for Miller indices.
     * @param h Miller index h.
     * @param k Miller index k.
     * @param l Miller index l.
     * @return The complex geometrical structure factor.
     */
    G4complex ComputeStructureFactorGeometrical(G4int h, G4int k, G4int l);

    /**
     * Adds an atomic bond to the crystal structure.
     * @param aBond Pointer to the G4AtomicBond to add.
     */
    void AddAtomicBond(G4AtomicBond* aBond) { theAtomicBondVector.push_back(aBond); };
    /**
     * Gets an atomic bond at the specified index.
     * @param idx The index of the atomic bond.
     * @return Pointer to the G4AtomicBond at the given index.
     */
    G4AtomicBond* GetAtomicBond(G4int idx) { return theAtomicBondVector[idx]; };
    /**
     * Gets the vector of all atomic bonds.
     * @return Vector of G4AtomicBond pointers.
     */
    std::vector<G4AtomicBond*> GetAtomicBondVector() { return theAtomicBondVector; };

  protected:

    Elasticity fElasticity;  // Full 4D elasticity tensor
    ReducedElasticity fElReduced;  // Reduced 2D elasticity tensor

  private:

    G4Material* fMaterial;

    // Crystal cell description, i.e. space group
    // and cell parameters
    G4CrystalUnitCell* theUnitCell{nullptr};

    // Map of atom positions for each element
    // The coordinate system is the unit cell
    std::map<const G4Element*, G4CrystalAtomBase*> theCrystalAtomBaseMap;

    // Bond between atoms. Each bond is mapped with two Elements
    // and the number of the atoms in the corresponding G4CrystalBaseAtomPos
    std::vector<G4AtomicBond*> theAtomicBondVector;
};

#endif
