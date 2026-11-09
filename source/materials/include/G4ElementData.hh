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
// GEANT4 Class file
//
// Description: Data structure for cross sections, shell cross sections,
//              isotope cross sections. Data access via integer variable
//              Z (atomic number in majority of applications), which may
//              be in the interval 0 <= Z < length. For isotope like
//              data a second parameter idx or data ID code are used.
//              In most cases ID = A - atomic weight number.
//              There are run time const methods, in which input is not checked
//              assuming responsibility of consumer code. Another run time
//              check input and may throwgh a fatal exception
//
// Author:      V.Ivanchenko 10.03.2011
//
// Modifications: 30.09.2023 Extended functionality, data size defined in constructor
//
//----------------------------------------------------------------------------
//

#ifndef G4ELEMENTDATA_HH
#define G4ELEMENTDATA_HH

#include "G4Physics2DVector.hh"
#include "G4PhysicsVector.hh"
#include "globals.hh"

#include <vector>

class G4ElementDataRegistry;

/**
 * @class G4ElementData
 * Data structure for managing cross sections and isotope/shell data.
 *
 * Provides storage and access for cross sections, shell cross sections, and
 * isotope cross sections indexed by atomic number Z and optional component IDs.
 */
class G4ElementData
{
  public:

    /**
     * Constructor for G4ElementData.
     * @param length Maximum number of elements (default 99 for all periodic table elements).
     */
    explicit G4ElementData(G4int length = 99);

    /** Destructor. */
    ~G4ElementData();

    /** Assignment operator (deleted). */
    G4ElementData& operator=(const G4ElementData& right) = delete;
    /** Copy constructor (deleted). */
    G4ElementData(const G4ElementData&) = delete;

    /**
     * Reserve memory for 1D data storage.
     * @param size Number of elements to reserve.
     */
    void Reserve1D(std::size_t size);
    /**
     * Reserve memory for 2D data storage.
     * @param size Number of elements to reserve.
     */
    void Reserve2D(std::size_t size);

    /**
     * Initialize cross section data for an element.
     * @param Z Atomic number.
     * @param v Pointer to the physics vector containing cross section data.
     */
    void InitialiseForElement(G4int Z, G4PhysicsVector* v);

    /**
     * Initialize 2D cross section data for an element.
     * @param Z Atomic number.
     * @param v Pointer to the 2D physics vector containing cross section data.
     */
    void InitialiseForElement(G4int Z, G4Physics2DVector* v);

    /**
     * Reserve vector of components for an element.
     * @param Z Atomic number.
     * @param nComponents Number of components to reserve.
     */
    void InitialiseForComponent(G4int Z, G4int nComponents = 0);

    /**
     * Reserve vector of 2D components for an element.
     * @param Z Atomic number.
     * @param nComponents Number of 2D components to reserve.
     */
    void InitialiseFor2DComponent(G4int Z, G4int nComponents = 0);

    /**
     * Add a component to the component vector for an element.
     * @param Z Atomic number.
     * @param id Component ID, used to look it up
     * @param v Pointer to the physics vector containing component data.
     */
    void AddComponent(G4int Z, G4int id, G4PhysicsVector* v);

    /**
     * Add a 2D component to the component vector for an element.
     * @param Z Atomic number.
     * @param id Component ID, used to look it up
     * @param v Pointer to the 2D physics vector containing component data.
     */
    void Add2DComponent(G4int Z, G4int id, G4Physics2DVector* v);

    /**
     * Create a new 1D physics vector.
     * @param Z Atomic number.
     * @param ne Number of energy points.
     * @return Pointer to the newly created 1D physics vector.
     */
    G4PhysicsVector* New1DVector(G4int Z, G4int ne);

    /**
     * Create a new 2D physics vector.
     * @param Z Atomic number.
     * @param ny Number of Y points.
     * @param ne Number of energy points.
     * @return Pointer to the newly created 2D physics vector.
     */
    G4Physics2DVector* New2DVector(G4int Z, G4int ny, G4int ne);

    /**
     * Set the name of the dataset.
     * @param nam Name string for the dataset.
     */
    inline void SetName(const G4String& nam);

    //--------------------------------------------------------------
    // run time const methods - no check on validity of input
    // it is a responsibility of the consume code to check the input
    //--------------------------------------------------------------

    /**
     * Get the name of the dataset.
     * @return Const reference to the dataset name string.
     */
    inline const G4String& GetName() const;

    /**
     * Get the 1D physics vector for an element.
     * @param Z Atomic number.
     * @return Pointer to the physics vector (no bounds check).
     */
    inline G4PhysicsVector* GetElementData(G4int Z) const;

    /**
     * Get the 2D physics vector for an element.
     * @param Z Atomic number.
     * @return Pointer to the 2D physics vector (no bounds check).
     */
    inline G4Physics2DVector* GetElement2DData(G4int Z) const;

    /**
     * Get a component vector by component ID.
     * @param Z Atomic number.
     * @param id Component ID to search for.
     * @return Pointer to the physics vector, or nullptr if not found (no bounds check).
     */
    inline G4PhysicsVector* GetComponentDataByID(G4int Z, G4int id) const;

    /**
     * Get a 2D component vector by component ID.
     * @param Z Atomic number.
     * @param id Component ID to search for.
     * @return Pointer to the 2D physics vector, or nullptr if not found (no bounds check).
     */
    inline G4Physics2DVector* Get2DComponentDataByID(G4int Z, G4int id) const;

    /**
     * Get the cross section value for an element at a given kinetic energy.
     * @param Z Atomic number.
     * @param kinEnergy Kinetic energy value.
     * @return Cross section value (no bounds check).
     */
    inline G4double GetValueForElement(G4int Z, G4double kinEnergy) const;

    //--------------------------------------------------------------
    // run time const methods with input parameters control
    //--------------------------------------------------------------

    /**
     * Get the number of components for an element.
     * @param Z Atomic number.
     * @return Number of components, or 0 if element has no components.
     */
    inline std::size_t GetNumberOfComponents(G4int Z) const;

    /**
     * Get the number of 2D components for an element.
     * @param Z Atomic number.
     * @return Number of 2D components, or 0 if element has no 2D components.
     */
    inline std::size_t GetNumberOf2DComponents(G4int Z) const;

    /**
     * Get the component ID at a given index.
     * @param Z Atomic number.
     * @param idx Index of the component.
     * @return Component ID (e.g., nucleon number or shell number), or 0 if index is out of bounds.
     */
    inline G4int GetComponentID(G4int Z, std::size_t idx) const;

    /**
     * Get a component vector by index.
     * @param Z Atomic number.
     * @param idx Component index.
     * @return Pointer to the physics vector, or nullptr if index is out of bounds.
     */
    inline G4PhysicsVector* GetComponentDataByIndex(G4int Z, std::size_t idx) const;

    /**
     * Get a 2D component vector by index.
     * @param Z Atomic number.
     * @param idx Component index.
     * @return Pointer to the 2D physics vector, or nullptr if index is out of bounds.
     */
    inline G4Physics2DVector* Get2DComponentDataByIndex(G4int Z, std::size_t idx) const;

    /**
     * Get the cross section value for a component at a given kinetic energy.
     * @param Z Atomic number.
     * @param idx Component index.
     * @param kinEnergy Kinetic energy value.
     * @return Cross section value, or 0.0 if index is out of bounds.
     */
    inline G4double GetValueForComponent(G4int Z, std::size_t idx, G4double kinEnergy) const;

  private:

    /**
     * Internal error handling method.
     * @param Z Atomic number.
     * @param text Error message text.
     */
    void DataError(G4int Z, const G4String&);

    /** Maximum number of elements in the periodic table. */
    const G4int maxNumElm;

    /** Pointer to the element data registry. */
    G4ElementDataRegistry* fRegistry{nullptr};

    /** 1D element cross section data indexed by atomic number. */
    std::vector<G4PhysicsVector*> elmData;
    /** Component-wise 1D cross section data indexed by atomic number. */
    std::vector<std::vector<std::pair<G4int, G4PhysicsVector*>>*> compData;

    /** 2D element cross section data indexed by atomic number. */
    std::vector<G4Physics2DVector*> elm2Data;
    /** Component-wise 2D cross section data indexed by atomic number. */
    std::vector<std::vector<std::pair<G4int, G4Physics2DVector*>>*> comp2D;

    /** Name of the dataset. */
    G4String name{""};
};

//--------------------------------------------------------------
// run time const methods without check on validity of input
//--------------------------------------------------------------

inline void G4ElementData::SetName(const G4String& nam)
{
  name = nam;
}

inline const G4String& G4ElementData::GetName() const
{
  return name;
}

inline G4PhysicsVector* G4ElementData::GetElementData(G4int Z) const
{
  return elmData[Z];
}

inline G4Physics2DVector* G4ElementData::GetElement2DData(G4int Z) const
{
  return elm2Data[Z];
}

inline G4PhysicsVector* G4ElementData::GetComponentDataByID(G4int Z, G4int id) const
{
  G4PhysicsVector* v = nullptr;
  for (auto const& p : *(compData[Z]))
  {
    if (id == p.first)
    {
      v = p.second;
      break;
    }
  }
  return v;
}

inline G4Physics2DVector* G4ElementData::Get2DComponentDataByID(G4int Z, G4int id) const
{
  G4Physics2DVector* v = nullptr;
  for (auto const& p : *(comp2D[Z]))
  {
    if (id == p.first)
    {
      v = p.second;
      break;
    }
  }
  return v;
}

inline G4double G4ElementData::GetValueForElement(G4int Z, G4double kinEnergy) const
{
  return elmData[Z]->Value(kinEnergy);
}

//--------------------------------------------------------------
// run time const methods with check on validity of input
//--------------------------------------------------------------

inline std::size_t G4ElementData::GetNumberOfComponents(G4int Z) const
{
  return (nullptr != compData[Z]) ? compData[Z]->size() : 0;
}

inline std::size_t G4ElementData::GetNumberOf2DComponents(G4int Z) const
{
  return (nullptr != comp2D[Z]) ? comp2D[Z]->size() : 0;
}

inline G4int G4ElementData::GetComponentID(G4int Z, std::size_t idx) const
{
  return (idx < GetNumberOfComponents(Z)) ? (*(compData[Z]))[idx].first : 0;
}

inline G4PhysicsVector* G4ElementData::GetComponentDataByIndex(G4int Z, std::size_t idx) const
{
  return (idx < GetNumberOfComponents(Z)) ? (*(compData[Z]))[idx].second : nullptr;
}

inline G4Physics2DVector* G4ElementData::Get2DComponentDataByIndex(G4int Z, std::size_t idx) const
{
  return (idx < GetNumberOf2DComponents(Z)) ? (*(comp2D[Z]))[idx].second : nullptr;
}

inline G4double G4ElementData::GetValueForComponent(G4int Z, std::size_t idx, G4double e) const
{
  return (idx < GetNumberOfComponents(Z)) ? (*(compData[Z]))[idx].second->Value(e) : 0.0;
}

#endif
