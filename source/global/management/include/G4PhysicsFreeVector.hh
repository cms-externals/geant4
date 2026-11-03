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
// G4PhysicsFreeVector
//
// Class description:
//
// A physics vector which has values of energy-loss, cross-section,
// and other physics values of a particle in matter in a given
// range of the energy, momentum, etc. The scale of energy/momentum
// bins is in free, i.e. it is NOT need to be linear or log. Only
// restriction is that bin values alway have to increase from
// a lower bin to a higher bin. This is necessary for the binary
// search to work correctly.

// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 06 Jun. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige: Use STL vector for dataVector and binVector
// - 04 Feb. 2021, V.Ivanchenko moved implementation of all free vectors
//                 to this class
// --------------------------------------------------------------------
#ifndef G4PHYSICSFREEVECTOR_HH
#define G4PHYSICSFREEVECTOR_HH

#include "G4PhysicsVector.hh"
#include "globals.hh"

#include <vector>

/** Physics vector specialised for general bins
 * @ingroup global_management
 * Supports bin edges with any spacing, as long as these are increasing (consecutive edges may be
 * equal)
 * @check Interpolation will always select the lower-numbered bin in case consecutive energies are
 * set equal
 */

class G4PhysicsFreeVector : public G4PhysicsVector
{
  public:

    /** @brief Construct empty vector */
    explicit G4PhysicsFreeVector(G4bool spline = false);
    /**  Construct vector of given length initialised to 0 */
    explicit G4PhysicsFreeVector(G4int length);
    /**  Construct vector of given length initialised to 0 */
    explicit G4PhysicsFreeVector(std::size_t length, G4bool spline = false);

    /** @deprecated Obsolete constructor */
    explicit G4PhysicsFreeVector(std::size_t length, G4double emin, G4double emax,
                                 G4bool spline = false);

    /** Construct and fill a vector
     *
     * @param energies Vector of energy bin values
     * @param values Vector of corresponding values
     * @param spline Whether to use spline interpolation
     * @pre The energies and values vectors are the same length
     * @pre Energies vector is (non-strictly) increasing
     * @post Energy and value are filled OR an exception is raised if the lengths do not match
     */
    explicit G4PhysicsFreeVector(const std::vector<G4double>& energies,
                                 const std::vector<G4double>& values, G4bool spline = false);
    /** Construct and fill a vector
     *
     * @param energies Array of energy bin values
     * @param values Array of corresponding values
     * @param spline Whether to use spline interpolation
     * @pre The energies and values pointers are valid arrays of length 'length'
     * @pre Energies array is (non-strictly) increasing
     * @post Energy and value are filled
     */
    explicit G4PhysicsFreeVector(const G4double* energies, const G4double* values,
                                 std::size_t length, G4bool spline = false);

    ~G4PhysicsFreeVector() override = default;

    /** Put data for energy bin and value at a given index
     * @param index Index in the vector
     * @param energy Energy value
     * @param value Data value
     * @pre index is valid
     * @post Energy bin and value at the given index are set
     * @post If the bin is first or last, then edgeMin or edgeMax are updated accordingly
     */
    void PutValues(const std::size_t index, const G4double energy, const G4double value);

    /** Insert an extra pair of energy and value.
     *
     * Inserts energy and value at the appropriate location. If energy is the
     * same as an existing energy bin value, the new value is added after.
     * @param energy Energy value
     * @param value Data value
     */
    void InsertValues(const G4double energy, const G4double value);

    /** Enable logarithmic bin search.
     * @param n Number of bins (default 1)
     */
    void EnableLogBinSearch(const G4int n = 1);

    /** Obsolete. Fill the vector with a value at a given index and energy.
     * @deprecated
     * @param index Index in the vector
     * @param e Energy value
     * @param value Data value
     */
    void PutValue(const std::size_t index, const G4double e, const G4double value);
};

#endif
