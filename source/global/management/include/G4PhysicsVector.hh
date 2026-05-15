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
// G4PhysicsVector
//
// Class description:
//
// A physics vector which has values of energy-loss, cross-section,
// and other physics values of a particle in matter in a given
// range of energy, momentum, etc.
// This class serves as the base class for a vector having various
// energy scale, for example like 'log', 'linear', 'free', etc.

// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 03 Mar. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige: Use STL vector for dataVector and binVector
// --------------------------------------------------------------------
#ifndef G4PHYSICSVECTOR_HH
#define G4PHYSICSVECTOR_HH

#include <fstream>
#include <iostream>
#include <vector>

#include "G4Log.hh"
#include "G4PhysicsVectorType.hh"
#include "G4ios.hh"
#include "globals.hh"

/** @brief Physics vector
* @ingroup global_management 
*
* A physics vector which has values of energy-loss, cross-section,
* and other physics values of a particle in matter in a given
* range of energy, momentum, etc.
* This class serves as the base class for a vector having various
* energy scale, for example like 'log', 'linear', 'free', etc.
* @author G.cosmo @author K.Amako @author H.Kurashige
*/
class G4PhysicsVector
{
public:
  /**
   * @brief Default constructor
   *
   * Vector should be filled via G4PhysicsVector::Retrieve() method. Free vector may be filled via G4PhysicsVector::InsertValue() instead
   * @param spline Use spline interpolation
   */
  explicit G4PhysicsVector(G4bool spline = false);

  /**
   * @brief Copy constructor
   * @param other Source vector
   */
  G4PhysicsVector(const G4PhysicsVector&) = default;

  /**
   * @brief Assignment operator
   * @param other Source vector
   * @return Reference to this vector
   */
  G4PhysicsVector& operator=(const G4PhysicsVector&) = default;

  /**
   * @brief Deleted move constructor
   */
  G4PhysicsVector(const G4PhysicsVector&&) = delete;
  /**
   * @brief Deleted move assignment operator
   */
  G4PhysicsVector& operator=(const G4PhysicsVector&&) = delete;
  /**
   * @brief Deleted equality operator
   */
  G4bool operator==(const G4PhysicsVector& right) const = delete;
  /**
   * @brief Deleted inequality operator
   */
  G4bool operator!=(const G4PhysicsVector& right) const = delete;

  /**
   * @brief Destructor
   */
  virtual ~G4PhysicsVector() = default;

  /**
   * @brief Get the cross-section/energy-loss value
   *
   * Get the cross-section/energy-loss value corresponding to the given energy.
   * An appropriate interpolation is used to calculate the value. The parameter lastidx is updated to contain the bin location where the value was located, which can be reused in a subsequent call - this may avoid having to recompute the bin
   * @param energy Energy value
   * @param[inout] lastidx Index (updated for next call)
   * @return Interpolated value
   */
  inline G4double Value(const G4double energy, std::size_t& lastidx) const;

  /**
   * @brief Get the cross-section/energy-loss value
   *
   * Get the cross-section/energy-loss value corresponding to the given energy.
   * An appropriate interpolation is used to calculate the value.
   * This method should be used if bin location cannot be kept in the user code
   * @param energy Energy value
   * @return Interpolated value
   */

  inline G4double Value(const G4double energy) const;

  /**
   * @brief Obsolete method to get value
   * @deprecated Kept for compatibility. isOutRange is not used anymore
   * @param energy Energy value
   * @param isOutRange Set if energy is outside the binned range
   * @return Value
   */
  inline G4double GetValue(const G4double energy, G4bool& isOutRange) const;

  /**
   * @brief Get value
   *
   * Same as the Value() method above but specialised for log-vector type.
   * Note, unlike the general Value() method above, this method will work properly only for G4PhysicsLogVector.
   * @param energy Energy value
   * @param theLogEnergy Log of the energy value
   * @pre theLogEnergy == log(energy)
   * @pre this is a G4PhysicsLogVector
   * @return Value
   */
  inline G4double LogVectorValue(const G4double energy,
                                 const G4double theLogEnergy) const;
  /**
   * @brief Get value
   *
   * Same as the Value() method above but specialised for free log-vector type.
   * Note, unlike the general Value() method above, this method will work properly only for G4PhysicsLogVector. 
   * @param energy Energy value
   * @param theLogEnergy Log of the energy value
   * @pre theLogEnergy == log(energy)
   * @pre this is a G4PhysicsLogVector
   * @return Value
   */
   inline G4double LogFreeVectorValue(const G4double energy,
                                     const G4double theLogEnergy) const;

  /**
   * @brief Locate bin
   *
   * Locates the bin for the given energy. lastidx is used as an initial guess for the location. If the energy is out of range, lastidx is set to the min or max edge and the return value is false
   * @param energy Energy value
   * @param lastidx Bin location
   * @return Whether given energy value can/should be interpolated
   */
  inline G4bool CheckIndex(const G4double energy, std::size_t& lastidx) const;

  /**
   * @brief Get value
   *
   * NOTE: does not check if index is in range
   * @param index Data index
   * @return Value
   */
  inline G4double operator[](const std::size_t index) const;
  /**
   * @brief Get value
   * NOTE: does not check if index is in range
   * @param index Data index
   * @return Value
   */
  inline G4double operator()(const std::size_t index) const;

  /**
   * @brief Put data into the vector at 'index' position
   *
   * If index is out of range, a G4 exception is raised
   * @param index Data index
   * @param value Value to set
   * @pre index is in range[0, size)
   * @pre the energy values have been set
   * @post If index is in range, the value at idx is set, else a G4 exception is raised
   */
  inline void PutValue(const std::size_t index, const G4double value);

  /**
   * @brief Returns the energy bin at 'index'
   *
   * No bounds checking is performed. Use this when compute cross-section, dEdx, or other value before filling the vector using PutValue().
   * @param index Energy index
   * @return Energy value
   */
  inline G4double Energy(const std::size_t index) const;
  /**
   * @brief Returns the low edge energy for the specified index
   * @param index Energy index
   * @return Low edge energy
   */
  inline G4double GetLowEdgeEnergy(const std::size_t index) const;

  /**
   * @brief Returns the energy of the first point of the vector
   * @return Minimum energy
   */
  inline G4double GetMinEnergy() const;
  /**
   * @brief Returns the energy of the last point of the vector
   * @return Maximum energy
   */
  inline G4double GetMaxEnergy() const;

  /**
   * @brief Returns the data of the first point of the vector
   * @return Minimum value
   */
  inline G4double GetMinValue() const;
  /**
   * @brief Returns the data of the last point of the vector
   * @return Maximum value
   */
  inline G4double GetMaxValue() const;

  /**
   * @brief Get the total length of the vector
   * @return Vector length
   */
  inline std::size_t GetVectorLength() const;

  /**
   * @brief Computes energy bin
   * @param logenergy Log of energy
   * @pre this is a G4PhysicsLogVector
   * @return Bin index
   */
  inline std::size_t ComputeLogVectorBin(const G4double logenergy) const;

  /**
   * @brief Get physics vector type tag
   * @return Physics vector type
   */
  inline G4PhysicsVectorType GetType() const;

  /**
   * @brief True if using spline interpolation
   * @return Spline flag
   */
  inline G4bool GetSpline() const;

  /**
   * @brief Define verbosity level
   * @param value Verbosity level
   */
  inline void SetVerboseLevel(G4int value);

  /**
   * @brief Find energy using linear interpolation for vector filled by cumulative probability function
   * @param rand Random value
   * @deprecated Possibly deprecated, used once
   * @check Consider deprecation
   * @return Interpolated energy
   */
  inline G4double FindLinearEnergy(const G4double rand) const;

  /**
   * @brief Find low edge index of a bin for given energy (obsolete)
   * @param energy Energy value
   * @param idx Starting index
   * @deprecated 
   * @return Bin index
   */
  std::size_t FindBin(const G4double energy, std::size_t idx) const;

  /**
   * @brief Scale all values of the vector by factorV, energies by vectorE
   * @param factorE Energy scale factor
   * @param factorV Value scale factor
   */
  void ScaleVector(const G4double factorE, const G4double factorV);

  /**
   * @brief Fill second derivatives for spline interpolation
   *
   * There are 3 types of second derivative computations:
   *  fSplineSimple -     2d derivative continues
   *  fSplineBase -       3d derivative continues (the default)
   *  fSplineFixedEdges - 3d derivatives continues, 1st and last derivatives are fixed
   * @param type Spline type
   * @param dir1 First point derivative
   * @param dir2 End point derivative
   * @pre The vector is fully filled
   * @post The second derivatives are filled
   */
  void FillSecondDerivatives(const G4SplineType type = G4SplineType::Base,
                             const G4double dir1 = 0.0,
                             const G4double dir2 = 0.0);

// This method may be applied only once.
// Force length of data using std::vector::resize() with the
// the default value 0; partial cross section vector is resized
// only if the number of partial x-sections is above zero. 
  void SetDataLength(G4int dlength);

  /**
   * @brief Get energy for a given value
   * @param value Value to search for
   * @pre energy and data for this vector are strictly monotonic increasing (e.g. a cumulative PDF)
   * @return Energy
   */
  G4double GetEnergy(const G4double value) const;

  /**
   * @brief Store persistent data to file stream
   * @param fOut Output file stream
   * @param ascii Store as ASCII if true
   * @return Success flag
   */
  G4bool Store(std::ofstream& fOut, G4bool ascii = false) const;
  /**
   * @brief Retrieve persistent data from file stream
   * @param fIn Input file stream
   * @param ascii Retrieve as ASCII if true
   * @return Success flag
   */
  G4bool Retrieve(std::ifstream& fIn, G4bool ascii = false);

  /**
   * @brief Print vector to stream
   * @param os Output stream
   * @param vec Vector to print
   * @return Output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const G4PhysicsVector& vec);

  /**
   * @brief Print vector values
   * @param unitE Energy unit
   * @param unitV Value unit
   */
  void DumpValues(G4double unitE = 1.0, G4double unitV = 1.0) const;

protected:

  // The default implements a free vector initialisation.
  virtual void Initialise();

  void PrintPutValueError(std::size_t index, G4double value, 
                          const G4String& text);

private:

  // Internal methods for computing of spline coeffitients
  void ComputeSecDerivative0();
  void ComputeSecDerivative1();
  void ComputeSecDerivative2(const G4double firstPointDerivative,
                             const G4double endPointDerivative);

  // Linear or spline interpolation.
  inline G4double Interpolation(const std::size_t idx,
                                const G4double energy) const;
  
  // Assuming (edgeMin <= energy <= edgeMax).
  inline std::size_t LogBin(const G4double energy, const G4double loge) const;
  inline std::size_t BinaryBin(const G4double energy) const;
  inline std::size_t GetBin(const G4double energy) const;

protected:

  G4double edgeMin = 0.0;  /**< Energy of first point*/
  G4double edgeMax = 0.0;  /**< Energy of the last point*/

  G4double invdBin = 0.0;  /**< 1/Bin width for linear and log vectors*/
  G4double logemin = 0.0;  /**< used only for log vector*/

  G4double iBin1 = 0.0;  /**< 1/Bin width for scale log vector */
  G4double lmin1 = 0.0;  /**< used for log search of free vector */

  G4int verboseLevel = 0;
  std::size_t idxmax = 0; /**< Maximum index value (@check for spline interpolation is this reduced? */
  std::size_t imax1 = 0;
  std::size_t numberOfNodes = 0;
  std::size_t nLogNodes = 0;

  G4PhysicsVectorType type = T_G4PhysicsFreeVector;
  // The type of PhysicsVector (enumerator)

  std::vector<G4double> binVector;      /**< energy bins */
  std::vector<G4double> dataVector;     /**< crossection/energyloss values */
  std::vector<G4double> secDerivative;  /**< second derivatives */
  std::vector<std::size_t> scale;       // log seach

private:

  G4bool useSpline = false;
};

#include "G4PhysicsVector.icc"

#endif
