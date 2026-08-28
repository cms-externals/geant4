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
// G4Physics2DVector
//
// Class description:
//
// A 2-dimentional vector with linear interpolation.

// Author: Vladimir Ivanchenko, 25.09.2011
// --------------------------------------------------------------------
#ifndef G4PHYSICS2DVECTOR_HH
#define G4PHYSICS2DVECTOR_HH

#include "G4PhysicsVectorType.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>
#include <iostream>
#include <vector>

using G4PV2DDataVector = std::vector<G4double>;

/**
 * @brief A 2-dimentional vector with linear interpolation.
 *
 * @author Vladimir Ivanchenko
 * @date 25.09.2011
 * @ingroup global_management
 */
class G4Physics2DVector
{
  public:

    /**
     * @brief Default constructor.
     *
     * @see G4Physics2DVector::Retrieve
     */
    G4Physics2DVector();
    // Vector will be filled via Retrieve method

    /**
     * @brief Constructor with dimensions.
     * @param nx Number of X bins
     * @param ny Number of Y bins
     */
    explicit G4Physics2DVector(std::size_t nx, std::size_t ny);
    // Vector will be filled via Put methods

    /**
     * @brief Copy constructor.
     * @param other Source vector
     */
    G4Physics2DVector(const G4Physics2DVector&);

    /**
     * @brief Assignment operator.
     * @param other Source vector
     * @return Reference to this vector
     */
    G4Physics2DVector& operator=(const G4Physics2DVector&);

    G4bool operator==(const G4Physics2DVector& right) const = delete;

    G4bool operator!=(const G4Physics2DVector& right) const = delete;

    /**
     * @brief Destructor.
     */
    ~G4Physics2DVector();

    /**
     * @brief Interpolate 2D vector.
     * @param x X value
     * @param y Y value
     * @param lastidx Last X index
     * @param lastidy Last Y index
     * @return Interpolated value
     */
    G4double Value(G4double x, G4double y, std::size_t& lastidx, std::size_t& lastidy) const;
    /**
     * @brief Interpolate 2D vector.
     * @param x X value
     * @param y Y value
     * @return Interpolated value
     */
    G4double Value(G4double x, G4double y) const;
    // Main method to interpolate 2D vector
    // Consumer class should provide initial values of lastidx and lastidy

    /**
     * @brief Set X value at index.
     * @param idx Index (0-based)
     * @param value Value to set
     */
    inline void PutX(std::size_t idx, G4double value);
    /**
     * @brief Set Y value at index.
     * @param idy Index (0-based)
     * @param value Value to set
     */
    inline void PutY(std::size_t idy, G4double value);
    /**
     * @brief Set value at (idx, idy).
     * @param idx X index (0-based)
     * @param idy Y index (0-based)
     * @param value Value to set
     */
    inline void PutValue(std::size_t idx, std::size_t idy, G4double value);
    /**
     * @brief Set X and Y vectors.
     * @param vecX X vector
     * @param vecY Y vector
     */
    void PutVectors(const std::vector<G4double>& vecX, const std::vector<G4double>& vecY);

    /**
     * @brief Scale all values of the vector by factor.
     * This method may be applied for example after Retrieve a vector from an external file to
     * convert values into Geant4 units
     *
     * @param factor Scale factor
     */
    void ScaleVector(G4double factor);

    /**
     * @brief Find X using linear interpolation for Y-vector filled by cumulative probability.
     * @param rand Random value [0,1]
     * @param y Y value
     * @param lastidy Last Y index
     * @return Interpolated X value
     */
    G4double FindLinearX(G4double rand, G4double y, std::size_t& lastidy) const;
    /**
     * @brief Find X using linear interpolation for Y-vector filled by cumulative probability.
     * @param rand Random value [0,1]
     * @param y Y value
     * @return Interpolated X value
     */
    inline G4double FindLinearX(G4double rand, G4double y) const;
    // Find Y using linear interpolation for Y-vector filled by cumulative
    // probability function value of rand should be between 0 and 1

    /**
     * @brief Get X value at index.
     * @param index Index
     * @return X value
     */
    inline G4double GetX(std::size_t index) const;
    /**
     * @brief Get Y value at index.
     * @param index Index
     * @return Y value
     */
    inline G4double GetY(std::size_t index) const;
    /**
     * @brief Get value at (idx, idy).
     * @param idx X index
     * @param idy Y index
     * @return Value
     */
    inline G4double GetValue(std::size_t idx, std::size_t idy) const;
    // Returns simply the values of the vector by index
    // of the energy vector. The boundary check will not be done

    /**
     * @brief Find the bin# in which theEnergy belongs. Starting from 0
     * @param x X value
     * @param lastidx Last X index
     * @return Bin index
     */
    inline std::size_t FindBinLocationX(const G4double x, const std::size_t lastidx) const;
    /**
     * @brief Find the bin# in which theEnergy belongs. Starting from 0
     * @param y Y value
     * @param lastidy Last Y index
     * @return Bin index
     */
    inline std::size_t FindBinLocationY(const G4double y, const std::size_t lastidy) const;
    // Find the bin# in which theEnergy belongs. Starting from 0

    /**
     * @brief Get the lengths of the vector (X dimension).
     * @return Length X
     */
    inline std::size_t GetLengthX() const;
    /**
     * @brief Get the lengths of the vector (Y dimension).
     * @return Length Y
     */
    inline std::size_t GetLengthY() const;
    // Get the lengths of the vector

    /**
     * @brief Get physics vector type.
     * @return Physics vector type
     */
    inline G4PhysicsVectorType GetType() const;
    // Get physics vector type

    /**
     * @brief Activate/deactivate bicubic interpolation.
     * @param flag Enable or disable
     */
    inline void SetBicubicInterpolation(G4bool flag);
    // Activate/deactivate bicubic interpolation

    /**
     * @brief Store persistent data to file stream.
     * @param fOut Output file stream
     */
    void Store(std::ofstream& fOut) const;
    /**
     * @brief Retrieve persistent data from file stream.
     * @param fIn Input file stream
     * @return Success flag
     */
    G4bool Retrieve(std::ifstream& fIn);
    // To store/retrieve persistent data to/from file streams

    /**
     * @brief Set verbose level.
     * @param value Verbosity level
     */
    inline void SetVerboseLevel(G4int value);

  protected:

    void PrepareVectors();

    void ClearVectors();

    void CopyData(const G4Physics2DVector& vec);

    G4double BicubicInterpolation(const G4double x, const G4double y, const std::size_t idx,
                                  const std::size_t idy) const;
    // Bicubic interpolation of 2D vector

    inline std::size_t FindBin(const G4double z, const G4PV2DDataVector&, const std::size_t idz,
                               const std::size_t idzmax) const;

  private:

    G4double InterpolateLinearX(G4PV2DDataVector& v, G4double rand) const;

    inline G4double DerivativeX(std::size_t idx, std::size_t idy, G4double fac) const;
    inline G4double DerivativeY(std::size_t idx, std::size_t idy, G4double fac) const;
    inline G4double DerivativeXY(std::size_t idx, std::size_t idy, G4double fac) const;
    // computation of derivatives

    G4PhysicsVectorType type = T_G4PhysicsFreeVector;
    // The type of PhysicsVector (enumerator)

    std::size_t numberOfXNodes = 0;
    std::size_t numberOfYNodes = 0;

    G4PV2DDataVector xVector;
    G4PV2DDataVector yVector;
    std::vector<G4PV2DDataVector*> value;

    G4int verboseLevel = 0;
    G4bool useBicubic = false;
};

#include "G4Physics2DVector.icc"

#endif
