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
// G4VCSGfaceted
//
// Class description:
//
// Virtual class defining CSG-like type shape that is built entirely
// of G4CSGface faces.

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4VCSGFACETED_HH
#define G4VCSGFACETED_HH

#include "G4VSolid.hh"
#include "G4VCSGface.hh"

#include <CLHEP/Units/PhysicalConstants.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

class G4VisExtent;

/**
 * @brief G4VCSGfaceted is a virtual class defining a CSG-like type shape
 * that is built entirely of G4CSGface faces.
 * @ingroup geometry_solids_specific
 */

class G4VCSGfaceted : public G4VSolid
{
  public:

    /**
     * Constructor taking a 'name'.
     */
    G4VCSGfaceted(const G4String& name);

    /**
     * Destructor.
     */
    ~G4VCSGfaceted() override;

    /**
     * Copy constructor and assignment operator.
     */
    G4VCSGfaceted(const G4VCSGfaceted& source);
    G4VCSGfaceted& operator=(const G4VCSGfaceted& source);

    /**
     * Calculates the minimum and maximum extent of the solid, when under the
     * specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pTransform The internal transformation applied to the solid.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the solid is intersected by the extent region.
     */
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform, G4double& pmin,
                           G4double& pmax) const override;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in G4VSolid.
     */
    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false, G4bool* validNorm = nullptr,
                           G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4CSGfaceted" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Returns a pointer to a generated polyhedron used for visualisation.
     */
    G4Polyhedron* CreatePolyhedron() const override = 0;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent() const override;
    G4Polyhedron* GetPolyhedron() const override;

    /**
     * Accessors and modifiers for capacity and area computation.
     */
    G4int GetCubVolStatistics() const;
    G4double GetCubVolEpsilon() const;
    void SetCubVolStatistics(G4int st);
    void SetCubVolEpsilon(G4double ep);
    G4int GetAreaStatistics() const;
    G4double GetAreaAccuracy() const;
    void SetAreaStatistics(G4int st);
    void SetAreaAccuracy(G4double ep);

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units. Caches the computed value
     * once computed the first time.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4VCSGfaceted(__void__&);

  protected:

    /**
     * Protected method used in DistanceToIn() and DistanceToOut().
     */
    virtual G4double DistanceTo(const G4ThreeVector& p, const G4bool outgoing) const;

    /**
     * Builds the Z-section optimisation structure for faceted solids whose
     * faces are arranged along Z. The optimisation is conservative: it only
     * rejects faces whose cached Z extent cannot improve the current query.
     */
    void BuildZSectionIndex();

    /**
     * Returns a random point located on the surface of the solid
     * in case of generic Polycone or generic Polyhedra.
     */
    G4ThreeVector GetPointOnSurfaceGeneric() const;

    /**
     * Copy parameters from other solid or reset them.
     * Used in copy constructor and assignment operator.
     */
    void CopyStuff(const G4VCSGfaceted& source);
    void DeleteStuff();

    /**
     * Z section auxiliary methods for fast retrieval.
     */
    inline G4bool HasZSectionIndex() const;
    inline EInside InsideNoZ(const G4ThreeVector& p) const;
    inline EInside InsideZ(const G4ThreeVector& p) const;
    inline std::size_t GetZSlabIndex(G4double z) const;
    static std::size_t FindZSlabIndex(const std::vector<G4double>& boundaries,
                                      G4double z);
    static std::size_t FindZSlabIndex(const std::vector<G4double>& boundaries,
                                      G4double z, G4double edgeTolerance);
    static void BuildZVoxelizer(const std::vector<G4double>& sectionZMin,
                                const std::vector<G4double>& sectionZMax,
                                G4double minHalfZ, std::vector<G4double>& slabBoundary,
                                std::vector<std::vector<G4int>>& slabCandidate);
    template <typename Section>
    static void BuildZVoxelizer(const std::vector<Section>& sections,
                                G4double minHalfZ, std::vector<G4double>& slabBoundary,
                                std::vector<std::vector<G4int>>& slabCandidate);
    const std::vector<G4int>& FindZCandidates(const std::vector<G4double>& boundaries,
                                const std::vector<std::vector<G4int>>& candidates,
                                G4double z) const;
    EInside PhiState(G4bool phiIsOpen, G4double startPhi, G4double endPhi,
                     G4double sStart, G4double cStart,
                     G4double sEnd, G4double cEnd,
                     const G4ThreeVector& p) const;
    G4double DistanceToPhiBoundary(G4bool phiIsOpen,
                                   G4double sStart, G4double cStart,
                                   G4double sEnd, G4double cEnd,
                                   const G4ThreeVector& p,
                                   G4ThreeVector* normal = nullptr) const;
    const std::vector<G4int>& GetZCandidates(G4double z) const;
    inline G4bool FaceCanImprovePoint(G4int i, G4double z,
                                      G4double best, G4double tolerance) const;
    inline G4bool RayIntersectsFaceZ(G4int i, const G4ThreeVector& p,
                                     const G4ThreeVector& v, G4double best,
                                     G4double tolerance,
                                     G4bool radialCull = true) const;
    static inline G4bool RayMissesRadialLimit(const G4ThreeVector& p,
                                              const G4ThreeVector& v,
                                              G4double t0, G4double t1,
                                              G4double radialLimit);
    static inline EInside CombineInside(EInside rz, EInside phi);

  protected:

    /**
     * R/Z section handling data and functions.
     */
    struct G4ZCandidateData
    {
      std::vector<G4double> boundary;
      std::vector<std::vector<G4int>> candidate;
      std::vector<G4int> global;
      // Conservative radial bounds for each Z slab.  They allow point
      // distance queries to reject a complete slab using an exact R/Z
      // rectangle lower bound instead of visiting every intervening slab.
      std::vector<G4double> rMin;
      std::vector<G4double> rMax;
    };

    struct G4RZSection
    {
      G4double r0 = 0., z0 = 0., r1 = 0., z1 = 0.;
      G4double rS = 0., zS = 0., length = 0.;
      G4double rNorm = 0., zNorm = 0.;
      G4double zMin = 0., zMax = 0.;
    };

    static std::unique_ptr<G4ZCandidateData>
           CopyZCandidateData(const std::unique_ptr<G4ZCandidateData>& source);
    static std::unique_ptr<G4ZCandidateData>
           MakeZCandidateData(G4int nFace);

    template <typename Corner>
    static G4bool HasRZEdge(const Corner& a, const Corner& b);

    template <typename Section, typename Corner>
    static G4bool InitRZSection(Section& section,
                                const Corner& a, const Corner& b);

    template <typename Section>
    void AddRZSectionCandidate(const Section& section,
                               std::vector<Section>& sections,
                               G4ZCandidateData* data) const;

    template <typename Section>
    void BuildRZVoxel(const std::vector<Section>& sections,
                      std::unique_ptr<G4ZCandidateData>& data) const;

    template <typename Section>
    inline G4double DistanceToRZSections(G4double radius, G4double z,
                                 const std::vector<Section>& sections,
                                 const std::unique_ptr<G4ZCandidateData>& data,
                                 G4bool outgoing,
                                 G4double metricScale = 1.) const;

    template <typename Section>
    inline EInside InsideRZSections(G4double radius, G4double z,
                             const std::vector<Section>& sections,
                             const std::unique_ptr<G4ZCandidateData>& data) const;

    template <typename Visit>
    static G4bool VisitZSlabsBySafety(const std::vector<G4double>& boundaries,
                                      G4double z, const G4double& best,
                                      G4double tolerance, Visit&& visit);

    template <typename Visit>
    static G4bool VisitRZSlabsBySafety(const G4ZCandidateData& data,
                                      G4double radius, G4double z,
                                      const G4double& best,
                                      G4double tolerance, Visit&& visit);

  protected:

    /** Stack-backed visited tracking for facets, to avoid per-query heap
        allocation in voxel traversal. */
    class G4StackVisited
    {
      public:

        inline explicit G4StackVisited(G4int entries);
        inline G4bool TestAndSet(G4int i);

      private:

        std::array<unsigned long long, 16> fStack = {};
        std::vector<unsigned long long> fHeap;
    };

    G4int numFace = 0;
    G4VCSGface** faces = nullptr;
    G4double fCubicVolume = 0.0;
    G4double fSurfaceArea = 0.0;
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
    std::vector<G4double> fZSectionMin;
    std::vector<G4double> fZSectionMax;
    std::vector<G4double> fRSectionMax;
    std::vector<G4double> fZSlabBoundary;
    std::vector<std::vector<G4int>> fZSlabCandidate;
    std::vector<G4int> fZGlobalCandidate;

    static constexpr G4int kMinZSectionIndexFaces = 17;

  private:

    /** Statistics, error accuracy for volume estimation. */
    G4int fStatistics;
    G4double fCubVolEpsilon;
    G4double fAreaAccuracy;
};

// Inlined methods

#include "G4VCSGfaceted.icc"

#endif
