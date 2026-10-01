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
// G4LogicalVolume
//
// Class description:
//
// Represents a leaf node or unpositioned subtree in the geometry hierarchy.
// Logical volumes are named, and may have daughters ascribed to them.
// They are responsible for retrieval of the physical and tracking attributes
// of the physical volume that it represents: solid, material, magnetic field,
// and optionally, user limits, sensitive detectors, regions, biasing weights.
//
// Get and Set functionality is provided for all attributes, but note that
// most set functions should not be used  when the geometry is `closed'.
// As a  further development, `Guard' checks can be added to ensure
// only legal operations at tracking time.
//
// On construction, solid, material and name must be specified.
//
// Daughters are ascribed and managed by means of a simple
// GetNoDaughters,Get/SetDaughter(n),AddDaughter interface.
//
// Smart voxels as used for tracking optimisation. They're also an attribute.
//
// Logical volumes self register to the logical volume Store on construction,
// and deregister on destruction.
//
// NOTE: This class is NOT meant to act as base class, except for exceptional
//       circumstances of extended types used in the kernel.
//
// Data members:
//
//    std::vector<G4VPhysicalVolume*> fDaughters
//    - Vector of daughters. Given initial size of 0.
//    G4FieldManager* fFieldManager
//    - Pointer (possibly 0) to (magnetic or other) field manager object.
//    G4Material* fMaterial
//    - Pointer to material at this node.
//    G4String fName
//    - Name of logical volume.
//    G4VSensitiveDetector *fSensitiveDetector
//    - Pointer (possibly 0) to `Hit' object.
//    G4VSolid* fSolid
//    - Pointer to solid.
//    G4UserLimits* fUserLimits
//    - Pointer (possibly 0) to user Step limit object for this node.
//    G4SmartVoxelHeader* fVoxel
//    - Pointer (possibly 0) to optimisation info objects.
//    G4bool fOptimise
//    - Flag to identify if optimisation should be applied or not.
//    G4bool fRootRegion
//    - Flag to identify if the logical volume is a root region.
//    G4double fSmartless
//    - Quality for optimisation, average number of voxels to be spent
//      per content.
//    const G4VisAttributes* fVisAttributes
//    - Pointer (possibly 0) to visualization attributes.
//    G4Region* fRegion
//    - Pointer to the cuts region (if any)
//    G4MaterialCutsCouple* fCutsCouple
//    - Pointer (possibly 0) to associated production cuts.
//    G4double fBiasWeight
//    - Weight used in the event biasing technique.
//
// Following data members has been moved to G4Region - M.Asai (Aug/18/2005)
//    G4FastSimulationManager* fFastSimulationManager
//    - Pointer (possibly 0) to G4FastSimulationManager object.
//    G4bool fIsEnvelope
//    - Flags if the Logical Volume is an envelope for a FastSimulationManager.

// Author: Paul Kent (CERN), 11.07.1995 - Initial version
// ------------------------------------------------------------------------
#ifndef G4LOGICALVOLUME_HH
#define G4LOGICALVOLUME_HH

#include <atomic>
#include <memory>
#include <vector>

#include "G4GeomSplitter.hh"  // Needed for MT RW data splitting
#include "G4LVData.hh"        // Split-class MT storage for SD/field manager
#include "G4LVMutableData.hh" // Split-class MT storage for solid/material/mass/cuts
#include "G4Region.hh"        // Required by inline methods
#include "G4Threading.hh"
#include "G4Types.hh"
#include "G4VPhysicalVolume.hh"  // Need operator == for vector fdaughters

// Forward declarations
//
class G4FieldManager;
class G4Material;
class G4VSensitiveDetector;
class G4VSolid;
class G4UserLimits;
class G4SmartVoxelHeader;
class G4FastSimulationManager;
class G4MaterialCutsCouple;
class G4VisAttributes;

/**
 * @brief G4LogicalVolume represents a leaf node or unpositioned subtree in the
 * geometry hierarchy.
 * @ingroup geometry_management
 *
 * Logical volumes are named, and may have daughters
 * ascribed to them. They are responsible for retrieval of the physical and
 * tracking attributes of the physical volume that it represents: solid,
 * material, magnetic field, and optionally, user limits, sensitive detectors,
 * regions, biasing weights.
 */

class G4LogicalVolume
{
  public:

    /**
     * Constructor for G4LogicalVolume. The solid and material pointer must be
     * non null. The parameters for field, detector and user limits are optional.
     * The volume also enters itself into the logical volume Store.
     * Optimisation of the geometry (voxelisation) for the volume hierarchy is
     * applied by default. For parameterised volumes in the hierarchy,
     * optimisation is -always- applied.
     *  @param[in] pSolid Pointer to the associated solid primitive.
     *  @param[in] pMaterial Pointer to the associated material.
     *  @param[in] name The volume name.
     *  @param[in] pFieldMgr Pointer to optional magnetic field manager.
     *  @param[in] pSDetector Pointer to optional associated sensitive detector.
     *  @param[in] pULimits Pointer to optional user limits.
     *  @param[in] optimise Flag to enable/disable optimisation structure.
     */
    G4LogicalVolume(G4VSolid* pSolid, G4Material* pMaterial, const G4String& name,
                    G4FieldManager* pFieldMgr = nullptr, G4VSensitiveDetector* pSDetector = nullptr,
                    G4UserLimits* pULimits = nullptr, G4bool optimise = true);

    /**
     * Destructor. Removes the logical volume from the logical volume Store.
     * This class is NOT meant to act as base class, except for exceptional
     * circumstances of extended types used in the kernel.
     */
    virtual ~G4LogicalVolume();

    /**
     * Copy-constructor and assignment operator not allowed.
     */
    G4LogicalVolume(const G4LogicalVolume&) = delete;
    G4LogicalVolume& operator=(const G4LogicalVolume&) = delete;

    /**
     * Returns and sets the name of the logical volume.
     */
    inline const G4String& GetName() const;
    void SetName(const G4String& pName);

    /**
     * Returns the number of daughters (0 to n).
     */
    inline std::size_t GetNoDaughters() const;

    /**
     * Returns the ith daughter. Note numbering starts from 0,
     * and no bounds checking is performed.
     */
    inline G4VPhysicalVolume* GetDaughter(const std::size_t i) const;

    /**
     * Adds the volume 'p' as a daughter of the current logical volume.
     */
    void AddDaughter(G4VPhysicalVolume* p);

    /**
     * Returns true if the volume 'p' is a daughter of the current
     * logical volume.
     */
    inline G4bool IsDaughter(const G4VPhysicalVolume* p) const;

    /**
     * Returns true if the volume 'p' is part of the hierarchy of volumes
     * established by the current logical volume. Scans recursively the volume
     * tree.
     */
    G4bool IsAncestor(const G4VPhysicalVolume* p) const;

    /**
     * Removes the volume 'p' from the list of daughters of the current
     * logical volume.
     */
    void RemoveDaughter(const G4VPhysicalVolume* p);

    /**
     * Clears the list of daughters. Used by the physical volume store when
     * the geometry tree is cleared, since modified at run-time.
     */
    void ClearDaughters();

    /**
     * Returns the total number of physical volumes (replicated or placed)
     * in the tree represented by the current logical volume.
     */
    G4int TotalVolumeEntities() const;

    /**
     * Characterises the daughters of this logical volume.
     */
    inline EVolume CharacteriseDaughters() const;

    /**
     * Utility method used by CharacteriseDaughters().
     */
    inline EVolume DeduceDaughtersType() const;

    /**
     * Gets and sets the current solid.
     */
    inline G4VSolid* GetSolid() const;
    void SetSolid(G4VSolid* pSolid);

    /**
     * Gets and sets the current material.
     */
    inline G4Material* GetMaterial() const;
    void SetMaterial(G4Material* pMaterial);

    /**
     * Sets the material and corresponding Material-Cuts-Couple.
     * This method is invoked by G4Navigator while it is navigating through
     * material parameterisation.
     */
    void UpdateMaterial(G4Material* pMaterial);

    /**
     * Returns the mass of the logical volume tree computed from the
     * estimated geometrical volume of each solid and material associated
     * to the logical volume and (by default) to its daughters.
     *  @note The computation may require a considerable amount of time,
     *        depending from the complexity of the geometry tree.
     *        The returned value is cached and can be used for successive
     *        calls (default), unless recomputation is forced by providing
     *        'true' for the Boolean argument in input. Computation should
     *        be forced if the geometry setup has changed after the previous
     *        call. By setting the 'propagate' Boolean flag to 'false' the
     *        method returns the mass of the present logical volume only
     *        (subtracted for the volume occupied by the daughter volumes).
     *        An optional argument to specify a material is also provided.
     *  @param[in] forced Flag to force recomputation or use cached value.
     *  @param[in] propagate Flag to limit or not computation to daughters.
     *  @param[in] parMaterial Optional pointer to a custom material, usually
     *             used in parameterisations.
     */
    G4double GetMass(G4bool forced = false, G4bool propagate = true,
                     G4Material* parMaterial = nullptr);

    /**
     * Resets the cached value of mass. Ensures that cached value of mass is
     * invalidated due to change in  state, e.g. change of the size of the
     * solid, change of the type of solid, or the addition/deletion of a
     * daughter volume.
     */
    inline void ResetMass();

    /**
     * Gets current Field Manager pointer.
     */
    G4FieldManager* GetFieldManager() const;

    /**
     * Sets the Field Manager and propagates it.
     *  @param[in] pFieldMgr Pointer to the field manager.
     *  @param[in] forceToAllDaughters Flag to force propagation to all
     *             daughters (true), or only to daughters having null pointer
     *             to the field manager (false).
     */
    void SetFieldManager(G4FieldManager* pFieldMgr, G4bool forceToAllDaughters);

    /**
     * Gets and sets the current sensitive detector (can be a null pointer).
     * Out-of-line (defined in G4LogicalVolume.cc, not the .icc): it always
     * dereferences the per-thread split-class storage, and that dereference
     * must stay inside the geometry management DLL/library that owns it
     * (see the "no split-class offset access in headers" note above
     * SolidSlot() below).
     */
    G4VSensitiveDetector* GetSensitiveDetector() const;
    void SetSensitiveDetector(G4VSensitiveDetector* pSDetector);

    /**
     * Gets and sets the current User Limits.
     */
    inline G4UserLimits* GetUserLimits() const;
    inline void SetUserLimits(G4UserLimits* pULimits);

    /**
     * Gets and sets the current Voxel Header.
     */
    inline G4SmartVoxelHeader* GetVoxelHeader() const;
    inline void SetVoxelHeader(G4SmartVoxelHeader* pVoxel);

    /**
     * Gets and sets the user defined optimisation quality associated to
     * the volume.
     */
    inline G4double GetSmartless() const;
    inline void SetSmartless(G4double s);

    /**
     * Replies if the geometry optimisation (voxelisation) is to be applied
     * for this volume hierarchy.
     */
    inline G4bool IsToOptimise() const;

    /**
     * Specifies if to apply or not geometry optimisation to the volume
     * hierarchy. For parameterised volumes in the hierarchy, optimisation is
     * always applied.
     */
    inline void SetOptimisation(G4bool optim);

    /**
     * Replies if the logical volume represents a root region or not.
     */
    inline G4bool IsRootRegion() const;

    /**
     * Sets/unsets the volume as a root region for cuts.
     */
    inline void SetRegionRootFlag(G4bool rreg);

    /**
     * Replies if the logical volume is part of a cuts region or not.
     */
    inline G4bool IsRegion() const;

    /**
     * Sets/unsets the volume as cuts region.
     */
    inline void SetRegion(G4Region* reg);

    /**
     * Returns the region to which the volume belongs, if any.
     */
    inline G4Region* GetRegion() const;

    /**
     * Propagates region pointer to daughters.
     */
    inline void PropagateRegion();

    /**
     * Accessor and modifier for production cuts.
     */
    inline const G4MaterialCutsCouple* GetMaterialCutsCouple() const;
    void SetMaterialCutsCouple(G4MaterialCutsCouple* cuts);

    /**
     * Equality defined by address only.
     * Returns true if objects are at same address, else false.
     */
    G4bool operator==(const G4LogicalVolume& lv) const;

    /**
     * Accessor and modifiers for visualization attributes.
     * Arguments are converted to shared_ptr.
     */
    const G4VisAttributes* GetVisAttributes() const;
    void SetVisAttributes(const G4VisAttributes* pVA);
    void SetVisAttributes(const G4VisAttributes& VA);

    /**
     * Gets the current FastSimulationManager pointer if existing,
     * otherwise null.
     */
    inline G4FastSimulationManager* GetFastSimulationManager() const;

    /**
     * Sets and gets the bias weight.
     */
    inline void SetBiasWeight(G4double w);
    inline G4double GetBiasWeight() const;

    /**
     * Returns true if it is not a base-class object.
     */
    virtual G4bool IsExtended() const;

    /**
     * Returns current Field Manager for the master thread.
     */
    inline G4FieldManager* GetMasterFieldManager() const;

    /**
     * Returns current Sensitive Detector for the master thread.
     */
    inline G4VSensitiveDetector* GetMasterSensitiveDetector() const;

    /**
     * Returns current Solid for the master thread.
     */
    inline G4VSolid* GetMasterSolid() const;

    /**
     * Returns the instance ID.
     */
    inline G4int GetInstanceID() const;

    /**
     * Returns the private data instance manager.
     */
    static const G4LVManager& GetSubInstanceManager();

    /**
     * Returns the private data instance manager for the mutable
     * (solid/material/mass/cuts-couple) per-thread data.
     */
    static const G4LVMutableManager& GetMutableSubInstanceManager();

    /**
     * Marks this logical volume's solid, material, mass and cuts-couple as
     * needing a per-thread copy, instead of being shared across all
     * threads. Called automatically when the volume is used as the
     * daughter of a G4PVParameterised, G4PVDivision or G4ReplicatedSlice
     * (i.e. VolumeType() == kParameterised); not needed, and not to be
     * called, for plain placements or replicas. Idempotent. Must be called
     * only during initial detector construction, before the geometry is
     * first closed.
     */
    void MarkMaterialMutable();

    /**
     * Returns true if this logical volume's solid/material/mass/cuts-couple
     * are stored per-thread (mutable), false if they are shared across all
     * threads. See MarkMaterialMutable().
     */
    inline G4bool IsMaterialMutable() const;

    /**
     * Clears memory allocated by sub-instance managers.
     */
    static void Clean();

    /**
     * Sets lock identifier for final deletion of entity.
     */
    inline void Lock();

    /**
     * This method is similar to the constructor. It is used by each worker
     * thread to achieve the partial effect as that of the master thread.
     */
    void InitialiseWorker(G4LogicalVolume* ptrMasterObject, G4VSolid* pSolid,
                          G4VSensitiveDetector* pSDetector);

    /**
     * This method is similar to the destructor. It is used by each worker
     * thread to achieve the partial effect as that of the master thread.
     */
    void TerminateWorker(G4LogicalVolume* ptrMasterObject);

    /**
     * Sets the Field Manager only at this level (does not push down hierarchy)
     */
    void AssignFieldManager(G4FieldManager* fldMgr);

    /**
     * Changes the type of the daughters volume to be of type 'atype'.
     * Meant for the user adopting an external navigator for the contents
     * of a volume.
     *  @returns Success (true) or failure (false).
     */
    G4bool ChangeDaughtersType(EVolume atype);

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4LogicalVolume(__void__&);

  private:

    /**
     * Out-of-line: dereferences subInstanceManager.offset, a dllimport
     * thread_local on Windows outside this class's owning DLL. Must stay
     * out-of-line (defined in G4LogicalVolume.cc) so that dereference never
     * appears in a header seen by other modules. See the note above
     * SolidSlot() below for the general rule this follows.
     */
    G4LVData& GetLVData() const;

    /**
     * Returns the mass of a daughter logical volume tree computed from the
     * estimated geometrical volume of each associated solid and material.
     *  @param[in] physDaughter Pointer to the daughter physical volume.
     *  @param[in] motherDensity Density value of the mother volume.
     *  @param[in] propagate Flag to limit or not computation to daughters.
     */
    G4double ComputeDaughterMassContribution(G4VPhysicalVolume* physDaughter,
                                             G4double motherDensity,
                                             G4bool propagate);

    /**
     * Computes mass without reading or publishing the cached value. Used for
     * recursive daughter contributions, including parameterised material cases.
     */
    G4double ComputeMassNoCache(G4bool propagate, G4Material* parMaterial);

    /**
     * Accessors which isolate, in one place, the storage difference between
     * - a "mutable" logical volume (per-thread G4LVMutableData, chosen for
     *   volumes used inside a parameterised/division-type placement) and
     * - a "shared" one (fields of this object, identical on every thread)
     *   for the solid, material and cuts-couple, decided at runtime by
     * fInstanceIDifMutable. See MarkMaterialMutable().
     *
     * Used by Get/Set Solid, Get/Set Material, UpdateMaterial,
     * and Get/Set MaterialCutsCouple, so that each method has a single
     * implementation valid regardless of whether the volume is mutable or
     * shared.
     */

    G4VSolid*& SolidSlot() const;
    G4Material*& MaterialSlot() const;
    G4MaterialCutsCouple*& CutsCoupleSlot() const;

    /** Out-of-line mutable-branch bodies for GetSolid/GetMaterial/
     *  GetMaterialCutsCouple/GetCachedMass/SetCachedMass. See the note
     *  above SolidSlot() for why these must not be inlined. */
    G4VSolid* GetSolidMutable() const;
    G4Material* GetMaterialMutable() const;
    G4MaterialCutsCouple* GetCutsCoupleMutable() const;
    G4double GetCachedMassMutable() const;
    void SetCachedMassMutable(G4double m);

    /**
     * Accessor/modifier for the cached mass value, isolating the storage
     * difference between a mutable (per-thread plain double, in
     * G4LVMutableData) and a shared (this object's atomic) logical volume.
     * Inline only for the shared (fImmutable) fast path; see the note
     * above SolidSlot().
     */
    inline G4double GetCachedMass() const;
    inline void SetCachedMass(G4double m);

    /**
     * A logical volume's solid/material/cuts-couple are either "shared" (this
     * object's own fields, identical on every thread — the common case) or
     * "mutable" (a per-thread G4LVMutableData slot, for volumes used inside a
     * parameterised/division-type placement; see MarkMaterialMutable()).
     * Setters on a *shared* volume are forbidden from worker threads (they
     * would silently corrupt the single field read by every other thread);
     * setters on a *mutable* volume are expected to be called from worker
     * threads routinely (that's what the navigator does for every
     * parameterised-volume entry).
     */
    void CheckGeometryConfiguration(const char* methodName,
                                    const G4String& volName,
                                    G4bool isMutable);
   
  private:

    using G4PhysicalVolumeList = std::vector<G4VPhysicalVolume*>;

    /** This field helps in the use of the class G4LVManager. */
    G4GEOM_DLL static G4LVManager subInstanceManager;

    /** This field helps in the use of the class G4LVMutableManager, used
     *  only for logical volumes marked mutable (see below.)
     */
    G4GEOM_DLL static G4LVMutableManager subInstanceManagerMutable;

    /** Vector of daughters. Given initial size of 0. */
    G4PhysicalVolumeList fDaughters;

    /** Name of logical volume. */
    G4String fName;

    /** Pointer (possibly nullptr) to user step limit object for this node. */
    G4UserLimits* fUserLimits = nullptr;

    /** Pointer (possibly nullptr) to optimisation info objects. */
    G4SmartVoxelHeader* fVoxel = nullptr;

    /** Optimisation quality, average number of voxels to be spent per content. */
    G4double fSmartless = 2.0;

    /** Pointer to the cuts region (if any). */
    G4Region* fRegion = nullptr;

    /** Weight used in the event biasing technique. */
    G4double fBiasWeight = 1.0;

    /** Pointer to visualization attributes. */
    std::shared_ptr<const G4VisAttributes> fVisAttributes;

    /** Originally shadow of master pointers.
     *  * Each worker thread can access this field from the master thread
     *    through these pointers.
     *  If attributes are immutable, they also used in threads!
     */
    G4VSolid* fSolid = nullptr;
    G4VSensitiveDetector* fSensitiveDetector = nullptr;
    G4FieldManager* fFieldManager = nullptr;

   /**
    * 'Master' values for attributes: for solid/material/cuts-couple/mass.
    * Also used as common, shared storage whenever all its uses are
    * unchanging - ie typically placements or replica
    *
    * For a *mutable* volume, fSolid/fMaterial/fCutsCouple instead act as a
    * master-thread shadow: 
    */
    G4Material* fMaterial = nullptr;
    G4MaterialCutsCouple* fCutsCouple = nullptr;
    std::atomic<G4double> fMass{0.0};

    /** This new field is used as instance ID. */
    G4int instanceID;

    /**
     *  Index for per-thread G4LVMutableData slots, or -1 if immutable.
     */
    G4int fInstanceIDifMutable = -1;

    /**
     *  Mark volume's attributes as 'unchangeable' - ie Placement / Replica
     *  Mirrors (fInstanceIDifMutable == -1) -- it is kept up-to-date. 
     *  Only an optimisation for 'hot' accessors.
     */
    G4bool fImmutable = true;

    /** Are contents of volume placements, replica, parameterised or external? */
    EVolume fDaughtersVolumeType;

    /** Flag to identify if optimisation should be applied or not. */
    G4bool fOptimise = true;

    /** Flag to identify if the logical volume is a root region. */
    G4bool fRootRegion = false;

    /** Flag to identify if entity is locked for final deletion. */
    G4bool fLock = false;
};

#include "G4LogicalVolume.icc"

// NOTE: "instanceID" is the subscript, for this instance, into the array of
// G4LVData objects (one per thread) managed by "subInstanceManager". See
// G4LVData.hh for the full split-class mechanism.

#endif
