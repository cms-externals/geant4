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
// G4GeometryManager
//
// Class description:
//
// A class responsible for high level geometrical functions, and for
// high level objects in the geometry subdomain.
// The class is a 'singleton', with access via the static method
// G4GeometryManager::GetInstance().

// Author: Paul Kent (CERN), 26.07.1995 - Initial version
//         John Apostolakis (CERN), 12.06.2024 - Added parallel optimisation
// --------------------------------------------------------------------
#ifndef G4GEOMETRYMANAGER_HH
#define G4GEOMETRYMANAGER_HH

#include "G4SmartVoxelStat.hh"
#include "G4Types.hh"
#include "G4ios.hh"

#include <memory>
#include <vector>

class G4VPhysicalVolume;
class G4Timer;
class G4VoxelisationHelper;

/**
 * @brief G4GeometryManager is a singleton class responsible for high level
 * geometrical functions, and for high level objects in the geometry subdomain.
 * @ingroup geometry_management
 */

class G4GeometryManager
{
  public:

    /**
     * Destructor; called by G4RunManagerKernel.
     */
    ~G4GeometryManager();

    /**
     * Returns a pointer to the singleton instance of the class, creating
     * it if not existing.
     */
    static G4GeometryManager* GetInstance();

    /**
     * Simply returns a pointer to the singleton instance.
     */
    static G4GeometryManager* GetInstanceIfExist();

    /**
     * Closes ('locks') the geometry: performs sanity and 'completion' checks
     * and optionally [default=yes] builds the optimisation structure.
     * Applies to just a specific subtree if a physical volume is specified.
     *  @param[in] pOptimise Flag to enabling/disabling optimisation structure.
     *  @param[in] verbose Flag for verbosity.
     *  @param[in] vol Optional pointer to a physical volume (subtree) for
     *             optimisation.
     *  @returns true if process succeeds.
     */
    G4bool CloseGeometry(G4bool pOptimise = true, G4bool verbose = false,
                         G4VPhysicalVolume* vol = nullptr);

    /**
     * Opens ('unlocks') the geometry and removes the optimisation structure if
     * present. Applies to just a specific subtree if a physical volume is
     * specified.
     *  @param[in] vol Optional pointer to a physical volume (subtree) for
     *             optimisation.
     */
    void OpenGeometry(G4VPhysicalVolume* vol = nullptr);

    /**
     * Returns true/false according to the state of the optimised geometry.
     */
    G4bool IsGeometryClosed();

    /**
     * Sets the maximum extent of the world volume.
     * The operation is allowed only if NO solids have been created already.
     */
    void SetWorldMaximumExtent(G4double worldExtent);

    // Methods to choose or undertake Parallel Optimisation
    // ----------------------------------------------------

    /**
     * Requests optimisation using threads (if MT is enabled & used ).
     */
    void OptimiseInParallel(G4bool val = true);

    /**
     * Method that contributes to (voxel) optimisation until all work is done.
     * To be called by a worker thread initialisation - not by the user.
     */
    void UndertakeOptimisation();

    /**
     * Detailed method for user to request parallel optimisation (if verbosity
     * is required). Calling this method is enough to ask for it.
     * It will be used if Geant4 is built with MT/tasks.
     * Parallelism will be used if Geant4 is built with MT/tasks.
     */
    void RequestParallelOptimisation(G4bool val = true, G4bool verbose = true);

    /**
     * Checks whether parallel optimisation was requested and configured.
     */
    G4bool IsParallelOptimisationConfigured();

    /**
     * Reports whether parallel optimisation was completed.
     */
    G4bool IsParallelOptimisationFinished();

    /**
     * Reports verbose parallel voxel statistics once optimisation is complete.
     * To be called by the master after synchronising with the workers.
     */
    void ReportParallelOptimisationStats();

    /**
     * Forces G4LogicalVolume::GetMass() daughter traversal to run serially.
     */
    void ForceSerialGetMass(G4bool val = true);

    /**
     * Reports whether serial G4LogicalVolume::GetMass() traversal is forced.
     */
    G4bool IsSerialGetMassForced();

    /**
     * Forces GetPointOnSurface()-related on-demand solid estimates to run
     * serially.
     */
    void ForceSerialGetPointOnSurface(G4bool val = true);

    /**
     * Reports whether serial GetPointOnSurface() processing is forced.
     */
    G4bool IsSerialGetPointOnSurfaceForced();

    /**
     * Forces geometry overlap checks to run serially.
     */
    void ForceSerialOverlapsCheck(G4bool val = true);

    /**
     * Reports whether serial overlap checking is forced.
     */
    G4bool IsSerialOverlapsCheckForced();

    /**
     * Controls the consistency pass at the start of CloseGeometry() that
     * detects (and repairs) 'unorthodox' parameterisations which reassign
     * the logical volume of their physical volume, via
     * G4VPhysicalVolume::SetLogicalVolume(), to logical volumes that were not
     * marked mutable at construction. On by default (as of the introduction
     * of UseNestedParameterisationCheckAtClose(), below); call with val=false
     * to disable, e.g. for performance-sensitive geometries known not to use
     * such parameterisations. See G4ParameterisationChecker.
     */
    void CheckParameterisedReassignmentAtClose(G4bool val = true);

    /**
     * Reports whether the parameterisation-reassignment consistency pass has
     * been requested.
     */
    G4bool IsParameterisedReassignmentCheckRequested() const;

    /**
     * Selects the stronger, nested-aware consistency pass
     * (G4ParameterisationCheckVisitor, in the navigation module) in place of
     * the flat G4ParameterisationChecker used by default. The nested-aware
     * pass visits every leaf of the geometry tree with a live parent
     * touchable, so it also properly exercises parameterisations nested
     * inside one or two levels of replica (G4VNestedParameterisation),
     * which the flat checker cannot. Off by default, since it is more
     * costly (a full tree walk rather than a flat pass over declared
     * multiplicities). Has no effect unless the navigation module has
     * registered itself (see RegisterNestedParameterisationChecker()) and
     * the basic check above is also enabled; otherwise CloseGeometry() falls
     * back to the flat checker and issues a warning.
     */
    void UseNestedParameterisationCheckAtClose(G4bool val = true);

    /**
     * Reports whether the nested-aware parameterisation check has been
     * requested.
     */
    G4bool IsNestedParameterisationCheckRequested() const;

    /**
     * Signature of the nested-aware parameterisation checker entry point
     * (G4ParameterisationCheckVisitor::CheckGeometry), matched by
     * RegisterNestedParameterisationChecker() below.
     */
    using NestedParamCheckCallback = G4int (*)(G4VPhysicalVolume*, G4bool, G4bool);

    /**
     * Registers the nested-aware parameterisation checker implemented in the
     * navigation module. The management module must not depend on
     * navigation headers, so navigation registers its entry point here
     * instead, mirroring the pattern used by
     * G4GeomTaskHelper::RegisterWorkspaceCallbacks().
     */
    static void RegisterNestedParameterisationChecker(NestedParamCheckCallback checker);

  private:

    /**
     * Private Constructor. Set the geometry to be open.
     */
    G4GeometryManager();

    // Create the voxel optimisation information. Return whether work is completed
    // (In case of parallel optimisation, it will be done later by the workers)

    /**
     * Optimises all volumes or just multi-volumes (parameterisations, etc. ).
     */
    G4bool BuildOptimisations(G4bool allOpt, G4bool verbose = false);

    /**
     * Optimises one volume or subtree only.
     */
    G4bool BuildOptimisations(G4bool allOpt, G4VPhysicalVolume* vol);

    /**
     * Builds the optimisation structure in sequential mode.
     */
    void BuildOptimisationsSequential(G4bool allOpts, G4bool verbose = true);

    /**
     * Methods to delete the optimisation structure.
     */
    void DeleteOptimisations();
    void DeleteOptimisations(G4VPhysicalVolume* vol);

  private:

    /** The static instance. */
    static G4ThreadLocal G4GeometryManager* fgInstance;

    /** Shared ownership of voxelisation helper across thread-local instances. */
    static std::shared_ptr<G4VoxelisationHelper> fParallelVoxeliser;

    /** Flag indicating if geometry is closed or not. */
    G4bool fIsClosed = false;

    // Flags for parallel initialization
    // ---------------------------------

    /** Request parallel optimisation. */
    static G4bool fParallelVoxelOptimisationRequested;

    /** Not just requested, but adopted (i.e. also in MT/tasking mode). */
    static G4bool fOptimiseInParallelConfigured;

    // Flags for forced serial processing tasks
    // ----------------------------------------

    /** Force serial processing of mass computation. */
    static G4bool fForceSerialGetMass;

    /** Force serial processing of point-on-surface related estimates. */
    static G4bool fForceSerialGetPointOnSurface;

    /** Force serial processing of overlaps checking. */
    static G4bool fForceSerialOverlapsCheck;

    /** Request the parameterisation-reassignment consistency pass at close. */
    static G4bool fCheckParamReassignment;

    /** Request the stronger, nested-aware variant of that pass. */
    static G4bool fUseNestedParamChecker;

    /** Nested-aware checker entry point, registered by the navigation module. */
    static NestedParamCheckCallback fNestedParamChecker;
};

#endif
