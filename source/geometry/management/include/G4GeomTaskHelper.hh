//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the  terms  and *
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
// * This  code  implementation is the result of the scientific and   *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4GeomTaskHelper
//
// Class description:
//
// A collection of helper utilities to setup tasking for geometry
// calculations. It centralises use of the PTL thread pool, task-local
// thread data, registered geometry workspaces and deterministic pseudo-random
// sampling helpers used by geometry management and overlap checks.
//
// Author: Gabriele Cosmo (CERN), June 2026
// --------------------------------------------------------------------
#ifndef G4GEOMTASKHELPER_HH
#define G4GEOMTASKHELPER_HH

#include "G4QuickRand.hh"
#include "G4TaskGroup.hh"
#include "G4ThreadPool.hh"
#include "globals.hh"

#include "PTL/ThreadData.hh"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <thread>
#include <utility>

/**
 * @brief G4GeomTaskHelper provides common tasking utilities for geometry
 * calculations.
 *
 * The helper owns no persistent global state. It can reuse an existing
 * PTL/G4 tasking thread pool on the calling thread, create a temporary local
 * pool when appropriate, and restore the previous PTL thread data afterwards.
 *
 * The class is intended for internal geometry code that parallelises
 * independent chunks of work on demand, such as solid estimates, mass
 * traversal and overlap checks.
 *
 * @ingroup geometry_management
 */

class G4GeomTaskHelper
{
  public:

    using WorkspaceCreateCallback = G4bool (*)();
    using WorkspaceCleanupCallback = void (*)();

  private:

    struct WorkspaceCallbacks
    {
        WorkspaceCreateCallback create = nullptr;
        WorkspaceCleanupCallback cleanup = nullptr;
    };

    static constexpr std::size_t kMaxWorkspaceCallbacks = 4;

    static std::array<WorkspaceCallbacks, kMaxWorkspaceCallbacks>&
    GetWorkspaceCallbacks();

  public:

    /**
     * Registers callbacks for task-local workspaces owned by dependent
     * geometry modules.
     *
     * The management module must not depend on modules owning concrete
     * workspaces. Those modules can register their setup and cleanup hooks
     * here while the task helper keeps only opaque function pointers.
     */
    static void RegisterWorkspaceCallbacks(std::size_t slot,
                                           WorkspaceCreateCallback create,
                                           WorkspaceCleanupCallback cleanup);

    /**
     * @brief RAII holder for the tasking pool used by geometry operations.
     *
     * The context reuses the thread pool stored in the current PTL thread
     * data, when present. Otherwise it may create a temporary local
     * G4ThreadPool sized from the requested number of work items. On
     * destruction, any locally created pool is destroyed and the previous PTL
     * thread data is restored.
     *
     * If @p allowWithinTask is false, construction refuses to reuse a context
     * already executing inside a task worker, avoiding nested task submission.
     */
    class PoolContext
    {
      public:

        /**
         * Creates a context able to process up to @p workItems independent
         * chunks. In serial builds the context remains invalid.
         */
        explicit PoolContext(std::size_t workItems,
                             G4bool allowWithinTask = false)
        {
#ifdef G4MULTITHREADED
          fPreviousTaskThreadData = PTL::ThreadData::GetInstance();
          if (fPreviousTaskThreadData != nullptr)
          {
            if (fPreviousTaskThreadData->within_task && !allowWithinTask)
            {
              return;
            }
            fPool = fPreviousTaskThreadData->thread_pool;
          }

          if (fPool == nullptr)
          {
            auto threadCount = std::min<std::size_t>(
              workItems, std::max(1u, std::thread::hardware_concurrency()));
            if (threadCount < 2)
            {
              return;
            }

            G4ThreadPool::Config config;
            config.pool_size = threadCount;
            fLocalPool = std::make_unique<G4ThreadPool>(config);
            fPool = fLocalPool.get();
          }
#else
          (void)workItems;
          (void)allowWithinTask;
#endif
        }

        ~PoolContext() { Reset(); }

        /**
         * Copying is disabled because the context may own a local thread pool
         * and manages PTL thread-local state.
         */
        PoolContext(const PoolContext&) = delete;
        PoolContext& operator=(const PoolContext&) = delete;

        /**
         * Returns the selected thread pool, or nullptr if no pool is usable.
         */
        G4ThreadPool* GetPool() const { return fPool; }

        /**
         * Returns true if a tasking pool is available.
         */
        G4bool IsValid() const { return fPool != nullptr && fPool->size() > 0; }

        /**
         * Returns the number of worker threads in the selected pool.
         */
        std::size_t Size() const
        {
          return (fPool != nullptr) ? fPool->size() : 0;
        }

        /**
         * Destroys any locally owned pool and restores PTL thread data.
         * It is safe to call explicitly; the destructor calls it as well.
         */
        void Reset()
        {
#ifdef G4MULTITHREADED
          if (fLocalPool != nullptr)
          {
            auto* localTaskThreadData = PTL::ThreadData::GetInstance();
            fLocalPool->destroy_threadpool();
            fLocalPool.reset();
            if (localTaskThreadData != fPreviousTaskThreadData)
            {
              delete localTaskThreadData;
            }
            PTL::ThreadData::GetInstance() = fPreviousTaskThreadData;
          }
#else
          (void)fPreviousTaskThreadData;
#endif
          fPool = nullptr;
        }

      private:

        G4ThreadPool* fPool = nullptr;
        std::unique_ptr<G4ThreadPool> fLocalPool;
        PTL::ThreadData* fPreviousTaskThreadData = nullptr;
    };

    /**
     * @brief Ensures PTL thread data exists while scheduling geometry tasks.
     *
     * If the calling thread has no PTL thread data, this guard installs a
     * temporary one bound to the supplied pool and removes it on destruction.
     */
    class ScopedTaskThreadData
    {
      public:

        /**
         * Installs temporary PTL thread data for @p pool when needed.
         */
        explicit ScopedTaskThreadData(G4ThreadPool* pool)
          : fData(pool), fPrevious(PTL::ThreadData::GetInstance())
        {
          if (fPrevious == nullptr)
          {
            PTL::ThreadData::GetInstance() = &fData;
          }
        }

        /**
         * Restores the previous PTL thread-data state.
         */
        ~ScopedTaskThreadData()
        {
          if (fPrevious == nullptr)
          {
            PTL::ThreadData::GetInstance() = nullptr;
          }
        }

        ScopedTaskThreadData(const ScopedTaskThreadData&) = delete;
        ScopedTaskThreadData& operator=(const ScopedTaskThreadData&) = delete;

      private:

        PTL::ThreadData fData;
        PTL::ThreadData* fPrevious = nullptr;
    };

    /**
     * @brief Provides task-local geometry and solids workspaces.
     *
     * Geometry tasks may run on threads that do not yet have the workspaces
     * normally created by worker initialisation. This guard invokes the
     * workspace callbacks registered by dependent modules and destroys only
     * the workspaces it created on destruction.
     */
    class ScopedGeometryWorkspaces
    {
      public:

        /**
         * Creates registered workspaces for the current thread if they are
         * missing.
         */
        ScopedGeometryWorkspaces()
        {
          const auto& callbacks = GetWorkspaceCallbacks();
          for (std::size_t i = 0; i < callbacks.size(); ++i)
          {
            if (callbacks[i].create != nullptr)
            {
              fOwnWorkspace[i] = callbacks[i].create();
            }
          }
        }

        /**
         * Destroys workspaces created by this guard.
         */
        ~ScopedGeometryWorkspaces()
        {
          const auto& callbacks = GetWorkspaceCallbacks();
          for (std::size_t i = 0; i < callbacks.size(); ++i)
          {
            if (fOwnWorkspace[i] && callbacks[i].cleanup != nullptr)
            {
              callbacks[i].cleanup();
            }
          }
        }

        ScopedGeometryWorkspaces(const ScopedGeometryWorkspaces&) = delete;
        ScopedGeometryWorkspaces& operator=(const ScopedGeometryWorkspaces&) = delete;

      private:

        std::array<G4bool, kMaxWorkspaceCallbacks> fOwnWorkspace{};
    };

    /**
     * Mixes @p value into a deterministic non-zero 64-bit seed.
     */
    static std::uint64_t MixSeed(std::uint64_t value)
    {
      value += 0x9e3779b97f4a7c15ULL;
      value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
      value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
      value ^= value >> 31;
      return (value == 0) ? 1 : value;
    }

    /**
     * Returns the next seed from an independent geometry stream.
     *
     * This stream is used for on-demand geometry computations and diagnostics
     * so that they do not advance the simulation random engine.
     */
    static std::uint64_t NextSeed()
    {
      static std::atomic<std::uint64_t> seed{0x6a09e667f3bcc909ULL};
      return MixSeed(seed.fetch_add(0x9e3779b97f4a7c15ULL));
    }

    /**
     * Seeds the thread-local G4QuickRand stream from a 64-bit geometry seed.
     * A zero seed leaves the current quick-random stream unchanged.
     */
    static void SeedQuickRand(std::uint64_t seed)
    {
      if (seed == 0)
      {
        return;
      }

      auto quickSeed = static_cast<std::uint32_t>(MixSeed(seed));
      G4QuickRand((quickSeed == 0) ? 1 : quickSeed);
    }

    /**
     * Returns a deterministic uniform random value in [0,1), updating
     * @p state with MixSeed().
     */
    static G4double UnitRand(std::uint64_t& state)
    {
      state = MixSeed(state);
      return static_cast<G4double>(state >> 11) / 9007199254740992.0;
    }

    /**
     * Executes @p task over [0,@p count) in chunks of @p grain.
     *
     * The task is called as task(begin,end,worker), where begin/end delimit
     * the chunk and worker is the scheduled worker index. The function creates
     * or reuses a tasking pool internally and returns false if parallel
     * execution is not available or not appropriate.
     */
    template <typename F>
    static G4bool ParallelFor(std::size_t count, std::size_t grain,
                              F&& task, G4bool allowWithinTask = false)
    {
#ifdef G4MULTITHREADED
      if (count == 0 || grain == 0)
      {
        return false;
      }

      auto chunkCount = (count + grain - 1) / grain;
      PoolContext context(chunkCount, allowWithinTask);
      if (context.Size() < 2)
      {
        return false;
      }
      return ParallelFor(context, count, grain, std::forward<F>(task));
#else
      (void)count;
      (void)grain;
      (void)task;
      (void)allowWithinTask;
      return false;
#endif
    }

    /**
     * Executes chunked work with a per-worker context.
     *
     * The @p workerContext callable is invoked once on each scheduled worker
     * as workerContext(worker,run), and must call run() to process that
     * worker's chunk stream. This allows callers to keep task-local guards,
     * such as geometry workspaces, alive for the whole worker task.
     */
    template <typename F, typename W>
    static G4bool ParallelForWithWorkerContext(std::size_t count,
                                               std::size_t grain,
                                               F&& task,
                                               W&& workerContext,
                                               G4bool allowWithinTask = false)
    {
#ifdef G4MULTITHREADED
      if (count == 0 || grain == 0)
      {
        return false;
      }

      auto chunkCount = (count + grain - 1) / grain;
      PoolContext context(chunkCount, allowWithinTask);
      if (context.Size() < 2)
      {
        return false;
      }
      return ParallelForWithWorkerContext(context, count, grain,
                                          std::forward<F>(task),
                                          std::forward<W>(workerContext));
#else
      (void)count;
      (void)grain;
      (void)task;
      (void)workerContext;
      (void)allowWithinTask;
      return false;
#endif
    }

    /**
     * Executes chunked work using an existing PoolContext.
     *
     * This overload is useful when the caller needs to know whether a pool is
     * available before preparing the work, for example to select deterministic
     * seeds for parallel execution.
     */
    template <typename F>
    static G4bool ParallelFor(PoolContext& context, std::size_t count,
                              std::size_t grain, F&& task)
    {
      return ParallelForWithWorkerContext(
        context, count, grain, std::forward<F>(task),
        [](std::size_t, const auto& run) { run(); });
    }

    /**
     * Executes chunked work with a per-worker context using an existing
     * PoolContext.
     */
    template <typename F, typename W>
    static G4bool ParallelForWithWorkerContext(PoolContext& context,
                                               std::size_t count,
                                               std::size_t grain,
                                               F&& task,
                                               W&& workerContext)
    {
#ifdef G4MULTITHREADED
      if (count == 0 || grain == 0 || !context.IsValid())
      {
        return false;
      }

      auto chunkCount = (count + grain - 1) / grain;
      std::atomic<std::size_t> next = 0;
      auto workerCount = std::min<std::size_t>(chunkCount, context.Size());
      ScopedTaskThreadData threadData(context.GetPool());
      G4TaskGroup<void> tasks(context.GetPool());
      for (std::size_t worker = 0; worker < workerCount; ++worker)
      {
        tasks.exec([&, worker]()
        {
          workerContext(worker, [&]()
          {
            for (std::size_t begin = next.fetch_add(grain);
                             begin < count; begin = next.fetch_add(grain))
            {
              task(begin, std::min(begin + grain, count), worker);
            }
          });
        });
      }
      tasks.wait();

      return true;
#else
      (void)context;
      (void)count;
      (void)grain;
      (void)task;
      (void)workerContext;
      return false;
#endif
    }

};

#endif
