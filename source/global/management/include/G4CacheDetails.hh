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
// G4CacheDetails
//
// Class description:
//
//    The classes contained in this header files are used by
//    G4Cache to store a TLS instance of the cached object.
//    These classes should not be used by client code, but
//    are used by one of the G4Cache classes.
//
//    G4Cache is a container of the cached value.
//    Not memory efficient, but CPU efficient (constant time access).
//    A different version with a map instead of a vector should be
//    memory efficient and less CPU efficient (log-time access).
//    If really a lot of these objects are used
//    we may want to consider the map version to save some memory.
//
//    These are simplified "split-classes" without any
//    copy-from-master logic. Each cached object is associated
//    a unique identified (an integer), that references an instance
//    of the cached value (of template type VALTYPE) in a
//    static TLS data structure.
//
//    In case the cache is used for a cached object the object class
//    has to provide a default constructor. Alternatively pointers to
//    objects can be stored in the cache and this limitation is removed
//    but explicit handling of memory (new/delete) of cached object becomes
//    client responsibility.
//
// Todos: - Understand if map based class can be more efficent than
//          vector one.
//        - Evaluate use of specialized allocator for TLS "new".

// Author: A.Dotti, 21 October 2013 - First implementation
// --------------------------------------------------------------------
#ifndef G4CACHEDETAILS_HH
#define G4CACHEDETAILS_HH

#include "G4Threading.hh"
#include "globals.hh"

#include <memory>
#include <vector>

// A TLS storage for a cache of type T
//
template<class T>
class G4CacheReference
{
  public:

    /**
     * @brief Initialize TLS storage
     * @param id Instance identifier
     */
    inline void Initialize(unsigned int id);

    /**
     * @brief Cleanup TLS storage for instance id. If last==true destroy and cleanup object
     * container
     * @param id Instance identifier
     * @param last If true, destroy and cleanup object container
     */
    inline void Destroy(unsigned int id, G4bool last);

    /**
     * @brief Returns cached value for instance id
     * @param id Instance identifier
     * @return Reference to cached value
     */
    inline T& GetCache(unsigned int id) const;

  private:

    // Implementation detail: the cached object is stored as a
    // pointer. Object is stored as a pointer to avoid too large
    // std::vector in case of stored objects and allow use of
    // specialized allocators
    using cache_container = std::vector<std::unique_ptr<T>>;
    static cache_container*& cache();
};

// Template specialization for pointers
// Note: Objects are not owned by cache, for this version of the cache
//       the explicit new/delete of the cached object must be handled
//       by the client.
//
template<class T>
class G4CacheReference<T*>
{
  public:

    /**
     * @brief Initliaze TLS storage
     * @param id Instance identifier
     */
    inline void Initialize(unsigned int id);

    /**
     * @brief Cleanup TLS storage for instance id. If last==true destroy and cleanup object
     * container
     * @param id Instance identifier
     * @param last If true, destroy and cleanup object container
     */
    inline void Destroy(unsigned int id, G4bool last);

    /**
     * @brief Returns cached pointer for instance id
     * @param id Instance identifier
     * @return Reference to cached pointer
     */
    inline T*& GetCache(unsigned int id) const;

  private:

    using cache_container = std::vector<T*>;
    static cache_container*& cache();
};

// Template specialization for probably the most used case: double
// Be more efficient avoiding unnecessary "new/delete"
//
template<>
class G4CacheReference<G4double>
{
  public:

    /**
     * @brief Initliaze TLS storage
     * @param id Instance identifier
     */
    inline void Initialize(unsigned int id);

    /**
     * @brief Cleanup TLS storage for instance id. If last==true destroy and cleanup object
     * container
     * @param id Instance identifier
     * @param last If true, destroy and cleanup object container
     */
    inline void Destroy(unsigned int id, G4bool last);

    /**
     * @brief Returns cached double for instance id
     * @param id Instance identifier
     * @return Reference to cached double
     */
    inline G4double& GetCache(unsigned int id) const;

  private:

    using cache_container = std::vector<G4double>;
    static G4GLOB_DLL cache_container*& cache();
};

//======= Implementation: G4CacheReference<T>
//===========================================

template<class T>
void G4CacheReference<T>::Initialize(unsigned int id)
{
  // Create cache container
  if (cache() == nullptr)
  {
    cache() = new cache_container;
  }
  if (cache()->size() <= id)
  {
    cache()->resize(id + 1);
  }
  if (!(*cache())[id])
  {
    (*cache())[id] = std::make_unique<T>();
  }
}

template<class T>
void G4CacheReference<T>::Destroy(unsigned int id, G4bool last)
{
  // We DO NOT raise an exception for calling Destroy on an empty cache
  // - We don't do anything as this is not neccesarily indicative of a memory leak
  // - IF client code has used "static G4Cache<T>" then that may be destructed
  //   AFTER the static backing store is. This is nominally undefined behaviour, though
  //   access should not cause actual errors unless the value is relied on after
  //   program exit (which would be extraordinarily weird).
  if (cache() != nullptr)
  {
    if (cache()->size() < id)
    {
      G4ExceptionDescription msg;
      msg << "Internal fatal error. Invalid G4Cache size (requested id: " << id
          << " but cache has size: " << cache()->size();
      msg << " Possibly client created G4Cache object in a thread and"
          << " tried to delete it from another thread!";
      G4Exception("G4CacheReference<T>::Destroy", "Cache001", FatalException, msg);
      return;
    }
    if (cache()->size() > id)
    {
      ((*cache())[id]).reset();
    }
    if (last)
    {
      delete cache();
      cache() = nullptr;
    }
  }
}

template<class T>
T& G4CacheReference<T>::GetCache(unsigned int id) const
{
  return *((*cache())[id]);
}

template<class T>
typename G4CacheReference<T>::cache_container*& G4CacheReference<T>::cache()
{
  G4ThreadLocalStatic cache_container* _instance = nullptr;
  return _instance;
}

//======= Implementation: G4CacheReference<T*>
//============================================

template<class T>
void G4CacheReference<T*>::Initialize(unsigned int id)
{
  if (cache() == nullptr)
  {
    cache() = new cache_container;
  }
  if (cache()->size() <= id)
  {
    cache()->resize(id + 1, nullptr);
  }
}

template<class T>
inline void G4CacheReference<T*>::Destroy(unsigned int id, G4bool last)
{
  // We DO NOT raise an exception for calling Destroy on an empty cache
  // See explanation in G4CacheReference<T>::Destroy
  if (cache() != nullptr)
  {
    if (cache()->size() < id)
    {
      G4ExceptionDescription msg;
      msg << "Internal fatal error. Invalid G4Cache size (requested id: " << id
          << " but cache has size: " << cache()->size();
      msg << " Possibly client created G4Cache object in a thread and"
          << " tried to delete it from another thread!";
      G4Exception("G4CacheReference<T*>::Destroy", "Cache001", FatalException, msg);
      return;
    }
    if (cache()->size() > id)
    {
      // Client **must** delete the object themselves. We just null the one we have
      (*cache())[id] = nullptr;
    }
    if (last)
    {
      delete cache();
      cache() = nullptr;
    }
  }
}

template<class T>
T*& G4CacheReference<T*>::GetCache(unsigned int id) const
{
  return (*cache())[id];
}

template<class T>
typename G4CacheReference<T*>::cache_container*& G4CacheReference<T*>::cache()
{
  G4ThreadLocalStatic cache_container* _instance = nullptr;
  return _instance;
}

//======= Implementation: G4CacheReference<double>
//============================================

void G4CacheReference<G4double>::Initialize(unsigned int id)
{
  if (cache() == nullptr)
  {
    cache() = new cache_container;
  }
  if (cache()->size() <= id)
  {
    cache()->resize(id + 1, 0.0);
  }
}

void G4CacheReference<G4double>::Destroy(unsigned int /*id*/, G4bool last)
{
  if (last)
  {
    delete cache();
    cache() = nullptr;
  }
}

G4double& G4CacheReference<G4double>::GetCache(unsigned int id) const
{
  return (*cache())[id];
}

#endif
