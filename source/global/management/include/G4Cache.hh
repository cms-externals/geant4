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
// G4Cache
//
// Class Description:
//
//      Helper classes for Geant4 Multi-Threaded.
//      The classes defined in this header file provide a thread-private
//      cache to store a thread-local variable V in a class instance
//      shared among threads.
//      These are templated classes on the to-be-stored object.
//
// Example:
//      Let's assume an instance myObject of class G4Shared is sharead between
//      threads. Still a data member of this class needs to be thread-private.
//      A typical example of this being a "cache" for a local calculation.
//      The helper defined here can be used to guarantee thread-safe operations
//      on the thread-private object.
//      Example:
//         class G4Shared
//         {
//           G4double sharedData;
//           G4Cache<G4double> threadPrivate;
//           void foo()
//           {
//             G4double priv = threadPrivate.Get();
//             if ( priv < 10 ) priv += sharedData;
//             threadPrivate.Put( priv );
//           }
//         }
//
//      Two variants of the base G4Cache exist. The first one being
//      G4VectorCache similar to std::vector.
//      Example:
//         G4VectorCache<G4double> aVect;
//         aVect.Push_back( 3.2 );
//         aVect.Push_back( 4.1 );
//         std::cout << aVect[0] << std::endl;
//      The second one being:
//      G4MapCache, similar to std::map.
//      Example:
//         G4MapCache<G4int, G4double> aMap;
//         aMap[320]=1.234;
//
//      See classes definition for details.

// Author: A.Dotti, 21 October 2013 - First implementation
// --------------------------------------------------------------------
#ifndef G4CACHE_HH
#define G4CACHE_HH

// Debug this code
// #define g4cdebug 1

#include "G4AutoLock.hh"
#include "G4CacheDetails.hh"  // Thread Local storage details are here

#include <atomic>
#include <map>
#include <system_error>

/**
 * @brief A templated cache to store thread-private data of type VALTYPE
 */
template<class T>
class G4Cache
{
  public:

    using value_type = T;

  public:

    /**
     * @brief Default constructor.
     */
    G4Cache();

    /**
     * @brief Construct cache object with initial value.
     * @param v Initial value.
     */
    G4Cache(const value_type& v);

    /**
     * @brief Default destructor.
     */
    virtual ~G4Cache();

    /**
     * @brief Copy constructor.
     *
     * Deep-copy - Copies the content of the cache
     * @param rhs Source cache.
     */
    G4Cache(const G4Cache& rhs);

    /**
     * @brief Assignment operator.
     * @param rhs Source cache.
     * @return Reference to this cache.
     */
    G4Cache& operator=(const G4Cache& rhs);

    /**
     * @brief Gets reference to cached value of this thread.
     * @return Reference to cached value.
     */
    inline value_type& Get() const;

    /**
     * @brief Sets this thread's cached value.
     * @param val Value to set.
     */
    inline void Put(const value_type& val) const;

    /**
     * @brief Gets copy of cached value.
     * @return Copy of cached value.
     */
    inline value_type Pop();

  protected:

    /**
     * @brief Gets cache ID.
     * @return Cache ID.
     */
    const G4int& GetId() const { return id; }

  private:

    G4int id;
    mutable G4CacheReference<value_type> theCache;
    inline static std::atomic<unsigned int> instancesctr{0};
    inline static std::atomic<unsigned int> dstrctr{0};

    /**
     * @brief Gets thread-local cache reference.
     * @return Reference to thread-local cache.
     */
    inline value_type& GetCache() const
    {
      theCache.Initialize(id);
      return theCache.GetCache(id);
    }
};

/**
 * @class G4VectorCache
 * @brief Vector version of thread-private cache.
 * @tparam T Type of vector element.
 */
template<class T>
class G4VectorCache : public G4Cache<std::vector<T>>
{
  public:

    using value_type = T;
    using vector_type = typename std::vector<value_type>;
    using size_type = typename vector_type::size_type;
    using iterator = typename vector_type::iterator;
    using const_iterator = typename vector_type::const_iterator;

  public:

    /**
     * @brief Default constructor.
     */
    G4VectorCache() = default;

    /**
     * @brief Creates a vector cache of nElems elements.
     * @param nElems Number of elements.
     */
    G4VectorCache(G4int nElems);

    /**
     * @brief Creates a vector cache with elements from an array.
     * @param nElems Number of elements.
     * @param vals Array of values.
     */
    G4VectorCache(G4int nElems, value_type* vals);

    /**
     * @brief Default destructor.
     */
    virtual ~G4VectorCache() = default;

    /**
     * @brief Adds an element to the vector cache.
     * @param val Value to add.
     */
    inline void Push_back(const value_type& val);

    /**
     * @brief Removes and returns the last element.
     * @return Last element.
     */
    inline value_type Pop_back();

    /**
     * @brief Access element by index.
     * @param idx Index.
     * @return Reference to element.
     */
    inline value_type& operator[](const G4int& idx);

    /**
     * @brief Returns iterator to beginning.
     * @return Iterator.
     */
    inline iterator Begin();

    /**
     * @brief Returns iterator to end.
     * @return Iterator.
     */
    inline iterator End();

    /**
     * @brief Clears the vector cache.
     */
    inline void Clear();
    /**
     * @brief Returns size of vector cache.
     * @return Size.
     */
    inline size_type Size() { return G4Cache<vector_type>::Get().size(); }
    //  Needs to be here for a VC9 compilation problem
};

/**
 * @class G4MapCache
 * @brief Map version of thread-private cache.
 *
 * Implements std::map interface so can be used directly as a std::map would be used.
 *
 * @tparam Key Type of key.
 * @tparam T Type of value.
 */
template<class Key, class T>
class G4MapCache : public G4Cache<std::map<Key, T>>
{
  public:

    using key_type = Key;
    using value_type = T;
    using map_type = typename std::map<key_type, value_type>;
    using size_type = typename map_type::size_type;
    using iterator = typename map_type::iterator;
    using const_iterator = typename map_type::const_iterator;

  public:

    /**
     * @brief Default destructor.
     */
    virtual ~G4MapCache() = default;

    /**
     * @brief Returns true if map contains element corresponding to key.
     * @param k Key.
     * @return True if key exists.
     */
    inline G4bool Has(const key_type& k);

    /**
     * @brief Inserts a key-value pair.
     * @param k Key.
     * @param v Value.
     * @return Pair of iterator and bool.
     */
    inline std::pair<iterator, G4bool> Insert(const key_type& k, const value_type& v);

    /**
     * @brief Returns iterator to beginning.
     * @return Iterator.
     */
    inline iterator Begin();

    /**
     * @brief Returns iterator to end.
     * @return Iterator.
     */
    inline iterator End();

    /**
     * @brief Finds element by key.
     * @param k Key.
     * @return Iterator.
     */
    inline iterator Find(const key_type& k);

    /**
     * @brief Gets value by key.
     * @param k Key.
     * @return Reference to value.
     */
    inline value_type& Get(const key_type& k);

    /**
     * @brief Erases element by key.
     * @param k Key.
     * @return Number of elements erased.
     */
    inline size_type Erase(const key_type& k);

    /**
     * @brief Access value by key.
     * @param k Key.
     * @return Reference to value.
     */
    inline value_type& operator[](const key_type& k);
    /**
     * @brief Returns size of map cache.
     * @return Size.
     */
    inline size_type Size() { return G4Cache<map_type>::Get().size(); }
    //  Needs to be here for a VC9 compilation problem
};

/**
 * @class G4Cache
 * @brief Thread-private cache for storing a thread-local variable.
 * @tparam T Type of the cached value.
 */
//========= Implementation: G4Cache<T> ====================================

template<class T>
G4Cache<T>::G4Cache()
{
  G4AutoLock l(G4TypeMutex<G4Cache<T>>());
  id = instancesctr++;
  theCache.Initialize(id);
}

template<class T>
G4Cache<T>::G4Cache(const G4Cache<T>& rhs)
{
  // Copy is special, we need to copy the content
  // of the cache, not the cache object

  if (this == &rhs) return;

  G4AutoLock l(G4TypeMutex<G4Cache<T>>());
  id = instancesctr++;
  T aCopy = rhs.GetCache();
  Put(aCopy);
}

template<class T>
G4Cache<T>& G4Cache<T>::operator=(const G4Cache<T>& rhs)
{
  if (this == &rhs) return *this;

  // Force copy of cached data
  T aCopy = rhs.GetCache();
  Put(aCopy);
  return *this;
}

template<class T>
G4Cache<T>::G4Cache(const T& v)
{
  G4AutoLock l(G4TypeMutex<G4Cache<T>>());
  id = instancesctr++;
  Put(v);
}

template<class T>
G4Cache<T>::~G4Cache()
{
  // don't automatically lock --> wait until we can catch an error
  // without scoping the G4AutoLock
  G4AutoLock l(G4TypeMutex<G4Cache<T>>(), std::defer_lock);

  // sometimes the mutex is unavailable in destructors so
  // try to lock the associated mutex, but catch if it fails
  try
  {
    // a system_error in lock means that the mutex is unavailable
    // we want to throw the error that comes from locking an unavailable
    // mutex so that we know there is a memory leak
    // if the mutex is valid, this will hold until the other thread finishes
    //
    l.lock();
  }
  catch (std::system_error& e)
  {
    // the error that comes from locking an unavailable mutex
#ifdef g4cdebug
    std::cout << "Non-critical error: mutex lock failure in ~G4Cache<" << typeid(V).name() << ">. "
              << std::endl
              << "If the RunManagerKernel has been deleted, it failed to "
              << "delete an allocated resource" << std::endl
              << "and this destructor is being called after the statics "
              << "were destroyed." << std::endl;
    std::cout << "Exception: [code: " << e.code() << "] caught: " << e.what() << std::endl;
#endif
  }
  ++dstrctr;
  G4bool last = (dstrctr == instancesctr);
  theCache.Destroy(id, last);
  if (last)
  {
    instancesctr.store(0);
    dstrctr.store(0);
  }
}

template<class T>
T& G4Cache<T>::Get() const
{
  return GetCache();
}

template<class T>
void G4Cache<T>::Put(const T& val) const
{
  GetCache() = val;
}

// Should here remove from cache element?
template<class T>
T G4Cache<T>::Pop()
{
  return GetCache();
}

//========== Implementation: G4VectorCache<V> ===========================

template<class T>
G4VectorCache<T>::G4VectorCache(G4int nElems)
{
  vector_type& cc = G4Cache<vector_type>::Get();
  cc.resize(nElems);
}

template<class T>
G4VectorCache<T>::G4VectorCache(G4int nElems, T* vals)
{
  vector_type& cc = G4Cache<vector_type>::Get();
  cc.resize(nElems);
  for (G4int idx = 0; idx < nElems; ++idx)
    cc[idx] = vals[idx];
}

template<class T>
void G4VectorCache<T>::Push_back(const T& val)
{
  G4Cache<vector_type>::Get().push_back(val);
}

template<class T>
T G4VectorCache<T>::Pop_back()
{
  vector_type& cc = G4Cache<vector_type>::Get();
  value_type val = cc[cc.size() - 1];
  cc.pop_back();
  return val;
}

template<class T>
T& G4VectorCache<T>::operator[](const G4int& idx)
{
  vector_type& cc = G4Cache<vector_type>::Get();
  return cc[idx];
}

template<class T>
typename G4VectorCache<T>::iterator G4VectorCache<T>::Begin()
{
  return G4Cache<vector_type>::Get().begin();
}

template<class T>
typename G4VectorCache<T>::iterator G4VectorCache<T>::End()
{
  return G4Cache<vector_type>::Get().end();
}

template<class T>
void G4VectorCache<T>::Clear()
{
  G4Cache<vector_type>::Get().clear();
}

// template<class V>
// typename G4VectorCache<V>::size_type G4VectorCache<V>::Size()
//{
//    return G4Cache<vector_type>::Get().size();
//}

//======== Implementation: G4MapType<K,V> ===========================

template<class Key, class T>
std::pair<typename G4MapCache<Key, T>::iterator, G4bool> G4MapCache<Key, T>::Insert(const Key& k,
                                                                                    const T& v)
{
  return G4Cache<map_type>::Get().insert(std::pair<key_type, value_type>(k, v));
}

// template<class K, class V>
// typename G4MapCache<K,V>::size_type G4MapCache<K,V>::Size()
//{
//    return G4Cache<map_type>::Get().size();
//}

template<class Key, class T>
typename G4MapCache<Key, T>::iterator G4MapCache<Key, T>::Begin()
{
  return G4Cache<map_type>::Get().begin();
}

template<class Key, class T>
typename G4MapCache<Key, T>::iterator G4MapCache<Key, T>::End()
{
  return G4Cache<map_type>::Get().end();
}

template<class Key, class T>
typename G4MapCache<Key, T>::iterator G4MapCache<Key, T>::Find(const Key& k)
{
  return G4Cache<map_type>::Get().find(k);
}

template<class Key, class T>
G4bool G4MapCache<Key, T>::Has(const Key& k)
{
  return (Find(k) != End());
}

template<class Key, class T>
T& G4MapCache<Key, T>::Get(const Key& k)
{
  return Find(k)->second;
}

template<class Key, class T>
typename G4MapCache<Key, T>::size_type G4MapCache<Key, T>::Erase(const Key& k)
{
  return G4Cache<map_type>::Get().erase(k);
}

template<class Key, class T>
T& G4MapCache<Key, T>::operator[](const Key& k)
{
  return (G4Cache<map_type>::Get())[k];
}

#endif
