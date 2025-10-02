// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_sys_mutex
#define tools_sys_mutex

#ifdef _WIN32
#  include <windows.h>
#else
#  include <pthread.h>
#endif

namespace tools
{

class mutex
{
  public:
#ifdef _WIN32
    bool lock()
    {
      if (!m_imp) return false;
      CRITICAL_SECTION* cs = (CRITICAL_SECTION*)m_imp;
      ::EnterCriticalSection(cs);
      return true;
    }

    bool unlock()
    {
      if (!m_imp) return false;
      CRITICAL_SECTION* cs = (CRITICAL_SECTION*)m_imp;
      ::LeaveCriticalSection(cs);
      return true;
    }

    bool try_lock()
    {
      if (!m_imp) return false;
      CRITICAL_SECTION* cs = (CRITICAL_SECTION*)m_imp;
      if (::TryEnterCriticalSection(cs)) return true;
      return false;
    }
#else

    bool lock()
    {
      if (!m_imp) return false;
      pthread_mutex_t* _mutex = (pthread_mutex_t*)m_imp;
      int status = ::pthread_mutex_lock(_mutex);
      return status ? false : true;
    }

    bool unlock()
    {
      if (!m_imp) return false;
      pthread_mutex_t* _mutex = (pthread_mutex_t*)m_imp;
      int status = ::pthread_mutex_unlock(_mutex);
      return status ? false : true;
    }

    bool try_lock()
    {
      if (!m_imp) return false;
      pthread_mutex_t* _mutex = (pthread_mutex_t*)m_imp;
      int status = ::pthread_mutex_trylock(_mutex);
      return status ? false : true;
    }

#endif

  public:
    mutex() : m_imp(0) {}
#ifdef _WIN32
    virtual ~mutex()
    {
      if (m_imp) {
        CRITICAL_SECTION* cs = (CRITICAL_SECTION*)m_imp;
        ::DeleteCriticalSection(cs);
        delete cs;
      }
    }
#else
    virtual ~mutex()
    {
      if (m_imp) {
        pthread_mutex_t* _mutex = (pthread_mutex_t*)m_imp;
        ::pthread_mutex_destroy(_mutex);
        delete _mutex;
      }
    }
#endif

  public:
#ifdef _WIN32
    bool initialize()
    {
      if (m_imp) return true;  // done.
      CRITICAL_SECTION* cs = new CRITICAL_SECTION;
      ::InitializeCriticalSection(cs);
      m_imp = cs;
      return true;
    }
#else
    bool initialize()
    {
      if (m_imp) return true;  // done.
      pthread_mutex_t* pm = new pthread_mutex_t;
      int status = ::pthread_mutex_init(pm, 0);
      if (status) {
        delete pm;
        return false;
      }
      m_imp = pm;
      return true;
    }
#endif

  protected:
    mutex(const mutex&) : m_imp(0) {}
    mutex& operator=(const mutex&) { return *this; }

  protected:
    void* m_imp;
};

}  // namespace tools

#endif
