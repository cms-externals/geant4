// Copyright (C) 2010 Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_vec
#define tools_vec

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <type_traits>
#include <vector>

namespace tools
{

template<class T, std::size_t N>
class vec
{
    static_assert(N >= 2 && N <= 4, "vec dimension must be 2, 3, or 4");

  public:

    using size_type = std::size_t;
    using value_type = T;

  public:

    vec() = default;

    // Variadic constructor for exactly N arguments
    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr vec(Args&&... args) : m_data{{static_cast<T>(args)...}}
    {}

    ~vec() = default;
    vec(const vec&) = default;
    vec& operator=(const vec&) = default;

    // Component accessors - always available for v0, v1
    const T& v0() const { return m_data[0]; }
    const T& v1() const { return m_data[1]; }

    void v0(const T& a_value) { m_data[0] = a_value; }
    void v1(const T& a_value) { m_data[1] = a_value; }

    // v2 - only for N >= 3
    template<std::size_t M = N, typename = std::enable_if_t<(M >= 3)>>
    const T& v2() const
    {
      return m_data[2];
    }

    template<std::size_t M = N, typename = std::enable_if_t<(M >= 3)>>
    void v2(const T& a_value)
    {
      m_data[2] = a_value;
    }

    // v3 - only for N == 4
    template<std::size_t M = N, typename = std::enable_if_t<(M >= 4)>>
    const T& v3() const
    {
      return m_data[3];
    }

    template<std::size_t M = N, typename = std::enable_if_t<(M >= 4)>>
    void v3(const T& a_value)
    {
      m_data[3] = a_value;
    }

    // Legacy x, y, z accessors
    const T& x() const { return v0(); }
    const T& y() const { return v1(); }

    template<std::size_t M = N, typename = std::enable_if_t<(M >= 3)>>
    const T& z() const
    {
      return v2();
    }

    // Set value - variadic template
    template<typename... Args, typename = std::enable_if_t<sizeof...(Args) == N>>
    void set_value(Args&&... args)
    {
      m_data = {{static_cast<T>(args)...}};
    }

    // Mathematical operations
    T length() const
    {
      T sum = T();
      for (std::size_t i = 0; i < N; ++i)
      {
        sum += m_data[i] * m_data[i];
      }
      return std::sqrt(sum);
    }

    T normalize()
    {
      T norme = length();
      if (norme != T())
      {
        divide(norme);
      }
      return norme;
    }

    T dot(const vec& aV) const
    {
      T sum = T();
      for (std::size_t i = 0; i < N; ++i)
      {
        sum += m_data[i] * aV.m_data[i];
      }
      return sum;
    }

    // Cross product for vec3 - returns via reference
    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    void cross(const vec& aV, vec& a_value) const
    {
      a_value.m_data[0] = m_data[1] * aV.m_data[2] - m_data[2] * aV.m_data[1];
      a_value.m_data[1] = m_data[2] * aV.m_data[0] - m_data[0] * aV.m_data[2];
      a_value.m_data[2] = m_data[0] * aV.m_data[1] - m_data[1] * aV.m_data[0];
    }

    // 2D cross product for vec2 - returns scalar
    template<std::size_t M = N, typename = std::enable_if_t<M == 2>>
    T cross(const vec& aV) const
    {
      return m_data[0] * aV.m_data[1] - m_data[1] * aV.m_data[0];
    }

    bool divide(const T& a_T)
    {
      if (a_T == T()) return false;
      for (auto& val : m_data)
      {
        val /= a_T;
      }
      return true;
    }

    void multiply(const T& a_T)
    {
      for (auto& val : m_data)
      {
        val *= a_T;
      }
    }

    void add(const vec& a_v)
    {
      for (std::size_t i = 0; i < N; ++i)
      {
        m_data[i] += a_v.m_data[i];
      }
    }

    // Only vec2 and vec3 have subtract in original
    template<std::size_t M = N, typename = std::enable_if_t<(M <= 3)>>
    void subtract(const vec& a_v)
    {
      for (std::size_t i = 0; i < N; ++i)
      {
        m_data[i] -= a_v.m_data[i];
      }
    }

    // equals with epsilon - only vec3 in original
    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    bool equals(const vec& a_v, const T a_epsil) const
    {
      T dist_sq = T();
      for (std::size_t i = 0; i < N; ++i)
      {
        T d = m_data[i] - a_v.m_data[i];
        dist_sq += d * d;
      }
      return dist_sq <= a_epsil;
    }

    // Operators
    T& operator[](std::size_t a_index) { return m_data[a_index]; }

    const T& operator[](std::size_t a_index) const { return m_data[a_index]; }

    vec operator+(const vec& a_v) const
    {
      vec result;
      for (std::size_t i = 0; i < N; ++i)
      {
        result.m_data[i] = m_data[i] + a_v.m_data[i];
      }
      return result;
    }

    vec operator-(const vec& a_v) const
    {
      vec result;
      for (std::size_t i = 0; i < N; ++i)
      {
        result.m_data[i] = m_data[i] - a_v.m_data[i];
      }
      return result;
    }

    // Unary minus - only vec3 in original
    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    vec operator-() const
    {
      vec result;
      for (std::size_t i = 0; i < N; ++i)
      {
        result.m_data[i] = -m_data[i];
      }
      return result;
    }

    vec operator*(const T& a_v) const
    {
      vec result;
      for (std::size_t i = 0; i < N; ++i)
      {
        result.m_data[i] = m_data[i] * a_v;
      }
      return result;
    }

    // Division operator - only vec3 in original
    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    vec operator/(const T& a_v) const
    {
      if (a_v == T()) return vec();
      vec result;
      for (std::size_t i = 0; i < N; ++i)
      {
        result.m_data[i] = m_data[i] / a_v;
      }
      return result;
    }

    vec& operator+=(const vec& a_v)
    {
      add(a_v);
      return *this;
    }

    // operator*= - only vec3 in original
    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    vec& operator*=(const T& a_v)
    {
      multiply(a_v);
      return *this;
    }

    bool operator==(const vec& a_v) const { return m_data == a_v.m_data; }

    bool operator!=(const vec& a_v) const { return !operator==(a_v); }

    // Utility
    size_type size() const { return N; }
    const T* data() const { return m_data.data(); }
    constexpr size_type dimension() const { return N; }

    // Static factory methods for vec3
    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    static const vec& s_x()
    {
      static const vec v(1, 0, 0);
      return v;
    }

    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    static const vec& s_y()
    {
      static const vec v(0, 1, 0);
      return v;
    }

    template<std::size_t M = N, typename = std::enable_if_t<M == 3>>
    static const vec& s_z()
    {
      static const vec v(0, 0, 1);
      return v;
    }

  private:

    std::array<T, N> m_data{};
};

// Non-member operator for scalar * vec (commutative multiplication)
template<class T, std::size_t N>
inline vec<T, N> operator*(const T& scalar, const vec<T, N>& v)
{
  return v * scalar;
}

// partial specialization so that sf<float>*vec<float,N> compiles
template<std::size_t N>
inline vec<float, N> operator*(const float& scalar, const vec<float, N>& v)
{
  return v * scalar;
}

// Type aliases for backward compatibility
template<class T>
using vec2 = vec<T, 2>;
template<class T>
using vec3 = vec<T, 3>;
template<class T>
using vec4 = vec<T, 4>;

using vec2f = vec2<float>;
inline const std::string& stype(const vec2f&)
{
  static const std::string s_v("tools::vec2f");
  return s_v;
}

using vec3f = vec3<float>;
inline const std::string& stype(const vec3f&)
{
  static const std::string s_v("tools::vec3f");
  return s_v;
}

using vec4f = vec4<float>;
inline const std::string& stype(const vec4f&)
{
  static const std::string s_v("tools::vec4f");
  return s_v;
}

// Helper functions
template<class T, std::size_t N>
inline const T* get_data(const vec<T, N>& a_v)
{
  return a_v.data();
}

template<class T, std::size_t N>
inline bool set_from_vec(vec<T, N>& a_v, const std::vector<T>& a_sv)
{
  if (a_sv.size() != N) return false;
  for (std::size_t i = 0; i < N; ++i)
  {
    a_v[i] = a_sv[i];
  }
  return true;
}

// Stream output
template<class T, std::size_t N>
inline std::ostream& operator<<(std::ostream& a_out, const vec<T, N>& a_this)
{
  a_out << "x = " << a_this.v0() << ",y = " << a_this.v1();
  if constexpr (N >= 3) a_out << ",z = " << a_this.v2();
  if constexpr (N >= 4) a_out << ",t = " << a_this.v3();
  return a_out;
}

// vec3-specific utility functions
template<class T>
inline vec3<T> direction(const vec3<T>& a_p0, const vec3<T>& a_p1, const vec3<T>& a_p2)
{
  // Orientation is computed by taking (p1 - p0) x (p2 - p0)
  vec3<T> P = a_p1 - a_p0;
  vec3<T> P2 = a_p2 - a_p0;
  vec3<T> result;
  P.cross(P2, result);
  return result;
}

// vec2-specific geometric functions
template<class T>
inline double is_left(const vec2<T>& P0, const vec2<T>& P1, const vec2<T>& P2)
{
  return ((P1.v0() - P0.v0()) * (P2.v1() - P0.v1()) - (P2.v0() - P0.v0()) * (P1.v1() - P0.v1()));
}

template<class T>
inline bool is_inside(const vec2<T>& a_P, const std::vector<vec2<T>>& a_V)
{
  std::size_t n = a_V.size() - 1;
  int wn = 0;

  for (std::size_t i = 0; i < n; i++)
  {
    if (a_V[i].v1() <= a_P.v1())
    {
      if (a_V[i + 1].v1() > a_P.v1())
        if (is_left(a_V[i], a_V[i + 1], a_P) > 0) ++wn;
    }
    else
    {
      if (a_V[i + 1].v1() <= a_P.v1())
        if (is_left(a_V[i], a_V[i + 1], a_P) < 0) --wn;
    }
  }

  return (wn != 0);
}

template<class T>
inline bool intersect(const vec2<T>& P1, const vec2<T>& Q1, const vec2<T>& P2, const vec2<T>& Q2,
                      vec2<T>& a_out)
{
  vec2<T> A = Q1 - P1;
  vec2<T> B = P2 - Q2;
  vec2<T> C = P2 - P1;

  T a = A.v0();
  T b = B.v0();
  T c = C.v0();
  T d = A.v1();
  T e = B.v1();
  T f = C.v1();

  T det = a * e - d * b;
  if (det == T()) return false;

  T r = (c * e - f * b) / det;

  a_out = P1 + A * r;
  return true;
}

}  // namespace tools

#endif
