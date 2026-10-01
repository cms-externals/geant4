// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_sprintf
#define tools_sprintf

#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

namespace tools
{

inline int vsnpf(char* a_s, size_t a_n, const char* a_fmt, va_list args)
{
  return std::vsnprintf(a_s, a_n, a_fmt, args);
}

inline int snpf(char* a_s, size_t a_n, const char* a_fmt, ...)
{
  va_list args;
  va_start(args, a_fmt);
  int n = vsnpf(a_s, a_n, a_fmt, args);
  va_end(args);
  return n;
}

inline bool vsprintf(std::string& a_string, int a_length, const char* a_format, va_list a_args)
{
  a_string.clear();
  if (a_length < 0 || !a_format) return false;

  std::vector<char> buf(static_cast<size_t>(a_length) + 1);
  int n = vsnpf(buf.data(), buf.size(), a_format, a_args);

  if (n < 0 || n > a_length) return false;
  a_string.assign(buf.data(), static_cast<size_t>(n));
  return true;
}

inline bool sprintf(std::string& a_string, int a_length, const char* a_format, ...)
{
  va_list args;
  va_start(args, a_format);
  bool result = vsprintf(a_string, a_length, a_format, args);
  va_end(args);
  return result;
}

}  // namespace tools

#endif
