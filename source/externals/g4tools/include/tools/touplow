// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_touplow
#define tools_touplow

#include <algorithm>
#include <string>

namespace tools
{

inline void tolowercase(std::string& a_string)
{
  // Though we could use std::tolower, we need to understand use cases
  // better before migrating
  std::transform(a_string.begin(), a_string.end(), a_string.begin(),
                 [](unsigned char c) { return (c >= 'A' && c <= 'Z') ? c - 'A' + 'a' : c; });
}

inline void touppercase(std::string& a_string)
{
  // Though we could use std::toupper, we need to understand use cases
  // better before migrating
  std::transform(a_string.begin(), a_string.end(), a_string.begin(),
                 [](unsigned char c) { return (c >= 'a' && c <= 'z') ? c - 'a' + 'A' : c; });
}

}  // namespace tools

#endif
