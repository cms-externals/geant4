// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_schar
#define tools_schar

#include <algorithm>
#include <string>

namespace tools
{

inline bool is_spaces(const std::string& a_s)
{
  return std::all_of(a_s.begin(), a_s.end(), [](char c) { return c == ' '; });
}

}  // namespace tools

#endif
