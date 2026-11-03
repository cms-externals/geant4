// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_rcmp
#define tools_rcmp

// used in safe cast.
#include <algorithm>
#include <string_view>

namespace tools
{

inline bool rcmp(std::string_view a_1, std::string_view a_2)
{
  return (a_1.size() == a_2.size()) && std::equal(a_1.rbegin(), a_1.rend(), a_2.rbegin());
}

}  // namespace tools

#endif
