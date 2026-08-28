// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_srep
#define tools_srep

#include <string>
#include <utility>
#include <vector>

namespace tools
{

inline bool replace(std::string& a_string, const std::string& a_old, const std::string& a_new)
{
  // return true : some replacement done.
  // return false : nothing replaced.
  if (a_old.empty()) return false;
  bool status = false;
  size_t pos = 0;
  while ((pos = a_string.find(a_old, pos)) != std::string::npos)
  {
    a_string.replace(pos, a_old.length(), a_new);
    pos += a_new.length();
    status = true;
  }
  return status;
}

inline bool replace(std::vector<std::string>& a_strings, const std::string& a_old,
                    const std::string& a_new)
{
  for (auto& it : a_strings)
  {
    if (!replace(it, a_old, a_new)) return false;
  }
  return true;
}

inline std::string to_xml(const std::string& a_string)
{
  std::string _s = a_string;
  // > : &lt;
  // < : &gt;
  // & : &amp;
  // " : &quot;
  // ' : &apos;
  replace(_s, "&", "&amp;");  // must be first.
  replace(_s, "<", "&lt;");
  replace(_s, ">", "&gt;");
  replace(_s, "\"", "&quot;");
  replace(_s, "'", "&apos;");
  return _s;
}

}  // namespace tools

#endif
