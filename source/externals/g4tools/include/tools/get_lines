// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_get_lines
#define tools_get_lines

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

namespace tools
{

inline void get_lines(const std::string& a_string, std::vector<std::string>& a_lines)
{
  // a_string is a list separated by "\n" or "\\n".
  // For "xxx\n\nxxx", {"xxx","","xxx"} will be created.
  // WARNING : if a_string is a Windows file name, it may
  //  contains a \n which is not a delimiter ; like ..\data\ntuples.hbook.
  a_lines.clear();
  std::string_view sv(a_string);
  size_t pos = 0;
  while (pos < sv.size())
  {
    size_t next = sv.find('\n', pos);
    size_t esc = sv.find("\\n", pos);
    // Use explicit template parameter to workaround Windows macros
    // See, e.g.:
    // https://stackoverflow.com/questions/5004858/why-is-stdmin-failing-when-windows-h-is-included
    size_t delim = std::min<size_t>(next, esc);

    if (delim == std::string_view::npos)
    {
      a_lines.emplace_back(sv.substr(pos));
      break;
    }

    if (delim == esc)
    {
      a_lines.emplace_back(sv.substr(pos, esc - pos));
      pos = esc + 2;
    }
    else
    {
      a_lines.emplace_back(sv.substr(pos, next - pos));
      pos = next + 1;
    }

    // If delimiter is at the end, add an empty string
    if (pos == sv.size())
    {
      a_lines.emplace_back("");
      break;
    }
  }
}

}  // namespace tools

#endif
