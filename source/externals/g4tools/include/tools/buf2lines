// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_buf2lines
#define tools_buf2lines

#include <algorithm>
#include <cstring>
#include <numeric>
#include <string>
#include <vector>

namespace tools
{

inline bool strings2buf(const std::vector<std::string>& a_strings, size_t& a_size, char*& a_buffer)
{
  // For {"aa","bbb"}, it produces "aa0bbb00".
  // For {""}, it produces "00".
  // For {}, it produces "0".
  // Calculate total size: sum of string sizes + 1 per string + 1 final '\0'
  a_size = std::accumulate(a_strings.begin(), a_strings.end(), size_t{1},
                           [](size_t sum, const std::string& s) { return sum + s.size() + 1; });

  a_buffer = new char[a_size];
  if (!a_buffer)
  {
    a_size = 0;
    return false;
  }

  char* pos = a_buffer;
  for (const auto& s : a_strings)
  {
    std::copy_n(s.c_str(), s.size() + 1, pos);  // copy string + null
    pos += s.size() + 1;
  }
  *pos = '\0';  // final null terminator
  return true;
}

inline bool buf2strings(size_t a_size, char* a_buffer, std::vector<std::string>& a_strings)
{
  if (a_size <= 1) return false;
  a_strings.clear();
  char* begin = a_buffer;
  char* end = a_buffer + a_size - 1;  // exclude final '\0'
  while (begin < end)
  {
    size_t len = std::strlen(begin);  // length up to next '\0'
    a_strings.emplace_back(begin, len);
    begin += len + 1;
  }
  return true;
}

}  // namespace tools

#endif
