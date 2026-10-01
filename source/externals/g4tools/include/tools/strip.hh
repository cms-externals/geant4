// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_strip
#define tools_strip

#include <string>
#include <vector>

namespace tools
{

enum what
{
  leading,
  trailing,
  both
};

inline void leading_strip(std::string& str, char ch)
{
  auto startIndex = str.find_first_not_of(ch);
  str.erase(0, startIndex);
}

inline void trailing_strip(std::string& str, char ch)
{
  auto endIndex = str.find_last_not_of(ch);
  if (endIndex == std::string::npos)
  {
    str = "";
  }
  else
  {
    str.erase(endIndex + 1);
  }
}

inline bool strip(std::string& a_string, what a_type = both, char a_char = ' ')
{
  // return true = some stripping had been done.
  if (a_string.empty())
  {
    return false;
  }

  std::string original = a_string;

  switch (a_type)
  {
    case leading:
      leading_strip(a_string, a_char);
      break;
    case trailing:
      trailing_strip(a_string, a_char);
      break;
    case both:
      leading_strip(a_string, a_char);
      trailing_strip(a_string, a_char);
      break;
  }
  return a_string != original;
}

}  // namespace tools

#endif
