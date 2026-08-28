#include <string>
#include <vector>

#ifndef G4VAL_SPLITSTRING_C
#define G4VAL_SPLITSTRING_C

void SplitString( const std::string& str, const std::string del, std::vector<std::string>& tokens )
{
         
   string::size_type last = str.find_first_not_of(del, 0);
   string::size_type pos  = str.find_first_of(del, last);
   while (std::string::npos != pos || std::string::npos != last)
   {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(last, pos-last));
      // Skip delimiters.  Note the "not_of"
      last = str.find_first_not_of(del, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(del, last);
   }

   return;

}

#endif
