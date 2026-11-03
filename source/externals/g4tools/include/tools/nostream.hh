// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#ifndef tools_nostream
#define tools_nostream

// dummy ostream class in case we need
// to pass a std::ostream as argument
// but having no concrete implementation
// at hand.

#include <ostream>

namespace tools
{

class nostream final : public std::ostream
{
  public:

    nostream() : std::ostream(nullptr) {}
    ~nostream() override = default;
    nostream(const nostream&) = delete;
    nostream& operator=(const nostream&) = delete;
};

}  // namespace tools

#endif
