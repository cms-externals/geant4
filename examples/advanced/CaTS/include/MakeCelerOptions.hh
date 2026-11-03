//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file include/MakeCelerOptions.hh
//---------------------------------------------------------------------------//
#pragma once

#if defined(WITH_CELERITAS)
#  include <accel/SetupOptions.hh>
#else
namespace celeritas
{
struct SetupOptions;
}
#endif

namespace celeritas
{
// Build options to set up Celeritas
SetupOptions MakeCelerOptions();
}  // namespace celeritas
