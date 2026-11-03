//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/MakeCelerOptions.cc
//---------------------------------------------------------------------------//
#include "MakeCelerOptions.hh"

#if defined(WITH_CELERITAS)
#  include <accel/AlongStepFactory.hh>
#  include <accel/SetupOptions.hh>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Build options to set up Celeritas.
 */
SetupOptions MakeCelerOptions()
{
  SetupOptions opts;

  opts.ignore_processes = {"CoulombScat"};

  // Set along-step factory with zero field
  opts.make_along_step = UniformAlongStepFactory();

  //// Configure optical tracking loop for GPU offload of optical photons
  OpticalSetupOptions optical;
  optical.capacity.tracks = 50650;
  const auto num_tracks = optical.capacity.tracks.value();
  optical.capacity.primaries = 8 * num_tracks;
  optical.capacity.generators = 2 * num_tracks;
  optical.generator = inp::OpticalOffloadGenerator{};
  opts.optical = std::move(optical);

  // Clear offload particles: only optical distribution data is pushed
  opts.offload_particles = {};

  // Enable per-action timing for GPU diagnostics
  opts.action_times = true;

  return opts;
}
}  // namespace celeritas
#else
//---------------------------------------------------------------------------//
/*!
 * Stub for CaTSClassesDict (without Celeritas include path).
 * Never called at runtime — CaTS executable provides the real impl.
 */
#  include <cstdint>
#  include <string>
#  include <vector>

namespace celeritas
{
struct SetupOptions
{
    using size_type = unsigned int;
    using VecString = std::vector<std::string>;

    std::string output_file;
    size_type max_num_tracks{};
    size_type initializer_capacity{};
    VecString ignore_processes;
};

//---------------------------------------------------------------------------//
/*!
 * Build options to set up Celeritas.
 */
SetupOptions MakeCelerOptions()
{
  return {};
}
}  // namespace celeritas
#endif
