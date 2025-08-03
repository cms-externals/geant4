//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// created by jrmadsen on Sun Oct 21 21:31:53 2018
//
//
//

#include <cstdint>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
#include <limits>
#include <thread>

#include <map>
#include <unordered_map>
#include <vector>
#include <deque>
#include <set>

#include "G4ConvergenceTester.hh"
#include "G4StatAnalysis.hh"
#include "G4THitsMap.hh"
#include "G4THitsVector.hh"
#include "Randomize.hh"
#include "testing.hh"

//----------------------------------------------------------------------------//
// constants
static intmax_t nevent = 1000000;
static intmax_t nvoxel = 100;
static intmax_t nexclude = 5;
static const G4double epsilon = std::numeric_limits<G4float>::epsilon();

//----------------------------------------------------------------------------//
// typedefs
template <typename _Tp, typename _Up>
using map_type = G4VTHitsMap<_Tp, std::map<G4int, _Up>>;

template <typename _Tp, typename _Up>
using mmap_type = G4VTHitsMap<_Tp, std::multimap<G4int, _Up>>;

template <typename _Tp, typename _Up>
using uomap_type = G4VTHitsMap<_Tp, std::unordered_map<G4int, _Up>>;

template <typename _Tp, typename _Up>
using uommap_type = G4VTHitsMap<_Tp, std::unordered_multimap<G4int, _Up>>;

template <typename _Tp, typename _Up>
using vector_type = G4VTHitsVector<_Tp, std::vector<_Up>>;

template <typename _Tp, typename _Up>
using deque_type = G4VTHitsVector<_Tp, std::deque<_Up>>;

//----------------------------------------------------------------------------//

// this is the standard structure used at the event level
static G4THitsMap<G4double> event_map;

//----------------------------------------------------------------------------//

// these are the variants for G4Run
static map_type     <G4double, G4double*> dbl_ptr_run_map;
static mmap_type    <G4double, G4double*> dbl_ptr_run_mmap;
static uomap_type   <G4double, G4double*> dbl_ptr_run_uomap;
static uommap_type  <G4double, G4double*> dbl_ptr_run_uommap;
static deque_type   <G4double, G4double*> dbl_ptr_run_deque;
static vector_type  <G4double, G4double*> dbl_ptr_run_vector;
static deque_type   <G4double, G4double>  dbl_pod_run_deque;
static vector_type  <G4double, G4double>  dbl_pod_run_vector;

static map_type     <G4StatAnalysis, G4StatAnalysis*> stat_ptr_run_map;
static mmap_type    <G4StatAnalysis, G4StatAnalysis*> stat_ptr_run_mmap;
static uomap_type   <G4StatAnalysis, G4StatAnalysis*> stat_ptr_run_uomap;
static uommap_type  <G4StatAnalysis, G4StatAnalysis*> stat_ptr_run_uommap;
static deque_type   <G4StatAnalysis, G4StatAnalysis*> stat_ptr_run_deque;
static vector_type  <G4StatAnalysis, G4StatAnalysis*> stat_ptr_run_vector;
static deque_type   <G4StatAnalysis, G4StatAnalysis>  stat_pod_run_deque;
static vector_type  <G4StatAnalysis, G4StatAnalysis>  stat_pod_run_vector;

//----------------------------------------------------------------------------//

template <typename _Tp, typename _Container>
void add(const G4THitsMap<_Tp>& evt_map, _Container& run_container)
{
    run_container += evt_map;
}

//----------------------------------------------------------------------------//

template <typename _Container>
void merge(const char* type, const _Container& run_local)
{
    _Container run_total;
    run_total += run_local;

    for(auto local_itr = run_local.begin(); local_itr != run_local.end(); ++local_itr)
    {
        auto  local_idx = run_local.GetIndex(local_itr);
        auto* local_obj = run_local.GetObject(local_itr);
        auto* total_obj = run_total.GetObject(local_idx);

        if(local_obj || total_obj)
        {
            // different pointers
            EXPECT_NE(local_obj, total_obj);

            std::cout << "[" << type << "] Testing merge "
                      << std::setw(3) << local_idx << " ... "
                      << std::flush;
            if(!local_obj)
            {
                std::cout << "warning! nullptr to local_obj. This can occur with "
                          << "non-pointer types ... " << std::flush;
                ASSERT_TRUE(total_obj != nullptr);
                G4double _local = 0.0;
                G4double _total = *total_obj;
                ASSERT_CLOSE(_local, _total, epsilon);
            }
            else
            {
                ASSERT_TRUE(local_obj != nullptr);
                ASSERT_TRUE(total_obj != nullptr);
                G4double _local = *local_obj;
                G4double _total = *total_obj;
                ASSERT_CLOSE(_local, _total, epsilon);
            }
        }
        else
        {
            std::cout << "[" << type << "] Testing merge "
                      << std::setw(3) << local_idx << " ... "
                      << "passed.  null pointers" << std::endl;
        }
    }
}

//----------------------------------------------------------------------------//

template <typename _Container>
void verify(const char* type, const _Container& run_container)
{
    for(auto itr = run_container.begin(); itr != run_container.end(); ++itr)
    {
        auto idx = run_container.GetIndex(itr);
        auto* obj = run_container.GetObject(itr);
        if(obj)
        {
            std::cout << "[" << type << "] Testing index "
                      << std::setw(3) << idx << " ... "
                      << std::flush;
            auto* evt = event_map[idx];
            if(!evt)
            {
                std::cout << "warning! nullptr to evt. This can occur with "
                          << "non-pointer types ... " << std::flush;
                ASSERT_TRUE(obj != nullptr);
                G4double _evt = 0.0;
                G4double _obj = *obj;
                ASSERT_CLOSE(_evt, _obj, epsilon);

            }
            else
            {
                ASSERT_TRUE(evt != nullptr);
                ASSERT_TRUE(obj != nullptr);
                G4double _evt = *evt;
                G4double _obj = *obj;
                ASSERT_CLOSE(_evt, _obj, epsilon);
            }
        }
        else
        {
            std::cout << "[" << type << "] Testing index "
                      << std::setw(3) << idx << " ... "
                      << "passed.  null pointers" << std::endl;
        }
    }
}

//----------------------------------------------------------------------------//

void check_map_dbl_ptr()        { verify(__FUNCTION__, dbl_ptr_run_map);       }
void check_mmap_dbl_ptr()       { verify(__FUNCTION__, dbl_ptr_run_mmap);      }
void check_uomap_dbl_ptr()      { verify(__FUNCTION__, dbl_ptr_run_uomap);     }
void check_uommap_dbl_ptr()     { verify(__FUNCTION__, dbl_ptr_run_uommap);    }
void check_deque_dbl_ptr()      { verify(__FUNCTION__, dbl_ptr_run_deque);     }
void check_vector_dbl_ptr()     { verify(__FUNCTION__, dbl_ptr_run_vector);    }
void check_deque_dbl_pod()      { verify(__FUNCTION__, dbl_pod_run_deque);     }
void check_vector_dbl_pod()     { verify(__FUNCTION__, dbl_ptr_run_vector);    }

void check_map_stat_ptr()       { verify(__FUNCTION__, stat_ptr_run_map);      }
void check_mmap_stat_ptr()      { verify(__FUNCTION__, stat_ptr_run_mmap);     }
void check_uomap_stat_ptr()     { verify(__FUNCTION__, stat_ptr_run_uomap);    }
void check_uommap_stat_ptr()    { verify(__FUNCTION__, stat_ptr_run_uommap);   }
void check_deque_stat_ptr()     { verify(__FUNCTION__, stat_ptr_run_deque);    }
void check_vector_stat_ptr()    { verify(__FUNCTION__, stat_ptr_run_vector);   }
void check_deque_stat_pod()     { verify(__FUNCTION__, stat_pod_run_deque);    }
void check_vector_stat_pod()    { verify(__FUNCTION__, stat_ptr_run_vector);   }

//----------------------------------------------------------------------------//

void merge_map_dbl_ptr()        { merge(__FUNCTION__, dbl_ptr_run_map);        }
void merge_mmap_dbl_ptr()       { merge(__FUNCTION__, dbl_ptr_run_mmap);       }
void merge_uomap_dbl_ptr()      { merge(__FUNCTION__, dbl_ptr_run_uomap);      }
void merge_uommap_dbl_ptr()     { merge(__FUNCTION__, dbl_ptr_run_uommap);     }
void merge_deque_dbl_ptr()      { merge(__FUNCTION__, dbl_ptr_run_deque);      }
void merge_vector_dbl_ptr()     { merge(__FUNCTION__, dbl_ptr_run_vector);     }
void merge_deque_dbl_pod()      { merge(__FUNCTION__, dbl_pod_run_deque);      }
void merge_vector_dbl_pod()     { merge(__FUNCTION__, dbl_ptr_run_vector);     }

void merge_map_stat_ptr()       { merge(__FUNCTION__, stat_ptr_run_map);       }
void merge_mmap_stat_ptr()      { merge(__FUNCTION__, stat_ptr_run_mmap);      }
void merge_uomap_stat_ptr()     { merge(__FUNCTION__, stat_ptr_run_uomap);     }
void merge_uommap_stat_ptr()    { merge(__FUNCTION__, stat_ptr_run_uommap);    }
void merge_deque_stat_ptr()     { merge(__FUNCTION__, stat_ptr_run_deque);     }
void merge_vector_stat_ptr()    { merge(__FUNCTION__, stat_ptr_run_vector);    }
void merge_deque_stat_pod()     { merge(__FUNCTION__, stat_pod_run_deque);     }
void merge_vector_stat_pod()    { merge(__FUNCTION__, stat_ptr_run_vector);    }

//============================================================================//

int main(int argc, char** argv)
{
    G4StatAnalysis::ResetCpuClock();

    if(argc > 1)
        nevent = atol(argv[1]);
    if(argc > 2)
        nvoxel = atol(argv[2]);

    std::cout << "\nExecuting " << argv[0] << ". # events = "
              << nevent << ", # voxels = " << nvoxel << "..." << std::endl;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values
    // from the mean
    std::normal_distribution<> score_dist{5,2};
    std::uniform_int_distribution<G4int> voxel_dist(0, nvoxel);

    std::set<intmax_t> exclude_voxels;
    while(exclude_voxels.size() < static_cast<size_t>(nexclude))
    {
        G4int voxel = voxel_dist(gen);
        exclude_voxels.insert(voxel);
    }

    {
        std::stringstream ss;
        ss << "Excluded voxels: ";
        for(auto itr = exclude_voxels.begin(); itr != exclude_voxels.end(); ++itr)
        {
            ss << *itr;
            auto _dist = std::distance(exclude_voxels.begin(), itr);
            if(static_cast<size_t>(_dist+1) != exclude_voxels.size())
                ss << ", ";
        }
        std::cout << "\n" << ss.str() << "\n" << std::endl;
    }

    for(intmax_t i = 0; i < nevent; ++i)
    {
        G4int voxel = voxel_dist(gen);
        // make sure not in "exclude" list
        while(exclude_voxels.count(voxel) != 0)
            voxel = voxel_dist(gen);

        G4double score = score_dist(gen);
        // G4ConvergenceTester does not accept negative values
        if(score < 0.0)
            score -= 2.0*score;
        event_map.add(voxel, score);
    }

    // here we test operator+=
    add(event_map, dbl_ptr_run_map);
    add(event_map, dbl_ptr_run_mmap);
    add(event_map, dbl_ptr_run_uomap);
    add(event_map, dbl_ptr_run_uommap);
    add(event_map, dbl_ptr_run_deque);
    add(event_map, dbl_ptr_run_vector);
    add(event_map, dbl_pod_run_deque);
    add(event_map, dbl_pod_run_vector);

    add(event_map, stat_ptr_run_map);
    add(event_map, stat_ptr_run_mmap);
    add(event_map, stat_ptr_run_uomap);
    add(event_map, stat_ptr_run_uommap);
    add(event_map, stat_ptr_run_deque);
    add(event_map, stat_ptr_run_vector);
    add(event_map, stat_pod_run_deque);
    add(event_map, stat_pod_run_vector);

    int num_test = 0;
    int num_fail = 0;

    // run the tests
    RUN_TEST(check_map_dbl_ptr, num_test, num_fail);
    RUN_TEST(check_mmap_dbl_ptr, num_test, num_fail);
    RUN_TEST(check_uomap_dbl_ptr, num_test, num_fail);
    RUN_TEST(check_uommap_dbl_ptr, num_test, num_fail);
    RUN_TEST(check_deque_dbl_ptr, num_test, num_fail);
    RUN_TEST(check_vector_dbl_ptr, num_test, num_fail);
    RUN_TEST(check_deque_dbl_pod, num_test, num_fail);
    RUN_TEST(check_vector_dbl_pod, num_test, num_fail);

    RUN_TEST(check_map_stat_ptr, num_test, num_fail);
    RUN_TEST(check_mmap_stat_ptr, num_test, num_fail);
    RUN_TEST(check_uomap_stat_ptr, num_test, num_fail);
    RUN_TEST(check_uommap_stat_ptr, num_test, num_fail);
    RUN_TEST(check_deque_stat_ptr, num_test, num_fail);
    RUN_TEST(check_vector_stat_ptr, num_test, num_fail);
    RUN_TEST(check_deque_stat_pod, num_test, num_fail);
    RUN_TEST(check_vector_stat_pod, num_test, num_fail);

    RUN_TEST(merge_map_dbl_ptr, num_test, num_fail);
    RUN_TEST(merge_mmap_dbl_ptr, num_test, num_fail);
    RUN_TEST(merge_uomap_dbl_ptr, num_test, num_fail);
    RUN_TEST(merge_uommap_dbl_ptr, num_test, num_fail);
    RUN_TEST(merge_deque_dbl_ptr, num_test, num_fail);
    RUN_TEST(merge_vector_dbl_ptr, num_test, num_fail);
    RUN_TEST(merge_deque_dbl_pod, num_test, num_fail);
    RUN_TEST(merge_vector_dbl_pod, num_test, num_fail);

    RUN_TEST(merge_map_stat_ptr, num_test, num_fail);
    RUN_TEST(merge_mmap_stat_ptr, num_test, num_fail);
    RUN_TEST(merge_uomap_stat_ptr, num_test, num_fail);
    RUN_TEST(merge_uommap_stat_ptr, num_test, num_fail);
    RUN_TEST(merge_deque_stat_ptr, num_test, num_fail);
    RUN_TEST(merge_vector_stat_ptr, num_test, num_fail);
    RUN_TEST(merge_deque_stat_pod, num_test, num_fail);
    RUN_TEST(merge_vector_stat_pod, num_test, num_fail);

    // print how many tests passed or failed
    TEST_SUMMARY(argv[0], num_test, num_fail);

    return (num_fail > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}

//============================================================================//


//============================================================================//
