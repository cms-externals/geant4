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

#ifndef testing_hh_
#define testing_hh_

// C headers
#include <cassert>
#include <cstdint>
#include <cstdio>

// C++ headers
#include "G4ios.hh"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

//----------------------------------------------------------------------------//

// EXPECT_EQ
// EXPECT_NE
// ASSERT_FALSE
// ASSERT_TRUE
// ASSERT_CLOSE

//----------------------------------------------------------------------------//

#define EXPECT_EQ(lhs, rhs)                                                                     \
  if (lhs != rhs)                                                                               \
  {                                                                                             \
    std::stringstream ss;                                                                       \
    ss << "Expectation failed : " << #lhs << " != " << #rhs << " @ line " << __LINE__ << " of " \
       << __FILE__ << " (" << lhs << " != " << rhs << ")";                                      \
    std::cerr << ss.str() << std::endl;                                                         \
    throw std::runtime_error(ss.str());                                                         \
  }

//----------------------------------------------------------------------------//

#define EXPECT_NE(lhs, rhs)                                                                     \
  if (lhs == rhs)                                                                               \
  {                                                                                             \
    std::stringstream ss;                                                                       \
    ss << "Expectation failed : " << #lhs << " == " << #rhs << " @ line " << __LINE__ << " of " \
       << __FILE__ << " (" << lhs << " == " << rhs << ")";                                      \
    std::cerr << ss.str() << std::endl;                                                         \
    throw std::runtime_error(ss.str());                                                         \
  }

//----------------------------------------------------------------------------//

#define ASSERT_FALSE(expr)                                    \
  if (expr)                                                   \
  {                                                           \
    std::stringstream ss;                                     \
    ss << "Assertion failed : "                               \
       << "Expression: ( " << #expr << " ) "                  \
       << "failed @ line " << __LINE__ << " of " << __FILE__; \
    std::cerr << ss.str() << std::endl;                       \
    throw std::runtime_error(ss.str());                       \
  }

//----------------------------------------------------------------------------//

#define ASSERT_TRUE(expr)                                     \
  if (!(expr))                                                \
  {                                                           \
    std::stringstream ss;                                     \
    ss << "Assertion failed : "                               \
       << "Expression: !( " << #expr << " ) "                 \
       << "failed @ line " << __LINE__ << " of " << __FILE__; \
    std::cerr << ss.str() << std::endl;                       \
    throw std::runtime_error(ss.str());                       \
  }

//----------------------------------------------------------------------------//

#define ASSERT_CLOSE(lhs, rhs, epsilon)                                                         \
  if (fabs(lhs - rhs) > epsilon)                                                                \
  {                                                                                             \
    std::stringstream ss;                                                                       \
    double _diff = fabs(lhs - rhs);                                                             \
    ss << "failed. \n\tExpression: fabs( " << lhs << " - " << rhs << " ) < " << std::scientific \
       << std::setprecision(6) << epsilon << " [ diff = " << _diff << " ] "                     \
       << "@ line " << __LINE__ << " of " << __FILE__;                                          \
    std::cerr << ss.str() << std::endl;                                                         \
    throw std::runtime_error(ss.str());                                                         \
  }                                                                                             \
  else                                                                                          \
  {                                                                                             \
    std::stringstream ss;                                                                       \
    double _diff = fabs(lhs - rhs);                                                             \
    std::cout << "passed.  " << std::flush;                                                     \
    ss << std::scientific << std::setprecision(6) << "lhs = " << lhs << ", rhs = " << rhs       \
       << ", diff = " << _diff << " < epsilon = " << epsilon;                                   \
    std::cout << ss.str() << std::endl;                                                         \
  }

//----------------------------------------------------------------------------//

#define PRINT_HERE printf(" [%s@'%s':%i]\n", __FUNCTION__, __FILE__, __LINE__)

//----------------------------------------------------------------------------//

#define TEST_SUMMARY(argv_0, ntest_counter, nfail_counter)                              \
  {                                                                                     \
    std::stringstream ss;                                                               \
    ss << "\nTesting completed.\n" << std::endl;                                        \
    ss << "[" << argv_0 << "] ";                                                        \
    if (num_fail > 0)                                                                   \
      ss << "Tests failed: " << nfail_counter << "/" << ntest_counter << std::endl;     \
    else                                                                                \
      ss << "Tests passed: " << (ntest_counter - nfail_counter) << "/" << ntest_counter \
         << std::endl;                                                                  \
    std::cout << ss.str();                                                              \
  }

//----------------------------------------------------------------------------//
// Usage:
//  try
//  {
//      RUN_TEST(test_serialize, num_test, num_fail);
//  }
//  catch(std::exception& e)
//  {
//      std::cerr << e.what() << std::endl;
//  }
//
#define RUN_TEST(func, ntest_counter, nfail_counter)                     \
  {                                                                      \
    try                                                                  \
    {                                                                    \
      std::cout << "\nRunning test " << #func << " ... \n" << std::endl; \
      ntest_counter += 1;                                                \
      func();                                                            \
    }                                                                    \
    catch (std::exception & e)                                           \
    {                                                                    \
      std::cerr << e.what() << std::endl;                                \
      nfail_counter += 1;                                                \
    }                                                                    \
  }

//----------------------------------------------------------------------------//

#endif
