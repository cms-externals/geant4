# checkCXXFeatures
Simple checks for C++17 (and beyond) support by the compiler and standard
library for Geant4 development.

Whilst located in the Geant4 development repository and integrated with
the CTest/CDash system as "test300", it may also be built standalone
to test outside of this system.

## Building and running
Requirements: CMake >= 3.12 and a C/C++ compiler.

```
$ cmake .
$ cmake --build .
$ ./checkCXXFeatures
```

The `cmake` stage will perform a normal CMake project configuration, requiring
that the compiler support C++17 at minimum. Additional checks of the language and standard
libraries are performed, and do *not* cause configuration failures. Success/Fail
of these checks are reported in the CMake output. As checks are implemented using
cmake's `try_compile` command, debugging them locally can be done with the cmake's
`--debug-trycompile` command line option. For debugging checks on CI systems, a CMake
option `CCF_BUILD_FAILING_CHECKS` may be set to `ON` to add failing checks to the main
build so that their failure may be seen in CI output.

The results of these checks are used to configure the `checkCXXFeatures` program that can be run
to report information on the system, compiler, and set of features tested.
It returns `0` if all tested features are supported, and `1`
otherwise.

## Adding new checks
Checking for a C++ header is usually a simple case of using `check_include_file_cxx`:

```cmake
# -- std::optional
check_include_file_cxx(optional CXXSTDLIB_HAS_OPTIONAL)
if(NOT CXXSTDLIB_HAS_OPTIONAL)
  check_include_file_cxx(experimental/optional CXXSTDLIB_HAS_EXPERIMENTAL_OPTIONAL)
endif()
```

Checking a C++ language/standard library feature requires

- Addition of a C++ source file in the `cmake` directory implementing the test
  - The file _must_ have a `main()` function, even if empty.
  - The file _must_ fail to compile if the C++ feature being checked is not supported.
- Addition of a call to the `check_cxx_feature` function in the top level `CMakeLists.txt`:

  ```
  check_cxx_feature(<var> <bindir> <src>)
  ```

  where `<var>` is the name of a CMake variable in which to store the result of
  the check (true/false), `<bindir>` is the directory to use to compile the test,
  and `<src>` is the full path to C++ file implementing the check.

  The project sets up a common base directory for checks, stored in the `testCxxFeatures_CHECK_DIR`
  CMake variable, so a typical invocation looks like:

  ```cmake
  check_cxx_feature(CXX_HAS_AUTO_RETURN_TYPE
    ${testCxxFeatures_CHECK_DIR}/cxx_has_auto_return_type
    ${PROJECT_SOURCE_DIR}/cmake/check_cxx_auto_return_type.cc)
  ```
- Addition of a reporting function in `checkCXXFeatures.cc.in` for the check. These are
  typically implemented, using the preceeding example, as::

  ```cpp
  // ...

  // Before the main()
  #cmakedefine01 CXX_HAS_AUTO_RETURN_TYPE
  bool reportAutoReturnType(std::ostream& os) {
    os << "- auto Return Type Deduction : " << AS_YES_OR_NO(CXX_HAS_AUTO_RETURN_TYPE) << std::endl;
    return CXX_HAS_AUTO_RETURN_TYPE == 1;
  }

  // ...

  // in the main()
  bool cxxStdIsSupported = true;

  std::cout << "## C++ Language Features" << std::endl;
  cxxStdIsSupported = reportAutoReturnType(std::cout) && cxxStdIsSupported;
  // ...
  ```

## Known Issues/Todos
- Additional of further checks and documentation to provide a developer guide to
  C++17 and beyond (allowed/recommended usage)
- Only C++17 support is checked (to add C++20 and C++23, but these require newer
  CMake versions).
