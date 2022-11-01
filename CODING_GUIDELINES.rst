Geant4 Coding Guidelines 
=====

**This document is a work in progress!** Additions and comments, through Merge
Requests and Issues are welcome.

This document only covers recommended use of C++/CMake in Geant4. For instructions on how
to develop and contribute changes to Geant4, please see the `Contribution Guide <./CONTRIBUTING.rst>`_.

Guidelines on C++ Standard and Feature Use
=====
ISO Standard and Platform Support
-----
Geant4 is developed using C++17, with the following language and standard library
features required (full support for preceeding standards is assumed, along with all
removals in C++17):

- Language:

  - `Auto return type deduction <https://en.cppreference.com/w/cpp/language/function#Return_type_deduction_.28since_C.2B.2B14.29>`_
  - `[[attribute-list]] <https://en.cppreference.com/w/cpp/language/attributes) (but no specific attributes are checked>`_
  - `Constexpr if <https://en.cppreference.com/w/cpp/language/if>`_
  - `Class template argument deduction <https://en.cppreference.com/w/cpp/language/class_template_argument_deduction>`_
  - `__has_include <https://en.cppreference.com/w/cpp/preprocessor/include>`_
  - `Initializer in if/switch <https://en.cppreference.com/w/cpp/language/if>`_
  - `Inline variables <https://en.cppreference.com/w/cpp/language/inline>`_
  - `Relaxed constexpr <https://en.cppreference.com/w/cpp/language/constant_expression>`_
  - `Generic lambdas <https://en.cppreference.com/w/cpp/language/lambda>`_
  - `Structured Bindings <https://en.cppreference.com/w/cpp/language/structured_binding>`_

- Standard Library (only presence checked for now):

  - `any <https://en.cppreference.com/w/cpp/header/any>`_
  - `make_unique <https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique>`_
  - `optional <https://en.cppreference.com/w/cpp/header/optional>`_
  - `string_view <https://en.cppreference.com/w/cpp/header/string_view>`_
  - `variant <https://en.cppreference.com/w/cpp/header/variant>`_
  - `filesystem <https://en.cppreference.com/w/cpp/header/filesystem>`_

The currently supported compilers/platforms are:

- `GNU GCC <https://gcc.gnu.org>`_ 7 or newer
- `Clang <https://clang.llvm.org>`_ 8 or newer, with suitable standard library support
- `Xcode/AppleClang <https://developer.apple.com/xcode/>`_ 11 or newer on macOS 10.15 or newer
- `MSVC <https://visualstudio.microsoft.com/vs/>`_ 2019 or newer
- `Intel <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/dpc-compiler.html#gs.25ps3h>`_ 19 or newer, with suitable standard library support


Recommended Usage of C++ features
-----
Use ``G4Filesystem.hh`` in place of ``<filesystem>``
^^^^^
Some older compilers we support supply ``<filesystem>`` in the ``std::experimental`` namespace,
and the implementation may be in a separate library from the main standard library. To avoid
repition of the code to work around namespace/library differences, Geant4 supplies this in the
``G4Filesystem.hh`` header of the ``G4globman`` module. It simply detects the appropriate header
and namespace to use, creating a ``G4fs`` namespace alias to ``std::filesystem`` or ``std::experimental::filesystem``
as required. Geant4 code requiring use of filesystem can then be written as, e.g.

.. code-block:: cpp

  #include "G4Filesystem.hh"

  // Instead of this...
  // std::filesystem::path p;

  // .. use this
  G4fs::path p;

All types in ``std::filesystem`` may be used this way.


C++ Coding Style and Formatting Guidelines
=====
The following general style guidelines should be followed to ensure long term maintainability and readability of Geant4's C++ code

- Readability

  - The ``public``, ``protected`` and ``private`` keywords must be used explicitly in the class declaration, and must appear in that order.
  - English and self-explaining names for constants, variables and functions should be used.
  - Avoid the use of underscore "_" characters within variables or function names (i.e. prefer ``theTotalEnergy``, ``SetEnergyTable()`` to ``the_Total_Energy`` or ``Set_Energy_Table()``).
  - The code must be properly indented with 2 spaces (Tabs must be replaced with spaces)

- Consistency

  - Each class name must begin with ``G4`` (ex. ``G4Particle``)
  - Each header file must contain only one or related class declarations, and must use a filename of the form ``G4<name>.hh``
  - Each class implementation's code must go into a single source file which must use a filename of the form ``G4<name>.cc``
  - Template, or inline, class and function implementations should follow their declarations in the same header file

- Maintainability

  - Each header file must be protected from multiple inclusions to avoid multiple declarations and circular dependences. Ex.:

    .. code-block:: cpp

      #ifndef NAME_HH
      #define NAME_HH
      // ...
      #endif


Code Formatting and use of ``clang-format``
-----
A ``.clang-format`` style file is provided in the root of the Geant4 repository that defines the lower level
layout of code. The `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ program may then be used
to automatically format code files with this style, either manually or `through integration with your IDE/Editor of choice <https://clang.llvm.org/docs/ClangFormat.html>`_. Packages supplying ``clang-format`` are available on most platforms, either directly
or as part of an ``llvm...`` or ``clang...`` package and you should consult the database of your package manager for details.
Visual Studio on Windows supplies it `directly with the IDE <https://devblogs.microsoft.com/cppblog/clangformat-support-in-visual-studio-2017-15-7-preview-1/>`_. If you have Linux/macOS and CVMFS, it is also available
via any LCG view based on clang.

If you need to explicitly disable formatting from being applied to a block of code (e.g. numeric tables), then it may be wrapped
using the `special comment blocks <https://clang.llvm.org/docs/ClangFormatStyleOptions.html#disabling-formatting-on-a-piece-of-code>`_

.. code-block:: cpp

   int formatted_code;
   // clang-format off
   void    unformatted_code  ;
   // clang-format on
   void formatted_code_again;

Formatting for an entire file may be switched off by having the comment ``// clang-format off`` at the
top of the file.

At present, application of formatting is optional but recommended in Geant4. However, automatic checks for, and
application of, formatting will be gradually rolled out as part of the Merge Request process, so you should
familiarize yourself with the process.


Organization of C++ Code into Source Code Modules and Compilation into Libraries
=====
Source Code Modules in Geant4
-----
The lowest level grouping of C++ code in Geant4 is into so-called *source code modules* (or just *modules*), with each
organized on disk as:

.. code-block:: console

   AModule/
   - include/
     - G4AModuleInterface.hh
     - ...
   - src/
     - G4AModuleInterface.cc
     - ...
   - History 
   - sources.cmake

The `include/` subdirectory should only contain header files for the interfaces, i.e. classes and functions, of the module,
with implementations in in the `src/` subdirectory. The ``History`` file is a high level changelog for the module, and is
described in detail in `the Contribution Guide <CONTRIBUTING.rst#making-a-merge-request>`_. 

The ``sources.cmake`` file is a CMake script declaring the module to Geant4's build system. Geant4's final libraries are
each composed and built from 1-N source code modules, the composition of each library being based on toolkit functionality such
as physics, geometry. These "categories", and thus the source code modules under it, are managed by specfic Working Groups within
the collaboration. A module developer only needs to declare their source code module to the build system in terms of the code it 
provides and what other module interfaces it uses without needing to know which library any module will end up composed into.
This is the task of the source code module's ``sources.cmake`` file, described in the following sections.

Writing and Maintaining ``sources.cmake`` for a Module
-----
``sources.cmake`` is a `CMake <https://cmake.org>`_ script, and thus is written using the `CMake Language <https://cmake.org/cmake/help/v3.16/manual/cmake-language.7.html>`_ and `commands <https://cmake.org/cmake/help/v3.16/manual/cmake-commands.7.html>`_. Geant4 has
a minimum CMake version of 3.16, and so language features and commands from newer versions should not be used to ensure backward compatibility. To declare source code modules to Geant4's build, a set of CMake commands are provided for developers. These
largely follow CMake's commands such as ``add_library`` (``geant4_add_module``) and the various ``target_xxx`` (``geant4_module_xxx``)
commands for declaring targets and their `usage requirements <https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#build-specification-and-usage-requirements>`_ though with some distinctions due to there not being a one-to-one correspondance
of source code modules to final libraries. The minimal example of a ``sources.cmake`` for a source code module looks like:

.. code-block:: cmake

   # Declare module and inputs
   geant4_add_module(G4foo
     PUBLIC_HEADERS
       G4Foo.hh
       ...
     SOURCES
       G4Foo.cc
   )

   # Declare modules and/or external libraries consumed by the module's code 
   geant4_module_link_libraries(G4foo
     PUBLIC
       G4bar
     PRIVATE
       G4baz
   )

Here, ``geant4_add_module`` declares a module named ``G4Foo`` to the build and lists its headers and sources in ``PUBLIC_HEADERS`` and
``SOURCES`` respectively. Note that subdirectories are not required as it is assumed (as required by the module layout rules)
that headers/sources are present in the module's ``include``/``src`` subdirectories. CMake will check this at configuration time and
raise an error if a file cannot be found. A common question here is "why can't wildcards/globbing be used so I don't have to
explicitly specify the sources"? There are several technical and policy reasons for this:

- CMake is a build system *generator* and `explicitly does not recommend the use of globbing <https://cmake.org/cmake/help/v3.16/command/file.html?highlight=file#filesystem>`_ as it can lead to inconsistent builds or build time costs.
- Geant4 has source code modules with optional sources dependent on configuration arguments, so globbing/wildcarding would require special casing/filtering 
  - As a historical note, the old Configure/GNUmake used this globbing, and had to workaround this issue by an awkward and error
    prone system of preprocessor flags both at Geant4 and application build time 
- Explicitly listed sources are clearer for ongoing and long term development/maintenance, especially when optional sources are involved
  - The build/test can more easily check for inconsistent configurations
- It forces developers to think about the content and build of their modules

The majority of Geant4 modules use interfaces from other modules for their implementations, and these dependencies are declared
using the ``geant4_module_link_libraries`` command. In our example above, the command takes the name of the module whose dependencies we want to declare as its first argument, followed by *usage requirements* stating which other modules are used, and how they are consumed.
These are derived from how other module's code is used locally, with the command:

.. code-block:: cmake

   geant4_module_link_libraries(G4foo
     PUBLIC
       G4bar
     PRIVATE
       G4baz
   )

following from this use of ``G4bar`` and ``G4baz`` interfaces in ``G4foo``:

.. code-block:: c++

   // G4Foo.hh
   #include "G4Bar.hh"

   class G4Foo : public G4Bar
   {
      ...
      void DoSomething();
   };

   // G4Foo.cc
   #include "G4Foo.hh"

   #include "G4Baz.hh"

   void DoSomething()
   {
      G4Baz x;
      G4double theAnswer = x.Calculate();
      ...
   }

Thus, ``G4bar`` is a ``PUBLIC`` dependency of ``G4foo`` because the latter exposes an interface of the
former in its own public interface. Correspondingly, ``G4baz`` is a ``PRIVATE`` dependency of ``G4foo`` 
because the latter only consumes interfaces of the former in its implementation details (or, no interfaces
of ``G4baz`` are exposed to users of ``G4foo``). In the case that a module consumes another in both ``PUBLIC``
and ``PRIVATE`` contexts, declare the dependency as ``PUBLIC`` as this has higher precedence.

Note that you do *not* need to know which library a given source code module is eventually compiled into.
The CMake scripts will determine this and resolve the final library dependencies and linking appropriately.
``geant4_module_link_libraries`` can also take external libraries as usage requirements, for example

.. code-block:: cmake

   geant4_module_link_libraries(G4foo
     PUBLIC
       G4bar
     PRIVATE
       G4baz
       ${ZLIB_LIBRARIES}
       XercesC::XercesC
   )

with the same rules as modules for declaring them as ``PUBLIC/PRIVATE``. However, before using an external
library in Geant4, you **must** consult with the Software Management Working Group to check that it can
be supported and is compatible with the Geant4 License. TBD: Document supported libs/how to use?

Most modules should only need the above two commands to integrate them in the Geant4 build, but a few extra
commands are available for more advanced use cases. First, if a module has sources that are only added if
a particular confiuration option is set, they may be added after module creation with ``geant4_module_sources``, e.g.:

.. code-block:: cmake

   geant4_add_module(G4foo
     PUBLIC_HEADERS
       G4Foo.hh
     SOURCES
       G4Foo.cc
   )

   if(GEANT4_USE_CORGE)
     geant4_module_sources(G4foo
       PUBLIC_HEADERS 
         G4UseCorge.hh
      SOURCES
         G4UseCorge.cc
     )
   endif()

that takes the same arguments as the ``geant4_add_module`` command. An explicit conditional must be used at
present, as the ``geant4_xxx`` commands do not yet support CMake `generator expressions <https://cmake.org/cmake/help/v3.16/manual/cmake-generator-expressions.7.html>`_.

Additional compile definitions (i.e. preprocessor defines) may be added to the compilation flags with ``geant4_module_compile_definitions``, e.g.

.. code-block:: cmake

   geant4_module_compile_definitions(G4foo
     PUBLIC AFLAG
     PRIVATE BFLAG
   )

This would add (on UNIX) ``-DAFLAG`` and ``-DBFLAG`` to the compiler flags for building ``G4foo``, and ``-DAFLAG`` to the compiler 
flags for building any module using the ``G4foo`` module (and thus library it is composed into). This command **must not** be
used to add arbitrary developer only debugging flags. Please consult with the Software Management Working Group if you need this
functionality.
     

Tools for Checking Source Code Module Interfaces and Dependencies
-----
Source code modules in Geant4 **must** be organised and designed so that

- There are no circular dependencies between modules, direct or indirect.
- If a module includes headers (i.e. uses interfaces) from a module or external package, it
  must declare a dependency on this via ``geant4_module_link_libraries``.
- A module must not declare a dependency on a module or external package that it does 
  not use.

To help developers work with modules and identify issues with module dependencies, a Python 3 script 
``geant4_module_check.py`` is generated in the build directory. This may be run manually,
or as a dedicated test via ``ctest``. Note that in both cases only the source code modules
in the current build will be considered, e.g. optional modules such as ``G4gdml`` must
be enabled to analyse them. Full help and a list of command line argument for the script
may be printed via

.. code-block:: console

   $ ./geant4_module_check.py --help

The following sections will walk through some of the more common use cases and options.
We will always assume that the script is being run from the build directory.


Querying Modules and their Interfaces
^^^^^
A list of all source code modules enabled in the current build may be printed using:

.. code-block:: console

   $ ./geant4_module_check.py --list

The public interface of a given module, i.e. the public headers that it provides,
may be printed using:

.. code-block:: console

   $ ./geant4_module_check.py -i <modulename>

To determine which module provides a given header, the ``--provides/-p`` argument
can be used, e.g.

.. code-block:: console

   $ ./geant4_module_check.py -p G4String.hh
   G4globman

The directory in the source tree where the module code is located may also be printed:

.. code-block:: console

   $ ./geant4_module_check.py -s <modulename>
   ... system dependent path ...

Querying Library/Module Composition
^^^^^
The library into which a given source code module is compiled may be printed using:

.. code-block:: console

   $ ./geant4_module_check.py --library G4globman
   G4global

A list of all defined libraries and which modules they are composed from can also be printed:

.. code-block:: console

   $ ./geant4_module_check.py --libraries

Checking for Circular Dependencies
^^^^^
At the global level, circular dependencies between source code modules may be detected
using the ``--find-cycles`` argument:

.. code-block:: console

   $ ./geant4_module_check.py --find-cycles
   No cycles detected in module dependency graph

If a cycle is detected, it will print the chain of dependencies leading to the cycle
and exit with a non-zero code, e.g.:

.. code-block:: console

   $ ./geant4_module_check.py --find-cycles
   Cycles detected in module dependency graph:
   G4partman -> G4leptons -> G4partman

Here, the cycle is printed as the sequence of modules in the cycle, and should
be read from left to right, with ``->`` meaning "depends on". The first and
last modules should always be the same. Note that if there is more than one
cycle in the module dependencies, only one will be printed. It is up to the
developer to fix the identified cycle first before trying to detect/fix any further
issues. Cycle detection is also added as a direct test in ``ctest`` and will be
run in any invocation of this. It may also be run in isolation via

.. code-block:: console

   $ ctest -R validate-no-module-cycles


Checking for Inconsistent Dependencies in Modules
^^^^^
Here we define a *consistent* source code module as follows:

- If module ``A`` ``#include`` s a header from module ``B`` in any of its own header (``.hh``)
  files, then ``A`` is consistent only if it declares ``B`` as a ``PUBLIC`` dependency
  in ``geant4_module_link_libraries``.
- If module ``A`` ``#include`` s a header from module ``B`` in any of its own implementation (``.cc``)
  files *only*, then ``A`` is consistent only if it declares ``B`` as a ``PRIVATE`` dependency
  in ``geant4_module_link_libraries``.
- If any header from ``B`` is included by *both* ``A`` s header and implementation files, then ``B``
  must be a ``PUBLIC`` dependency of ``A``.

Inconsistencies between the headers ``#include`` ed by module sources and the dependencies declared 
to CMake may be detected for a given module using the ``--check-consistency`` argument or its short ``-c``
form:

.. code-block:: console

   $ ./geant4_module_check.py -c G4globman
   G4globman appears consistent

If any inconsistencies are found, they will be printed to standard error, e.g.

.. code-block:: console

   $ ./geant4_module_check.py -c G4phys_ctor_em
   G4phys_ctor_em has inconsistent dependencies:
     + may require PUBLIC or INTERFACE dependencies: {'G4emdna-processes', 'G4emlowenergy'}
     + may require PRIVATE dependencies: {'G4transportation', 'G4procman', 'G4materials', 'G4partman'}
     - may not require PUBLIC dependencies: {'G4decay', 'G4procman', 'G4materials', 'G4partman'}
     - may not require PRIVATE dependencies: {'G4hadronic_xsect', 'G4emdna-processes', 'G4emlowenergy', 'G4hadronic_proc'}

The same check may be run over all modules in the build at the same same using:

.. code-block:: console

   $ ./geant4_module_check.py --find-inconsistencies

If any inconsistencies are found, they will be printed to standard error module by module in the same
format as shown above. The reported inconsistencies describe the following cases for a
module "A":

- *+ may require PUBLIC or INTERFACE dependencies: {<module1>, ... <moduleN>}*
  
  Module "A" ``#include`` s headers from the listed modules in its own header files,
  but has not declared those modules as ``PUBLIC`` dependencies. An ``INTERFACE``
  dependency is *only* needed if "A" is header only.

  The listed modules should be added to "A"s ``PUBLIC`` dependencies.

- *+ may require PRIVATE dependencies: {<module1>, ... <moduleN>}*

  Module "A" ``#include`` s headers from the listed modules only in its source files (```.cc``),
  but has not declared those modules as ``PRIVATE`` dependencies.

  The listed modules should be added to "A"s ``PRIVATE`` dependencies.

- *- may not require [PUBLIC|PRIVATE] dependencies: {<module1>, ... <moduleN>}*

  Here, module "A" has declared that the listed modules are dependencies, but
  none of its header/source files ``#include`` headers from them.

  The listed modules should be removed from "A"s ``PUBLIC/PRIVATE`` dependencies
  as appropriate.

It's important to note upfront that these checks are only high level. No code is
actually compiled, only roughly parsed for ``#include`` statements, and so cannot
detect problems caused by reliance on transitive includes. This is the case that
the module's code uses, e.g. ``G4String``, but has no explicit ``#include "G4String.hh"``,
relying instead on another header including that file. The checks also cannot pick
up the case that a header is included but no interfaces in that header are used, though
even compilation cannot pick up this problem and it is the developer's responsibility to
only include what is required. Nevertheless, it provides
a handy set of tools that should be used to detect the most common module dependency
related issues.

To pick up these lower level issues, a dedicated module-by-module build is performed 
as part of Continuous testing for every Merge Request. This "GranularBuild" check
*must* pass for the Merge Request to progress, and by using a full isolated compilation
of each module provides a highly reliable check and detailed information on resolving 
any remaining dependency issues.


Guidelines for use of Static Analysis and Sanitizer Tools
=====
Static Analysis with Coverity
-----
Detailed static analysis is performed on each monthly reference tag of the ``master`` branch using the 
`Coverity <https://coverity.cern.ch>`_ tool. A report is emailed to developers each month, and the 
`detailed reports <https://coverity.cern.ch>`_ for code you are responsible should be reviewed at this 
point to triage and fix issues through Merge Requests.


Static Analysis and Maintenance using Clang Tidy
-----
A ``.clang-tidy`` check file is provided in the root of the Geant4 repository that defines a minimal set of
checks for code clarity and robustness including:

- Modernization
- Performance (not a substitute for detailed profiling/benchmarking)
- Readability

The `clang-tidy <https://clang.llvm.org/extra/clang-tidy/>`_ program may then be used
to check code for issues. It will warn about these issues, suggest the recommended fix, and
optionally apply this automatically including reformatting if required using ``clang-format``.
Packages supplying ``clang-tidy`` are available on most platforms, either directly
or as part of an ``llvm...`` or ``clang...`` package and you should consult the database of your package manager for details.
Visual Studio on Windows supplies it `when installing C++ support <https://docs.microsoft.com/en-us/cpp/code-quality/clang-tidy?view=msvc-160>`_. If you have Linux/macOS and CVMFS, it is also available
via any LCG view based on clang. Geant4's ``.clang-tidy`` file has been tested on LLVM 7 upwards, with
LLVM 8 the recommended minimum.

``clang-tidy`` may be incorporated into your development workflow using the ``run-clang-tidy`` (``run-clang-tidy.py`` in some installs) program provided by ``clang-tidy`` installs (recommended) or CMake's native support. 
`Integrations with many IDEs <https://clang.llvm.org/extra/clang-tidy/Integrations.html>`_ are also available, though these
are currently untested. Feedback and documentation on the  use of these is welcome.

To run ``clang-tidy`` directly, the ``run-clang-tidy`` program may be used. Depending on how ``clang-tidy`` was packaged,
this may be present alongside ``clang-tidy``, or present under the ``share/llvm/clang`` directory under the main LLVM 
install prefix. It may also be named ``run-clang-tidy.py`` in older LLVM versions, so substitute that command in the examples
below. It's recommended that you have the directories holding ``clang-tidy`` and ``run-clang-tidy`` appended to your ``PATH`` 
so that they and the tools they run are located easily. To use ``run-clang-tidy``, first configure your build
of Geant4 as normal, adding the CMake argument:

.. code-block:: console

   $ cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON <otherargs>

to generate a compile command database for use by ``clang-tidy``. The ``run-clang-tidy`` tool may then
be run from the build directory as

.. code-block:: console

   $ run-clang-tidy

or

.. code-block:: console

   $ run-clang-tidy -jN

to run ``N`` parallel jobs to run (just like ``make -jN``). These will however run ``clang-tidy`` over every
compiled source file in Geant4. It is better to pass a subset of the code to run the tool over
by giving it one or more paths to directories relative to the root of the source tree:

.. code-block:: console

   $ run-clang-tidy source/global

or individual files:

.. code-block:: console

   $ run-clang-tidy source/global/management/src/G4UnitsTable.cc

These paths pattern match, so the above could be reduced to:

.. code-block:: console

   $ run-clang-tidy global/
   ...
   $ run-clang-tidy G4UnitsTable.cc

By default, ``run-clang-tidy`` will use the configuration from the ``.clang-tidy`` file in the
root of the Geant4 repository. These checks may be overriden using the ``-checks`` argument, e.g.
to run only the `performance-inefficient-vector-operation <https://clang.llvm.org/extra/clang-tidy/checks/performance/inefficient-vector-operation.html>`_ check:

.. code-block:: console

   $ run-clang-tidy -checks="-*,performance-inefficient-vector-operation"

However, it is strongly recommended that you consult with the Software Management working group
during any Merge Request that adds fixes suggested by checks outside those in the main ``.clang-tidy``
file or those described in `Additional Suggested Checks`_. Not all may be suitable, or may clash with other 
requirements (e.g. some `checks are stylistic in nature <https://clang.llvm.org/extra/clang-tidy/checks/modernize/use-trailing-return-type.html>`_)

Automatic fixing of the detected issues can be enabled using the ``-fix`` option, e.g.

.. code-block:: console

   $ run-clang-tidy -fix global/

and to additionally reformat with ``clang-format``:

.. code-block:: console

   $ run-clang-tidy -fix -format global/

It is strongly recommended that you apply fixes in one or more separate commits on your Topic Branches. This assists
in the Merge Request review process, as well as enabling easy revert/correction if needed.

As currently implemented Geant4's ``.clang-tidy`` file does not check Geant4 header files and so will only report issues
in implementation files. This is done at present to avoid a cascade of fixes being applied across
multiple categories other than the one being developed. Headers for a given category can be included into the checks
and automatic fixes using the ``-header-filter`` argument to ``run-clang-tidy``. This can only be used to
accept, not reject, header patterns, so is best used when running checks over individual categories/sub-categories:

.. code-block:: console

  # Check and fix the G4intercoms category headers and sources only
  $ run-clang-tidy --header-filter="intercoms/" intercoms/

It is strongly recommended to recompile after each invocation of ``run-clang-tidy`` as it can
introduce syntax errors when fixing some corner cases. See the section below for details on suppressing
checks on particular lines using ``// NOLINT`` comments.

As with ``clang-format``, linting by ``clang-tidy`` can be disabled on specific lines if required.
If you need to disable checks, then `special comments <https://clang.llvm.org/extra/clang-tidy/index.html#suppressing-undesired-diagnostics>`_
may be used

.. code-block:: cpp

  // Silence all diagnostics on this line
  Foo(int param); // NOLINT

  // Silence all diagnostics on the line following the special comment
  // NOLINTNEXTLINE
  Foo(std::string param)

However, note that ``clang-tidy`` can only be disabled for single lines, rather than blocks, of code.


``clang-tidy`` may also be run automatially as part of the build generated by CMake. To enable this, 
simply configure the build as:

.. code-block:: console
   $ cmake -DCMAKE_CXX_CLANG_TIDY="/path/to/clang-tidy" <otherargs>

Building Geant4 will then run ``clang-tidy`` alongside the full compilation and report detected
issues, and suggested fixes, as warnings, e.g.

.. code-block:: console
   
   $ ninja
   ...
   [50/51] Building CXX object source/CMakeFiles/G4global.dir/global/management/src/G4UnitsTable.cc.o
   ...
   /src/geant4-dev.git/source/global/management/src/G4UnitsTable.cc:586:13: warning: use '= default' to define a trivial destructor [modernize-use-equals-default]
   G4BestUnit::~G4BestUnit() {}
               ^             ~~
                             = default;

However, it is not recommended to use the CMake integration unless you know what you are doing, as it will slow
down compilation times as it runs ``clang-tidy`` over every compile file in Geant4. It will also not, in general, apply
fixes to files cleanly, leading to compile errors and corrupted source files. The CMake integration
is therefore best used for incremental development to act like an additional set of compiler warnings.


Additional Suggested Checks
^^^^^
The default set of checks in Geant4's ``.clang-tidy`` file have been selected on the basis of providing most benefit
and with cleanly applyable fixes. A range of additional fixes are listed below which developers should consider
for application on a priority and case-by-case basis, reviewing the applied fixes for correctness and applicability.
These may be added to the default set of checks, or their priority changed, as the code evolves. 

Recommended for application. They are not in default purely as they may not apply automatic fixes cleanly
or be appropriate for all use cases.

- `modernize-avoid-c-arrays <https://clang.llvm.org/extra/clang-tidy/checks/modernize/avoid-c-arrays.html>`_
- `modernize-loop-convert <https://clang.llvm.org/extra/clang-tidy/checks/modernize/loop-convert.html>`_

  - Whilst an obvious modernization, it is not yet in the default set of checks as its fixes should be reviewed by the
    developer to correct the `range declaration <https://en.cppreference.com/w/cpp/language/range-for>`_ if needed. For
    certain types of collections, especially those of pointers, it may not use the correct `const` qualifier.

- `modernize-use-default-member-init <https://clang.llvm.org/extra/clang-tidy/checks/modernize/use-default-member-init.html>`_, `readability-redundant-member-init <https://clang.llvm.org/extra/clang-tidy/checks/readability/redundant-member-init.html>`_

  - It is recommended to apply these together as they provide a coherent check on member initialization.

- `modernize-use-equals-delete <https://clang.llvm.org/extra/clang-tidy/checks/modernize/use-equals-delete.html>`_
- `modernize-use-emplace <https://clang.llvm.org/extra/clang-tidy/checks/modernize/use-emplace.html>`_
- `modernize-use-override <https://clang.llvm.org/extra/clang-tidy/checks/modernize/use-override.html>`_
- `readability-uniqueptr-delete-release <https://clang.llvm.org/extra/clang-tidy/checks/readability/uniqueptr-delete-release.html>`_

Suggested for application at developer discretion

- `readability-container-size-empty <https://clang.llvm.org/extra/clang-tidy/checks/readability/container-size-empty.html>`_
- `readability-implicit-bool-conversion <https://clang.llvm.org/extra/clang-tidy/checks/readability/implicit-bool-conversion.html>`_
- `performance-move-const-arg <https://clang.llvm.org/extra/clang-tidy/checks/performance/move-const-arg.html>`_
- `performance-unnecessary-value-param <https://clang.llvm.org/extra/clang-tidy/checks/performance/unnecessary-value-param.html>`_

  - The suggested fix is typically to change a function argument from ``T`` to ``const T&``. This does lead to a change in
    function signature that can affect user code, so should generally not be considered if the function is a public interface.

- `performance-move-constructor-init <https://clang.llvm.org/extra/clang-tidy/checks/performance/move-constructor-init.html>`_
- `performance-no-int-to-ptr <https://clang.llvm.org/extra/clang-tidy/checks/performance/no-int-to-ptr.html>`_
- `performance-noexcept-move-constructor <https://clang.llvm.org/extra/clang-tidy/checks/performance/noexcept-move-constructor.html>`_

Optional for application. Primarily checks for style/clarity

- `readability-braces-around-statements <https://clang.llvm.org/extra/clang-tidy/checks/readability/braces-around-statements.html>`_
- `readability-else-after-return <https://clang.llvm.org/extra/clang-tidy/checks/readability/else-after-return.html>`_
- `readability-isolate-declaration <https://clang.llvm.org/extra/clang-tidy/checks/readability/isolate-declaration.html>`_
- `readability-simplify-boolean-expr <https://clang.llvm.org/extra/clang-tidy/checks/readability/simplify-boolean-expr.html>`_
- `modernize-return-braced-init-list <https://clang.llvm.org/extra/clang-tidy/checks/modernize/return-braced-init-list.html>`_
- `modernize-pass-by-value <https://clang.llvm.org/extra/clang-tidy/checks/modernize/pass-by-value.html>`_

  - It's important to consider the type being passed here, as this will only provide a benefit if the type is always
    expensive to copy. The canonical case here is ``G4String/std::string`` which seems an obvious candidate, but only benefits
    if the string is always going to be longer than the small string optimization threshold.


Use of `G4Backtrace` for backtracing
-----
Whilst you should use full debuggers and santizers for in depth debugging,
a simple signal handler is available in Geant4 for tracing the source of simple
errors such as segmentation faults. This may be enabled automatically by building
Geant4 with the CMake option `GEANT4_BUILD_BUILTIN_BACKTRACE` set to `ON`.

It can also be enabled on demand in applications by including the relevant header
and calling the `Enable` static member function:

.. code-block:: cpp

   #include "G4Backtrace.hh"

   ...

   G4Backtrace::Enable();

before calling any other Geant4 functionality. See the `G4Backtrace.hh` header
for further documentation on its functionality.
  
`G4Backtrace` is always enabled in Continuous and Nightly CI builds to help report
and triage any immediate issues on all platforms, especially those you may not have access to.

Development and Maintenance with Sanitizers
-----
Sanitizers instrument code to pick up issues such as memory errors and thread
Geant4's build system currently supports memory, thread, and undefined behaviour sanitizers when
using the GNU or Clang compilers. Only one sanitizer may be enabled in a given build, and may
be chosen using the ``GEANT4_BUILD_SANITIZER`` option:

- ``address``: enable `Address Sanitizer <https://github.com/google/sanitizers/wiki/AddressSanitizer>`_
- ``thread``: enable `Thread Sanitizer <https://github.com/google/sanitizers/wiki/ThreadSanitizerCppManual>`_
- ``undefined``: enable `Undefined Behaviour Sanitizer <https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html>`_

Both GCC and Clang implement these sanitizers in similar ways and so may be used in the same way. See the
links above for examples of the types of problem that can be picked up, and options for running the instrumented
programs.



