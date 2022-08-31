Coding Guidelines for Geant4
=====

**This document is a work in progress!** Additions and comments, through Merge
Requests and Issues are welcome.

This document only covers recommended use of C++/CMake in Geant4. For instructions on how
to develop and contribute changes to Geant4, please see the `Contribution Guide <./CONTRIBUTING.rst>`_.

Geant4 C++ Guidelines
=====
Supported Standard and Platforms
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


Geant4 Code Documentation and Formatting
=====
Coding Style and Formatting Guidelines
-----
The following coding-style guidelines should be followed to ensure long term maintainability and readability

- Readability
  - The ``public``, ``protected`` and ``private`` keywords must be used explicitly in the class declaration.
  - English and self-explaining names for constants, variables and functions should be used.
  - Avoid the use of underscore "_" characters within variables or function names (prefer ``theTotalEnergy``, ``SetEnergyTable()``).
  - The code must be properly indented with 2 spaces (Tabs must be replaced with spaces)
- Consistency
  - Each class name must begin with ``G4`` (ex. ``G4Particle``)
  - Each header file must contain only one or related class declarations.
  - Each class implementation source code must go into a single source file.
- Maintainability
  - Each header file must be protected from multiple inclusions to avoid multiple declarations and circular dependences. Ex.:

    .. code-block:: cpp

      #ifndef NAME_HH
      #define NAME_HH
      // ...
      #endif

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



Geant4 Static Analysis/Sanitizer Tools and Guidelines
=====
Development and Maintenance using Clang Tidy
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


Analysis of Reference Tags with Coverity
-----
Detailed static analysis is performed on each monthly reference tag using the `Coverity <https://coverity.cern.ch>`_ tool.
A report is emailed to developers each month, and the `detailed reports <https://coverity.cern.ch>`_ for code
you are responsible should be reviewed at this point to identify and fix issues.


Geant4 CMake Use and Guidelines
=====
The minimum supported CMake version for building Geant4 is 3.16, with versions up to 3.22 checked.
