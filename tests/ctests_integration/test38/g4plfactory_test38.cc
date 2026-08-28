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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysListFactory.hh"  // original version of factory
#include "G4PhysListFactoryAlt.hh"  // new version

#ifndef _WIN32
// dynamic library loading ...
#  include <dlfcn.h>

#  ifdef USE_GETOPT
extern "C"
{
#    include <getopt.h>
}
#    ifndef GETOPTDONE  // Some getopt.h's have this, some don't...
#      define GETOPTDONE (-1)
#    endif
#  endif
#else
// dynamic library loading ...
#  include <windows.h>
// this defines max() and min() macros that interfere with
// std::max() etc in <algorithms>
#  undef max
#  undef min
#endif

#include <algorithm>
#include <algorithm>  // std::max(), std::min by G4PhysicsVector.icc etc
#include <cassert>
#include <cstdlib>  // abort(), atoi()
#include <ostream>
#include <sstream>
#include <typeinfo>
#include <vector>

#if defined(__clang__) || defined(__GNUC__)
// c++ name mangling
#  include <cxxabi.h>
#endif

#undef USE_GETOPT

#ifndef DYNAMIC_LOADER_PREFIX
#  define DYNAMIC_LOADER_PREFIX "cmake-failed-PREFIX"
#endif
#ifndef DYNAMIC_LOADER_SUFFIX
#  define DYNAMIC_LOADER_SUFFIX "cmake-failed-SUFFIX"
#endif

#include "G4PhysicsConstructorRegistry.hh"
#include "G4RunManager.hh"

//
// pull in a physics list definition into the main program
//
#include "G4PhysListStamper.hh"

#include "MyPL0.hh"
G4_DECLARE_PHYSLIST_FACTORY(MyPL0);

#if (defined G4TEST38_TRY_EXTRAS)

#  if (defined G4TEST38_SHARED)
#    if (defined _WIN32) && (!defined G4TEST38_WIN_NO_PRELINK1)
// if we pre-link to a library in Windows it is like static linking to
// it; unlike other cases where it still treats it like a shared library
G4_REFERENCE_PHYSLIST_FACTORY(MyPL1);
#    endif
#  endif

#  if (defined G4TEST38_STATIC)
// factories built into external libraries ... linking against those
// libraries isn't sufficient to pull them into self-register with
// the registry, so we need to tickle them with the following:

G4_REFERENCE_PHYSLIST_FACTORY(MyPL1);
// should expand to:
/*
class G4VModularPhysicsList;
template <class T> class TMyPL1;
typedef TMyPL1<G4VModularPhysicsList> MyPL1;
extern const G4PhysListStamper<MyPL1>& MyPL1Factory;
const G4PhysListStamper<MyPL1>& MyPL1FactoryRef = MyPL1Factory;
*/

G4_REFERENCE_PHYSLIST_FACTORY(MyPL2);

G4_REFERENCE_PHYSLIST_FACTORY_NS(myns::MyNSPL3, myns, MyNSPL3);
// should expand to:
/*
class G4VModularPhysicsList;
namespace myns {
  template <class T> class TMyNSPL3;
  typedef TMyNSPL3<G4VModularPhysicsList> MyNSPL3;
  extern const G4PhysListStamper<MyNSPL3>& MyNSPL3Factory;
  const G4PhysListStamper<MyNSPL3>& MyNSPL3FactoryRef = MyNSPL3Factory;
}
*/

// and similarly for physics constructors we want to use
#    include "G4PhysicsConstructorFactory.hh"
G4_REFERENCE_PHYSCONSTR_FACTORY(G4NewDecayPhysics);
G4_REFERENCE_PHYSCONSTR_FACTORY_NS(myns::G4NewExoticPhysics, myns, G4NewExoticPhysics);
// note that REFERENCE (vs. DECLARE) doesn't require header to be included
#  endif
#endif

// CommandLineParser copied from test70
// (namespace changed from DNAPARSER to TEST38)
#include "CommandLineParser.hh"
using namespace TEST38;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// some flags to indicate what behaviour we expect
const int STATUS_NONE = 0;
const int STATUS_NEW = 0x01;
const int STATUS_OLD = 0x02;
const int STATUS_BOTH = (STATUS_NEW | STATUS_OLD);

#if (!defined _WIN32) && (defined USE_GETOPT)
// nothing special ... using 'getopt' on non-WIN32 platform
#else
CommandLineParser* parser = 0;
#endif

class G4plfTest;  // forward declaration
std::ostream& operator<<(std::ostream& os, const G4plfTest& atest);

// class G4plfTest holds the results of a single test of the factory/factories
class G4plfTest
{
  public:

    G4plfTest(G4String plname, bool isEnv, int verbose, int factoryVerbose, int expectation);
    ~G4plfTest() { ; }  // delete _physicsNew; delete _physicsOld; }

    int RunTest(bool doOld, bool doNew = true);
    int GetSuccess() const { return _success; }
    int GetExpectation() const { return _expectation; }
    G4String StateAsString() const;

    void UpdateTypeInfo(G4VModularPhysicsList* physics, G4String& mangled, G4String& demangled);
    void PrintTestStart(G4String msg);

  public:

    G4String _plname;  ///< physics list name to be tested
    bool _isEnv;  ///< true _plname is to be set as $PHYLIST
                  ///<    if true & plname="" then unset env variable
    int _verbose;  ///< verbosity when printing
    int _factoryVerbose;  ///< verbosity of the factory
    int _expectation;  ///< expectation new / old bits
    int _success;  ///< success status new / old bits
    G4VModularPhysicsList* _physicsNew;
    G4VModularPhysicsList* _physicsOld;
    G4String _mangledNew;
    G4String _mangledOld;
    G4String _demangledNew;
    G4String _demangledOld;

    static const G4String _smpl;
    static const G4String _snullptr;
    static const G4String _snodemangle;
};
const G4String G4plfTest::_smpl = "G4VModularPhysicsList";
const G4String G4plfTest::_snullptr = "<<null-ptr>>";
const G4String G4plfTest::_snodemangle = "<<no-cxa-demangle>>";
G4plfTest::G4plfTest(G4String plname, bool isEnv, int verbose, int factoryVerbose, int expectation)
  : _plname(plname),
    _isEnv(isEnv),
    _verbose(verbose),
    _factoryVerbose(factoryVerbose),
    _expectation(expectation),
    _success(0),
    _physicsNew(0),
    _physicsOld(0)
{
  ;
}
int G4plfTest::RunTest(bool doOld, bool doNew)
{
  g4alt::G4PhysListFactory factoryNew;
  G4PhysListFactory factoryOld;
  // set this test's factory verbosity
  // new is a singleton, but old isn't so needs initialization

  if (_factoryVerbose != -999)
  {
    factoryNew.SetVerbose(_factoryVerbose);
    factoryOld.SetVerbose(_factoryVerbose);
  }

  if (!doNew) _expectation &= ~STATUS_NEW;  // clear bit if we didn't try
  if (!doOld) _expectation &= ~STATUS_OLD;  // clear bit if we didn't try

  if (!_isEnv)
  {
    if (doNew)
    {
      if (_verbose > 0) PrintTestStart("---NEW GetReferencePhysList");
      _physicsNew = factoryNew.GetReferencePhysList(_plname);
    }
    if (doOld)
    {
      if (_verbose > 0) PrintTestStart("---OLD GetReferencePhysList");
      _physicsOld = factoryOld.GetReferencePhysList(_plname);
    }
  }
  else
  {
    // for env variables
    // grab a copy of what is currently set, so we can restore it
#ifndef _WIN32
    const char* prevEnvVal = std::getenv("PHYSLIST");
    if (_plname == "")
      unsetenv("PHYSLIST");  // force UNsetting ...
    else
      setenv("PHYSLIST", _plname, 1);  // force user value
#else
    G4cout << "_WIN32 alternative to getenv/unsetenv/setenv ..." << G4endl;
    // can't find equivalent to unsetenv, _putenv to "" isn't quite the same
    const char* prevEnvVal = std::getenv("PHYSLIST");
    if (_plname == "")
      _putenv_s("PHYSLIST", "");  // force UNsetting ...
    else
      _putenv_s("PHYSLIST", _plname);  // force user value
#endif
    if (doNew)
    {
      if (_verbose > 0) PrintTestStart("---NEW ReferencePhysList $PHYLIST");
      _physicsNew = factoryNew.ReferencePhysList();
    }
    if (doOld)
    {
      if (_verbose > 0) PrintTestStart("---OLD ReferencePhysList $PHYLIST");
      _physicsOld = factoryOld.ReferencePhysList();
    }
    // restore old value
#ifndef _WIN32
    if (prevEnvVal) setenv("PHYSLIST", prevEnvVal, 1);
#else
    if (prevEnvVal) _putenv_s("PHYSLIST", prevEnvVal);
#endif
  }
  UpdateTypeInfo(_physicsNew, _mangledNew, _demangledNew);
  UpdateTypeInfo(_physicsOld, _mangledOld, _demangledOld);
  if (_physicsNew) _success |= STATUS_NEW;
  if (_physicsOld) _success |= STATUS_OLD;

  int _failure = 0;
  if (_success != _expectation)
  {
    // someone failed to meet expectations
    _failure = (_expectation - _success);
  }
  if (_verbose || _failure)
  {
    G4cout << "--Results begin:" << G4endl;
    if (doNew)
      G4cout << "NEW Requested \"" << _plname << "\" got " << _mangledNew << " (" << _demangledNew
             << ")" << G4endl;
    if (doOld)
      G4cout << "OLD Requested \"" << _plname << "\" got " << _mangledOld << " (" << _demangledOld
             << ")" << G4endl;
    G4cout << "   _success " << _success << "   _expectation " << _expectation << "   _failure "
           << _failure << G4endl;
    if (_failure)
    {
#ifndef _WIN32
      char esc = '\033';  // 0x1B;
      G4String makeRedBack = G4String(1, esc) + G4String("[41m");
      G4String makeWhite = G4String(1, esc) + G4String("[1;37m");
      G4String makeCyan = G4String(1, esc) + G4String("[0;36m");
      G4String makeNormal = G4String(1, esc) + G4String("[0m");
#else
      G4String makeRedBack = "";
      G4String makeWhite = "";
      G4String makeCyan = "";
      G4String makeNormal = "";
#endif
      G4cout << makeCyan << makeRedBack << G4endl << "                                      "
             << G4endl << "   unexpected result:  failure code " << _failure << "  " << G4endl
             << "                                      " << G4endl << makeNormal << G4endl;
    }
    G4cout << "--Results end\n" << G4endl;
  }
  return _failure;
}
G4String G4plfTest::StateAsString() const
{
  std::stringstream s;
  G4String rname = _mangledNew;
  if (_demangledNew != _snodemangle) rname = _demangledNew;
  size_t loc = rname.find(_smpl);
  if (loc != std::string::npos) rname.erase(loc, _smpl.size());
  G4String req = _plname;
  if (_isEnv)
  {
    if (_plname != "")
      req = "$PHYLIST=" + _plname;
    else
      req = "$PHYLIST=<<unset>>";
  }
  s << std::setw(22) << std::left << req
    << " |"
#ifndef _WIN32
    // \033=<esc>  [41m = red background [1,37m = white text
    // [0m = restore no color
    << ((_success == _expectation) ? " okay " : "\033[41m\033[1;37m fail \033[0m")
#else
    // windows doesn't use same escape code for coloring terminal
    << ((_success == _expectation) ? " okay " : " fail ")
#endif
    //<< " " << success[i] << "^" << expect[i]
    << "| " << std::setw(40) << std::left << rname << " | ";

  return s.str();
}
void G4plfTest::UpdateTypeInfo(G4VModularPhysicsList* physics, G4String& mangled,
                               G4String& demangled)
{
  if (!physics)
  {
    mangled = _snullptr;
    demangled = _snullptr;
    return;
  }
  const std::type_info& ti = typeid(*physics);
  mangled = ti.name();
#if defined(HAVE_CXA_DEMANGLE) || defined(__clang__) || defined(__GNUC__)
  int status;
  char* realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
  demangled = realname;
  free(realname);
#else
  demangled = _snodemangle;
#endif
}
void G4plfTest::PrintTestStart(G4String msg)
{
  G4cout << msg << " \"" << _plname << "\"" << G4endl;
}
std::ostream& operator<<(std::ostream& os, const G4plfTest& atest)
{
  if (os.good())
  {
    if (os.tie()) os.tie()->flush();
    os << atest.StateAsString();
  }
  if (os.flags() & std::ios::unitbuf) os.flush();
  return os;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Forward declarations
void PrintUsage();
void ParseArgs(int argc, char** argv, std::vector<G4plfTest>& testlist);
void AddDefaults(std::vector<G4plfTest>& testlist);
void AddEnvTests(std::vector<G4plfTest>& testlist);
int String2Expect(std::string expects);
bool TryManySharedLibraries(const G4String& libbase);
bool LoadSharedLibrary(const G4String& libname, G4String& msg);

std::vector<std::string> tokenizeString(std::string values, std::string sepchar,
                                        bool gVerbose = false);

// Globals (because I'm too lazy to pass these around)
static std::string gPrgname = "g4plftest";
static std::string gPrgnameFull = "g4plftest";
// use this to test env variable (unless overridden by --env flag)
// set to "skip" to not do these two tests
static G4String gPlenv = "QGSP_BERT";
static int gFactoryVerbosity = -999;  // -V flag, if -999 then don't set
static int gPrintFactory = 0;  // -f
static int gPrintFactoryOld = 0;  // -F
static int gPrintCtorList = 0;  // -c
static int gPrintPLRegList = 0;  // -r
static int gAddDefaults = 0;  // -D
static int gLEND = 0;  // --lend
static int gXYZZY = 0;  // --xyzzy
static int gUnknownFatal = 0;  // --fatal
static int gVerbosity = 0;
static int gTestOld = 0;
static int gDefaultExpect = STATUS_BOTH;

static int gExitCode = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  gPrgnameFull = std::string(argv[0]);
  std::size_t lastSlash = gPrgnameFull.find_last_of("/\\");
  if (lastSlash != std::string::npos)
  {
    gPrgname = gPrgnameFull.substr(lastSlash + 1);
  }
  else
  {
    gPrgname = gPrgnameFull;
  }

  G4String sepline =
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    "%%%%%\n";

  // list of PhysLists and whether they should work (0=no, 1=yes)
  std::vector<G4plfTest> testlist;

  ParseArgs(argc, argv, testlist);
  if (testlist.empty() || gAddDefaults) AddDefaults(testlist);
  if (gPlenv != "skip") AddEnvTests(testlist);

  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // Access to registries and factories
  //
  G4PhysicsConstructorRegistry* g4pcr = G4PhysicsConstructorRegistry::Instance();

  G4PhysListRegistry* g4plr = G4PhysListRegistry::Instance();

  g4alt::G4PhysListFactory factoryNew;
  G4PhysListFactory factoryOld;

  if (gFactoryVerbosity != -999)
  {
    factoryNew.SetVerbose(gFactoryVerbosity);
    // old factory isn't a singleton so this doesn't affect the
    // one used by the individual tests
    factoryOld.SetVerbose(gFactoryVerbosity);
  }
  factoryNew.SetUnknownFatal(gUnknownFatal);

  // if requested print the state of the _OLD_ factory
  if (gPrintFactoryOld > 0)
  {
    // old factory doesn't have a print state method
    G4cout << sepline << gPrgname << ": state of the _old_ factory:-" << G4endl << sepline
           << G4endl;
    std::vector<G4String> availBase = factoryOld.AvailablePhysLists();
    std::vector<G4String> availEM = factoryOld.AvailablePhysListsEM();
    G4cout << "old factory base physlists:" << G4endl;
    for (size_t i = 0; i < availBase.size(); ++i)
    {
      G4cout << "  [" << std::setw(2) << i << "] = " << availBase[i] << G4endl;
    }
    G4cout << "old factory EM extensions:" << G4endl;
    for (size_t i = 0; i < availEM.size(); ++i)
    {
      G4cout << "  [" << std::setw(2) << i << "] = " << availEM[i] << G4endl;
    }
    G4cout << G4endl << sepline << G4endl;
  }

  // print the state of the factory (or registry) before loading 2nd library
  if (gPrintCtorList > 1 || gPrintPLRegList > 1)
  {
#ifdef G4TEST38_SHARED
    G4cout << sepline << gPrgname << ": is pre-linked against g4plft1 library, "
           << "  list what is available before loading g4plft2 library" << G4endl << sepline
           << G4endl;
#endif
    if (gPrintCtorList > 1) g4pcr->PrintAvailablePhysicsConstructors();
    if (gPrintPLRegList > 1) g4plr->PrintAvailablePhysLists();
    G4cout << sepline << G4endl;
  }
  else
  {
    size_t npc1 = g4pcr->AvailablePhysicsConstructors().size();
    size_t npl1base = g4plr->AvailablePhysLists().size();
    size_t npl1ext = g4plr->AvailablePhysicsExtensions().size();
    size_t npl1em = g4plr->AvailablePhysListsEM().size();
    G4cout << gPrgname << ": (before) phys ctors " << npc1 << ", phys lists " << npl1base
           << ", ext " << npl1ext << ", em " << npl1em << G4endl;
  }

#if (defined G4TEST38_SHARED) && (defined G4TEST38_TRY_EXTRAS)
#  if (defined _WIN32) && (defined G4TEST38_WIN_NO_PRELINK1)
  // for windows don't prelink library when shared
  TryManySharedLibraries("g4plftest38_lib1");
#  endif

  // load a 2nd library to show run-time extension
  TryManySharedLibraries("g4plftest38_lib2");
#else
  G4cout << "skip attempt to load g4plftest38_lib[1|2]" << G4endl;
#endif

  G4cout << gPrgname << ": now add extensions mappings ALTDK, NEWPHY and XYZZY " << G4endl;

  g4plr->AddPhysicsExtension("ALTDK", "G4NewDecayPhysics");
  g4plr->AddPhysicsExtension("NEWPHY", "myns::G4NewExoticPhysics");

  // this one has no corresponding physics constructor
  // testing whether a mis-register works as expected
  g4plr->AddPhysicsExtension("XYZZY", "NoSuchPhysics");

  // print state of the factory after loading 2nd library
  size_t npc2 = g4pcr->AvailablePhysicsConstructors().size();
  size_t npl2base = g4plr->AvailablePhysLists().size();
  size_t npl2ext = g4plr->AvailablePhysicsExtensions().size();
  size_t npl2em = g4plr->AvailablePhysListsEM().size();
  G4cout << gPrgname << ": (after)  phys ctors " << npc2 << ", phys lists " << npl2base << ", ext "
         << npl2ext << ", em " << npl2em << G4endl;
  G4cout << G4endl;

  if (gPrintCtorList > 0 || gPrintPLRegList > 0)
  {
#ifdef G4TEST38_SHARED
    G4cout << sepline << gPrgname << ": after dynamically loading g4plft2 library" << G4endl
           << sepline << G4endl;
#endif
    if (gPrintCtorList > 0) g4pcr->PrintAvailablePhysicsConstructors();
    if (gPrintPLRegList > 0) g4plr->PrintAvailablePhysLists();
    G4cout << sepline << G4endl;
  }

  if (gPrintFactory > 0)
  {
    G4cout << sepline << gPrgname << ": state of the new factory:-" << G4endl << sepline << G4endl;
    factoryNew.PrintAvailablePhysLists();
    G4cout << sepline << G4endl;
  }

  // do all the tests
  for (size_t i = 0; i < testlist.size(); ++i)
  {
    // accumulate status codes for succes/failure
    if (gVerbosity > 1)
      G4cout << "--------------------------------------------------------------------------"
             << G4endl;
    gExitCode |= testlist[i].RunTest(gTestOld);
  }

  // summary table
  G4cout << G4endl;
  G4cout << "===============================================================================\n"
         << " Results Summary (removed \"" << G4plfTest::_smpl << "\" from names) \n"
         << "===============================================================================\n"
         << "     Requested              |Result| Returned                                 |\n"
         << "----------------------------+------+------------------------------------------+"
         << G4endl;
  for (size_t i = 0; i < testlist.size(); ++i)
  {
    G4cout << "[" << std::setw(2) << i << std::right << "] " << testlist[i] << G4endl;
  }
  G4cout << "----------------------------+------+------------------------------------------+"
         << G4endl;

  // if ( physics ) runManager->SetUserInitialization(physics);

  G4cout << G4endl;

  delete runManager;

  /*
  // Summary of results

  for (size_t i=0; i < pllist.size(); ++i ) {
    G4String rname = plreturned[i];
    size_t loc = rname.find(mpl);
    if ( loc != std::string::npos ) rname.erase(loc,mpl.size());
    G4cout << "[" << std::setw(2) << std::right << i << "] "
           << std::setw(22) << std::left << pllist[i] << " | "
           << ((success[i]==expect[i])?"okay":"fail")
      //<< " " << success[i] << "^" << expect[i]
           << " | "
           << std::setw(40) << std::left << rname << " | "
           << G4endl;
  }
  */

  G4cout << G4endl << "Returning exit code " << gExitCode << G4endl;
  return gExitCode;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void PrintUsage()
{
  G4cout << gPrgname << ":  G4PLFactoryTest - "
         << "a simplified Geant4 app for testing " << G4endl << "the G4PhysListFactory" << G4endl;
  G4cout << "   " << gPrgname << " [options] [physList1 physList2[=N]]" << G4endl;
  G4cout << "  -h --help       this output" << G4endl;
  G4cout << "  -f              print phylist factory status" << G4endl;
  G4cout << "  -F              print old phylist factory availability" << G4endl;
  G4cout << "  -c              print physics ctor list" << G4endl;
  G4cout << "  -r              print physics list registry list" << G4endl
         << "                    repeat to print before adding 2nd library" << G4endl;
  G4cout << "  -v --verbose    increase program verbosity" << G4endl;
  G4cout << "  -V <n>          set factory verbosity" << G4endl;
  G4cout << "  -D --defaults   add default tests even if user supplied tests" << G4endl;
  G4cout << "  -o --old        test old factory" << G4endl;
  G4cout << "  -e --env=PNAME  PhysicsList to use as env variable [" << gPlenv << "]" << G4endl;
  G4cout << "                       use \"skip\" to skip these 2 tests" << G4endl;
  G4cout << "     --lend       try ShieldingLEND (needs special data) in default list" << G4endl;
  G4cout << "     --xyzzy      try to add non-existent physics ctor in default list" << G4endl
         << "                   (will though throw G4Exception w/ --fatal)" << G4endl;
  G4cout << "     --fatal      throw exception if new factory can't satisfy request" << G4endl;
  G4cout << " " << G4endl;
  G4cout << "  If given, the list of physics lists to try override the default set.\n"
         << "  User can specify if they expect each to work with the \n"
         << "    new (" << STATUS_NEW << "), old (" << STATUS_OLD << "), both (" << STATUS_BOTH
         << ") or neither (" << STATUS_NONE << ") factory;\n"
         << "    if unspecified, assumes " << gDefaultExpect << ".\n";
  G4cout << " " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#if (!defined _WIN32) && (defined USE_GETOPT)
// parse the arguments using 'getopt' on non-Windows machines
void ParseArgs(int argc, char** argv, std::vector<G4plfTest>& testlist)
{
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"old", no_argument, 0, 'o'},  // exercise old factory
    {"env", required_argument, 0, 'e'},
    {"defaults", no_argument, 0, 'D'},
    {"Verbose", required_argument, 0, 'V'},  // factory verbosity
    {"status", no_argument, 0, 'f'},  // factory status
    {"oldstatus", no_argument, 0, 'F'},  // old physlist factory list
    {"pctor", no_argument, 0, 'c'},  // physics ctor registry list
    {"reglist", no_argument, 0, 'r'},  // pl registry list (up to 2x)
    {"lend", no_argument, 0, 'l'},
    {"xyzzy", no_argument, 0, 'x'},
    {"fatal", no_argument, 0, 'E'},
    // dont' forget to add short char to getopt_long list below
    {0, 0, 0, 0}  // signal end
  };
// getopt_long stored the option index here
#  ifndef MACOSX
  optind = 0;  // getopt.h: Reset getopt to start of arguments
#  else
  optind = 1;  // skip 0th argument (executable) for MACOSX
#  endif

  int option_index = 0;
  int c;
  while (true)
  {
    c = getopt_long(argc, argv, "hvoe:DV:fFcrlxE", long_options, &option_index);

    if (c == -1) break;  // detect end of options

    switch (c)
    {
      case 0: {
        // if option set a flag, do nothing else
        if (long_options[option_index].flag != 0) break;

        G4String oname = long_options[option_index].name;
        if (gVerbosity > 3)
        {
          G4cout << " oname=\"" << oname << "\" optarg=\"" << optarg << "\"" << G4endl;
        }
        else if (oname == "env")
          gPlenv = optarg;
        break;
      }
      case 'h':
        PrintUsage();
        exit(0);
      case 'v':
        gVerbosity++;
        break;
      case 'o':
        gTestOld++;
        break;
      case 'f':
        gPrintFactory++;
        break;
      case 'F':
        gPrintFactoryOld++;
        break;
      case 'c':
        gPrintCtorList++;
        break;
      case 'r':
        gPrintPLRegList++;
        break;
      case 'e':
        // G4cout << " -e optarg=\"" << optarg << "\"" << G4endl;
        gPlenv = optarg;
        break;
      case 'D':
        gAddDefaults++;
        break;
      case 'V':
        // G4cout << " -V optarg=\"" << optarg << "\"" << G4endl;
        gFactoryVerbosity = atoi(optarg);
        break;
      case 'l':
        gLEND++;
        break;
      case 'x':
        gXYZZY++;
        break;
      case 'E':
        gUnknownFatal++;
        break;
      case '?':
        // getopt_long already printed an error message
        PrintUsage();
        break;
      default:
        G4cout << "unknown opt=\"" << char(c) << "\""
               << " optarg=\"" << optarg << "\"" << G4endl;
        PrintUsage();
        abort();
    }
  }

  for (int index = optind; index < argc; ++index)
  {
    std::vector<std::string> v = tokenizeString(argv[index], "=");
    std::string plname = v[0];
    int expecti = gDefaultExpect;
    if (v.size() > 1) expecti = String2Expect(v[1]);
    if (gVerbosity > 3)
    {
      G4cout << " [" << index << "]=\"" << argv[index] << "\" : plname=\"" << plname
             << "\" expecti " << expecti << G4endl;
    }
    testlist.push_back(G4plfTest(plname, false, gVerbosity, gFactoryVerbosity, expecti));
  }
}
#else
// using CommandLineParser for Windows (or if not 'getopt')
void ParseArgs(int argc, char** argv, std::vector<G4plfTest>& testlist)
{
  G4cout << gPrgname << " using CommandLineParser" << G4endl;
  parser = CommandLineParser::GetParser();

  parser->AddCommand("--verbose", Command::OptionNotCompulsory, "increase verbosity [N]", "1");
  parser->AddCommand("-v", Command::OptionNotCompulsory, "increase verbosity [N]", "1");
  // parser->AddCommand("&",Command::OptionNotCompulsory,"increase verbosity [N]","1");

  parser->AddCommand("--old", Command::WithoutOption, "exercise old factory");
  parser->AddCommand("-o", Command::WithoutOption, "exercise old factory");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--env", Command::WithOption, "exercise env variable usage", "QGSP");
  parser->AddCommand("-e", Command::WithOption, "exercise env variable usage", "QGSP");
  // parser->AddCommand("&",Command::WithOption,"exercise env variable usage","QGSP");

  parser->AddCommand("--defaults", Command::WithoutOption, "add default set");
  parser->AddCommand("-D", Command::WithoutOption, "add default set");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--Verbose", Command::WithOption, "factory verbosity", "1");
  parser->AddCommand("-V", Command::WithOption, "factory verbosity", "1");
  // parser->AddCommand("&",Command::WithOption,"factory verbosity","1");

  parser->AddCommand("--status", Command::WithoutOption, "factory status");
  parser->AddCommand("-f", Command::WithoutOption, "factory status");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--oldstatus", Command::WithoutOption, "old physlist factory list");
  parser->AddCommand("-F", Command::WithoutOption, "old physlist factory list");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--pctor", Command::WithoutOption, "physics ctor registry list");
  parser->AddCommand("-c", Command::WithoutOption, "physics ctor registry list");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--reglist", Command::OptionNotCompulsory, "pl registry list [N]", "1");
  parser->AddCommand("-r", Command::OptionNotCompulsory, "pl registry list [N]", "1");
  // parser->AddCommand("&",Command::OptionNotCompulsory,"pl registry list [N]","1");

  parser->AddCommand("--lend", Command::WithoutOption, "try LEND physlist");
  parser->AddCommand("-l", Command::WithoutOption, "try LEND physlist");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--xyzzy", Command::WithoutOption, "try non-existent list");
  parser->AddCommand("-x", Command::WithoutOption, "try non-existent list");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--fatal", Command::WithoutOption, "make error fatal");
  parser->AddCommand("-E", Command::WithoutOption, "make error fatal");
  // parser->AddCommand("&",Command::WithoutOption);

  parser->AddCommand("--usage", Command::WithoutOption, "get usage text");
  parser->AddCommand("-u", Command::WithoutOption, "get usage text");
  // parser->AddCommand("&",Command::WithoutOption);

  // int pstat =
  parser->Parse(argc, argv);
  // G4cout << gPrgname << " CommandLineParser returned " << pstat << G4endl;

  Command* cmdLong(0);
  Command* cmdShort(0);

  if ((cmdLong = parser->GetCommandIfActive("--verbose"))
      || (cmdShort = parser->GetCommandIfActive("-v")))
  {
    if (cmdLong)
    {
      if (cmdLong->GetOption() == "")
        gVerbosity++;
      else
      {
        gVerbosity += atoi(cmdLong->GetOption().c_str());
      }
    }
    if (cmdShort)
    {
      if (cmdShort->GetOption() == "")
        gVerbosity++;
      else
      {
        gVerbosity += atoi(cmdShort->GetOption().c_str());
      }
    }
  }

  if ((cmdLong = parser->GetCommandIfActive("--old"))
      || (cmdShort = parser->GetCommandIfActive("-o")))
  {
    gTestOld++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--env"))
      || (cmdShort = parser->GetCommandIfActive("-e")))
  {
    if (cmdLong) gPlenv = cmdLong->GetOption();
    if (cmdShort) gPlenv = cmdShort->GetOption();
  }

  if ((cmdLong = parser->GetCommandIfActive("--defaults"))
      || (cmdShort = parser->GetCommandIfActive("-D")))
  {
    gAddDefaults++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--Verbose"))
      || (cmdShort = parser->GetCommandIfActive("-V")))
  {
    gFactoryVerbosity = 0;
    if (cmdLong)
    {
      if (cmdLong->GetOption() == "")
        gFactoryVerbosity++;
      else
      {
        gFactoryVerbosity = atoi(cmdLong->GetOption().c_str());
      }
    }
    if (cmdShort)
    {
      if (cmdShort->GetOption() == "")
        gFactoryVerbosity++;
      else
      {
        gFactoryVerbosity = atoi(cmdShort->GetOption().c_str());
      }
    }
  }

  if ((cmdLong = parser->GetCommandIfActive("--status"))
      || (cmdShort = parser->GetCommandIfActive("-f")))
  {
    gPrintFactory++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--oldstatus"))
      || (cmdShort = parser->GetCommandIfActive("-F")))
  {
    gPrintFactoryOld++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--pctor"))
      || (cmdShort = parser->GetCommandIfActive("-c")))
  {
    gPrintCtorList++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--reglist"))
      || (cmdShort = parser->GetCommandIfActive("-r")))
  {
    if (cmdLong)
    {
      if (cmdLong->GetOption() == "")
        gPrintPLRegList++;
      else
      {
        gPrintPLRegList += atoi(cmdLong->GetOption().c_str());
      }
    }
    if (cmdShort)
    {
      if (cmdShort->GetOption() == "")
        gPrintPLRegList++;
      else
      {
        gPrintPLRegList += atoi(cmdShort->GetOption().c_str());
      }
    }
  }

  if ((cmdLong = parser->GetCommandIfActive("--lend"))
      || (cmdShort = parser->GetCommandIfActive("-l")))
  {
    gLEND++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--xyzzy"))
      || (cmdShort = parser->GetCommandIfActive("-x")))
  {
    gXYZZY++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--fatal"))
      || (cmdShort = parser->GetCommandIfActive("-E")))
  {
    gUnknownFatal++;
  }

  if ((cmdLong = parser->GetCommandIfActive("--usage"))
      || (cmdShort = parser->GetCommandIfActive("-u")))
  {
    PrintUsage();
    exit(0);
  }

  // G4cout << gPrgname << " CommandLineParser left " << argc
  //        << " options" << G4endl;

  int optindCLP = 0;
  for (int index = optindCLP; index < argc; ++index)
  {
    // program name is not a physics list to try ... skip it
    if (argv[index] == gPrgname) continue;
    if (argv[index] == gPrgnameFull) continue;

    std::vector<std::string> v = tokenizeString(argv[index], "=");
    std::string plname = v[0];
    int expecti = gDefaultExpect;
    if (v.size() > 1) expecti = String2Expect(v[1]);
    if (gVerbosity > 3)
    {
      G4cout << " [" << index << "]=\"" << argv[index] << "\" : plname=\"" << plname
             << "\" expecti " << expecti << G4endl;
    }
    testlist.push_back(G4plfTest(plname, false, gVerbosity, gFactoryVerbosity, expecti));
  }

  CommandLineParser::DeleteInstance();
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AddDefaults(std::vector<G4plfTest>& testlist)
{
  int v1 = gVerbosity;
  int v2 = gFactoryVerbosity;

  // some standard physics lists
  testlist.push_back(G4plfTest("FTFP_BERT", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("FTFP_BERT_HP", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("NuBeam", false, v1, v2, STATUS_BOTH));
  if (gLEND)
  {
    testlist.push_back(G4plfTest("ShieldingLEND", false, v1, v2, STATUS_BOTH));
  }
  testlist.push_back(G4plfTest("ShieldingM", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("QGSP_INCLXX", false, v1, v2, STATUS_BOTH));
  testlist.push_back(
    G4plfTest("G4GenericPhysicsList", false, v1, v2, STATUS_NEW));  //! NOTE old doesn't know

  // test the standard EM extensions
  testlist.push_back(G4plfTest("FTFP_BERT_EMV", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("FTFP_BERT_EMX", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("FTFP_BERT_EMY", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("FTFP_BERT_EMZ", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("FTFP_BERT_LIV", false, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest("FTFP_BERT_PEN", false, v1, v2, STATUS_BOTH));
  // yes .. double underscores in original factory ...
  testlist.push_back(G4plfTest("QGSP_BIC_AllHP__GS", false, v1, v2, STATUS_BOTH));

  // lists created for this test
  // MyPL0 = compiled w/ main program
  // MyPL1 = in separate library, linked at build time
  // MyPL2 = in separate library loaded at run time
  //    currently don't try run time loading for _WIN32
  // myns::MyNSPL3 = like MyPL2 but within a namespace
  testlist.push_back(G4plfTest("MyPL0", false, v1, v2, STATUS_NEW));

#ifdef G4TEST38_STATIC
  G4cout << "G4TEST38_STATIC enabled ... " << G4endl;
  // LoadSharedLibrary(), if run, might report success at loading ...
  // but new physlists and physic ctor's don't actually show up
  // these should have been pulled in via G4_REFERENCE_PHYSLIST_FACTORY
  // and G4_REFERENCE_PHYSCONSTR_FACTORY
#endif
#ifdef G4TEST38_SHARED
  G4cout << "G4TEST38_SHARED enabled ... " << G4endl;
#endif

#ifdef G4TEST38_TRY_EXTRAS
  G4cout << "G4TEST38_TRY_EXTRAS enabled ... " << G4endl;

  testlist.push_back(G4plfTest("MyPL1", false, v1, v2, STATUS_NEW));
  testlist.push_back(G4plfTest("MyPL2", false, v1, v2, STATUS_NEW));
  testlist.push_back(G4plfTest("myns::MyNSPL3", false, v1, v2, STATUS_NEW));

  // test the extensions added in 'g4plft2'
  testlist.push_back(G4plfTest("MyPL0+NEWPHY+ALTDK_EMV", false, v1, v2, STATUS_NEW));
#else
  G4cout << "G4TEST38_TRY_EXTRAS disabled ... " << G4endl;
#endif

  if (gXYZZY)
  {
    // test case where mapping extension defined,
    // but no underlying physics process exists
    // we _expect_ this to fail ... so failing is okay (i.e. STATUS_NONE)
    testlist.push_back(G4plfTest("MyPL0+XYZZY", false, v1, v2, STATUS_NONE));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AddEnvTests(std::vector<G4plfTest>& testlist)
{
  int v1 = gVerbosity;
  int v2 = gFactoryVerbosity;

  testlist.push_back(G4plfTest("", true, v1, v2, STATUS_BOTH));
  testlist.push_back(G4plfTest(gPlenv, true, v1, v2, STATUS_BOTH));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int String2Expect(std::string expects)
{
  std::transform(expects.begin(), expects.end(), expects.begin(), ::tolower);

  if (expects == "none") return STATUS_NONE;
  if (expects == "neither") return STATUS_NONE;
  if (expects == "new") return STATUS_NEW;
  if (expects == "old") return STATUS_OLD;
  if (expects == "both") return STATUS_BOTH;

  char c1 = expects[0];
  char s2[2] = {c1, '\0'};
  int statnum = atoi(s2);
  if (statnum == STATUS_NONE) return STATUS_NONE;
  if (statnum == STATUS_NEW) return STATUS_NEW;
  if (statnum == STATUS_OLD) return STATUS_OLD;
  if (statnum == STATUS_BOTH) return STATUS_BOTH;

  // everything else use default
  return gDefaultExpect;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4TEST38_SHARED
bool TryManySharedLibraries(const G4String& libbase)
{
  // compiler should be passed something like
  //    -DDYNAMIC_LOADER_PREFIX=\"lib\"
  //    -DDYNAMIC_LOADER_SUFFIX=\".so\"

  // on Mac CMAKE_SHARED_MODULE_SUFFIX _reports_ .so .... but
  // running it makes a .dylib ... so, just try all possible combinations
  //
  //  std::string libname = "";
  //  libname += DYNAMIC_LOADER_PREFIX;
  //  libname += libbase;
  //  libname += DYNAMIC_LOADER_SUFFIX;

  bool success = false;

  const char* prefixes[] = {"lib", "l", ""};
  const char* suffixes[] = {".so", ".dylib", ".dll", ".lib", ".a"};

  const int npre = sizeof(prefixes) / sizeof(const char*);
  const int nsuf = sizeof(suffixes) / sizeof(const char*);

  std::ostringstream errorMsgs;

  for (int ipre = 0; ipre < npre; ++ipre)
  {
    for (int isuf = 0; isuf < nsuf; ++isuf)
    {
      G4String libname = "";
      libname += prefixes[ipre];
      libname += libbase;
      libname += suffixes[isuf];

      G4String msg = "";

      success = LoadSharedLibrary(libname, msg);
      if (success)
      {
        G4cout << gPrgname << ": dynamically loaded '" << libname << "'" << G4endl;
        return true;
      }
      else
      {
        errorMsgs << msg << std::endl;
      }
    }
  }

  G4cout << gPrgname << ": tried to dynamically load '" << libbase << "'"
         << " but found no SUFFIX/PREFIX combination that worked" << G4endl << G4endl
         << errorMsgs.str() << G4endl;

  return false;
}

bool LoadSharedLibrary(const G4String& libname, G4String& msg)
{
#  ifndef _WIN32
  // basic UNIX style loading
  void* handle = dlopen(libname.c_str(), RTLD_LAZY);
  if (!handle)
  {
    std::ostringstream ssmsg;
    ssmsg << gPrgname << ": cannot open library: '" << libname << "'"
          << " due to: " << dlerror() << std::endl;
    msg = ssmsg.str();
    return false;
  }

#  else
  // Windows ...
  G4cout << gPrgname << ": WIN32 dynamic library loading try '" << libname << "'" << G4endl;

  HMODULE lh;
  int length = MultiByteToWideChar(CP_UTF8, 0, libname.c_str(), -1, NULL, 0);
  wchar_t* wchars = new wchar_t[length + 1];
  wchars[0] = '\0';
  MultiByteToWideChar(CP_UTF8, 0, libname.c_str(), -1, wchars, length);
  lh = LoadLibraryW(wchars);
  delete[] wchars;
  std::ostringstream ssmsg;

  if (!lh)
  {
    ssmsg << gPrgname << ": cannot open library: '" << libname << "'" << std::endl;

    LPVOID lpMsgBuf = NULL;

    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM, NULL, GetLastError(),
                  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),  // Default language
                  (LPTSTR)&lpMsgBuf, 0, NULL);

    if (!lpMsgBuf)
    {
      ssmsg << gPrgname << ": no lpMsgBuf " << std::endl;
      return 0;
    }

    static char* str = 0;
    delete[] str;
    str = strcpy(new char[strlen((char*)lpMsgBuf) + 1], (char*)lpMsgBuf);
    // Free the buffer.
    LocalFree(lpMsgBuf);

    ssmsg << gPrgname << ": load message: '" << str << "'" << std::endl;

    msg = ssmsg.str();
    return 0;
  }
#  endif  // _WIN32

  G4cout << gPrgname << ": successfully loaded '" << libname << "'" << G4endl;
  return true;
}
#endif  // G4TEST38_SHARED

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<std::string> tokenizeString(std::string values, std::string sepchar, bool gVerbose)
{
  // Separate "values" string into elements under the assumption
  // that they are separated by any of the characters in "spechar".

  std::vector<std::string> rlist;

  if (gVerbose) G4cout << "values " << values << " separated by \"" << sepchar << "\"" << G4endl;

  size_t pos_beg = 0;
  size_t str_end = values.size();
  while (pos_beg != std::string::npos && pos_beg < str_end)
  {
    size_t pos_end = values.find_first_of(sepchar.c_str(), pos_beg);
    std::string onevalue = values.substr(pos_beg, pos_end - pos_beg);
    if (gVerbose)
      G4cout << " onevalue \"" << onevalue << "\" in [" << pos_beg << "," << pos_end << ")"
             << G4endl;
    pos_beg = pos_end + 1;
    if (pos_end == std::string::npos) pos_beg = str_end;
    if (onevalue != "")
    {
      rlist.push_back(onevalue);
    }
  }

  return rlist;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
