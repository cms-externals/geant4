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
// G4FPEDetection
//
// Description:
//
// This global method should be used on LINUX/gcc or MacOS/clang platforms
// for activating NaN detection and FPE signals, and forcing abortion of
// the application at the time these are detected.
// Meant to be used for debug purposes, can be activated by compiling the
// "run" module with the flag G4FPE_DEBUG set in the environment.

// Author: G.Cosmo, 14 July 2010 - First version
// --------------------------------------------------------------------
#ifndef G4FPEDETECTION_HH
#define G4FPEDETECTION_HH

#include <cfenv>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>

// Linux/glibc only appear to support FPE exception traps on non-ARM
// e.g.: https://sourceware.org/glibc/wiki/Release/2.27#AArch64
//     : https://gcc.gnu.org/pipermail/gcc-patches/2024-February/645667.html
// This does appear to be somewhat chip dependent, but we disable explicitly as
// it cannot be guaranteed.
// Workaround suggestions include using undefined behaviour sanitizers:
// - https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
#if (defined(__linux__) && !defined(__aarch64__))
#  if (defined(__GNUC__) && !defined(__clang__))
#    include <features.h>
// for G4StackBacktrace()
#    include <cxxabi.h>
#    include <execinfo.h>

#    define G4FPE_SIGNAL SIGFPE  // Linux/GCC FPE signal

/**
 * Print a demangled call stack to stderr.
 *
 * Uses @c backtrace() and @c backtrace_symbols() to capture the call stack
 * and @c abi::__cxa_demangle() to convert mangled C++ symbol names back to
 * human-readable form. Only available on Linux/GCC.
 */
static void G4StackBackTrace()
{
//   from http://linux.die.net/man/3/backtrace_symbols_fd
#    define BSIZE 50
  void* buffer[BSIZE];
  int nptrs = backtrace(buffer, BSIZE);
  char** strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL)
  {
    perror("backtrace_symbols");
    return;
  }
  std::cerr << std::endl << "Call Stack:" << std::endl;
  for (int j = 0; j < nptrs; j++)
  {
    std::cerr << nptrs - j - 1 << ": ";
    char* mangled_start = std::strchr(strings[j], '(') + 1;
    if (mangled_start) *(mangled_start - 1) = '\0';
    char* mangled_end = std::strchr(mangled_start, '+');
    if (mangled_end) *mangled_end = '\0';
    int status = 0;
    char* realname = 0;
    if (mangled_end && std::strlen(mangled_start))
      realname = abi::__cxa_demangle(mangled_start, 0, 0, &status);
    if (realname)
    {
      std::cerr << strings[j] << " : " << realname << std::endl;
      free(realname);
    }
    else
    {
      std::cerr << strings[j] << std::endl;
    }
  }
  free(strings);
  // c++filt can demangle:
  // http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
  //-------------------------------------------------------------------
}

/**
 * Return the @c FPE_* si_code for a Linux/GCC SIGFPE signal.
 *
 * On Linux, the kernel delivers SIGFPE with @c siginfo_t::si_code already set
 * to the appropriate @c FPE_FLT* / @c FPE_INT* constant, so no additional
 * decoding is needed.
 *
 * @param sinfo  Signal info structure provided by the signal handler.
 * @return @c sinfo->si_code, or @c -1 if @p sinfo is null.
 */
static int DecodeFPECode(int /*sig*/, siginfo_t* sinfo, void* /*context*/)
{
  return sinfo ? sinfo->si_code : -1;
}
#  endif /* GCC */

#elif defined(__MACH__) /* MacOS */
#  include <sys/ucontext.h>

#  include <cstdint>

#  if (defined(__i386__) || defined(__x86_64__))  // INTEL

#    define G4FPE_SIGNAL SIGFPE  //! Intel FPE signal

/**
 * Unmask (enable trapping for) the given floating-point exceptions on x86/x86-64 macOS.
 *
 * macOS does not expose @c feenableexcept() in its C library, so this
 * implementation manipulates the x87 control word and SSE MXCSR register
 * directly via @c fegetenv() / @c fesetenv().
 *
 * @param excepts  Bitwise OR of @c FE_* exception flags (e.g. @c FE_DIVBYZERO).
 * @return Previous exception mask, or @c -1 on error.
 */
static inline int feenableexcept(unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if (fegetenv(&fenv))
  {
    return -1;
  }
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  //
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr &= ~(new_excepts << 7);

  return (fesetenv(&fenv) ? -1 : old_excepts);
}

/**
 * Mask (disable trapping for) the given floating-point exceptions on x86/x86-64 macOS.
 *
 * Counterpart to the local @c feenableexcept(); masks exception traps by
 * setting the corresponding bits in the x87 control word and MXCSR register.
 *
 * @param excepts  Bitwise OR of @c FE_* exception flags to disable.
 * @return Previous exception mask, or @c -1 on error.
 */
static inline int fedisableexcept(unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // all previous masks

  if (fegetenv(&fenv))
  {
    return -1;
  }
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // mask
  //
  fenv.__control |= new_excepts;
  fenv.__mxcsr |= new_excepts << 7;

  return (fesetenv(&fenv) ? -1 : old_excepts);
}

/**
 * Return the @c FPE_* si_code for an Intel macOS SIGFPE signal.
 *
 * On x86/x86-64 macOS, the kernel delivers SIGFPE with @c siginfo_t::si_code
 * already set to the appropriate @c FPE_FLT* constant, so no additional
 * decoding is needed.
 *
 * @param sinfo  Signal info structure provided by the signal handler.
 * @return @c sinfo->si_code, or @c -1 if @p sinfo is null.
 */
static int DecodeFPECode(int /*sig*/, siginfo_t* sinfo, void* /*scp*/)
{
  return sinfo ? sinfo->si_code : -1;
}

#  elif (defined(__arm) || defined(__arm64) || defined(__aarch64__))  // Apple Silicon

#    define G4FPE_SIGNAL SIGILL  //! Apple Silicon FPE signal

/**
 * Maps @c FE_* exception flags to the corresponding @c __fpcr_trap_* enable bits.
 *
 * Used by @c feenableexcept() and @c fedisableexcept() to translate the
 * standard @c FE_* flags into the AArch64 FPCR trap-enable bit positions.
 *
 * References:
 * - https://stackoverflow.com/questions/69059981/how-to-trap-floating-point-exceptions-on-m1-macs
 * -
 * https://stackoverflow.com/questions/71821666/trapping-floating-point-exceptions-and-signal-handling-on-apple-silicon
 *
 * @param excepts  Bitwise OR of @c FE_* exception flags.
 * @return Corresponding bitwise OR of @c __fpcr_trap_* bits.
 */
static inline std::uint32_t G4FPEExceptToFPCRMask(unsigned int excepts)
{
  // clang-format off
  std::uint32_t mask = 0;
  if (excepts & FE_INVALID)   mask |= __fpcr_trap_invalid;
  if (excepts & FE_DIVBYZERO) mask |= __fpcr_trap_divbyzero;
  if (excepts & FE_OVERFLOW)  mask |= __fpcr_trap_overflow;
  if (excepts & FE_UNDERFLOW) mask |= __fpcr_trap_underflow;
  if (excepts & FE_INEXACT)   mask |= __fpcr_trap_inexact;
  // clang-format on
  return mask;
}

/**
 * Maps @c __fpcr_trap_* enable bits back to the corresponding @c FE_* exception flags.
 *
 * Inverse of @c G4FPEExceptToFPCRMask(); used to reconstruct the previous
 * exception mask for the glibc-compatible return value of @c feenableexcept()
 * and @c fedisableexcept().
 *
 * @param fpcr  Value of the @c __fpcr field from @c fenv_t.
 * @return Bitwise OR of @c FE_* flags for each trap bit that is set.
 */
static inline unsigned int G4FPCRMaskToFPEExcept(std::uint32_t fpcr)
{
  // clang-format off
  unsigned int excepts = 0;
  if (fpcr & __fpcr_trap_invalid)   excepts |= FE_INVALID;
  if (fpcr & __fpcr_trap_divbyzero) excepts |= FE_DIVBYZERO;
  if (fpcr & __fpcr_trap_overflow)  excepts |= FE_OVERFLOW;
  if (fpcr & __fpcr_trap_underflow) excepts |= FE_UNDERFLOW;
  if (fpcr & __fpcr_trap_inexact)   excepts |= FE_INEXACT;
  // clang-format on
  return excepts;
}

/**
 * Unmask (enable trapping for) the given floating-point exceptions on Apple Silicon.
 *
 * macOS does not expose @c feenableexcept() in its C library. On AArch64 the
 * trap-enable bits live in the AArch64 FPCR register, exposed via the private
 * @c __fpcr field of @c fenv_t. @c G4FPEExceptToFPCRMask() maps standard
 * @c FE_* flags to the corresponding @c __fpcr_trap_* bits.
 *
 * @param excepts  Bitwise OR of @c FE_* exception flags (e.g. @c FE_DIVBYZERO).
 * @return Previous set of enabled exceptions as a bitwise OR of @c FE_* flags,
 *         or @c -1 on error.
 */
static inline int feenableexcept(unsigned int excepts)
{
  std::fenv_t fenv;
  if (std::fegetenv(&fenv)) return -1;
  unsigned int old_excepts = G4FPCRMaskToFPEExcept(fenv.__fpcr);
  fenv.__fpcr |= G4FPEExceptToFPCRMask(excepts);
  return std::fesetenv(&fenv) ? -1 : static_cast<int>(old_excepts);
}

/**
 * Mask (disable trapping for) the given floating-point exceptions on Apple Silicon.
 *
 * Counterpart to the local @c feenableexcept(); clears the corresponding
 * @c __fpcr_trap_* bits in the AArch64 FPCR register so that the named
 * exceptions no longer trap.
 *
 * @param excepts  Bitwise OR of @c FE_* exception flags to disable.
 * @return Previous set of enabled exceptions as a bitwise OR of @c FE_* flags,
 *         or @c -1 on error.
 */
static inline int fedisableexcept(unsigned int excepts)
{
  std::fenv_t fenv;
  if (std::fegetenv(&fenv)) return -1;
  unsigned int old_excepts = G4FPCRMaskToFPEExcept(fenv.__fpcr);
  fenv.__fpcr &= ~G4FPEExceptToFPCRMask(excepts);
  return std::fesetenv(&fenv) ? -1 : static_cast<int>(old_excepts);
}

/**
 * Return FPE_FLT* code corresponding to AArch64 ESR
 *
 * On Apple Silicon, FPEs trap as SIGILL+ILL_ILLTRP, not SIGFPE.
 * Decode the AArch64 ESR and return the equivalent FPE_FLT* si_code.
 * ESR_EL1: EC[31:26]==0x2C = FPU trapping op; TFV(bit 23) gates flags[6:0].
 * References:
 * - https://stackoverflow.com/questions/69059981/how-to-trap-floating-point-exceptions-on-m1-macs
 * -
 * https://stackoverflow.com/questions/71821666/trapping-floating-point-exceptions-and-signal-handling-on-apple-silicon
 */
static int DecodeFPECode(int sig, siginfo_t* sinfo, void* scp)
{
  if (!sinfo || sig != SIGILL || sinfo->si_code != ILL_ILLTRP) return -1;
  auto* uctx = static_cast<ucontext_t*>(scp);
  if (!uctx) return -1;
  std::uint32_t esr = uctx->uc_mcontext->__es.__esr;
  if (((esr >> 26) & 0x3F) != 0x2C) return -1;  // not an FPU trap
  if (!((esr >> 23) & 1)) return -1;  // TFV not set
  std::uint32_t flags = esr & 0x7F;
  if (flags & (1 << 0)) return FPE_FLTINV;  // IOC
  if (flags & (1 << 1)) return FPE_FLTDIV;  // DZC
  if (flags & (1 << 2)) return FPE_FLTOVF;  // OFC
  if (flags & (1 << 3)) return FPE_FLTUND;  // UFC
  if (flags & (1 << 4)) return FPE_FLTRES;  // IXC
  return -1;
}

#  endif /* INTEL or Apple Silicon enabling */

/**
 * No-op stack backtrace on macOS.
 *
 * @c backtrace() / @c execinfo.h are not used on macOS. The
 * @c TerminationSignalHandler() calls this unconditionally, so an empty
 * stub is provided here to satisfy that call site.
 *
 * @todo Consider implementing via @c G4BackTrace if available.
 */
static void G4StackBackTrace() {}

#endif /* Linux or MacOS */

// Enable FPE handling if platform/compiler has defined G4FPE_SIGNAL
#if defined(G4FPE_SIGNAL)
/**
 * Signal handler that reports the FPE type and aborts the process.
 *
 * Installed by @c InvalidOperationDetection() on @c G4FPE_SIGNAL
 * (@c SIGFPE on Linux/Intel macOS, @c SIGILL on Apple Silicon).
 * Calls @c DecodeFPECode() to normalise the platform-specific signal
 * information to a standard @c FPE_FLT* code, prints a human-readable
 * message to stderr, optionally prints a stack backtrace, then calls
 * @c std::abort().
 *
 * @param sig      Signal number received.
 * @param sinfo    Pointer to @c siginfo_t provided by the OS.
 * @param context  Opaque pointer to the platform @c ucontext_t; forwarded
 *                 to @c DecodeFPECode() for register decoding on AArch64.
 */
static void TerminationSignalHandler(int sig, siginfo_t* sinfo, void* context)
{
  std::cerr << "ERROR: " << sig;
  std::string message = "Floating-point exception (FPE).";

  switch (DecodeFPECode(sig, sinfo, context))
  {
#  ifdef FPE_NOOP /* macOS uses this instead of FPE_INTDIV */
    case FPE_NOOP:
#  endif
    case FPE_INTDIV:
      message = "Integer divide by zero.";
      break;
    case FPE_INTOVF:
      message = "Integer overflow.";
      break;
    case FPE_FLTDIV:
      message = "Floating point divide by zero.";
      break;
    case FPE_FLTOVF:
      message = "Floating point overflow.";
      break;
    case FPE_FLTUND:
      message = "Floating point underflow.";
      break;
    case FPE_FLTRES:
      message = "Floating point inexact result.";
      break;
    case FPE_FLTINV:
      message = "Floating point invalid operation.";
      break;
    case FPE_FLTSUB:
      message = "Subscript out of range.";
      break;
    default:
      message = "Unknown error.";
      break;
  }

  std::cerr << " - " << message << std::endl;
  G4StackBackTrace();
  std::abort();
}

/**
 * Enable FPE trapping and install the termination signal handler.
 *
 * Enables trapping for @c FE_DIVBYZERO and @c FE_INVALID via
 * @c feenableexcept(), then installs @c TerminationSignalHandler on
 * @c G4FPE_SIGNAL using @c sigaction(). On Linux and Intel macOS the
 * trapped signal is @c SIGFPE; on Apple Silicon it is @c SIGILL.
 *
 * Overflow and underflow trapping are intentionally left disabled by
 * default but can be enabled by uncommenting the corresponding
 * @c feenableexcept() calls below.
 *
 * @note Intended for debug builds only. Activate by compiling the
 *       @c run module with @c G4FPE_DEBUG defined in the environment.
 */
static void InvalidOperationDetection()
{
  struct sigaction termaction, oldaction;

  std::cout << std::endl
            << "        "
            << "############################################" << std::endl
            << "        "
            << "!!! WARNING - FPE detection is activated !!!" << std::endl
            << "        "
            << "############################################" << std::endl;
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO);
  feenableexcept(FE_INVALID);
  // fedisableexcept( FE_OVERFLOW  );
  // fedisableexcept( FE_UNDERFLOW );

  sigfillset(&termaction.sa_mask);
  sigdelset(&termaction.sa_mask, G4FPE_SIGNAL);
  termaction.sa_sigaction = TerminationSignalHandler;
  termaction.sa_flags = SA_SIGINFO;
  sigaction(G4FPE_SIGNAL, &termaction, &oldaction);
}
#else
static void InvalidOperationDetection() {}
#endif /* G4FPE_SIGNAL */

#endif /* G4FPEDetection_hh */
