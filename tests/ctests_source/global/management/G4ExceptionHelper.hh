#ifndef __G4ExceptionHelper__
#define __G4ExceptionHelper__

#include <G4ExceptionSeverity.hh>
#include <G4StateManager.hh>
#include <G4VExceptionHandler.hh>

namespace G4ExceptionHelper
{

// Subclassing std::exception is not necessarily correct, but this allows us to
//  throw with a severity level and message easily
struct wrappedException : public std::runtime_error
{
    G4ExceptionSeverity severity;
    wrappedException(char const* exceptionMsg, G4ExceptionSeverity sev)
      : std::runtime_error(exceptionMsg)
    {
      severity = sev;
    }
};

template<G4ExceptionSeverity sev>
struct specifiedException : public std::runtime_error
{
    specifiedException(char const* exceptionMsg) : std::runtime_error(exceptionMsg) { ; }
};

class TestExceptionHandler final : public G4VExceptionHandler
{
  public:

    TestExceptionHandler()
    {
      // Forcibly insert this exception handler
      G4StateManager::GetStateManager()->SetExceptionHandler(this);
    }
    G4bool Notify(char const* originOfException, char const* exceptionCode,
                  G4ExceptionSeverity severity, char const* description) final;
};

G4bool TestExceptionHandler::Notify([[maybe_unused]] char const* originOfException,
                                    char const* exceptionCode, G4ExceptionSeverity severity,
                                    [[maybe_unused]] char const* description)
{
  switch (severity)
  {
    case FatalException:
      throw G4ExceptionHelper::specifiedException<FatalException>(exceptionCode);
      break;
    case FatalErrorInArgument:
      throw G4ExceptionHelper::specifiedException<FatalErrorInArgument>(exceptionCode);
      break;
    case RunMustBeAborted:
      throw G4ExceptionHelper::specifiedException<RunMustBeAborted>(exceptionCode);
      break;
    case EventMustBeAborted:
      throw G4ExceptionHelper::specifiedException<EventMustBeAborted>(exceptionCode);
      break;
    case JustWarning:
      throw G4ExceptionHelper::specifiedException<JustWarning>(exceptionCode);
      break;
    case IgnoreTheIssue:
      break;
  }
  return false;  // By default do not generate a core dump for tests
}
}  // namespace G4ExceptionHelper

#endif
