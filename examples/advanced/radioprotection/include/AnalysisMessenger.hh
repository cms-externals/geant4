#ifndef AnalysisMessenger_h
#define AnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4String.hh"

class G4UIdirectory;
class G4UIcmdWithABool;

class AnalysisMessenger: public G4UImessenger
{
public:
  AnalysisMessenger();
  ~AnalysisMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
  // functions to read parameters
  G4bool IsRootOutput() { return rootOutput; }
  
private:
  G4UIdirectory *analysisDir;      // UI directory
  
  G4UIcmdWithABool *useRootForOutput;
  
  // parameters to store
  G4bool rootOutput;
};
#endif

