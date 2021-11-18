#ifndef eRositaActionInitialization_h
#define eRositaActionInitialization_h 1

#include "globals.hh"

#include "G4VUserActionInitialization.hh"
#include "G4VSteppingVerbose.hh"

class eRositaActionInitialization : public G4VUserActionInitialization
{
public:
    eRositaActionInitialization();

    virtual ~eRositaActionInitialization();

    virtual void BuildForMaster() const;

    virtual void Build() const;
};
#endif
