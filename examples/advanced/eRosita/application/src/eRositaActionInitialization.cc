#include "eRositaActionInitialization.hh"
#include "eRositaPrimaryGeneratorAction.hh"
#include "eRositaEventAction.hh"
#include "eRositaRunAction.hh"
#include "eRositaSteppingAction.hh"

#include "G4GeneralParticleSource.hh"

eRositaActionInitialization::eRositaActionInitialization() : G4VUserActionInitialization()
{
}

eRositaActionInitialization::~eRositaActionInitialization()
{
}

void eRositaActionInitialization::BuildForMaster() const
{
    eRositaRunAction* runAction = new eRositaRunAction();
    SetUserAction(runAction);
}

void eRositaActionInitialization::Build() const
{
    eRositaPrimaryGeneratorAction* primaryGeneratorAction = new eRositaPrimaryGeneratorAction();
    SetUserAction(primaryGeneratorAction);

    eRositaRunAction* runAction = new eRositaRunAction();
    SetUserAction(runAction);
    
    //SetUserAction(new eRositaPrimaryGeneratorAction);
    //SetUserAction(new eRositaRunAction);
    SetUserAction(new eRositaEventAction);
    SetUserAction(new eRositaSteppingAction);
}
