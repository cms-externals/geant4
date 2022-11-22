#include "AnalysisMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

AnalysisMessenger::AnalysisMessenger()
{
	analysisDir = new G4UIdirectory("/analysis/");
	analysisDir -> SetGuidance("Customize data analysis output");
	
	useRootForOutput = new G4UIcmdWithABool("/analysis/useRoot", this);
	useRootForOutput -> SetGuidance("If true, a ROOT file is used for the output. If false, a plaintext csv file is used");
	useRootForOutput -> SetParameterName("Root", true);
	useRootForOutput -> AvailableForStates(G4State_PreInit);
	
	//default
	rootOutput = true;
}

AnalysisMessenger::~AnalysisMessenger()
{
	delete useRootForOutput;
	
	delete analysisDir;
}

void AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String commandContent)
{

	if( command == useRootForOutput )
	{
		rootOutput = G4UIcmdWithABool::GetNewBoolValue(commandContent);
		
		if( rootOutput == true ) G4cout << "Outputting to a ROOT file" << G4endl;
		else G4cout << "Outputting to a csv file" << G4endl;
	}
}
