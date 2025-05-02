// Author: Christian Velten (2025)

#include "G4MoleculeLocator.hh"

#include "G4GeometryManager.hh"
#include "G4Navigator.hh"
#include "G4ITTransportationManager.hh"

G4ThreadLocal G4MoleculeLocator* G4MoleculeLocator::fpInstance = nullptr;

G4MoleculeLocator::G4MoleculeLocator()
	: fIsInitialized(false)
{
	fNavigator = std::make_unique<G4Navigator>();
}

G4MoleculeLocator* G4MoleculeLocator::Instance()
{
	if (fpInstance == nullptr)
	{
		fpInstance = new G4MoleculeLocator;
	}

	if (!fpInstance->fIsInitialized)
	{
		fpInstance->Initialize();
	}

	return fpInstance;
}

void G4MoleculeLocator::Initialize()
{
	fNavigator->SetWorldVolume(G4ITTransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
	fIsInitialized = true;
}

G4TouchableHistory* G4MoleculeLocator::CreateTouchableHistory() const
{
	return fNavigator->CreateTouchableHistory();
}

G4VPhysicalVolume* G4MoleculeLocator::LocateGlobalPointAndSetup(const G4ThreeVector& point,
																const G4ThreeVector* direction)
{
	return fNavigator->LocateGlobalPointAndSetup(point, direction, false,
												 direction == nullptr);
}

G4VTouchable* G4MoleculeLocator::LocateGlobalPointAndReturnNewTouchable(const G4ThreeVector& position)
{
	fNavigator->ResetStackAndState();
	G4VTouchable *touchable = fNavigator->CreateTouchableHistory();
	fNavigator->LocateGlobalPointAndUpdateTouchable(position, touchable, false);
	return touchable;
}

G4VTouchable* G4MoleculeLocator::LocateGlobalPointAndReturnNewTouchable(const G4ThreeVector& position,
																		const G4ThreeVector& direction)
{
	fNavigator->ResetStackAndState();
	G4VTouchable *touchable = fNavigator->CreateTouchableHistory();
	fNavigator->LocateGlobalPointAndUpdateTouchable(position, direction, touchable, false);
	return touchable;
}

G4TouchableHandle G4MoleculeLocator::LocateGlobalPointAndReturnNewTouchableHandle(const G4ThreeVector& position)
{
	return G4TouchableHandle{LocateGlobalPointAndReturnNewTouchable(position)};
}

G4TouchableHandle G4MoleculeLocator::LocateGlobalPointAndReturnNewTouchableHandle(const G4ThreeVector& position,
																				  const G4ThreeVector& direction)
{
	return G4TouchableHandle{LocateGlobalPointAndReturnNewTouchable(position, direction)};
}

void G4MoleculeLocator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector& position,
															G4VTouchable* touchable)
{
	fNavigator->LocateGlobalPointAndUpdateTouchable(position, touchable, false);
}

void G4MoleculeLocator::LocateGlobalPointAndUpdateTouchable(const G4ThreeVector& position,
															const G4ThreeVector& direction,
															G4VTouchable* touchable)
{
	fNavigator->LocateGlobalPointAndUpdateTouchable(position, direction, touchable, false);
}
