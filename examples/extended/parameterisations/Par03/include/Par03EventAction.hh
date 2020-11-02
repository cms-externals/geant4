#ifndef PAR03EVENTACTION_HH
#define PAR03EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "G4Timer.hh"

class Par03DetectorConstruction;

/**
 * @brief Event action class for hits' analysis.
 *
 * Analysis of single-particle events and developed showers in the detector.
 * At the end of the event basic variables are calculated and saved in the
 * histograms.
 *
 */

class Par03EventAction : public G4UserEventAction
{
 public:
  Par03EventAction(Par03DetectorConstruction* aDetector);
  virtual ~Par03EventAction();

  /// Timer is started
  virtual void BeginOfEventAction(const G4Event* aEvent) final;
  /// Hits collection is retrieved, analysed, and histograms are filled.
  virtual void EndOfEventAction(const G4Event* aEvent) final;

 private:
  /// ID of a hit collection to analyse
  G4int fHitCollectionID;
  /// Timer measurement
  G4Timer fTimer;
  /// Pointer to detector construction to retrieve (once) the detector
  /// dimensions
  Par03DetectorConstruction* fDetector;
  /// Size of cell along Z axis
  G4double fCellSizeZ = 0;
  /// Size of cell along radius of cylinder
  G4double fCellSizeRho = 0;
};

#endif /* PAR03EVENTACTION_HH */
