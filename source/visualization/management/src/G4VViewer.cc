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
//
// $Id: G4VViewer.cc 111536 2018-08-13 07:43:41Z gcosmo $
//
// 
// John Allison  27th March 1996
// Abstract interface class for graphics views.

#include "G4VViewer.hh"

#include "G4ios.hh"
#include <sstream>

#include "G4VisManager.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VSceneHandler.hh"
#include "G4Scene.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4UImanager.hh"

G4VViewer::G4VViewer (G4VSceneHandler& sceneHandler,
		      G4int id, const G4String& name):
fSceneHandler (sceneHandler),
fViewId (id),
//fModified (true),
fNeedKernelVisit (true)
{
  if (name == "") {
    std::ostringstream ost;
    ost << fSceneHandler.GetName () << '-' << fViewId;
    fName = ost.str();
  }
  else {
    fName = name;
  }
  fShortName = fName (0, fName.find (' '));
  fShortName.strip ();

  fVP = G4VisManager::GetInstance()->GetDefaultViewParameters();
  fDefaultVP = fVP;
}

G4VViewer::~G4VViewer () {
  fSceneHandler.RemoveViewerFromList(this);
}

void G4VViewer::SetName (const G4String& name) {
  fName = name;
  fShortName = fName (0, fName.find (' '));
  fShortName.strip ();
}

void G4VViewer::NeedKernelVisit () {

  fNeedKernelVisit = true;

  // At one time I thought we'd better notify all viewers.  But I guess
  // each viewer can take care of itself, so the following code is
  // redundant (but keep it commented out for now).   (John Allison)
  // Notify all viewers that a kernel visit is required.
  // const G4ViewerList& viewerList = fSceneHandler.GetViewerList ();
  // G4ViewerListConstIterator i;
  // for (i = viewerList.begin(); i != viewerList.end(); i++) {
  //   (*i) -> SetNeedKernelVisit ();
  // }
  // ??...but, there's a problem in OpenGL Stored which seems to
  // require *all* viewers to revisit the kernel, so...
  //  const G4ViewerList& viewerList = fSceneHandler.GetViewerList ();
  //  G4ViewerListConstIterator i;
  //  for (i = viewerList.begin(); i != viewerList.end(); i++) {
  //    (*i) -> SetNeedKernelVisit (true);
  //  }
  // Feb 2005 - commented out.  Let's fix OpenGL if necessary.
}

void G4VViewer::FinishView () {}

void G4VViewer::ShowView () {}

void G4VViewer::ProcessView ()
{
  // If the scene has changed, or if the concrete viewer has decided
  // that it necessary to visit the kernel, perhaps because the view
  // parameters have changed significantly (this should be done in the
  // concrete viewer's DrawView)...
  if (fNeedKernelVisit) {
    // Reset flag.  This must be done before ProcessScene to prevent
    // recursive calls when recomputing transients...
    fNeedKernelVisit = false;
    fSceneHandler.ClearStore ();
    fSceneHandler.ProcessScene ();
  }
}

void G4VViewer::SetViewParameters (const G4ViewParameters& vp) {
  fVP = vp;
}

void G4VViewer::SetTouchable
(const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath)
{
  // Set the touchable for /vis/touchable/set/... commands.
  std::ostringstream oss;
  for (const auto& pvNodeId: fullPath) {
    oss
    << ' ' << pvNodeId.GetPhysicalVolume()->GetName()
    << ' ' << pvNodeId.GetCopyNo();
  }
  G4UImanager::GetUIpointer()->ApplyCommand("/vis/set/touchable" + oss.str());
}

void G4VViewer::TouchableSetVisibility
(const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath,
 G4bool visibiity)
{
  // Changes the Vis Attribute Modifiers WITHOUT triggering a rebuild.

  std::ostringstream oss;
  oss << "/vis/touchable/set/visibility ";
  if (visibiity) oss << "true"; else oss << "false";

  // The following is equivalent to
  //  G4UImanager::GetUIpointer()->ApplyCommand(oss.str());
  // (assuming the touchable has already been set), but avoids view rebuild.

  // Instantiate a working copy of a G4VisAttributes object...
  G4VisAttributes workingVisAtts;
  // and set the visibility.
  workingVisAtts.SetVisibility(visibiity);

  fVP.AddVisAttributesModifier
  (G4ModelingParameters::VisAttributesModifier
   (workingVisAtts,
    G4ModelingParameters::VASVisibility,
    fullPath));
  // G4ModelingParameters::VASVisibility (VAS = Vis Attribute Signifier)
  // signifies that it is the visibility that should be picked out
  // and merged with the touchable's normal vis attributes.

  // Record on G4cout (with #) for information.
  if (G4UImanager::GetUIpointer()->GetVerboseLevel() >= 2) {
    G4cout << "# " << oss.str() << G4endl;
  }
}

void G4VViewer::TouchableSetColour
(const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>& fullPath,
 const G4Colour& colour)
{
  // Changes the Vis Attribute Modifiers WITHOUT triggering a rebuild.

  std::ostringstream oss;
  oss << "/vis/touchable/set/colour "
  << colour.GetRed() << ' ' << colour.GetGreen()
  << ' ' << colour.GetBlue() << ' ' << colour.GetAlpha();

  // The following is equivalent to
  //  G4UImanager::GetUIpointer()->ApplyCommand(oss.str());
  // (assuming the touchable has already been set), but avoids view rebuild.

  // Instantiate a working copy of a G4VisAttributes object...
  G4VisAttributes workingVisAtts;
  // and set the colour.
  workingVisAtts.SetColour(colour);

  fVP.AddVisAttributesModifier
  (G4ModelingParameters::VisAttributesModifier
   (workingVisAtts,
    G4ModelingParameters::VASColour,
    fullPath));
  // G4ModelingParameters::VASColour (VAS = Vis Attribute Signifier)
  // signifies that it is the colour that should be picked out
  // and merged with the touchable's normal vis attributes.

  // Record on G4cout (with #) for information.
  if (G4UImanager::GetUIpointer()->GetVerboseLevel() >= 2) {
    G4cout << "# " << oss.str() << G4endl;
  }
}

#ifdef G4MULTITHREADED

void G4VViewer::DoneWithMasterThread () {
  // G4cout << "G4VViewer::DoneWithMasterThread" << G4endl;
}

void G4VViewer::MovingToMasterThread () {
  // G4cout << "G4VViewer::MovingToMasterThread" << G4endl;
}

void G4VViewer::SwitchToVisSubThread () {
  // G4cout << "G4VViewer::SwitchToVisSubThread" << G4endl;
}

void G4VViewer::DoneWithVisSubThread () {
  // G4cout << "G4VViewer::DoneWithVisSubThread" << G4endl;
}

void G4VViewer::MovingToVisSubThread () {
  // G4cout << "G4VViewer::MovingToVisSubThread" << G4endl;
}

void G4VViewer::SwitchToMasterThread () {
  // G4cout << "G4VViewer::SwitchToMasterThread" << G4endl;
}

#endif

std::ostream& operator << (std::ostream& os, const G4VViewer& v) {
  os << "View " << v.fName << ":\n";
  os << v.fVP;
  return os;
}
