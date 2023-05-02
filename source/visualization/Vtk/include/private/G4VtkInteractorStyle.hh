//
// Created by Stewart Boogert on 28/02/2023.
//

#ifndef G4VTKINTERACTORSTYLE_HH
#define G4VTKINTERACTORSTYLE_HH

#include "vtkInteractionStyleModule.h"  // For export macro
#include "vtkInteractorStyleTrackballCamera.h"
#include <vtkObjectFactory.h>

// Define interaction style
class VTKINTERACTIONSTYLE_EXPORT G4VtkInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static G4VtkInteractorStyle* New();
    vtkTypeMacro(G4VtkInteractorStyle, vtkInteractorStyleTrackballCamera)

      void OnLeftButtonDown() override
    {
      // Forward events
      vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    void OnMiddleButtonDown() override
    {
      // Forward events
      vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    }

    void OnRightButtonDown() override
    {
      // Forward events
      vtkInteractorStyleTrackballCamera::OnRightButtonDown();
    }
};

#endif  // G4VTKINTERACTORSTYLE_HH
