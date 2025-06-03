// The example class categories definitions for Doxygen

/// \file Doxymodules_visualization.h
/// \brief The page that defines the extended/visualization examples modules 


/** @defgroup extended_visualization visualization
 *  Extended visualization examples classes
 *  @{
 */

/** @defgroup extended_visualization_perspective perspective
 *  Visualization example perspective
 *  @ingroup extended_visualization
 *  @{
 */

  class PerspectiveVisAction {};
  class PerspectiveVisActionMessenger {};

/** @} */

/** @defgroup extended_visualization_standalone standalone
 *  Visualization example standalone
 *  @ingroup extended_visualization
 *  @{
 */

  class StandaloneVisAction {};

/** @} */

/** @defgroup extended_visualization_userVisAction userVisAction
 *  visualization example userVisAction
 *  @ingroup extended_visualization
 *  @{
 */

  class UVA_VisAction {};

/** @} */

/** @defgroup extended_visualization_vtk vtk
 *  visualization example vtk
 *  @ingroup extended_visualization
 *  @{
 */

  class VtkVis::ActionInitialization {};
  class VtkVis::DetectorConstruction {};
  class VtkVis::EventAction {};
  class VtkVis::PrimaryGeneratorAction {};
  class VtkVis::RunAction {};
  class VtkVis::SteppingAction {};

/** @} */

/** @} */
