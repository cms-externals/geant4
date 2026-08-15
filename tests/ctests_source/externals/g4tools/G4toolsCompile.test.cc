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

// This is a totally dumb, but simple way, to check that the tools headers compile
// correctly in case they are not included in a source file.
//
// It also provides a ways to run clang-tidy on individual tools headers as they must be
// #included in a source file for clang-tidy to see them. To run clang-tidy on an individual
// tools header, setup for Geant4 and clang-tidy as normal (see CODING_GUIDELINES.rst), then run
//
// $ run-clang-tidy tests/ctests_source/externals -header-filter="tools/HEADERTOTEST"
//
// We do leave out some .icc files and those which pull in externals like zlib, xml...

// TODO: Expand this list as needed
#include <tools/HEADER>
#include <tools/array>
#include <tools/buf2lines>
#include <tools/ccontour>
#include <tools/charmanip>
#include <tools/cid>
#include <tools/cids>
#include <tools/cmemT>
#include <tools/colorf>
#include <tools/colorfs>
#include <tools/colors>
#include <tools/columns>
#include <tools/curve>
#include <tools/eqT>
#include <tools/file>
#include <tools/fileis>
#include <tools/forit>
#include <tools/fpng>
#include <tools/gl2ps>
#include <tools/glprims>
#include <tools/handle>
#include <tools/hatcher>
#include <tools/hershey>
#include <tools/hls>
#include <tools/hplot>
#include <tools/img>
#include <tools/impi>
#include <tools/mapmanip>
#include <tools/mathd>
#include <tools/mathf>
#include <tools/mnmx>
#include <tools/nostream>
#include <tools/num2s>
#include <tools/path>
#include <tools/platform>
#include <tools/rcmp>
#include <tools/rntuple>
#include <tools/scast>
#include <tools/schar>
#include <tools/sep>
#include <tools/smath>
#include <tools/sout>
#include <tools/spline>
#include <tools/sprintf>
#include <tools/srep>
#include <tools/sto>
#include <tools/strip>
#include <tools/stype>
#include <tools/tokenize>
#include <tools/toojpeg>
#include <tools/tos>
#include <tools/touplow>
#include <tools/typedefs>
#include <tools/value>
#include <tools/version>
#include <tools/viewplot>
#include <tools/vmanip>
#include <tools/wps>

#include <tools/S_STRING>
#include <tools/aida_ntuple>
#include <tools/charp_out>
#include <tools/clist_contour>
#include <tools/data_axis>
#include <tools/file_reader>
#include <tools/get_env>
#include <tools/get_lines>
#include <tools/impi_world>
#include <tools/ntuple_binding>
#include <tools/ntuple_booking>
#include <tools/press_func>
#include <tools/raxml_out>
#include <tools/rcsv_histo>
#include <tools/rcsv_ntuple>
#include <tools/tess_contour>
#include <tools/wcsv_histo>
#include <tools/wcsv_ntuple>

// Histogram code
#include <tools/histo/axes>
#include <tools/histo/axis>
#include <tools/histo/b1>
#include <tools/histo/b2>
#include <tools/histo/b3>
#include <tools/histo/c1d>
#include <tools/histo/c2d>
#include <tools/histo/c3d>
#include <tools/histo/dps>
#include <tools/histo/h1>
#include <tools/histo/h1d>
#include <tools/histo/h1df>
#include <tools/histo/h2>
#include <tools/histo/h2d>
#include <tools/histo/h2df>
#include <tools/histo/h3>
#include <tools/histo/h3d>
#include <tools/histo/h3df>
#include <tools/histo/hd2mpi>
#include <tools/histo/hmpi>
#include <tools/histo/p1>
#include <tools/histo/p1d>
#include <tools/histo/p2>
#include <tools/histo/p2d>

#include <tools/histo/base_cloud>
#include <tools/histo/base_histo>
#include <tools/histo/histo_data>
#include <tools/histo/profile_data>

// IO code
#include <tools/io/irbuf>
#include <tools/io/iwbuf>

// Linear Algebra code
#include <tools/lina/box3>
#include <tools/lina/clip>
#include <tools/lina/line>
#include <tools/lina/mat4>
#include <tools/lina/mat>
#include <tools/lina/plane>
#include <tools/lina/qrot>
#include <tools/lina/vec>

// ROOT file reading
#include <tools/rroot/THistogram>
#include <tools/rroot/basket>
#include <tools/rroot/branch>
#include <tools/rroot/buffer>
#include <tools/rroot/cids>
#include <tools/rroot/clss>
#include <tools/rroot/date>
#include <tools/rroot/directory>
#include <tools/rroot/dummy>
#include <tools/rroot/fac>
#include <tools/rroot/file>
#include <tools/rroot/graph>
#include <tools/rroot/ifac>
#include <tools/rroot/ifile>
#include <tools/rroot/info>
#include <tools/rroot/iobject>
#include <tools/rroot/iro>
#include <tools/rroot/iros>
#include <tools/rroot/key>
#include <tools/rroot/leaf>
#include <tools/rroot/matrix>
#include <tools/rroot/named>
#include <tools/rroot/ntuple>
#include <tools/rroot/object>
#include <tools/rroot/rall>
#include <tools/rroot/rbuf>
#include <tools/rroot/seek>
#include <tools/rroot/streamers>
#include <tools/rroot/tree>
#include <tools/rroot/vector3>

#include <tools/rroot/base_leaf>
#include <tools/rroot/branch_element>
#include <tools/rroot/branch_object>
#include <tools/rroot/dummy_fac>
#include <tools/rroot/obj_array>
#include <tools/rroot/obj_list>
#include <tools/rroot/stl_vector>
#include <tools/rroot/streamer_fac>
#include <tools/rroot/tree_index>

// ROOT file writing
#include <tools/wroot/basket>
#include <tools/wroot/branch>
#include <tools/wroot/buffer>
#include <tools/wroot/bufobj>
#include <tools/wroot/cids>
#include <tools/wroot/date>
#include <tools/wroot/directory>
#include <tools/wroot/element>
#include <tools/wroot/file>
#include <tools/wroot/ibo>
#include <tools/wroot/icol>
#include <tools/wroot/idir>
#include <tools/wroot/ifile>
#include <tools/wroot/imutex>
#include <tools/wroot/info>
#include <tools/wroot/infos>
#include <tools/wroot/iobject>
#include <tools/wroot/itree>
#include <tools/wroot/key>
#include <tools/wroot/leaf>
#include <tools/wroot/named>
#include <tools/wroot/ntuple>
#include <tools/wroot/seek>
#include <tools/wroot/streamers>
#include <tools/wroot/to>
#include <tools/wroot/tree>
#include <tools/wroot/wbuf>

#include <tools/wroot/base_leaf>
#include <tools/wroot/base_pntuple>
#include <tools/wroot/base_pntuple_column_wise>
#include <tools/wroot/base_pntuple_row_wise>
#include <tools/wroot/branch_element>
#include <tools/wroot/branch_object>
#include <tools/wroot/free_seg>
#include <tools/wroot/impi_ntuple>
#include <tools/wroot/imt_ntuple>
#include <tools/wroot/mpi_basket_add>
#include <tools/wroot/mpi_create_basket>
#include <tools/wroot/mpi_ntuple_column_wise>
#include <tools/wroot/mpi_ntuple_row_wise>
#include <tools/wroot/mpi_protocol>
#include <tools/wroot/mpi_send_basket>
#include <tools/wroot/mt_basket_add>
#include <tools/wroot/mt_ntuple_column_wise>
#include <tools/wroot/mt_ntuple_row_wise>

// ??
#include <tools/waxml/begend>
#include <tools/waxml/histos>
#include <tools/waxml/ntuple>

// ??
#include <tools/xml/aidas>
#include <tools/xml/element>
#include <tools/xml/styles>
#include <tools/xml/tree>

#include <tools/xml/wrap_viewplot_fonts_google_style>

// SG
#include <tools/sg/action>
#include <tools/sg/axis>
#include <tools/sg/bcbk>
#include <tools/sg/blend>
#include <tools/sg/bmf>
#include <tools/sg/bsf>
#include <tools/sg/cbks>
#include <tools/sg/cloud2plot>
#include <tools/sg/colormap>
#include <tools/sg/cube>
#include <tools/sg/ecbk>
#include <tools/sg/ellipse>
#include <tools/sg/enums>
#include <tools/sg/event>
#include <tools/sg/field>
#include <tools/sg/group>
#include <tools/sg/gstos>
#include <tools/sg/h2plot>
#include <tools/sg/holder>
#include <tools/sg/keys>
#include <tools/sg/legend>
#include <tools/sg/lpat>
#include <tools/sg/markers>
#include <tools/sg/matrix>
#include <tools/sg/mf>
#include <tools/sg/mnmx>
#include <tools/sg/node>
#include <tools/sg/nodekit>
#include <tools/sg/noderef>
#include <tools/sg/normal>
#include <tools/sg/ortho>
#include <tools/sg/path>
#include <tools/sg/perspective>
#include <tools/sg/plots>
#include <tools/sg/plottable>
#include <tools/sg/plottables>
#include <tools/sg/plotter>
#include <tools/sg/rep>
#include <tools/sg/rgba>
#include <tools/sg/search>
#include <tools/sg/senum>
#include <tools/sg/senums>
#include <tools/sg/separator>
#include <tools/sg/sf>
#include <tools/sg/state>
#include <tools/sg/states>
#include <tools/sg/strings>
#include <tools/sg/style>
#include <tools/sg/text>
#include <tools/sg/tools>
#include <tools/sg/torche>
#include <tools/sg/vertices>
#include <tools/sg/viewer>

#include <tools/sg/_switch>
#include <tools/sg/atb_vertices>
#include <tools/sg/back_area>
#include <tools/sg/base_camera>
#include <tools/sg/base_freetype>
#include <tools/sg/base_tex>
#include <tools/sg/base_text>
#include <tools/sg/bbox_action>
#include <tools/sg/cloud2plot_cp>
#include <tools/sg/cursor_shape>
#include <tools/sg/device_interactor>
#include <tools/sg/draw_style>
#include <tools/sg/dummy_freetype>
#include <tools/sg/event_action>
#include <tools/sg/event_dispatcher>
#include <tools/sg/field_desc>
#include <tools/sg/get_matrix_action>
#include <tools/sg/gl2ps_action>
#include <tools/sg/gl2ps_manager>
#include <tools/sg/gstos_add>
#include <tools/sg/h2plot_cp>
#include <tools/sg/head_light>
#include <tools/sg/ifield_factory>
#include <tools/sg/infos_box>
#include <tools/sg/light_off>
#include <tools/sg/line_style>
#include <tools/sg/matrix_action>
#include <tools/sg/node_desc>
#include <tools/sg/pick_action>
#include <tools/sg/plots_viewer>
#include <tools/sg/plotter_some_styles>
#include <tools/sg/plotter_style>
#include <tools/sg/primitive_visitor>
#include <tools/sg/read_action>
#include <tools/sg/render_action>
#include <tools/sg/render_gstos>
#include <tools/sg/render_manager>
#include <tools/sg/search_action>
#include <tools/sg/sf_enum>
#include <tools/sg/sf_img>
#include <tools/sg/sf_mat4f>
#include <tools/sg/sf_rotf>
#include <tools/sg/sf_string>
#include <tools/sg/sf_vec3f>
#include <tools/sg/sf_vec4f>
#include <tools/sg/sf_vec>
#include <tools/sg/style_color>
#include <tools/sg/style_colormap>
#include <tools/sg/style_parser>
#include <tools/sg/tex_quadrilateral>
#include <tools/sg/tex_rect>
#include <tools/sg/text_hershey>
#include <tools/sg/text_hershey_marker>
#include <tools/sg/text_style>
#include <tools/sg/text_valop>
#include <tools/sg/visible_action>
#include <tools/sg/win_action>
#include <tools/sg/write_action>
#include <tools/sg/write_paper>
#include <tools/sg/zb_action>
#include <tools/sg/zb_manager>
#include <tools/sg/zb_viewer>

// This is a test suite as such, but a test to check that the tools headers compile.
// We put an empty case in here just so we have something to compile and link.
#include <gtest/gtest.h>

TEST(tools_clang_tidy, Basic) {}
