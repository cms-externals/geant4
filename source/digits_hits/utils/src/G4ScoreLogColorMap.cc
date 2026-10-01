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
//

#include "G4ScoreLogColorMap.hh"

#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4UIcommand.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

void G4ScoreLogColorMap::GetMapColor(G4double val, G4double color[4])
{
  if (fMinVal < 0. || fMaxVal < 0.)
  {
    G4ExceptionDescription ed;
    ed << "The minimum value (fMinVal=" << fMinVal << ") or the maximum value (fMaxVal=" << fMaxVal
       << ") is negative for this log-scale color map.";
    G4Exception("G4ScoreLogColorMap::GetMapColor()", "DigiHitsUtilsScoreLogColorMap000",
                JustWarning, ed);
    color[0] = 0.;
    color[1] = 0.;
    color[2] = 0.;
    color[3] = 0.;
    return;
  }

  if (val < 0.)
  {
    G4ExceptionDescription ed;
    ed << "Requested value " << val << " is negative.";
    G4Exception("G4ScoreLogColorMap::GetMapColor()", "DigiHitsUtilsScoreLogColorMap001",
                JustWarning, ed);
    color[0] = 0.;
    color[1] = 0.;
    color[2] = 0.;
    color[3] = -1.;
    return;
  }

  G4double logmin = fMinVal > 0. ? std::log10(fMinVal) : 0.;
  G4double logmax = fMaxVal > 0. ? std::log10(fMaxVal) : 0.;
  G4double value =
    (val > 0. && (logmax - logmin) > 0.) ? (std::log10(val) - logmin) / (logmax - logmin) : 0.;
  value = std::max(value, 0.);
  value = std::min(value, 1.);

  // color map
  const G4int NCOLOR = 6;
  struct ColorMap
  {
      G4double val;
      G4double rgb[4];
  } colormap[] = {{0.0, {1., 1., 1., 1.}},  // value, r, g, b, alpha
                  {0.2, {0., 0., 1., 1.}}, {0.4, {0., 1., 1., 1.}}, {0.6, {0., 1., 0., 1.}},
                  {0.8, {1., 1., 0., 1.}}, {1.0, {1., 0., 0., 1.}}};

  // search
  G4int during[2] = {0, 0};
  for (auto i = 1; i < NCOLOR; ++i)
  {
    if (colormap[i].val >= value)
    {
      during[0] = i - 1;
      during[1] = i;
      break;
    }
  }

  // interpolate
  G4double a = std::fabs(value - colormap[during[0]].val);
  G4double b = std::fabs(value - colormap[during[1]].val);
  for (auto i = 0; i < 4; ++i)
  {
    color[i] = (b * colormap[during[0]].rgb[i] + a * colormap[during[1]].rgb[i])
               / (colormap[during[1]].val - colormap[during[0]].val);
    if (color[i] > 1.) color[i] = 1.;
  }
}

void G4ScoreLogColorMap::DrawColorChartBar(G4int _nPoint)
{
  G4double min = fMinVal > 0 ? std::log10(fMinVal) : 0.;
  G4double max = fMaxVal > 0 ? std::log10(fMaxVal) : 0.;

  fVisManager->BeginDraw2D();

  G4double smin = -0.89, smax = smin + 0.05 * (_nPoint) * 0.83, step = 0.001;
  yoff = 1.7 * (smax - smin) * offset;
  G4double c[4];
  for (auto y = smin; y < smax; y += step)
  {
    G4double ra = (y - smin) / (smax - smin), rb = 1. - ra;
    G4Polyline line;
    line.push_back(G4Point3D(-0.96, y + yoff, 0.));
    line.push_back(G4Point3D(-0.91, y + yoff, 0.));
    G4double val = std::pow(10., (ra * max + rb * min) / (ra + rb));
    this->GetMapColor(val, c);
    if (c[0] == 0 && c[1] == 0 && c[2] == 0 && c[3] == 0) return;
    if (c[0] == 0 && c[1] == 0 && c[2] == 0 && c[3] == -1.) continue;
    G4Colour col(c[0], c[1], c[2]);
    G4VisAttributes att(col);
    line.SetVisAttributes(&att);
    fVisManager->Draw2D(line);
  }

  fVisManager->EndDraw2D();
}

void G4ScoreLogColorMap::DrawColorChartText(G4int _nPoint)
{
  G4double min = fMinVal > 0 ? std::log10(fMinVal) : 0.;
  G4double max = fMaxVal > 0 ? std::log10(fMaxVal) : 0.;

  fVisManager->BeginDraw2D();

  G4double c[4] = {1., 1., 1., 1.};
  G4Colour black(0., 0., 0.);
  for (auto n = 0; n < _nPoint; ++n)
  {
    G4double a = n / (_nPoint - 1.), b = 1. - a;
    G4double v = (a * max + b * min) / (a + b);

    this->GetMapColor(std::pow(10., v), c);
    if (c[0] == 0 && c[1] == 0 && c[2] == 0 && c[3] == 0) return;
    if (c[0] == 0 && c[1] == 0 && c[2] == 0 && c[3] == -1.) continue;

    // background color
    for (auto l = 0; l < 21; ++l)
    {
      G4Polyline line;
      line.push_back(G4Point3D(-0.90, -0.905 + 0.05 * n + 0.002 * l + yoff, 0.));
      line.push_back(G4Point3D(-0.80, -0.905 + 0.05 * n + 0.002 * l + yoff, 0.));
      G4VisAttributes attblack(black);
      line.SetVisAttributes(&attblack);
      fVisManager->Draw2D(line);
    }
    // text
    std::ostringstream oss;
    oss << std::setw(8) << std::setprecision(1) << std::scientific << std::pow(10., v);
    std::string str = oss.str();
    G4String value(str);
    G4Text text(value, G4Point3D(-0.9, -0.9 + 0.05 * n + yoff, 0));
    G4double size = 12.;
    text.SetScreenSize(size);
    G4Colour color(1., 1., 1.);
    G4VisAttributes att(color);
    text.SetVisAttributes(&att);

    fVisManager->Draw2D(text);
  }

  // draw ps name and unit
  // background
  G4String txt = fMSName + "/" + fPSName + " [" + fPSUnit + "]";
  G4double lpsname = txt.size();
  for (auto l = 0; l < 22; ++l)
  {
    G4Polyline line;
    line.push_back(G4Point3D(-0.96, -0.965 + 0.002 * l + yoff, 0.));
    line.push_back(G4Point3D(-0.96 + 0.01 * lpsname, -0.965 + 0.002 * l + yoff, 0.));
    G4VisAttributes attblack(black);
    line.SetVisAttributes(&attblack);
    fVisManager->Draw2D(line);
  }
  // ps name and unit
  G4Text txtpsname(txt, G4Point3D(-0.96, -0.96 + yoff, 0.));
  G4double size = 12.;
  txtpsname.SetScreenSize(size);
  G4Colour color(1., 1., 1.);
  G4VisAttributes att(color);
  txtpsname.SetVisAttributes(&att);
  fVisManager->Draw2D(txtpsname);

  fVisManager->EndDraw2D();
}
