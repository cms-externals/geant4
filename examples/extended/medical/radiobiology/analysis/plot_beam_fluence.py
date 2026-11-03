#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# ********************************************************************
# * License and Disclaimer                                           *
# *                                                                  *
# * The Geant4 software is copyright of the Copyright Holders of the  *
# * Geant4 Collaboration. It is provided under the terms and          *
# * conditions of the Geant4 Software License, included in the file   *
# * LICENSE and available at http://cern.ch/geant4/license.           *
# *                                                                  *
# * Neither the authors of this software system, nor their employing  *
# * institutes, nor the agencies providing financial support for     *
# * this work make any representation or warranty, express or         *
# * implied, regarding this software system or assume any liability   *
# * for its use.                                                      *
# * Please see the license in the file LICENSE and URL above for the  *
# * full disclaimer and limitation of liability.                      *
# *                                                                  *
# * This code implementation is the result of the scientific and      *
# * technical work of the Geant4 Collaboration. By using, copying,    *
# * modifying, or distributing the software you agree to acknowledge  *
# * its use in resulting scientific publications and accept all terms *
# * of the Geant4 Software License.                                   *
# ********************************************************************
"""Interactive plot of voxel-resolved particle-entry counts."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.widgets import CheckButtons


# =============================================================================
# User settings
# =============================================================================

fluence_path = Path("beamFluence.out")
selected_y_index = 0
selected_z_index = 0
use_logarithmic_y_axis = True
checkbox_size = 120


# =============================================================================
# Read beamFluence.out
# =============================================================================

if not fluence_path.is_file():
    raise FileNotFoundError(
        "The beam-fluence file was not found:\n"
        f"{fluence_path.resolve()}"
    )

data = pd.read_csv(
    fluence_path,
    sep=r"\s+",
    engine="python",
)

index_columns = (
    "x_index",
    "y_index",
    "z_index",
)

missing_columns = [
    column
    for column in index_columns
    if column not in data.columns
]

if missing_columns:
    raise ValueError(
        "The input file is missing the required column(s): "
        + ", ".join(missing_columns)
    )

for column in data.columns:
    data[column] = pd.to_numeric(
        data[column],
        errors="raise",
    )


# =============================================================================
# Detect all particle-entry columns automatically
# =============================================================================

fluence_columns = [
    column
    for column in data.columns
    if column not in index_columns
]

if not fluence_columns:
    raise ValueError(
        "No particle-entry columns were found in the input file."
    )

print(f"Loaded file: {fluence_path.resolve()}")
print(f"Number of rows: {len(data)}")
print(f"Number of particle-entry columns: {len(fluence_columns)}")
print("\nParticle-entry columns detected:")

for column in fluence_columns:
    print(f"  - {column}")


# =============================================================================
# Select the profile along x
# =============================================================================

profile = data[
    (data["y_index"] == selected_y_index)
    & (data["z_index"] == selected_z_index)
].copy()

if profile.empty:
    available_indices = (
        data[["y_index", "z_index"]]
        .drop_duplicates()
        .sort_values(["y_index", "z_index"])
    )

    available_pairs = ", ".join(
        f"({int(row.y_index)}, {int(row.z_index)})"
        for row in available_indices.itertuples(index=False)
    )

    raise ValueError(
        f"No data were found for y_index={selected_y_index} and "
        f"z_index={selected_z_index}.\n"
        f"Available pairs are: {available_pairs}"
    )

profile = profile.sort_values("x_index")


# =============================================================================
# Create the figure
# =============================================================================

figure_height = max(
    7,
    0.35 * len(fluence_columns) + 2,
)

fig, ax = plt.subplots(
    figsize=(14, figure_height)
)

plt.subplots_adjust(
    left=0.09,
    right=0.69,
    bottom=0.10,
    top=0.92,
)


# =============================================================================
# Plot all particle-entry columns
# =============================================================================

lines = {}

for column in fluence_columns:
    y_values = profile[column].copy()

    if use_logarithmic_y_axis:
        y_values = y_values.where(y_values > 0)

    line, = ax.plot(
        profile["x_index"],
        y_values,
        linestyle="--",
        marker="o",
        markersize=5,
        linewidth=1.5,
        label=column,
    )

    lines[column] = line


# =============================================================================
# Plot formatting
# =============================================================================

ax.set_xlabel("x index", fontsize=15)
ax.set_ylabel("Particles entering the voxel", fontsize=15)

if use_logarithmic_y_axis:
    ax.set_yscale("log")

ax.set_title(
    "Particle-entry profile along the beam-propagation direction",
    fontsize=16,
)

ax.tick_params(
    axis="both",
    labelsize=13,
)

ax.grid(
    True,
    which="both",
    alpha=0.3,
)


# =============================================================================
# Interactive checkbox panel
# =============================================================================

checkbox_axis = fig.add_axes(
    [0.71, 0.10, 0.28, 0.82]
)

checkbox_axis.set_title(
    "Particle-entry contributions",
    fontsize=12,
)

initial_visibility = [True for _ in fluence_columns]

try:
    checkboxes = CheckButtons(
        checkbox_axis,
        fluence_columns,
        initial_visibility,
        frame_props={"s": checkbox_size},
        check_props={"s": checkbox_size},
    )
except TypeError:
    checkboxes = CheckButtons(
        checkbox_axis,
        fluence_columns,
        initial_visibility,
    )

checkbox_symbols = {}

for label, column in zip(
    checkboxes.labels,
    fluence_columns,
):
    label.set_fontsize(9)

    _, y_position = label.get_position()
    label.set_position((0.43, y_position))

    curve = lines[column]

    symbol = Line2D(
        [0.20, 0.36],
        [y_position, y_position],
        transform=checkbox_axis.transAxes,
        color=curve.get_color(),
        linestyle=curve.get_linestyle(),
        marker=curve.get_marker(),
        markersize=curve.get_markersize(),
        linewidth=curve.get_linewidth(),
    )

    checkbox_axis.add_artist(symbol)
    checkbox_symbols[column] = symbol


# =============================================================================
# Checkbox callback
# =============================================================================

def toggle_curve(label):
    """Show or hide the selected particle-entry curve."""
    curve = lines[label]
    symbol = checkbox_symbols[label]

    new_visibility = not curve.get_visible()
    curve.set_visible(new_visibility)
    symbol.set_alpha(1.0 if new_visibility else 0.20)

    fig.canvas.draw_idle()


checkboxes.on_clicked(toggle_curve)


# =============================================================================
# Show the interactive plot
# =============================================================================

plt.show()
