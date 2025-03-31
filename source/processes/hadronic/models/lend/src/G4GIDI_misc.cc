/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <g4gidi_version.hh>
#include <G4GIDI.hh>

/* *********************************************************************************************************//**
 * Returns the version string for the version of G4GIDI.
 *
 * @return                      A std::string.
 ***********************************************************************************************************/

std::string G4GIDI_version( ) {

    std::string versionString( G4GIDI_VERSION );

    return versionString;
}

/* *********************************************************************************************************//**
 * Returns the major version number for the version of G4GIDI.
 *
 * @return                      An int.
 ***********************************************************************************************************/

int G4GIDI_versionMajor( ) {

    return G4GIDI_MAJOR;
}

/* *********************************************************************************************************//**
 * Returns the minor version number for the version of G4GIDI.
 *
 * @return                      An int.
 ***********************************************************************************************************/

int G4GIDI_versionMinor( ) {

    return G4GIDI_MINOR;
}

/* *********************************************************************************************************//**
 * Returns the patch level number for the version of G4GIDI.
 *
 * @return                      An int.
 ***********************************************************************************************************/

int G4GIDI_versionPatchLevel( ) {

    return G4GIDI_PATCHLEVEL;
}

/* *********************************************************************************************************//**
 * Returns the git repo hash for the G4GIDI repo for this version..
 *
 * @return                      A std::string.
 ***********************************************************************************************************/

std::string G4GIDI_GitHash( ) {

    return G4GIDI_GIT;
}

static std::vector<char const *> G4GIDI_chemicalElementSymbols { // Note, this is zero based. Ergo, "O" is at index 7 not 8.
         "H",  "He",  "Li",  "Be",   "B",   "C",   "N",   "O",   "F",  "Ne",
        "Na",  "Mg",  "Al",  "Si",   "P",   "S",  "Cl",  "Ar",   "K",  "Ca",
        "Sc",  "Ti",   "V",  "Cr",  "Mn",  "Fe",  "Co",  "Ni",  "Cu",  "Zn",
        "Ga",  "Ge",  "As",  "Se",  "Br",  "Kr",  "Rb",  "Sr",   "Y",  "Zr",
        "Nb",  "Mo",  "Tc",  "Ru",  "Rh",  "Pd",  "Ag",  "Cd",  "In",  "Sn",
        "Sb",  "Te",   "I",  "Xe",  "Cs",  "Ba",  "La",  "Ce",  "Pr",  "Nd",
        "Pm",  "Sm",  "Eu",  "Gd",  "Tb",  "Dy",  "Ho",  "Er",  "Tm",  "Yb",
        "Lu",  "Hf",  "Ta",   "W",  "Re",  "Os",  "Ir",  "Pt",  "Au",  "Hg",
        "Tl",  "Pb",  "Bi",  "Po",  "At",  "Rn",  "Fr",  "Ra",  "Ac",  "Th",
        "Pa",    "U", "Np",  "Pu",  "Am",  "Cm",  "Bk",  "Cf",  "Es",  "Fm",
        "Md",  "No",  "Lr",  "Rf",  "Db",  "Sg",  "Bh",  "Hs",  "Mt",  "Ds",
        "Rg",  "Cn",  "Nh",  "Fl",  "Mc",  "Lv",  "Ts",  "Og" };

/* *********************************************************************************************************//**
 * Returns the chemical elemental symbol for a specified atomic number *a_Z*.
 *
 * @param a_Z           [in]    The Z of the chemical element.
 *
 * @return                      The symbol as a std::string.
 ***********************************************************************************************************/

std::string G4GIDI_Misc_Z_toSymbol( int a_Z ) {

    int Zm1 = a_Z - 1;

    if( ( Zm1 < 0 ) || ( Zm1 >= static_cast<int>( G4GIDI_chemicalElementSymbols.size( ) ) ) ) {
        throw LUPI::Exception( "G4GIDI_Misc_Z_toSymbol: invalid Z = " + std::to_string( a_Z ) + "." );
    }

    return( std::string( G4GIDI_chemicalElementSymbols[Zm1] ) );
}

/* *********************************************************************************************************//**
 * Returns the PoPs for a specified atomic number *a_Z*, mass number *a_A* and meta-stable index *a_M*.
 *
 * @param a_Z           [in]    The Z of the chemical element.
 * @param a_A           [in]    The A of the isotope.
 * @param a_M           [in]    The meta-stable index of the nuclide.
 *
 * @return                      The symbol as a std::string.
 ***********************************************************************************************************/

std::string G4GIDI_Misc_Z_A_m_ToName( int a_Z, int a_A, int a_M ) {

    std::string id( G4GIDI_Misc_Z_toSymbol( a_Z ) + std::to_string( a_A ) );

    if( a_M > 0 ) id += "_" + std::to_string( a_M );

    return( id );
}
