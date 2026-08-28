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
// Check C++14's relaxed constexpr requirements work
// - https://en.cppreference.com/w/cpp/language/constexpr

// Can now declare local variables, use loop/conditionals
constexpr int my_strcmp( const char* str1, const char* str2 ) {
  int i = 0;
  for( ; str1[i] && str2[i] && str1[i] == str2[i]; ++i )
  {}
  if( str1[i] == str2[i] ) return 0;
  if( str1[i] < str2[i] ) return -1;
  return 1;
}

int main() {
  const char* a = "foobar";
  const char* b = "foobaz";

  return my_strcmp(a, b);
}
