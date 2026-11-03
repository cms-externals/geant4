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
// Author: Christian Velten (2025~26)

#ifndef G4MOLECULECOUNTERTEMPLATES_HH
#define G4MOLECULECOUNTERTEMPLATES_HH

#include "G4Types.hh"
#include "G4VMoleculeCounter.hh"

#include <map>
#include <set>
#include <vector>

namespace G4
{
namespace MoleculeCounter
{

//------------------------------------------------------------------------------

template<typename T>
G4String GetTemplateTypeName()
{
  return G4String(typeid(T).name());
}

//------------------------------------------------------------------------------

template<typename T>
G4bool Contains(const std::set<T>& _set, const T& _value)
{
  return std::binary_search(_set.cbegin(), _set.cend(), _value);
}

template<typename T>
G4bool Contains(const std::vector<T>& _vector, const T& _value)
{
  auto p = std::find(_vector.cbegin(), _vector.cend(), _value);
  return p != _vector.cend();
}

template<typename T, typename U>
G4bool ContainsKey(const std::map<T, U>& _map, const T& _key)
{
  auto keys = GetMapIndices(_map);
  return Contains(keys, _key);
}

//------------------------------------------------------------------------------

template<typename T, typename U, typename V>
const std::vector<T> GetMapIndices(const std::map<T, U, V>& _map)
{
  std::vector<T> output;
  for (auto const& it : _map)
    output.push_back(it.first);
  return output;
}

//------------------------------------------------------------------------------

template<typename T>
void DumpCounterMapIndices(const std::map<T, InnerCounterMapType>& map, G4bool includeEmpty = false)
{
  static_assert(
    std::is_base_of<G4VMoleculeCounterInternalBase::G4VMoleculeCounterIndexInterface, T>::value,
    "T must be derived from G4VMoleculeCounter::G4VMoleculeCounterIndex or "
    "from G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex! "
    "No forward declaration is allowed.");

  G4cout << "--- BEGIN COUNTER MAP INDEX DUMP ---" << G4endl;
  auto i = 0;
  for (auto const& it : map)
  {
    if (!includeEmpty && it.second.size() == 0) continue;
    G4cout << i++ << ": " << it.first.GetInfo() << G4endl;
  }
  G4cout << "---  END  COUNTER MAP INDEX DUMP ---" << G4endl;
}

//------------------------------------------------------------------------------

template<typename T>
void DumpCounterMapContents(const std::map<T, InnerCounterMapType>& map,
                            G4bool includeEmpty = false)
{
  static_assert(
    std::is_base_of<G4VMoleculeCounterInternalBase::G4VMoleculeCounterIndexInterface, T>::value,
    "T must be derived from G4VMoleculeCounter::G4VMoleculeCounterIndex or "
    "from G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex! "
    "No forward declaration is allowed.");

  G4cout << "--- BEGIN COUNTER MAP DUMP ---" << G4endl;
  for (auto const& it : map)
  {
    if (!includeEmpty && it.second.size() == 0) continue;
    G4cout << " :: " << it.first.GetInfo() << G4endl;
    for (auto const& it2 : it.second)
    {
      G4cout << std::setw(3) << std::setprecision(3) << " " << G4BestUnit(it2.first, "Time")
             << "    " << std::setw(5) << it2.second << G4endl;
    }
    G4cout << G4endl;
  }
  G4cout << "---  END  COUNTER MAP DUMP ---" << G4endl;
}

//------------------------------------------------------------------------------

template<typename... V>
void DumpCounterMapIndices(const std::map<std::tuple<V...>, InnerCounterMapType>& map,
                           G4bool includeEmpty = false)
{
  G4cout << "--- BEGIN COUNTER MAP INDEX(TUPLE) DUMP ---" << G4endl;
  auto i = 0;
  for (auto const& it : map)
  {
    if (!includeEmpty && it.second.size() == 0) continue;
    G4cout << i++ << ": ";
    std::apply([](auto&&... args) { ((G4cout << args << " "), ...); }, it.first);
    G4cout << G4endl;
  }
  G4cout << "---  END  COUNTER MAP INDEX(TUPLE) DUMP ---" << G4endl;
}

//------------------------------------------------------------------------------

template<typename... V>
void DumpCounterMapContents(const std::map<std::tuple<V...>, InnerCounterMapType>& map,
                            G4bool includeEmpty = false)
{
  G4cout << "--- BEGIN COUNTER MAP DUMP ---" << G4endl;
  for (auto const& it : map)
  {
    if (!includeEmpty && it.second.size() == 0) continue;
    G4cout << " :: ";
    std::apply([](auto&&... args) { ((G4cout << args << " "), ...); }, it.first);
    G4cout << G4endl;
    for (auto const& it2 : it.second)
    {
      G4cout << std::setw(3) << std::setprecision(3) << " " << G4BestUnit(it2.first, "Time")
             << "    " << std::setw(5) << it2.second << G4endl;
    }
    G4cout << G4endl;
  }
  G4cout << "---  END  COUNTER MAP DUMP ---" << G4endl;
}

//------------------------------------------------------------------------------

template<typename... V>
void DumpCounterMapContentsPretty(const std::map<std::tuple<V...>, InnerCounterMapType>& data,
                                  const std::vector<std::string>& keyHeaders = {},
                                  std::ostream& out = G4cout, G4bool flush = true)
{
  constexpr std::size_t N = sizeof...(V);
  constexpr int colW = 14;

  // --- Header row ---
  for (std::size_t i = 0; i < N; ++i)
  {
    std::string h = (i < keyHeaders.size()) ? keyHeaders[i] : ("Key" + std::to_string(i + 1));
    out << std::setw(colW) << h;
  }
  out << std::setw(colW) << "Time" << std::setw(colW) << "Count" << "\n";

  // --- Separator ---
  out << std::string(colW * (N + 2), '-') << "\n";

  // --- Data rows ---
  for (const auto& [outerKey, innerMap] : data)
  {
    for (const auto& [time, count] : innerMap)
    {
      std::apply([&](const auto&... args) { ((out << std::setw(colW) << args), ...); }, outerKey);
      out << std::setw(colW) << time << std::setw(colW) << count << "\n";
    }
  }
  if (flush) out << G4endl;
}

//------------------------------------------------------------------------------

template<typename T>
std::set<G4double> GetRecordedTimes(const std::map<T, InnerCounterMapType>& map)
{
  std::set<G4double> output{};
  for (const auto& it : map)
  {
    G4cerr << it.second.size() << G4endl;
    for (const auto& it2 : it.second)
    {
      output.insert(it2.first);
    }
  }
  return output;
}

//------------------------------------------------------------------------------

template<typename TKey, typename TValue, typename TComp>
typename std::map<TKey, TValue, TComp>::iterator
FindClosestEntryForKey(std::map<TKey, TValue, TComp>& map, TKey key)
{
  if (map.empty()) return map.end();

  auto it_lb = map.lower_bound(key);

  if (it_lb == map.begin()) return it_lb;
  if (it_lb == map.end()) return --it_lb;

  auto prev_it = std::prev(it_lb);

  if (std::abs(it_lb->first - key) < std::abs(prev_it->first - key))
    return it_lb;
  else
    return prev_it;
}

//------------------------------------------------------------------------------

}  // namespace MoleculeCounter
}  // namespace G4

#endif
