#!/bin/bash
# To be reviewed: running this on Windows with Linux Subsystem

# Try to ensure execution rather than sourcing
if [ ! "${BASH_SOURCE[0]}" -ef "$0" ]
then
  echo "error: SetupForDevelopment.sh should be executed, not sourced!" && return 1
fi

# Though semi-protected from sourcing above, try and implement
# with as few side effects as possible.
cd "${BASH_SOURCE[0]%/*}/.." &&
if [ ! -d "$PWD/.git" ] ; then
  echo "Not a git worktree of Geant4" && exit 1
fi &&
GitUtilities/setup-user &&
GitUtilities/setup-config &&
GitUtilities/setup-hooks &&

# Record script version
# This will be checked in pre-commit hook to require rerun
# if a new version is detected
Geant4_GitUtilities_VERSION=3
git config geant4.GitUtilitiesVersion ${Geant4_GitUtilities_VERSION}

