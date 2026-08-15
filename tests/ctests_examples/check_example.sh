#!/bin/bash
#
# This script should be called from the example directory
# Usage:
# check_example.sh [-silent|-violations]
#     -silent:      no outout is printed
#     -violations:  only violations are printed
#
# Script for checking the Geant4 examples coding guidelines:
# - Avoid using tabulators
# - Presence of the agreed separator
# - Avoid using long lines (> 100 characters)
# - Documentation of all macros in README
#
# By I. Hrivnacova, IJCLab Orsay

#set -x

# Help function
function print_help()
{
  echo "Usage:"
  echo "check_example.sh [-silent|-violations]"
  echo "    -silent:      no output is printed"
  echo "    -violations:  only violations are printed"
}

# Process script arguments
SILENT="0"
VIOLATIONS="0"
for arg in "${@}"
do
  #echo "got: $arg"
  case $arg in
    "-silent"        ) SILENT="1" ;;
    "-violations"    ) VIOLATIONS="1";;
    *                ) echo "Unsupported option $arg chosen."
                       print_help
                       exit 1
                       ;;
  esac
done

# Global parameters
INDENTION="   "
MAX=100
TABS_FILES=" "
SEPARATOR_FILES=" "
LONGLINE_FILES=" "
MACRO_FILES=" "
FINAL_RESULT=0

#The correct separator with 80 characters
#SEPARATOR="//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......"
#Still tolerated separator with 78 characters
SEPARATOR="//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...."

# Check functions

# Check tabulations
# {1}: source file to be processed
# {2}: detected files list
function check_tabs()
{
  local RESULT=`cat ${1} | sed 's/\t/TAB_TO_BE_REMOVED/g' | grep TAB_TO_BE_REMOVED` 
  if [ ! "$RESULT" = "" ]; 
  then 
    TABS_FILES="$TABS_FILES ${1}"
  fi
}

# Check separator
# {1}: source file to be processed
function check_separator()
{
  local LINES=`cat ${1} | grep $SEPARATOR | wc -l `  
  local RESULT=`echo "$LINES > 0" | bc`
  if [ $RESULT -ne 1 ]; 
  then 
    SEPARATOR_FILES="$SEPARATOR_FILES ${1}"
  fi
}

# Check long lines
# {1}: source file to be processed
# {2}: detected files list
function check_long_lines()
{
  local LINE_LENGTH=`awk '{ if (length($0) > max) {max = length($0)} } END { print max }' "${1}"`
  local RESULT=`echo "$LINE_LENGTH > $MAX" | bc`
  if [ $RESULT -ne 0 ]; 
  then 
    LONGLINE_FILES="$LONGLINE_FILES ${1}"
  fi
}

# Check if macro is documented in README
# {1}: macro name to be checked
# {2}: detected macros list
function check_macro()
{
  EXAMPLE_NAME="$(basename $PWD)" 
  TEST_NAME1=$EXAMPLE_NAME".in"
  TEST_NAME2=example$EXAMPLE_NAME".in"

  for DOC in README README.md
  do
    if [ -f $DOC ]; then 
      # echo "Checking documentation for $MACRO"
      local LINES=`cat $DOC | grep ${1} | wc -l `  
      local RESULT=`echo "$LINES > 0" | bc`
      if [ $RESULT -ne 1 ] && [ "${1}" != "init_vis.mac" ] && [ "${1}" != "vis.mac" ] &&  [ "${1}" != "gui.mac" ];
      then
        # skip tests (should we require that they are also documented ?)
        if [ ${1} != $TEST_NAME1 ] && [ ${1} != $TEST_NAME2 ]; then
          MACRO_FILES="$MACRO_FILES ${1}"
        fi
      fi
    fi
  done
}

# Check if macro is documented in README
# {1}: list to be printed
# {2}: message
function print_check_result()
{
  if [ "${1}" != " " ]; then
    if [ "$SILENT" == "0" ]; then
      echo "$INDENTION""${2}" >&2
      echo "$INDENTION""   ${1}" >&2
    fi
    FINAL_RESULT=$((FINAL_RESULT+1))
  fi
}

# Process all checks
for SOURCE in `find . -name '*.hh' -o  -name '*.cc' -type f 2> /dev/null`
do
  check_tabs "$SOURCE" "$TAB_FILES"
  check_long_lines "$SOURCE" "$LONGLINE_FILES"
done

for SOURCE in `ls src/*.cc */src/*.cc 2> /dev/null`
do
  check_separator "$SOURCE" "$SEPARATOR_FILES"
done

for MACRO in `ls *.mac *.in *.g4 2> /dev/null`
do
  check_macro "$MACRO" "$MACRO_FILES"
done

# Evaluate check results
#

print_check_result "$TABS_FILES" "TAB detected in:"
print_check_result "$SEPARATOR_FILES" "NO SEPARATOR found in:"
print_check_result "$LONGLINE_FILES" "LONG LINE (> 100 characters) detected in:"
print_check_result "$MACRO_FILES" "Not docummented macros:"

if [ $FINAL_RESULT -eq 0 ]; then
  if [ "$SILENT" == "0" ] && [ "$VIOLATIONS" == "0" ]; then
    echo "Great ! No coding guidelines violations found."
  fi
  exit 0
else
  exit 1
fi
