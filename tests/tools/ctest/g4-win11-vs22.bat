@echo off 

rem - replace computer name by generic
if %COMPUTERNAME:~0,6% == G4-WIN SET CTEST_SITE=g4-win

echo starting for g4-win10-vs22

echo running on %COMPUTERNAME%.

set Platform=x64


set COMNTOOLS=C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\
rem ---Compiler------------------------------------------------------
call "%COMNTOOLS%\..\..\VC\Auxiliary\Build\vcvarsall.bat" %Platform%


rem ---Xerces-C------------------------------------------------------


rem -- XercesC search first line has highest priority, as this ends up last  

set xerc_dirs=D:\sw\lcg\external\XercesC\3.2.5\%Platform%-windows-vc22
rem - set xerc_dirs=D:\sw\lcg\external\XercesC\3.2.3\%Platform%-windows-vc22,%xerc_dirs%

for %%a IN (%xerc_dirs%) DO (
        echo %%a
	if exist %%a   set XERCESC_ROOT_DIR=%%a
)         
echo Using XercesC from %XERCESC_ROOT_DIR%

set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

set CMAKE_PREFIX_PATH=%XERCESC_ROOT_DIR%

rem ---------------------------------------


rem - location of call  .bat file
set THIS=%~d0%~p0

rem ---Define basic config parameters--------------------------------
call %THIS%g4-win-common.bat

rem net use j: \\cernbox-drive\project\g\geant4
rem dir \\cernbox-drive\project\g\geant4\dev\data

echo Geant4 CMake options - 2 : %G4_XOPTS%

rem - echo %Path%
rem - python --version 
rem - set

set Path=%Path%;C:\Python27

set CdashTrack=Nightly
if %MODE% == patch (
  CdashTrack=PatchesNightly
  MODE=nightly
)

python --version
set 

echo starting ctest with CdashTrack %CdashTrack% using %THIS%g4%MODE%.cmake
rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake
