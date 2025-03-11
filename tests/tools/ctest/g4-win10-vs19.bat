@echo off 

rem - replace computer name by generic
if %COMPUTERNAME:~0,6% == G4-WIN SET CTEST_SITE=g4-win

echo starting for g4-win10-vs19

echo running on %COMPUTERNAME% using %Platform%.


set COMNTOOLS=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\Tools\
rem ---Compiler------------------------------------------------------
call "%COMNTOOLS%\..\..\VC\Auxiliary\Build\vcvarsall.bat" %Platform%

rem ---Xerces-C------------------------------------------------------

dir D:\sw\lcg\external\XercesC\3.2.2\%Platform%-windows-vc14

rem -- XercesC search first line has highest priority, as this ends up last  

set xerc_dirs=D:\sw\lcg\external\XercesC\3.2.3\%Platform%-windows-vc14
set xerc_dirs=D:\sw\lcg\external\XercesC\3.2.2\%Platform%-windows-vc14,%xerc_dirs%
set xerc_dirs=D:\sw\lcg\external\XercesC\3.2.2\x86-windows-vc14,%xerc_dirs%

for %%a IN (%xerc_dirs%) DO (
	if exist %%a   set XERCESC_ROOT_DIR=%%a
)         
echo Using XercesC from %XERCESC_ROOT_DIR%

set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

rem ---------------------------------------

rem - set CONFIG=%Platform%-win7-vc15

rem - location of this .bat file
set THIS=%~d0%~p0

rem ---Define basic config parameters--------------------------------
call %THIS%g4-win-common.bat

set G4_XOPTS=%G4_XOPTS%;-DGEANT4_INSTALL_DATASETS_TENDL=ON

echo Geant4 CMake options - 2 : %G4_XOPTS%


rem - echo %Path%
rem - python --version 
rem - set

set Path=%Path%;C:\Python27

rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake 
