@echo off 

rem - replace computer name by generic
if %COMPUTERNAME:~0,8% == G4STTWIN SET CTEST_SITE=g4sttwin

echo starting for g4-win10, compiler %COMPILER%

set Platform=x64

echo running on %COMPUTERNAME%

if "%COMPILER%"=="VC141" (
     call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat" %Platform%

) ELSE IF "%COMPILER%"=="VC14" (
     call "%VS140COMNTOOLS%\..\..\vc\vcvarsall" %Platform%
) ELSE (
     echo "Error: Not configured for %COMPILER%"
     exit
)

rem ---Xerces-C------------------------------------------------------

dir D:\sw\lcg\external\XercesC\3.1.3\%Platform%-windows-vc14

set xerc_dirs=D:\sw\lcg\external\XercesC\3.1.2\x86-windows-vc14
set xerc_dirs=D:\sw\lcg\external\XercesC\3.1.3\%Platform%-windows-vc14,%xerc_dirs%
set xerc_dirs=D:\sw\lcg\external\XercesC\3.2.3\%Platform%-windows-vc14,%xerc_dirs%

for %%a IN (%xerc_dirs%) DO (
	if exist %%a   set XERCESC_ROOT_DIR=%%a
)         
echo Using XercesC from %XERCESC_ROOT_DIR%

set PATH=%XERCESC_ROOT_DIR%\bin;%PATH%

rem ---------------------------------------

rem set CONFIG=%Platform%-win7-vc15

rem - location of this .bat file
set THIS=%~d0%~p0

rem ---Define basic config parameters--------------------------------
call %THIS%g4-win-common.bat

echo Geant4 CMake options - 2 : %G4_XOPTS%


rem - echo %Path%
rem - python --version 
rem - set
set Path=%Path%;C:\Python27

rem ---Run the CTest script-------------------------------------------
ctest -V -S %THIS%g4%MODE%.cmake 
