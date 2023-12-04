@ECHO ON
COLOR A

CALL "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022

SET MyRepository=%~dp0
SET MyProjectBin=%~dp0\..\FEM-build

SET CMAKE_EXE="C:\Program Files\CMake\bin\cmake.exe"
SET CMAKE_EXE_GUI="C:\Program Files\CMake\bin\cmake-gui.exe"
SET CMAKE_GENERATOR_NAME="Visual Studio 17 2022"

%CMAKE_EXE% -G%CMAKE_GENERATOR_NAME% -Ax64 -S%MyRepository% -B%MyProjectBin% 

PAUSE
