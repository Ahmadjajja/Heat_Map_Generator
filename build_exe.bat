@echo off
echo Building Montana Heat Map Generator...
echo.

REM Clean previous builds
if exist "build" rmdir /s /q "build"
if exist "dist" rmdir /s /q "dist"

REM Build the executable using the spec file
pyinstaller --clean montana_heatmap.spec

echo.
echo Build completed!
echo Executable location: dist\Montana Heat Map Generator.exe
echo.
pause 