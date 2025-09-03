@echo off
echo Copying models and transformations to bin...
copy models\*.obj bin\ >nul 2>&1
copy capacitor_viz\transformations_*.txt bin\ >nul 2>&1

echo Running dynamic capacitor calculation...
cd bin
dynamic_capacitor.exe

echo.
echo Dynamic capacitor analysis complete!
echo Results saved in bin folder:
dir capacitance_*.csv /b 2>nul
cd ..
pause