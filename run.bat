@echo off
echo Copying models to bin...
copy models\*.obj bin\ >nul 2>&1
echo Running capacitor calculation...
cd bin
ray_distance.exe -v
move ray_data.txt .. >nul 2>&1
cd ..
echo.
echo Opening visualization...
cd capacitor_viz
call capacitor_venv\Scripts\activate.bat
py capacitor_viewer.py
cd ..
pause