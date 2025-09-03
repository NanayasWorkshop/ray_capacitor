@echo off
cd ..\bin
ray_distance.exe -v
move ray_data.txt ..\output\
cd ..\capacitor_viz
call capacitor_venv\Scripts\activate.bat
py capacitor_viewer.py
cd ..
pause