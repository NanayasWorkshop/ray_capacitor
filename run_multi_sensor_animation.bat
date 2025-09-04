@echo off
echo Copying models and transformations to bin...
copy models\*.obj bin\ >nul 2>&1
copy capacitor_viz\transformations_*.txt bin\ >nul 2>&1

echo Expected model files:
echo   - A1_model.obj, A2_model.obj
echo   - B1_model.obj, B2_model.obj  
echo   - C1_model.obj, C2_model.obj
echo   - stationary_negative.obj

echo.
echo Checking for model files in bin...
if exist "bin\A1_model.obj" echo   ✓ A1_model.obj found
if not exist "bin\A1_model.obj" echo   ✗ A1_model.obj missing
if exist "bin\A2_model.obj" echo   ✓ A2_model.obj found
if not exist "bin\A2_model.obj" echo   ✗ A2_model.obj missing
if exist "bin\B1_model.obj" echo   ✓ B1_model.obj found
if not exist "bin\B1_model.obj" echo   ✗ B1_model.obj missing
if exist "bin\B2_model.obj" echo   ✓ B2_model.obj found
if not exist "bin\B2_model.obj" echo   ✗ B2_model.obj missing
if exist "bin\C1_model.obj" echo   ✓ C1_model.obj found
if not exist "bin\C1_model.obj" echo   ✗ C1_model.obj missing
if exist "bin\C2_model.obj" echo   ✓ C2_model.obj found
if not exist "bin\C2_model.obj" echo   ✗ C2_model.obj missing
if exist "bin\stationary_negative.obj" echo   ✓ stationary_negative.obj found
if not exist "bin\stationary_negative.obj" echo   ✗ stationary_negative.obj missing

echo.
echo Checking for transformation files in bin...
if exist "bin\transformations_A.txt" echo   ✓ transformations_A.txt found
if not exist "bin\transformations_A.txt" echo   ✗ transformations_A.txt missing
if exist "bin\transformations_B.txt" echo   ✓ transformations_B.txt found
if not exist "bin\transformations_B.txt" echo   ✗ transformations_B.txt missing
if exist "bin\transformations_C.txt" echo   ✓ transformations_C.txt found
if not exist "bin\transformations_C.txt" echo   ✗ transformations_C.txt missing

echo.
echo Running multi-sensor capacitor calculation with ANIMATION EXPORT...
cd bin
multi_sensor_capacitor.exe -anim

echo.
echo Multi-sensor capacitor analysis complete!
echo Results saved in bin folder:
dir capacitance_*.csv /b 2>nul
echo.
echo Animation files generated:
dir ray_data_*_step_*.txt /b 2>nul | find /c ".txt" >temp_count.txt
set /p file_count=<temp_count.txt
del temp_count.txt
echo   %file_count% ray data files created
echo.
echo Expected output files:
echo   - capacitance_A1.csv, capacitance_A2.csv
echo   - capacitance_B1.csv, capacitance_B2.csv  
echo   - capacitance_C1.csv, capacitance_C2.csv
echo   - ray_data_[sensor]_step_[XXX].txt files
echo.
echo To view animation:
echo   python multi_sensor_animated_viewer.py animate
echo To view specific step:
echo   python multi_sensor_animated_viewer.py [step_number]
cd ..
pause