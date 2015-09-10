:: Copy DLL to UCP program folders on S drive, so you can run the debug version from there
@ECHO OFF

xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\AOMs\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Example Program\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Experiment Monitor\" 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Fast Cam Capture\" 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\HV Matsusada\" 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\HV Ultravolt\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\MaximCam\" 		 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\MOT Correction Coils\" 		 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\OPA\" 		 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Picomotor controller\" 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Programmable Pattern Generator\" 
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\ScanController\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\DummyVariable\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\ShowStatusAllClients\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Scope\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Spectrometer\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Power Supply TTi\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\VacuumShift\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Webcam\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED
xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\UCP\Software\Wavelength Pulsed Dye\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED

xcopy /D /Y /F "S:\UCP\Software\General\General_DLL.dll" "S:\GENERAL EXPERIMENTAL\Mantis & Legend\1 - General\TDU\"
IF %errorlevel% NEQ 0 CALL :ERROROCCURED

GOTO:EOF

::===================================::
::   -   S u b r o u t i n e s   -   ::
::===================================::

:ERROROCCURED
ECHO.
ECHO Error in copying
PAUSE
GOTO:EOF