:: Copy library of Matlab Scripts to local drive or
:: Backup Matlab Scripts to S drive
@ECHO OFF

:: First argument: 	Path of program 		e.g. D:\UCP Programs\ScanController
:: Second argument:	Name of script directory	e.g. Matlab
:: Third argument:	Mode: 1 for copying LIB from S drive to local drive, 2 for copying local script files to S drive

:: /MIR Mirrors a directory tree 
:: /XD Excludes dir with the specified names
:: /XO Excludes files tagged as “Older”
:: /R:0 Specifies the number of retries on failed copies
:: /W:5 Specifies the wait time between retries
:: /XF bat.log Excludes files with the specified names
:: /LOG:bat.log Redirects output to the specified file, overwriting the file if it already exists.

:: Change current director, so robocopy exe can be found and log file opened on an error
cd /D %1\%2

IF %3%==1 (
	REM CALL robocopy.exe "S:\bunches\cvi\Matlab Scripts\Analyze Data\LIB" %1\%2\LIB /MIR /XD .svn /XO /R:0 /W:5 /XF bat.log /LOG:bat.log
) ELSE (
	REM CALL robocopy.exe %1\%2 S:\bunches\cvi\ScanController\%2 /MIR /XO /R:0 /W:5 /XD LIB /XF bat.log  robocopy.exe robocopy.doc CopyScripts.bat /LOG:bat.log
)
:: If an error has occured, open log, else delete log
IF errorlevel 16 CALL :ERROROCCURED
IF errorlevel 8  CALL :ERROROCCURED
IF errorlevel 4  CALL :ERROROCCURED
DEL bat.log

GOTO:EOF

::===================================::
::   -   S u b r o u t i n e s   -   ::
::===================================::

:ERROROCCURED
ECHO.
ECHO Error in copying
bat.log
GOTO:EOF