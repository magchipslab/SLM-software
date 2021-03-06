@echo off  
rem INTELF11MSVS2008SHELLENGMATOPTS.BAT
rem
rem    Compile and link options used for building stand-alone engine or MAT
rem    programs using the Intel� Visual Fortran Compiler 11.1 with the
rem    Microsoft� Visual Studio� 2008 Shell linker.
rem    
rem    $Revision: 1.1.6.1 $  $Date: 2009/10/24 16:28:18 $
rem
rem ********************************************************************
rem General parameters
rem ********************************************************************
set MATLAB=%MATLAB%
set IFORT_COMPILER11=%IFORT_COMPILER11%
set VS90COMNTOOLS=%VS90COMNTOOLS%
set LINKERDIR=%VS90COMNTOOLS%\..\..
set VSINSTALLDIR=%LINKERDIR%
set VCINSTALLDIR=%VSINSTALLDIR%\VC
set VSAENVDIR=%CommonProgramFiles%\Microsoft Shared\VSA\9.0\VsaEnv
set PATH=%IFORT_COMPILER11%\Bin\ia32;%VSAENVDIR%;%VSINSTALLDIR%\Common7\IDE;%VCINSTALLDIR%\BIN;%VSINSTALLDIR%\Common7\Tools;%VSINSTALLDIR%\Common7\Tools\bin;%VCINSTALLDIR%\PlatformSDK\bin;%VCINSTALLDIR%\BIN;%MATLAB_BIN%;%PATH%
set INCLUDE=%IFORT_COMPILER11%\Include;%VCINSTALLDIR%\atlmfc\include;%VCINSTALLDIR%\include;%VCINSTALLDIR%\PlatformSDK\include;%INCLUDE%
set LIB=%IFORT_COMPILER11%\Lib\ia32;%VCINSTALLDIR%\lib;%VCINSTALLDIR%\PlatformSDK\lib;%LIB%
set MW_TARGET_ARCH=win32

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=ifort
set COMPFLAGS=/fpp /Qprec /I"%MATLAB%/extern/include" /c /nologo /fixed /fp:source /MD /assume:bscc /Qvc9
set OPTIMFLAGS=/O2 /DNDEBUG
set DEBUGFLAGS=/Z7
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win32\microsoft
set LINKER=link
set LINKFLAGS=/LIBPATH:"%LIBLOC%" libmx.lib libmat.lib libeng.lib /nologo /subsystem:console
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/debug /PDB:"%OUTDIR%%MEX_NAME%.pdb" /INCREMENTAL:NO
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%.exe"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=
set RC_LINKER=
set POSTLINK_CMDS1=mt -outputresource:"%OUTDIR%%MEX_NAME%.exe";1 -manifest "%OUTDIR%%MEX_NAME%.exe.manifest" 
set POSTLINK_CMDS2=del "%OUTDIR%%MEX_NAME%.exe.manifest" 
