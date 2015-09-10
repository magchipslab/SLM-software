/*------------------------------------------------------------------------------------------------------------------------
 *  File 		MatlabFunctions.c
 *  Author		Wouter Engelen
 *	Version		October 2009
 *		
 *  Details		Functions for calling Matlab from Labwindows
 *
------------------------------------------------------------------------------------------------------------------------*/
#include "Matlab.h"
#include "toolbox.h"
#include <userint.h>
#include <ansi_c.h>
#include <utility.h>
#include <formatio.h>
#include "General_DLL.h"

#define MATLAB_SCRIPT_DIR "Matlab"
#define MATLAB_DISPLAY_RESULTS_COLUMNS 2
#define MATLAB_DISPLAY_RESULTS_ROWS 4

CAObjHandle MatlabHandle;
int panelHandle_GUI, ControlID_Script_Ring, ControlID_Results_Numeric;
int Matlab_script_number_of_variables;
int *Matlab_script_gui_controlID;
char Matlab_script_name[200];

/*------------------------------------------------------------------------------------------------------------------------
 *  How to create the MATLAB Automation Server Type Library
 *  
 *	Remove the current MATLAB instrument and rebuild the MATLAB Automation Server Type Library by going to Tools»Create ActiveX Controller.
 *	The name of the MATLAB instrument will be MATLAB Application Type Library. 
 *	
 *	Open the source code for the MATLAB instrument and search and replace ALL instances of &MLApp_IID_DIMLApp with &IID_IDispatch.
 *
------------------------------------------------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------------------------------------------------
 *  Know bugs
 *  
 *	When changing an m file after the Automation server is started, the old m file is executed. After restarting the
 *  automation server the new m file is used. Also, when after changing a script, the script is called from the
 *	command window manually, the new version is used, and thereafter the new version is also used in the program
 *	
------------------------------------------------------------------------------------------------------------------------*/




/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Functions for calling Matlab
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ErrFunc
	Arguments		errCode 	the Automation error code
					errMesg 	the error message that Matlab returns in some functions
	Return value	-
					
	Description		Reports the Automation error associated with the code.
------------------------------------------------------------------------------------------------------------------------*/
void ErrFunc(HRESULT errCode, char *errMesg)
{
	char errStr[500];
    
    CA_GetAutomationErrorString(errCode, errStr, 500);
    message_window_write(errStr);
	
    if (errMesg)
    {
        message_window_write(errMesg);
        CA_FreeMemory(errMesg);
    }
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabRunCommand
	Arguments		command		command
	Return value	-
					
	Description		Executes a Matlab command
------------------------------------------------------------------------------------------------------------------------*/
void MatlabRunCommand(char *command)
{
    HRESULT stat;
	char *matlabMesg = NULL;
	
	stat = MLApp_DIMLAppExecute(MatlabHandle, NULL, command, &matlabMesg);
    if (stat < 0)
	{
		message_window_write("Error in MatlabRunCommand:");
		ErrFunc(stat, NULL);
	}
	
	//if (matlabMesg)
	if (matlabMesg[0] != 0)
	{
		message_window_write(matlabMesg);
		CA_FreeMemory(matlabMesg);
	}
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabStart
	Arguments		-
	Return value	0	Matlab started
					-1	Matlab failed to start
					
	Description		Start Matlab
------------------------------------------------------------------------------------------------------------------------*/
char MatlabStart(void)
{
	HRESULT stat;
	
	stat = MLApp_NewDIMLApp(NULL, 1, LOCALE_NEUTRAL, 0, &MatlabHandle);
    if (stat < 0) 
	{
		message_window_write("Error in MatlabStart:");
		ErrFunc(stat, NULL);
		return -1;
	}
	MLApp_DIMLAppSetVisible(MatlabHandle, NULL, 1); //Make Matlab window visible
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabStop
	Arguments		-
	Return value	0	Success
					-1	Error
					
	Description		Stop Matlab
------------------------------------------------------------------------------------------------------------------------*/
char MatlabStop(void)
{
	HRESULT stat;
	stat = MLApp_DIMLAppQuit(MatlabHandle, NULL);
    if (stat < 0)
	{
		message_window_write("Error in MatlabStop:");
		ErrFunc(stat, NULL);
		return -1;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabGetVariable
	Arguments		variable_name	Variable name
			   Out	variable_value	Variable value
	Return value	0	Success
					-1	Error
					
	Description		Get the value of a variable
------------------------------------------------------------------------------------------------------------------------*/
char MatlabGetVariable(char *variable_name, double *variable_value)
{
	HRESULT stat;
	VARIANT var;
	
	stat = MLApp_DIMLAppGetWorkspaceData(MatlabHandle, NULL, variable_name, "base", &var);
	if (stat < 0)
	{
		message_window_write("Error in MatlabGetVariable:");
		ErrFunc(stat, NULL);
		return -1;
	}
	CA_VariantGetDouble(&var, variable_value);
	return 0;
}
	
/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabGetMatrix
	Arguments		variable_name	Variable name
			   Out	matrixReal		Real part of the matrix.                          
			   Out	matrixImag 		Imaginary part of the matrix.                     
			   Out	rows 			Number of rows in the matrixs.                         
			   Out	columns			Number of columns in the matrix.                      
	Return value	-
					
	Description		Get the value of a matrix
------------------------------------------------------------------------------------------------------------------------*/
void MatlabGetMatrix(char *variable_name, double **matrixReal, double **matrixImag, unsigned *rows, unsigned *columns)
{
    SAFEARRAY *saReal = NULL;
    SAFEARRAY *saImag = NULL;
    HRESULT   stat;

    // Get the matrix's real and imaginary parts from MATLAB
    stat = MLApp_DIMLAppGetFullMatrix(MatlabHandle, NULL, variable_name, "base", &saReal, &saImag);
    if (stat < 0)
    {
        message_window_write("Error in MatlabGetMatrix:");
		ErrFunc(stat, NULL);
		*rows = -1;
		*columns = -1;
		return;
    }
    
    // Convert the SAFEARRAYs to C datatypes
    if ((saReal == NULL) && (saImag == NULL))
    {
        message_window_write("Error in MatlabGetMatrix: Empty array");
		*rows = -1;
		*columns = -1;
		return;
    }
    
    if (saReal)
    {
        stat = CA_SafeArrayGet2DSize(saReal, rows, columns);
        if (stat < 0)
        {
            message_window_write("Error in MatlabGetMatrix:");
			ErrFunc(stat, NULL);
			*rows = -1;
			*columns = -1;
			return;
        }
        stat = CA_SafeArrayTo2DArray((LPSAFEARRAY *)&saReal, CAVT_DOUBLE, matrixReal, rows, columns);
        if (stat < 0)
        {
            message_window_write("Error in MatlabGetMatrix:");
			ErrFunc(stat, NULL);
			*rows = -1;
			*columns = -1;
			return;
        }
        
    }
    if (saImag)
    {												  
        stat = CA_SafeArrayGet2DSize(saImag, rows, columns);
        if (stat < 0)
        {
            message_window_write("Error in MatlabGetMatrix:");
			ErrFunc(stat, NULL);
			*rows = -1;
			*columns = -1;
			return;
        }
        stat = CA_SafeArrayTo2DArray((LPSAFEARRAY *)&saImag, CAVT_DOUBLE, matrixImag, rows, columns);
        if (stat < 0)
        {
            message_window_write("Error in MatlabGetMatrix:");
			ErrFunc(stat, NULL);
			*rows = -1;
			*columns = -1;
			return;
        }
    }
}

/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Functions setting gui and reading matlab script settings files
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_strlen+33
	Arguments		panelHandle		Panel Handle
					ControlID		Control ID of the Ring to set the script names in
	Return value	-
					
	Description		Fill filelist with names of matlab scripts in [MATLAB_SCRIPT_DIR] folder
------------------------------------------------------------------------------------------------------------------------*/
void Create_MATLAB_Script_List(int panelHandle, int ControlID)
{
	char filename[MAX_FILENAME_LEN +1], filename2[MAX_FILENAME_LEN +1];
	int result;
	int i=0, number_of_files = 0;
	
	ClearListCtrl(panelHandle, ControlID);
	result = GetFirstFile(".\\" MATLAB_SCRIPT_DIR "\\*.txt", 1, 0, 0, 0, 0, 0, filename);
	if (result == -1)
	{
		message_window_write("Warning: no script files found");
		return;
	}
	
	while (result==0)
	{
		filename[strlen(filename)-4] = 0; //remove extension
		for (i=0; i<number_of_files; i++) //List filenames alphabettically
		{
			GetLabelFromIndex(panelHandle, ControlID, i, filename2);
			if ( strcmp(filename, filename2) < 0) //Filename becomes before current Filename in alphabet
			{
				InsertListItem(panelHandle, ControlID, i, filename, i); //Insert before filename2 in list
				break;
			}
		}
		if (i==number_of_files) 
			InsertListItem(panelHandle, ControlID, -1, filename, i); //Add to end of list
		number_of_files++;
		result = GetNextFile(filename);
	}
	SetCtrlIndex(panelHandle, ControlID, 0);
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		ResetMatlabGui
	Arguments		-
	Return value	-
	
	Description		Remove controls for showing script return values from panel
------------------------------------------------------------------------------------------------------------------------*/
void ResetMatlabGui(void)
{
	int i;
	
	if (Matlab_script_number_of_variables == 0) return;
	SetCtrlAttribute(panelHandle_GUI, ControlID_Results_Numeric, ATTR_VISIBLE, 0);
	SetCtrlVal(panelHandle_GUI, ControlID_Results_Numeric, 0.0);
	
	for (i=1; i<Matlab_script_number_of_variables; i++)
		DiscardCtrl(panelHandle_GUI, Matlab_script_gui_controlID[i]);

	Matlab_script_number_of_variables = 0;
	free(Matlab_script_gui_controlID);
	Matlab_script_gui_controlID = NULL;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		ReadMatlabSettingsFromFile
	Arguments		FileName	File name
	Return value	-
	
	Description		Read matlab script settings from file
------------------------------------------------------------------------------------------------------------------------*/
void ReadMatlabSettingsFromFile(char *FileName)
{
	FILE *file; 
	char line[201];
	int i, top_position = -1, left_position = -1;
	
	ResetMatlabGui();
	
	file = fopen(FileName, "r"); //open the file for reading     
	fgets(line, 200, file);
	if ( strcmp(line, "MaximCam_MatlabFileDescriptor\n") != 0 )
	{
		message_window_write("Error reading script settings file");
		fclose(file); //close the file
		return;
	}
	
	if (fscanf(file, "%200s\t%i\n", line, &Matlab_script_number_of_variables) != 2)
	{
		message_window_write("Error reading script settings file");
		fclose(file); //close the file
		Matlab_script_number_of_variables = 0;
		return;
	}
	if ( strcmp(line, "Number_of_return_values") != 0 )
	{
		message_window_write("Error reading script settings file");
		fclose(file); //close the file
		Matlab_script_number_of_variables = 0;
		return;
	}
	if (Matlab_script_number_of_variables > MATLAB_DISPLAY_RESULTS_COLUMNS*MATLAB_DISPLAY_RESULTS_ROWS)
	{
		sprintf(line, "Warning: maximum %i variables allowed to be displayed in panel. Skipping the rest", MATLAB_DISPLAY_RESULTS_COLUMNS*MATLAB_DISPLAY_RESULTS_ROWS);
		message_window_write(line);
		Matlab_script_number_of_variables = MATLAB_DISPLAY_RESULTS_COLUMNS*MATLAB_DISPLAY_RESULTS_ROWS;
	}

	Matlab_script_gui_controlID = calloc(Matlab_script_number_of_variables, sizeof(int));
	
	fgets(line, 200, file); //read line Return_value_names
	if ( strcmp(line, "Return_value_names\n") != 0 )
	{
		message_window_write("Error reading script settings file");
		fclose(file); //close the file
		Matlab_script_number_of_variables = 0;
		return;
	}
	
	for (i=0; i<Matlab_script_number_of_variables; i++)
	{
		if (fscanf(file, "%200s\n", line) != 1)
		{
			message_window_write("Error reading script settings file");
			fclose(file); //close the file
			Matlab_script_number_of_variables = 0;
			return;
		}
		if (i==0)
		{
			SetCtrlAttribute(panelHandle_GUI, ControlID_Results_Numeric, ATTR_LABEL_TEXT, line);
			GetCtrlAttribute(panelHandle_GUI, ControlID_Results_Numeric, ATTR_TOP,  &top_position);
			GetCtrlAttribute(panelHandle_GUI, ControlID_Results_Numeric, ATTR_LEFT, &left_position);
			Matlab_script_gui_controlID[0] = ControlID_Results_Numeric;
			SetCtrlAttribute(panelHandle_GUI, ControlID_Results_Numeric, ATTR_VISIBLE, 1);
		}																
		else
			Matlab_script_gui_controlID[i] = DuplicateCtrl(panelHandle_GUI, ControlID_Results_Numeric, panelHandle_GUI, line, top_position + 25*(i%MATLAB_DISPLAY_RESULTS_ROWS), left_position + ((i/MATLAB_DISPLAY_RESULTS_ROWS)*120) );
		
	}
	
	if (fscanf(file, "%200s\n", line) != -1)
		message_window_write("Error reading script settings file: extra information found");
	
	fclose(file); //close the file 
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		SelectMatlabScript
	Arguments		-
	Return value	-
	
	Description		Select the active MATLAB Script in the ring index
------------------------------------------------------------------------------------------------------------------------*/
void SelectMatlabScript(void)
{
	char MATLAB_filename[300];
	int	script_index;
	
	GetCtrlIndex(panelHandle_GUI, ControlID_Script_Ring, &script_index);
	GetLabelFromIndex(panelHandle_GUI, ControlID_Script_Ring, script_index, Matlab_script_name);

	//Check if m file with same name exists
	sprintf(MATLAB_filename, ".\\%s\\%s.m", MATLAB_SCRIPT_DIR, Matlab_script_name);
	if ( FileExists(MATLAB_filename, 0) != 1)
	{
		char message[500];
		sprintf(message, "Error: Matlab script '%s' not found", Matlab_script_name);
		message_window_write(message);
		return;
	}
	
	//read settings from txt file
	sprintf(MATLAB_filename, ".\\%s\\%s.txt", MATLAB_SCRIPT_DIR, Matlab_script_name);
	ReadMatlabSettingsFromFile(MATLAB_filename);
	
	MatlabRunCommand("clear all");

	//check if matlab init file is present and run this
	sprintf(MATLAB_filename, ".\\%s\\%s_init.m", MATLAB_SCRIPT_DIR, Matlab_script_name);
	if (FileExists(MATLAB_filename, 0) == 1)
	{
		sprintf(MATLAB_filename, "%s_init", Matlab_script_name);
		MatlabRunCommand(MATLAB_filename);
	}
	
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		StartStopMatlabServer
	Arguments		start_matlab			1	Start matlab
											0	Stop matlab
					panelHandle				Panel Handle
					ControlIDScriptRing		Control ID of the Ring to set the script names in
					ControlIDResultsNumeric	Control ID of the Numeric control to show the script results in
	Return value	0	Success
					-1	Error
	
	Description		Start / Stop the Matlab server
					In the panel you need
					- A ring control (initially hidden) where the available matlab scripts are shown
					- A numeric control (initially hidden) to show the results of the script
					  There should be space in the panel to copy the control MATLAB_DISPLAY_RESULTS_ROWS-1
					  times below the current location, and MATLAB_DISPLAY_RESULTS_COLUMS-1 times right 
					  of the current location (+ x times below that)
------------------------------------------------------------------------------------------------------------------------*/
char StartStopMatlabServer(char start_matlab, int panelHandle, int ControlIDScriptRing, int ControlIDResultsNumeric)
{
	char return_val;
	char project_path[MAX_PATHNAME_LEN];
	char command[MAX_PATHNAME_LEN];
	
	GetDir(project_path);
	
	if (start_matlab)
	{
		SetWaitCursor(1);
		
		panelHandle_GUI = panelHandle;
		ControlID_Script_Ring = ControlIDScriptRing;
		ControlID_Results_Numeric = ControlIDResultsNumeric;
		
		return_val = MatlabStart();
		if (return_val == 0) //start succesfull
		{
			//Change matlab current dir
			sprintf(command, "cd \'%s\\%s\'", project_path, MATLAB_SCRIPT_DIR);
			MatlabRunCommand(command);

			Create_MATLAB_Script_List(panelHandle_GUI, ControlID_Script_Ring);
			SelectMatlabScript();
			ProcessDrawEvents();
			SetCtrlAttribute(panelHandle_GUI, ControlID_Script_Ring, ATTR_VISIBLE, 1);
		}
		SetWaitCursor(0);
	}
	else //stop matlab
	{
		SetCtrlAttribute(panelHandle_GUI, ControlID_Script_Ring, ATTR_VISIBLE, 0);  
		MatlabRunCommand("close all"); //Close all open figures (if any)
		ResetMatlabGui();
		
		//Copy matlab files to S drive
		sprintf(command, ".\\%s\\CopyScripts.bat \"%s\" %s", MATLAB_SCRIPT_DIR, project_path, MATLAB_SCRIPT_DIR);
		if ( LaunchExecutableEx(command, LE_HIDE, NULL) )
				message_window_write("Error in copying matlab files");
		
		return_val = MatlabStop();
	}
	return return_val;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabShowResults
	Arguments		-
	Return value	-
	
	Description		Show return value(s) of the script in the GUI
------------------------------------------------------------------------------------------------------------------------*/
void MatlabShowResults(void)
{
	double   *matrixReal = NULL;
    double   *matrixImag = NULL;
    unsigned rows        = 0;
    unsigned columns     = 0;
	int i;
	
	if (Matlab_script_number_of_variables==0) return; //nothing to show
	
	MatlabGetMatrix("ans", &matrixReal, &matrixImag, &rows, &columns);
	
	if (rows == 1 && columns == Matlab_script_number_of_variables) //Script succesfull
	{
		for (i=0; i<Matlab_script_number_of_variables ; i++)
			SetCtrlVal(panelHandle_GUI, Matlab_script_gui_controlID[i], matrixReal[i]);
	}
	else //Script failed
	{
		for (i=0; i<Matlab_script_number_of_variables ; i++)
			SetCtrlVal(panelHandle_GUI, Matlab_script_gui_controlID[i], 0.0);
	}
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabRunScript
	Arguments		arguments	Arguments to pass to script
	Return value	-
					
	Description		Call a script with (optional) arguments and show results.
------------------------------------------------------------------------------------------------------------------------*/
void MatlabRunScript(char *arguments)
{
	char command[704];
	if (arguments)	
	{
		sprintf(command, "%.200s(%.200s);", Matlab_script_name, arguments); // script_name([arguments]);
		MatlabRunCommand(command); 
	}
	else
		MatlabRunCommand(Matlab_script_name);
	
	MatlabShowResults();
}



/*
void MatlabExportArray(char *name, double *array, int len)
{
	LPSAFEARRAY SafeA; 
	CA_Array1DToSafeArray (array, CAVT_DOUBLE, len, &SafeA);            
	MLApp_DIMLAppPutFullMatrix (MatlabHandle, NULL, name, "base", SafeA, NULL);
	CA_SafeArrayDestroy(SafeA);
}

void MatlabExportArray2D(char *name, double **array, int len1, int len2)
{
	LPSAFEARRAY SafeA;  
	CA_Array2DToSafeArray (array, CAVT_DOUBLE, len1, len2, &SafeA); 
	MLApp_DIMLAppPutFullMatrix (MatlabHandle, NULL, name, "base", SafeA, NULL);
	CA_SafeArrayDestroy(SafeA);
}

void MatlabExportVar(char *name, double value)
{
	VARIANT var;
	CA_VariantSetDouble (&var, value);
	MLApp_DIMLAppPutWorkspaceData (MatlabHandle, NULL, name, "base", var);
	
}

void MatlabImportString(char *name, char *string)
{
	VARIANT var;
	MLApp_DIMLAppGetWorkspaceData (MatlabHandle, NULL, name, "base", &var);
	CA_VariantGetCStringBuf (&var, string,CMATLAB_MAXSTRLENGTH );
}
*/

/*---------------------------------------------------------------------------*/
/* This function minimizes/maximizes the MATLAB window.                      */
/* MatlabHandle is the handle to the MATLAB Application Object.                   */
/* Pass 0 in minmaxFlag to minimize, and 1 to maximize the window.           */
/*---------------------------------------------------------------------------*/
/*
int MinMaxMatlab(CAObjHandle MatlabHandle, int minmaxFlag)
{
    HRESULT     stat    =   0;
    
    if (minmaxFlag == 0)
    {
        // Minimize the MATLAB window
        stat = MLApp_DIMLAppMinimizeCommandWindow (MatlabHandle, NULL);
        if (stat < 0) return -1;             
    }
    else
    {
        // Maximize the MATLAB window
        stat = MLApp_DIMLAppMaximizeCommandWindow (MatlabHandle, NULL);
        if (stat < 0) return -1;  
    }
    return 0;

}
  */

/*---------------------------------------------------------------------------*/
/* This function executes a MATLAB .m file.                                  */
/* MatlabHandle is the handle to the MATLAB Application Object.                   */
/* mFilePath is the full path to the .m file to run.                         */
/*---------------------------------------------------------------------------*/

/* 
int RunMatlabScript(CAObjHandle MatlabHandle, char *mFilePath)
{
    int     index                           =   -1;
    char    fName[MAX_FILENAME_LEN];
    char    dirName[MAX_DIRNAME_LEN];
    char    driveName[MAX_DRIVENAME_LEN];
    char    command[300];
    int     result                          =   0;

    /* Split the path name into the drive, directory, and file names
    SplitPath (mFilePath, driveName, dirName, fName);
    if (driveName && dirName && fName)
    {
        /* Add the m-file directory to MATLAB's search path
        Fmt(command,"%s<path(path,'%s%s');",driveName,dirName); 
        RunMatlabCommand(MatlabHandle, command);   
        
        /* Copy the .m file name into the command string 
        Fmt(command,"%s<%s",fName);
        /* Remove the ".m" extension from the command 
        index=strlen(command)-2;
        command[index]='\0';
        
        /* Execute the .m file 
        result = RunMatlabCommand(MatlabHandle, command);
        return result;
    }
    else
    {
        return ERROR_FILE_PATH;
    }
                
}

/*---------------------------------------------------------------------------*/
/* This function sends a string to MATLAB       .                            */
/* MatlabHandle is the handle to the MATLAB Application Object.                   */
/* matStringName is the name of the MATLAB variable to store the  string in. */
/* CVIstring is the actual string to be transferred to MATLAB                */
/* MATLAB does not support BSTRs, and so we send the string as a double [].  */
/*---------------------------------------------------------------------------*/
/* 
int SendString(CAObjHandle MatlabHandle, char *matStringName, char *CVIString)
{
    double          *buffer_r         =   NULL;
    double          *buffer_i         =   NULL;
    char            command[300];
    int             index           =   0;
    int             i               =   0;
    int             result          =   0;
    
    if(CVIString == NULL)
    {
        return ERROR_EMPTY_STRING;
    }
    index=strlen(CVIString);
    /* Convert the string into a double array 
    buffer_r=(double *)malloc(index*sizeof(double));
    buffer_i=(double *)malloc(index*sizeof(double));
    for(i=0;i<index;++i)
    {
        buffer_r[i]=(double)CVIString[i];
        buffer_i[i]=0.0;
    }
        
    
    /* Send the double array to MATLAB 
    result = SendMatrix(MatlabHandle, matStringName, buffer_r, buffer_i, index, 1);
    free(buffer_r);
    free(buffer_i);
    if (result != SUCCESS) return result;
    
    /* Convert the double array into a MATLAB string 
    Fmt(command,"%s<%s=transpose(%s)",matStringName,matStringName);
    result = RunMatlabCommand(MatlabHandle, command);
    if (result != SUCCESS) return result;
    Fmt(command,"%s<%s=char(%s)",matStringName,matStringName);
    result = RunMatlabCommand(MatlabHandle, command);
    if (result != SUCCESS) return result;
    
    return SUCCESS;
}

/*---------------------------------------------------------------------------*/
/* This function gets a string from MATLAB      .                            */
/* MatlabHandle is the handle to the MATLAB Application Object.                   */
/* matStringName is the name of the MATLAB variable to get the  string from. */
/* The transferred string will be stored in *cString.                        */
/* MATLAB does not support BSTRs, and so we get the string as a double [].   */
/*---------------------------------------------------------------------------*/
/* 
int GetString(CAObjHandle MatlabHandle, char *matStringName, char **cString)
{
    double          *buffer     =   NULL;
    double          *dummy      =   NULL;
    unsigned        dim1        =   0;
    unsigned        dim2        =   0;
    int             i           =   0;
    char            command[300];
    int             result      =   0;
    
    /* convert the MATLAB string into a double matrix 
    Fmt(command,"%s<CVIString=transpose(double(%s));",matStringName);
    result = RunMatlabCommand(MatlabHandle,command);
    if (result != SUCCESS) return result;

    /* get the matrix into CVI and convert into string 
    result = GetMatrix(MatlabHandle, "CVIString", &buffer, &dummy, &dim1, &dim2);
    if (result != SUCCESS) return result;
    if (buffer == NULL)
    {
        return ERROR_NULL_POINTER;
    }
    
    *cString=(char *)malloc((dim1+1)*sizeof(char));
    for(i=0;i<dim1;++i)
        (*cString)[i]=(char)buffer[i];
    (*cString)[dim1]='\0';
    if (buffer != NULL) CA_FreeMemory(buffer);  
    if (dummy != NULL) CA_FreeMemory(dummy);
    
    /* clear the temporary MATLAB variable created 
    result = RunMatlabCommand(MatlabHandle, "clear CVIString;");
    if (result != SUCCESS) return result;
    return SUCCESS; 
}

/*---------------------------------------------------------------------------*/
/* This function sends a matrix to MATLAB.                                   */
/* MatlabHandle is the handle to the MATLAB Application Object.                   */
/* matlabName is the MATLAB variable that will contain the matrix.           */
/* matrixReal contains the real part of the matrix.                          */
/* matrixImag contains the imaginary part of the matrix.                     */
/* dim1 is the number of rows in the matrix.                                 */
/* dim2 is the number of columns in the matrix.                              */
/*---------------------------------------------------------------------------*/
/* 
int SendMatrix(CAObjHandle MatlabHandle, char *matlabName, double *matrixReal, 
                    double *matrixImag, unsigned dim1, unsigned dim2)
{
    LPSAFEARRAY         saReal  =   NULL;
    LPSAFEARRAY         saImag  =   NULL;
    HRESULT             stat;
    
    /* Check that valid pointers have been passed for the matrices
    if((!matrixReal) && (!matrixImag))
    {
        return ERROR_INVALID_MATRIX;
    }
    
    /* Convert the real and imaginary matrices into SAFEARRAYs 
    if(matrixReal)
    {
        stat = CA_Array2DToSafeArray (matrixReal, CAVT_DOUBLE, dim1, 
                dim2, &saReal);
        if (stat < 0)
        {
            return ERROR_ARRAY_CONVERSION;
        }
    }
    if(matrixImag)
    {
        stat = CA_Array2DToSafeArray (matrixImag, CAVT_DOUBLE, dim1, 
                dim2, &saImag);
        if (stat < 0)
        {
            return ERROR_ARRAY_CONVERSION;
        }
    }
    
    /* Send the matrices into MATLAB 
    stat = MLApp_DIMLAppPutFullMatrix (MatlabHandle, NULL, matlabName, "base",
                                (SAFEARRAY *)saReal, (SAFEARRAY *)saImag);
    if (stat < 0)
    {
        return ERROR_MATRIX_TRANSFER;
    }
    if (saReal != NULL) CA_SafeArrayDestroy(saReal);    
    if (saImag != NULL) CA_SafeArrayDestroy(saImag);
    return SUCCESS;
} */



