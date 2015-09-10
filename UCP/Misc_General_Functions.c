#include "toolbox.h"
#include <utility.h>
#include <userint.h>

/*------------------------------------------------------------------------------------------------------------------------
	Function		check_for_duplicate_program_running
	Arguments		0	Program is not running
					-1	Another instance of this program is already running
	Return value	-
					
	Description		Check if onother copy of the program is already running 
------------------------------------------------------------------------------------------------------------------------*/
char check_for_duplicate_program_running(void)
{
	int program_is_already_running;
	
	CheckForDuplicateAppInstance(ACTIVATE_OTHER_INSTANCE, &program_is_already_running);
	if (program_is_already_running)
	{
		MessagePopup("Error", "An other instance of this program is already running");
		return -1;
	}
	return 0;
}


/*------------------------------------------------------------------------------------------------------------------------
	Function		MainCallBack
					
	Description		-
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK MainCallBack (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	if (event == EVENT_END_TASK)
	{
		QuitUserInterface(0);
		return 1; //Swallow this event
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		install_main_callback
	Arguments		-
	Return value	-
					
	Description		Install a Main Callback function
					When the user closes the program by right-clicking -> close in the taskbar, this callback makes
					sure that the functions defined after RunUserInterface() are executed, which is not the case 
					when no Main Callback is installed
------------------------------------------------------------------------------------------------------------------------*/
void install_main_callback(void)
{
	InstallMainCallback(MainCallBack, 0, 0);
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		recall_panel_state
	Arguments		-
	Return value	-
					
	Description		Recall a previous panel state, saved in the file 'panelstate.dat'
------------------------------------------------------------------------------------------------------------------------*/
void recall_panel_state(void)
{
	if (FileExists("panelstate.dat",NULL)) 
		if (RecallPanelState(1, "panelstate.dat", 0) < 0) //Main panel has always a panelhandle of 1
		MessagePopup("Warning", "Error in restoring previous panel state");
}
	
/*------------------------------------------------------------------------------------------------------------------------
	Function		save_panel_state
	Arguments		-
	Return value	-
					
	Description		Save the current panel state to the file 'panelstate.dat'
------------------------------------------------------------------------------------------------------------------------*/
void save_panel_state(void)	
{
	SavePanelState(1, "panelstate.dat", 0); //Main panel has always a panelhandle of 1
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		save_panel_state
	Arguments		path_name	Path name of the directory to make
	Return value	-
					
	Description		Make a directory. The standard MakeDir function can only make one directory at a time
					E.g. if D:\ is empty it cannot make D:\Dir1\Subdir1 directly
------------------------------------------------------------------------------------------------------------------------*/
void make_directory(char *path_name)
{
	char *temp_path_string;
	int length_pathname, position_in_string=3;

	if ( !FileExists(path_name, NULL) ) 
	{
		//You can only create one folder at a time, so loop over directory tree
		length_pathname = strlen(path_name);
		temp_path_string = calloc(length_pathname+2, sizeof(char));
		while (1)
		{
			position_in_string += strcspn(path_name+position_in_string, "\\") + 1;
			strncpy(temp_path_string, path_name, position_in_string);
			temp_path_string[position_in_string] = 0; //terminate string
			if ( !FileExists(temp_path_string, NULL) ) //create folder if necessary
				if (MakeDir(temp_path_string) != 0)
				{
					MessagePopup("Error", "Error in creating directory");
					break;
				}		  
			if (position_in_string>=length_pathname)
				break;
		}
		free(temp_path_string);
	}				  
}
	



	
