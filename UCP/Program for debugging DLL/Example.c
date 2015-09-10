#include <userint.h>

/*------------------------------------------------------------------------------------------------------------------------
 *  File 		Example.c
 *  Author		Wouter Engelen
 *	Version		September 2010
 *		
 *  Details		Example program that shows how to connect to 
 				- Scancontroller: publish variables to controller and deal with changes made to these variables
				- Experiment monitor: publish variables to Experiment monitor
 *
------------------------------------------------------------------------------------------------------------------------*/
#define USE_MESSAGE_WINDOW 

// Include files
#include <ansi_c.h>
#include <utility.h>
#include "toolbox.h"
#include "Example.h"
//#include "ExpMonitor_Client.h"
//#include "ScanController_Client.h"

#include "General_DLL.h"

//#include "message_window_functions.h"

//Global variables
static int panelHandle;     

char Connected_To_ScanController;

//Variables that are published to Scancontroller
double ScanController_ExampleVar;

//Variable handle for publishing variables to the Experiment Monitor
int ExpMonitorVar;

char Analyze_In_Matlab = 0;

/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Code to control hardware
 *
-----------------------------------------------------------------------------------------------------------------------*/

void test(void)
{
	//make_directory("D:\\UCP\\Data\\12343\\");
	make_directory("D:\\UCP\\Data\\12343");
	
}

// Update hardware setpoints
char Update_value(void)
{
	SetCtrlVal(panelHandle, PANEL_SCANCONTR_NUMERIC, ScanController_ExampleVar); //write updated value to the panel
	return 1;
}


/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Experiment Monitor code
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ExpMonitor_Add_Variables
	Arguments		-
	Return value	-
					
	Description		Add variable(s) on to the Experiment monitor
------------------------------------------------------------------------------------------------------------------------*/
void ExpMonitor_Add_Variables(void)
{
	ExpMonitor_Client_Add_Var("Example_Prog", "Example_ExpMonVar", &ExpMonitorVar, panelHandle, PANEL_EXP_MON);
}


/*------------------------------------------------------------------------------------------------------------------------
	Function		ExpMonitor_Del_All_Variables
	Arguments		-
	Return value	-
					
	Description		Remove all variable(s) from the Experiment monitor
------------------------------------------------------------------------------------------------------------------------*/
void ExpMonitor_Del_All_Variables(void)
{
	ExpMonitor_Client_Del_All_Vars();
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		ExpMonitor_Del_Variable
	Arguments		-
	Return value	-
					
	Description		Remove variable from the Experiment monitor
					ExpMonitor_Del_Variable(ExpMonitorVar)
------------------------------------------------------------------------------------------------------------------------*/
void ExpMonitor_Del_Variable(int VarHandle)
{
	ExpMonitor_Client_Del_Var(VarHandle);
}



/*------------------------------------------------------------------------------------------------------------------------
 *
 *			ScanController
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Add_Variables
	Arguments		-
	Return value	-
					
	Description		Add variable(s) to the Scancontroller
					When the value of the variable is changed in the Scancontroller,
					the value in this program is also updated
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Add_Variables(void)
{
	ScanController_Struct_VarInfo varinfo={0};
	
	// delete all client variables on server
	ScanController_Client_ClearAllVars();
	
	// add vars
	varinfo.initvalue   = 2;
	varinfo.initmin     = 0;
	varinfo.initmax     = 5;
	strcpy( varinfo.name, "Example/1#23");
	strcpy( varinfo.unit, "[V]");
	ScanController_Client_AddVar(varinfo, &ScanController_ExampleVar, sizeof(ScanController_ExampleVar), NULL);
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Start_Callback
	Arguments		VarList
					startdata
	Return value	-
					
	Description		Called when the ScanController calls the start callback function
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Start_Callback(void) 
{
	BOOL update_succesfull;
	
	update_succesfull = Update_value();  
	if (!update_succesfull)
	{
		message_window_write("Error updating.");
		ScanController_Client_ClientError(); //signal server there is a problem
		return;
	}
	
	ScanController_Client_ClientReady(); // Signal server it is done 
}

//When acting as monitor client
/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Start_Callback
	Arguments		VarList
					startdata
	Return value	-
					
	Description		Called when the ScanController calls the start callback function
------------------------------------------------------------------------------------------------------------------------*/
/*
void ScanController_Start_Callback(void) 
{
	if (error_occurred_for_monitoring_variable)
	{
		ScanController_Client_ClientError(); //signal server there is a problem
	}
	else
	{
		ScanController_Client_ClientReady(); // Signal server it is done
	}
}
*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Start_Acq_Callback
	Arguments		VarList
					startdata
	Return value	-
					
	Description		Called when the ScanController calls the start acquisition callback function
------------------------------------------------------------------------------------------------------------------------*/
/*
void ScanController_Start_Acq_Callback(ScanController_Struct_AcqInfo startdata)
{
	char message[1000];
	
	message_window_write("Aquisition callback");
	
	
	sprintf(message,"Save data to %s\\%s_example.txt", startdata.path, startdata.filenameprefix); // Set save filename prefix
	message_window_write(message);
	/*
	save_succesfull = save_data(startdata.path, startdata.filenameprefix);  
	if (!save_succesfull)
	{
		message_window_write("Error saving.");
		ScanController_Client_ClientError(); //signal server there is a problem
		return;
	}
	* /
	ScanController_Client_ClientReady(); 
}
*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Disconnect
	Arguments		-
	Return value	-
					
	Description		Disconnect from ScanController
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Disconnect(void)
{   
	Connected_To_ScanController = 0; 
	
	// update GUI
	SetCtrlAttribute(panelHandle, PANEL_CLIENTMODECOMMAND, ATTR_LABEL_TEXT, "Start Clientmode");
	SetCtrlAttribute(panelHandle, PANEL_SCANCONTR_NUMERIC, ATTR_DIMMED, 0);
	
	// disconnect from scanning server
	ScanController_Client_Disconnect();
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Disconnect_Callback
	Arguments		-
	Return value	-
					
	Description		ScanController gives Disconnected command (i.e. ScanController program is closed)
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Disconnect_Callback(void) 
{
	ScanController_Disconnect();
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Connect
	Arguments		-
	Return value	-
					
	Description		Connect to Scancontroller
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Connect(void)
{
	int connection_succesfull;
	
	
	//Hardware Client
	connection_succesfull = ScanController_Client_Connect(ScanController_clienttype_hw,  "Example_HW_client", ScanController_Start_Callback, NULL, ScanController_Disconnect_Callback);
	
	//Acquisition Client (this client can also add variables)
	//connection_succesfull = ScanController_Client_Connect(ScanController_clienttype_acq,  "Example_Acq_client", ScanController_Variable_Update, ScanController_Start_Callback, ScanController_Start_Acq_Callback, ScanController_Disconnect_Callback);  
	
	//Monitor Client
	//connection_succesfull = ScanController_Client_Connect(ScanController_clienttype_mon,  "Example_Mon_client", NULL, ScanController_Start_Callback, NULL, ScanController_Disconnect_Callback);
	
	if (connection_succesfull)
	{
		Connected_To_ScanController = 1; 	
		// update GUI
		SetCtrlAttribute(panelHandle, PANEL_CLIENTMODECOMMAND, ATTR_LABEL_TEXT, "Stop Clientmode");
		SetCtrlAttribute(panelHandle, PANEL_SCANCONTR_NUMERIC, ATTR_DIMMED, 1); //make it impossible to manually enter value when in clientmode
		
		// register variables on server  
		ScanController_Add_Variables();	
	}
}





/*------------------------------------------------------------------------------------------------------------------------
	Function		matlab_script_callback
					
	Description		User selected a MATLAB script for image analysis
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK matlab_script_callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			SelectMatlabScript();
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		AnalyseInMatlab
					
	Description		User has checked / unchecked 'Analyse In Matlab'
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK analyze_in_matlab_callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	int checked;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelHandle, PANEL_ANALYSE_IN_MATLAB, &checked);
			if (StartStopMatlabServer(checked, panelHandle, PANEL_MATLAB_SCRIPT, PANEL_MATLAB_FIT_VALUE) == 0)
			{ 
				//start successfull
				Analyze_In_Matlab = checked;
				//if (Analyse_In_Matlab)
					//SetCtrlAttribute(panelHandle, PANEL_MATLAB_ANALYZE_LED, ATTR_VISIBLE, 1);  
				//else
					//SetCtrlAttribute(panelHandle, PANEL_MATLAB_ANALYZE_LED, ATTR_VISIBLE, 0);  
			}
			else
				SetCtrlVal(panelHandle, PANEL_ANALYSE_IN_MATLAB, !checked);
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		run_matlab_script_callback
					
	Description		
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK run_matlab_script_callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	double voltage;
	char argument[100];
	
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelHandle, PANEL_EXAMPLE_NUMERIC, &voltage);
			sprintf(argument, "%f", voltage);
			MatlabRunScript(argument);
			break;
	}
	return 0;
}


/*------------------------------------------------------------------------------------------------------------------------
 *
 *			UI Callbacks
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		clientmodecallback
					
	Description		Button 'Start / Stop Clientmode' on panel pressed
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK clientmodecallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			if (!Connected_To_ScanController) ScanController_Connect();
			else							  ScanController_Disconnect();
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		panelCB
					
	Description		User closes UI
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK panelCB (int panel, int event, void *callbackData,
        int eventData1, int eventData2)
{
    if (event == EVENT_CLOSE)
        QuitUserInterface (0);
    return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		changecallback
					
	Description		Value of constant changed
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK changecallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelHandle, PANEL_SCANCONTR_NUMERIC, &ScanController_ExampleVar);
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		main
					
	Description		Main function
------------------------------------------------------------------------------------------------------------------------*/
int main (int argc, char *argv[])
{
    int error = 0;
    int program_is_already_running; 
	
	/* initialize and load resources */
    nullChk (InitCVIRTE (0, argv, 0));
    errChk (panelHandle = LoadPanel (0, "Example.uir", PANEL));
	
	test();
	
	//Prevent starting the program if it is already running
	CheckForDuplicateAppInstance(ACTIVATE_OTHER_INSTANCE, &program_is_already_running);
	if (program_is_already_running)
	{
		MessagePopup("Error", "An other instance of this program is already running");
		return 0;
	}
	
	// init GUI
	
	ExpMonitor_Add_Variables();
	
	
    /* display the panel and run the user interface */
    errChk (DisplayPanel (panelHandle));
    
	
	errChk (RunUserInterface ());
	
	if (Analyze_In_Matlab)
	{
		SetCtrlVal(panelHandle, PANEL_ANALYSE_IN_MATLAB, 0);
		//SetCtrlAttribute(panelHandle, PANEL_MATLAB_ANALYZE_LED, ATTR_VISIBLE, 0);
		StartStopMatlabServer(0, panelHandle, PANEL_MATLAB_SCRIPT, PANEL_MATLAB_FIT_VALUE);
	}
	
	// Disconnect from scancontroller if needed
	if (Connected_To_ScanController) ScanController_Disconnect(); 
	ExpMonitor_Del_All_Variables();
	
	// save settings panel
	SavePanelState(panelHandle, "panelstate.dat", 0);
	
	
		
Error:
    /* clean up */
    DiscardPanel (panelHandle);
    return 0;
}
