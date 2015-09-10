//==============================================================================
//
// Title:       SLM Controller Client Mode
// Purpose:     Routines that allow the SLM Controller to operate in Client Mode,
//              in which it can be controlled by the Scan Controller application
//
// Created on:  9-6-2012 at 15:10:12 by Rick.
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

//==============================================================================
// Include files

#include <formatio.h>
#include <ansi_c.h>
#include <cvirte.h>     
#include <userint.h>
#include "SLM Control Panel.h"
#include "SLM Control Panel Internal.h"
#include "toolbox.h"
#include "SLM.h"
#include "SLM Camera Panel.h"
#include "UCP/General_DLL.h"

// indicator variable whether the program is in Client Mode or not
static int gClientMode = 0;

// the number of controls linked to the ScanController
static int gNumClientControls = 0;

// the list of controls
static int* gClientControl_Panel;
static int* gClientControl_Control;
static int* gClientControl_DataType;

// the list of values of the controls, to be set by the ScanController
void** gClientControl_Values;

// indicator variables used by the ScanController to indicate whether a variable was updated
char* gClientControl_Updated;

// indicator variables that indicate whether a change to this control would require a recalculation
// of the SLM base phase pattern
char* gClientControl_RequiresPatternCalculation;



char *replace(char *st, char *orig, char *repl) {
  static char buffer[4096];
  char *ch;
  if (!(ch = strstr(st, orig)))
   return st;
  strncpy(buffer, st, ch-st);  
  buffer[ch-st] = 0;
  sprintf(buffer+(ch-st), "%s%s", repl, ch+strlen(orig));
  return buffer;
  }



/// HIFN add a control to the list of controls to be operated in Client Mode
void addClientModeControl(int panel, int control, char pattern_recalc, char units[])
{
	// first, we add the control to the list of Client Mode Controls
	
	// make room for an extra control
	gNumClientControls++;
	gClientControl_Panel    = (int*) realloc(gClientControl_Panel,    gNumClientControls * sizeof(int));
	gClientControl_Control  = (int*) realloc(gClientControl_Control,  gNumClientControls * sizeof(int));
	gClientControl_DataType = (int*) realloc(gClientControl_DataType, gNumClientControls * sizeof(int));
	gClientControl_Values = (void**) realloc(gClientControl_Values,   gNumClientControls * sizeof(void*));
	gClientControl_Updated = (char*) realloc(gClientControl_Updated,  gNumClientControls * sizeof(char));
	gClientControl_RequiresPatternCalculation = (char*) realloc(gClientControl_RequiresPatternCalculation, gNumClientControls * sizeof(char));
	
	// obtain the control data type
	int datatype;
	GetCtrlAttribute(panel, control, ATTR_DATA_TYPE, &datatype);
	
	// fill out the data fields
	int index = gNumClientControls - 1;
	gClientControl_Panel   [index] = panel;
	gClientControl_Control [index] = control;
	gClientControl_DataType[index] = datatype;
	gClientControl_Updated [index] = 0;
	gClientControl_RequiresPatternCalculation[index] = pattern_recalc;

	// allocate memory for holding the variable value
	int datasize;
	switch (datatype)
	{
		case VAL_INTEGER:
			datasize = sizeof(int);
			break;
			
		case VAL_DOUBLE:
		default:
			datasize = sizeof(double);
			break;
	}
	gClientControl_Values[index] = (void*) malloc(1 * datasize);
	
	// structure for storing the variable information
	ScanController_Struct_VarInfo varinfo={0};	
	
	// get the control value and range as doubles
	switch (datatype)
	{
		case VAL_INTEGER:
		{			
			// get the value, min and max as integers
			int val, min, max;
			GetCtrlVal(panel, control, &val);
			GetCtrlAttribute(panel, control, ATTR_MIN_VALUE, &min);
			GetCtrlAttribute(panel, control, ATTR_MAX_VALUE, &max);
		
			// store them in the varinfor structure as doubles
			varinfo.initvalue = (double) val;
			varinfo.initmin   = (double) min;
			varinfo.initmax   = (double) max;
			
			break;
		}
		case VAL_DOUBLE:
		{
			// read out the value, min, and max directly into the varinfo struct
			GetCtrlVal(panel, control, &(varinfo.initvalue));
			GetCtrlAttribute(panel, control, ATTR_MIN_VALUE, &(varinfo.initmin));
			GetCtrlAttribute(panel, control, ATTR_MAX_VALUE, &(varinfo.initmax));
			
			break;
		}
	}
	
	// fill the remaining fields of the varinfo struct
	GetCtrlAttribute(panel, control, ATTR_LABEL_TEXT, (varinfo.name));
	
	for (int k = 0; k < strlen(varinfo.name); k++)
		if (varinfo.name[k] == ' ')
			varinfo.name[k] = '_';
	
	varinfo.vartype = datatype;
	strcpy(varinfo.unit, units);
	
	// add the variable to the ScanController
	ScanController_Client_AddVar(varinfo, gClientControl_Values[index], datasize, &gClientControl_Updated[index]);
	
	// dim the control such that it is no longer editable through the GUI
	SetCtrlAttribute(panel, control, ATTR_DIMMED, 1);
}


/// HIFN callback function called by the ScanController when it sets new variable values
void ScanController_Callback(void) 
{
	// check if the SLM pixel panel is in fullscreen mode, 
	// otherwise send an error message to ScanController
	int xsize, ysize;
	GetPanelAttribute(pnlSLMpixels, ATTR_WIDTH,  &xsize);
	GetPanelAttribute(pnlSLMpixels, ATTR_HEIGHT, &ysize);
	if ( (xsize != 1920) || (ysize != 1080) )
	{
		// signal the server that the SLM is not working properly, and abort
		ScanController_Client_ClientError();
		return;		
	}
	
	// variable to indicate whether the current pattern needs to be re-calculated
	int recalc_pattern = 0;
	
	// loop over all the Client Controls
	for (int n = 0; n < gNumClientControls; n++)
	{
		// check if the associated variable value was updated
		if (gClientControl_Updated[n])
		{
			// yes, update the control value
			switch (gClientControl_DataType[n])
			{
				case VAL_INTEGER:
					SetCtrlVal(gClientControl_Panel[n], gClientControl_Control[n], *((int*) gClientControl_Values[n]));
					break;
				case VAL_DOUBLE:
					SetCtrlVal(gClientControl_Panel[n], gClientControl_Control[n], *((double*) gClientControl_Values[n]));
					break;
				default:
					break;
			}
			
			// check if the pattern needs to be updated
			recalc_pattern |= gClientControl_RequiresPatternCalculation[n];
		}
	}
	
	// do we need to recalculate the SLM base phase pattern?
	if (recalc_pattern)
	{	
		switch (SLM_getCurrentPattern())
		{
			case SLM_BEAMSHAPE_STD:
		
				// call the standard beamshaping callback
				//BeamShape_Callback (int panel, int control, EVENT_COMMIT, NULL, SLM_NO_UPDATE, 0);
				
				break;
		
			case SLM_BEAMSHAPE_ARB:
		
				// call the beamshaping callback
				//LoadFile_Callback (int panel, int control, EVENT_COMMIT, NULL, SLM_NO_UPDATE, 0);
				
				break;
		
			case SLM_BEAMSPLIT_PROJECT:
		
				// call the spot pattern 'callback'
				setSpotPatternFromParameters();
				
				break;
		
			case SLM_BEAMSPLIT_IFTA:
		
				// call the beamsplitting callback
				//PhaseRetrieval_Callback (int panel, int control, EVENT_COMMIT, NULL, SLM_NO_UPDATE, 0);
				
				break;
		
			case SLM_BEAMSPLIT_PHASEGRATING:
				
				// set the phase grating
				setPhaseGratingFromParameters();
					
				break;
				
			case SLM_SHACKHARTMANN:
				
				ShackHartmann_Callback (TabPage_7, 0, EVENT_COMMIT, NULL, SLM_NO_UPDATE, 0);
				
				break;
				
			case SLM_BEAMSHAPE_CHAR:
				
				CharacterGo_Callback(TabPage_2_1, TABPANEL_8_CharacterGo, EVENT_COMMIT, NULL, SLM_NO_UPDATE, 0);
				
				break;
				
			default:
		
				// not supported
				
				break;	
		}
	}
	
	// update the SLM settings from the control panel values
	//setDefaultValuesFromPanel();
	Lens_Callback			(TabPage_0, TABPANEL_LensXphase,	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	Lens_Callback			(TabPage_0, TABPANEL_LensYphase,	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	HorizTrans_Callback		(TabPage_0, TABPANEL_HorizTrans, 	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	VertTrans_Callback		(TabPage_0, TABPANEL_VertTrans, 	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	
	// update the SLM
	SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, recalc_pattern);
	
	// wait a short while for the image to stabilise (?)
	Delay(0.2);
	
	// tell the ScanController that we are finished
	ScanController_Client_ClientReady(); 
}


/// HIFN disconnect from the ScanController and clear the ClientControl arrays
void disconnectScanController(void)
{
	// disconnect the Scan Controller
	ScanController_Client_Disconnect();
	
	// clear the Client Control arrays, and re-enable the controls
	for (int n = 0; n < gNumClientControls; n++)
	{
		// re-enable the control
		SetCtrlAttribute(gClientControl_Panel[n], gClientControl_Control[n], ATTR_DIMMED, 0);	
		
		// free the variable value memory
		free(gClientControl_Values[n]);
	}
	free(gClientControl_Panel);
	free(gClientControl_Control);
	free(gClientControl_Values);
	free(gClientControl_DataType);
	free(gClientControl_Updated);
	free(gClientControl_RequiresPatternCalculation);
	gClientControl_Panel 	= NULL;
	gClientControl_Control 	= NULL;
	gClientControl_Values 	= NULL;
	gClientControl_DataType = NULL;
	gClientControl_Updated 	= NULL;
	gClientControl_RequiresPatternCalculation = NULL;
	gNumClientControls = 0;
	
	// update the indicator variable for the Client Mode
	gClientMode = 0;
	
	// update the label on the Client Mode button
	SetCtrlAttribute(TabPage_0, TABPANEL_ClientMode, ATTR_LABEL_TEXT, "__Start Client Mode");
}


/// HIFN registers the controls that are to be operated in Client mode
void registerClientModeControls(void)
{
	// controls on the 'General' tab
	addClientModeControl(TabPage_0, TABPANEL_HorizTrans, 0, "[mm]");
	addClientModeControl(TabPage_0, TABPANEL_VertTrans,  0, "[mm]");
	addClientModeControl(TabPage_0, TABPANEL_LensXphase, 0, "[mm]");
	addClientModeControl(TabPage_0, TABPANEL_LensYphase, 0, "[mm]");
	addClientModeControl(TabPage_0, TABPANEL_Xpos,       0, "[mm]");
	addClientModeControl(TabPage_0, TABPANEL_Ypos,       0, "[mm]");
	addClientModeControl(TabPage_0, TABPANEL_SigmaX,     0, "[um]");
	addClientModeControl(TabPage_0, TABPANEL_SigmaY,     0, "[um]");
	addClientModeControl(TabPage_0, TABPANEL_Bias,       0, "[beard-seconds]");
	
	// add controls specific to the current pattern mode
	switch (SLM_getCurrentPattern())
	{
		case SLM_BEAMSHAPE_STD:
			
			// TODO: add indicator variable in the beam shaping .c file to distinguish between squares and Gaussians
		
			break;
		
		case SLM_BEAMSHAPE_ARB:
		
			break;
			
		case SLM_BEAMSHAPE_CHAR:
		
			addClientModeControl(TabPage_2_1, TABPANEL_8_CharacterNumber, 1, "[char nr. (1-80)]");
			
			break;			
		
		case SLM_BEAMSPLIT_PROJECT:
		
			break;
		
		case SLM_BEAMSPLIT_IFTA:
		
			break;
		
		case SLM_BEAMSPLIT_PHASEGRATING:
				
					
			break;
			
		case SLM_SHACKHARTMANN:
			
			addClientModeControl(TabPage_7, TABPANE_13_SHSpotSize, 1, "[um]");
			addClientModeControl(TabPage_7, TABPANE_13_SHSpotX, 1, "[um]");
			addClientModeControl(TabPage_7, TABPANE_13_SHSpotY, 1, "[um]");
			
			break;
	}
}


/// HIFN Toggle the client mode, i.e. transfer control between the GUI and the Scan Controller (remote program)
int CVICALLBACK ClientMode_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// are we in Client Mode?
			if (gClientMode)
			{
				// yes, stop client mode
				disconnectScanController();				
			}
			else
			{
				// no, start client mode
				if (ScanController_Client_Connect(ScanController_clienttype_hw, "SLM", ScanController_Callback , NULL, disconnectScanController, 1))
				{
					// register variables on server  
					registerClientModeControls();	
				
					// update the indicator variable
					gClientMode = 1;
				
					// update the label on the button
					SetCtrlAttribute(panel, control, ATTR_LABEL_TEXT, "__Stop Client Mode");
				}
			}

			break;
	}
	return 0;
}
