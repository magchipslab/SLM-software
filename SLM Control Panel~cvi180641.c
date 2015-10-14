//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Application for generating and displaying phase patterns for
//              digital holography with a Spatial Light Modulator (SLM).
//
// Created on:  9-9-2011 at 13:40:12 by Rick.
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

int pnlControl, pnlSLMpixels, pnlSimPanel, pnlCamera, pnlDebug;

int TabPage_0, TabPage_1, TabPage_2, TabPage_3, TabPage_4, TabPage_5, TabPage_6, TabPage_7, TabPage_1_0, TabPage_1_1, TabPage_1_2;
int TabPage_2, TabPage_2_0, TabPage_2_1;

//==============================================================================
// Constants

//==============================================================================
// Types

//==============================================================================
// Static global variables

static int panelHandle;

//==============================================================================
// Static functions

//==============================================================================
// Global variables

//==============================================================================
// Global functions

/// HIFN The main entry-point function.
int main (int argc, char *argv[])
{
	// load the panels from the user interface file
	if (InitCVIRTE (0, argv, 0) == 0)
		return -1;	/* out of memory */
	if ((pnlControl = LoadPanel (0, "SLM Control Panel.uir", Control)) < 0)
		return -1;
	if ((pnlSLMpixels = LoadPanel (0, "SLM Control Panel.uir", SLMpixels)) < 0)
		return -1;
	if ((pnlSimPanel = LoadPanel (0, "SLM Control Panel.uir", SimPanel)) < 0)
		return -1;
	if ((pnlCamera = LoadPanel (0, "SLM Control Panel.uir", PANEL)) < 0)
		return -1;
	if ((pnlDebug = LoadPanel (0, "SLM Control Panel.uir", DebugPanel)) < 0)
		return -1;
	
	// get the tab panel handles
	
	// the main tab control, page 0 (General Settings) page 1 (Spot Array) 
	// page 2 (Beam Shaping) page 3 (Feedback) page 4 (Phase Error),
	// page 5 (load and save), page 6 (Genetic Algorithms), page 7 (Shack-Hartmann)
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 0, &TabPage_0);
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 1, &TabPage_1);
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 2, &TabPage_2);	
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 3, &TabPage_3);		
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 4, &TabPage_4);
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 5, &TabPage_5); 
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 6, &TabPage_6);
	GetPanelHandleFromTabPage (pnlControl, Control_TAB, 7, &TabPage_7);  
	
	// the nested tab controls on page 1 (Spot Array) of the main tab, pages 0, 1, and 2
	// (Spot Projector, Phase Grating, and Beam Splitting respectively)
	GetPanelHandleFromTabPage (TabPage_1, TABPANEL_2_TAB, 0, &TabPage_1_0);	
	GetPanelHandleFromTabPage (TabPage_1, TABPANEL_2_TAB, 1, &TabPage_1_1);
	GetPanelHandleFromTabPage (TabPage_1, TABPANEL_2_TAB, 2, &TabPage_1_2);
	
	// the nested tab controls on page 2 (Beam Shaping), pages, 0 and 1
	// (Beam Shaping standard shapes, and Beam Shaping Arbitrary, respectively)
	GetPanelHandleFromTabPage (TabPage_2, TABPANEL_6_TAB, 0, &TabPage_2_0);
	GetPanelHandleFromTabPage (TabPage_2, TABPANEL_6_TAB, 1, &TabPage_2_1);
	
	// get the dimensions of the SLM pixel panel
	int xsize, ysize;
	GetPanelAttribute(pnlSLMpixels, ATTR_WIDTH,  &xsize);
	GetPanelAttribute(pnlSLMpixels, ATTR_HEIGHT, &ysize);
	
	// initialise the camera
	InitCamera(pnlCamera, PANEL_CameraCanvas, PANEL_HistoX, PANEL_HistoY, PANEL_CameraTimer); 
	
	// uncomment the following line to display the debug panel: a panel with several canvases and graphs which can be used for debug plots
	DisplayPanel (pnlDebug);
	
	// display the panels
	DisplayPanel (pnlSimPanel);
	DisplayPanel (pnlSLMpixels);	
	DisplayPanel (pnlCamera);                                        
	DisplayPanel (pnlControl);
	
	// get the size of the simulation canvas
	int SimX, SimY;
	GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_WIDTH,  &SimX);
	GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_HEIGHT, &SimY);

	// get the subsampling factor
	int subsample;
	GetCtrlVal(TabPage_0, TABPANEL_SubSample, &subsample);
	
	// initialise the SLM and simulation
	SLM_initialise(xsize, ysize, SimX, SimY, subsample);
	
	// set the default values of the SLM library global variables, to those of the controls on the panel
	setDefaultValuesFromPanel();
	
	// set the debug plotting controls
	SLM_setDebugPlottingControls(pnlDebug, DebugPanel_CANVAS, DebugPanel_CANVAS_2, DebugPanel_CANVAS_3,
		DebugPanel_CANVAS_4, DebugPanel_CANVAS_5, DebugPanel_CANVAS_6, DebugPanel_DBGGRAPH, DebugPanel_DBGGRAPH2, DebugPanel_DBGGRAPH3);
	
	// run the user interface
	RunUserInterface ();
	
	// clean up
	QuitCamera();
	DiscardPanel (pnlControl);
	DiscardPanel (pnlSLMpixels);
	DiscardPanel (pnlSimPanel);
	DiscardPanel (pnlCamera);
	
	return 0;
}

//==============================================================================
// UI callback function prototypes


/// HIFN Exit when the user dismisses the panel.
int CVICALLBACK Control_Callback (int panel, int event, void *callbackData,
		int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_CLOSE:
        	QuitUserInterface (0);
			break;
	}
	return 0;
}


/// HIFN Sets the variables of the SLM library to those of the control panel
void setDefaultValuesFromPanel()
{
	// call all the input control callbacks, but with the SLM_NO_UPDATE specifier
	
	// the lens and wavelength callbacks need to be executed first, they set the focal unit which is needed
	// in the calculations of the other controls
	Magnification_Callback  (TabPage_0, TABPANEL_Magnification,   EVENT_COMMIT,      NULL, SLM_NO_UPDATE, 0);
	LensFocalLength_Callback(TabPage_0, TABPANEL_LensFocalLength, EVENT_COMMIT, 	 NULL, SLM_NO_UPDATE, 0);
	Wavelength_Callback		(TabPage_0, TABPANEL_Wavelength, 	  EVENT_COMMIT, 	 NULL, SLM_NO_UPDATE, 0);
	
	
	// now we can continue with the rest of the callbacks
	InputIntensity_Callback	(TabPage_0, NULL, 					  EVENT_COMMIT, 	 NULL, SLM_NO_UPDATE, 0);
	Lens_Callback			(TabPage_0, TABPANEL_LensXphase,	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	Lens_Callback			(TabPage_0, TABPANEL_LensYphase,	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	SimSaturation_Callback	(pnlSimPanel, SimPanel_SimSaturation, EVENT_COMMIT, 	 NULL, SLM_NO_UPDATE, 0);
	SimToggle_Callback		(pnlSimPanel, SimPanel_SimToggle, 	  EVENT_COMMIT,		 NULL, SLM_NO_UPDATE, 0);
	HorizTrans_Callback		(TabPage_0, TABPANEL_HorizTrans, 	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	VertTrans_Callback		(TabPage_0, TABPANEL_VertTrans, 	  EVENT_VAL_CHANGED, NULL, SLM_NO_UPDATE, 0);
	Bias_Callback			(TabPage_0, TABPANEL_Bias, 			  EVENT_COMMIT,		 NULL, SLM_NO_UPDATE, 0);
	

	
	// update the information about the SLM
	UpdateSLMinfo();
	
	// update the SLM panel and the simulation Panel (if toggled)
	SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
	
}

void UpdateSLMinfo(void)
{

	
	// update the resolution and focal unit text box on the general panel
	char txtSLMres[1024];																								    
	char txtFocalUnit[1024];
	char txtEffectiveFocalLength[1024];
	sprintf(txtSLMres, "(%i x%i)", SLM_getXres(), SLM_getYres());
	sprintf(txtFocalUnit, "(%.2f x%.2f)", SLM_getFocalUnitX() / 1e-6, SLM_getFocalUnitY() / 1e-6);
	sprintf(txtEffectiveFocalLength, "(%.2f mm)", SLM_getLensFocalLength()*1.0e3);
	SetCtrlVal(TabPage_0, TABPANEL_SLMResolution, txtSLMres);
	SetCtrlVal(TabPage_0, TABPANEL_FocalUnit, txtFocalUnit);
	SetCtrlVal(TabPage_0, TABPANEL_EffectiveFocalLength, txtEffectiveFocalLength);
}



int CVICALLBACK ControlTab_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_ACTIVE_TAB_CHANGE:
        	
			// check if the new tab is the feedback tab
			if (eventData2 == 3)
			{
				// yes, update the signal canvas
				UpdateFeedbackCanvas(TabPage_3, TABPANE_10_FeedbackCanvas);	
			}
			
			/*// check if the new tab is the Genetic Aberration Correction tab
			if (eventData2 == 6)
			{
				// show the actual window dimensions
				char winlabel[256];
				int wx, wy, wxsize, wysize;
				SLM_getSignalWindowParameters(&wx, &wy, &wxsize, &wysize);
				sprintf(winlabel, "Window Size: %u x %u", wxsize, wysize);
				SetCtrlVal(TabPage_6, TABPANE_12_WindowSizeMsg, winlabel);		
			}*/
			
			// check if the new tab is the Shack-Hartmann aberration correction tab
			if (eventData2 == 7)
			{
				ShackHartmann_Callback(TabPage_7, NULL, EVENT_COMMIT, NULL, 0, 0);
			}
			
			
			break;
	}
	return 0;
}
