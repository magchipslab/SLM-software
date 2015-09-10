//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Callback functions for the simulation panel.
//
// Created on:  04-01-2012 by Rick.
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
//#include "Camera Controller.h"


// number of pixels of canvas and camera feed
static int gSimX, gSimY, gSimCanvasX, gSimCanvasY;

// the pixel position of the bitmap (camera feed) that should be centered in the canvas
static int gSimCenterX = 0, gSimCenterY = 0;

// zoom factor
static double gSimZoomFactor = 1.0;


/// HIFN Callback for the simulation panel resize event                           
int CVICALLBACK SimPanel_Callback (int panel, int event, void *callbackData,
		int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_PANEL_SIZE:
		{	
			// get the new panel size
			int SimPanelX, SimPanelY;
			GetPanelAttribute(panel, ATTR_WIDTH,  &SimPanelX);
			GetPanelAttribute(panel, ATTR_HEIGHT, &SimPanelY);
			
			// compute the new canvas size
			gSimX = SimPanelX - 6;
			gSimY = SimPanelY - 90;
	
			// adjust its canvas size to match as well
			SetCtrlAttribute(panel, SimPanel_CANVAS, ATTR_WIDTH,  gSimX);
			SetCtrlAttribute(panel, SimPanel_CANVAS, ATTR_HEIGHT, gSimY);
			
			// reinitialise the simulation for the new dimensions
			SLM_initialiseSimulation(gSimX, gSimY);
			
			// move the zoom, amplitude, and saturation slider
			SetCtrlAttribute(panel, SimPanel_SimZoom, 		ATTR_TOP, SimPanelY - 45);
			SetCtrlAttribute(panel, SimPanel_SimSaturation, ATTR_TOP, SimPanelY - 45);
			SetCtrlAttribute(panel, SimPanel_AmplMod, ATTR_TOP, SimPanelY - 45);
			
			// move the toggle buttons
			SetCtrlAttribute(panel, SimPanel_SimToggle,      ATTR_TOP, SimPanelY - 90);
			SetCtrlAttribute(panel, SimPanel_SimPhaseToggle, ATTR_TOP, SimPanelY - 70);
			SetCtrlAttribute(panel, SimPanel_SimSuperSample, ATTR_TOP, SimPanelY - 40);
			SetCtrlAttribute(panel, SimPanel_SimHelmholtz,   ATTR_TOP, SimPanelY - 30);
			
			// move the grid controls
			SetCtrlAttribute(panel, SimPanel_SimGridSpacing, ATTR_TOP, SimPanelY - 30);
			SetCtrlAttribute(panel, SimPanel_SimGrid,        ATTR_TOP, SimPanelY - 70);

			// move the save button
			SetCtrlAttribute(panel, SimPanel_SaveSim, ATTR_TOP, SimPanelY - 63);
			
			// update the simulation canvas (through the SLM canvas update routine)
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

			break;
		}
	}
	return 0;
}


/// HIFN Callback function for mouse events on the camera canvas control
int CVICALLBACK SimCanvas_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	// get the coordinates of the mouse relative to the canvas origin
	int MouseX;
	int MouseY;
	
	switch (event)
	{
		case EVENT_LEFT_CLICK:

			// get the mouse pixel position on the canvas
			GetRelativeMouseState(panel, control, &MouseX, &MouseY, 0, 0, 0);
			
			break;
			
		case EVENT_LEFT_DOUBLE_CLICK:
			
			// get the mouse pixel position on the canvas
			GetRelativeMouseState(panel, control, &MouseX, &MouseY, 0, 0, 0);

			// reset the zoom
			gSimZoomFactor = 1.0;
			
			break;			
	}
	
	return 0;
}


/// HIFN Callback function for the zoom sliderbar
int CVICALLBACK SimZoom_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:
		{	
			// get the sliderposition
			double sliderpos;
			GetCtrlVal (panel, control, &sliderpos);
			
			// set the zoom factor
			SLM_setSimZoom(sliderpos);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}

			break;
		}
	}
	return 0;
}


/// HIFN when the user changes the value of the simulation toggle button
int CVICALLBACK SimToggle_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// check if the simulation button is toggled
			unsigned int simtoggle;
			GetCtrlVal(panel, SimPanel_SimToggle, &simtoggle);
			
			// check if the PHASE simulation button is toggled
			unsigned int simphasetoggle;
			GetCtrlVal(panel, SimPanel_SimPhaseToggle, &simphasetoggle);
			
			// check if the SLM pixels simulation button is toggled
			unsigned int SLM_pixeltoggle;
			GetCtrlVal(panel, SimPanel_SLMpixeltoggle, &SLM_pixeltoggle); 
			
			// toggle the SLM simulation
			SLM_toggleSim(simtoggle, simphasetoggle, SLM_pixeltoggle);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			
			break;
		}
	}
	return 0;
}
				


/// HIFN Callback function for changing the simulation saturation threshold
int CVICALLBACK SimSaturation_Callback (int panel, int control, int event, void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			double sliderpos;
			GetCtrlVal(panel, control, &sliderpos);
			SLM_setSimSaturation(sliderpos);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
			}
			
			break;
		}
	}
	return 0;
}


/// HIFN Callback for the grid toggle button	  
int CVICALLBACK SimGrid_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// read out the setting of the grid toggle button
			int ShowGrid;
			GetCtrlVal(panel, control, &ShowGrid);
			
			// toggle the grid in the simulation
			SLM_setSimShowGrid(ShowGrid);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
			}
			
			break;
		}
	}
	return 0;
}


/// HIFN Callback for the grid-spacing numeric control
int CVICALLBACK SimGridSpacing_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// read out the grid spacing
			double GridSpacing;
			GetCtrlVal(panel, control, &GridSpacing);
			
			// toggle the grid in the simulation
			SLM_setSimGridSpacing(1e-6 * GridSpacing);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
			}

			break;
		}
	}
	return 0;
}

int CVICALLBACK SaveSim_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// generate a bitmap from the canvas			
			int TmpBitmap = -1;
			GetCtrlDisplayBitmap (panel, SimPanel_CANVAS, 0, &TmpBitmap);

			// let the user specify a filename for the bitmap
			char pngfilename[1024];
		    int SelectionStatus = FileSelectPopup("", "*.png", "*.png", "Save Screenshot", 
									VAL_OK_BUTTON, 0, 1, 1, 1, pngfilename);
			
			if (SelectionStatus > 0)
			{
				// save it to file
				SaveBitmapToPNGFile (TmpBitmap, pngfilename);
			}
			
			DiscardBitmap(TmpBitmap);

			break;
		}
	}
	return 0;
}

int CVICALLBACK AmplMod_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// read out the slider position
			double sliderpos;
			GetCtrlVal(panel, control, &sliderpos);
			
			// set the ampl.modulation of the SLM
			SLM_setAmplitudeModulation(sliderpos);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}

			break;
		}
	}
	return 0;
}

int CVICALLBACK SimSuperSample_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// read out the slider position
			int ssample;
			GetCtrlVal(panel, control, &ssample);
			
			// toggle the supersampling in the simulation
			SLM_setSimSuperSampling(ssample);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}

			break;
		}
	}
	return 0;
}

int CVICALLBACK SimHelmholtz_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// read out the slider position
			int helmholtz;
			GetCtrlVal(panel, control, &helmholtz);
			
			// toggle the supersampling in the simulation
			SLM_setSimHelmholtz(helmholtz);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			
			break;
	}
	return 0;
}

int CVICALLBACK SimZPos_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// read out the slider position
			double sliderpos;
			GetCtrlVal(panel, control, &sliderpos);
			
			// set the ampl.modulation of the SLM
			SLM_setSimZPos(sliderpos / 1000.0);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}

			break;
	}
	return 0;
}
