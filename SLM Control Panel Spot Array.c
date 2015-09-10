//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Callback functions for generating spot arrays.
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
#include "SLM_internal.h" 
//#include "Camera Controller.h"


//==============================================================================
// UI callback functions for the SLM spot pattern


/// HIFN gets the parameters for the phase grating pattern and sets the SLM
void setPhaseGratingFromParameters()
{
	// read out the parameters from the phase grating tab page
	double a1, a2, b1, b2, angle;	
	GetCtrlVal(TabPage_1_1, TABPANEL_4_Amplitude1, &a1);
	GetCtrlVal(TabPage_1_1, TABPANEL_4_Amplitude2, &a2);
	GetCtrlVal(TabPage_1_1, TABPANEL_4_Period1,    &b1);
	GetCtrlVal(TabPage_1_1, TABPANEL_4_Period2,   &b2);
	GetCtrlVal(TabPage_1_1, TABPANEL_4_Angle,      &angle);
	
	// set the SLM pattern to a phase grating
	SLM_setPhaseGrating(a1, b1, a2, b2, angle);
}

/// HIFN gets the parameters for the spot pattern and sets the SLM
void setSpotPatternFromParameters()
{
	// get spots input parameters
	unsigned int Nx, Ny;
	double xspacing, yspacing, randomampl, spotlens;
	GetCtrlVal(TabPage_1_0, TABPANEL_3_inputNx, &Nx);
	GetCtrlVal(TabPage_1_0, TABPANEL_3_inputNy, &Ny);
	GetCtrlVal(TabPage_1_0, TABPANEL_3_Xspacing, &xspacing);
	GetCtrlVal(TabPage_1_0, TABPANEL_3_Yspacing, &yspacing);
	GetCtrlVal(TabPage_1_0, TABPANEL_3_RandomAmpl, &randomampl);
	GetCtrlVal(TabPage_1_0, TABPANEL_3_Spotlens, &spotlens);
	//GetCtrlVal(Control, Control_LinSpotlens, &linspotlens);	
	unsigned int linspotlens = 0;
	
	// convert the spacings to m, and set the SLM spot pattern
	SLM_setSpotPattern(Nx, Ny, xspacing * 1e-6, yspacing * 1e-6, randomampl, spotlens, linspotlens);
}


/// HIFN when the user changes the value of the spot controls
int CVICALLBACK SpotControls_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:

			// get the spot pattern parameters and set the SLM pattern
			setSpotPatternFromParameters();
			
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


/// HIFN Callback function for changing the spotmask
int CVICALLBACK SpotSize_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:
		{
			double sliderpos;
			unsigned int Nx, Ny;
			GetCtrlVal(panel, TABPANEL_3_inputNx, &Nx);
			GetCtrlVal(panel, TABPANEL_3_inputNy, &Ny);
			
			
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


/// HIFN Callback function for the radio buttons for choosing the pattern generation mode
int CVICALLBACK rbPattern_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// check which radio button was toggled
			switch (control)
			{
				case TABPANEL_2_rbSpotProjector:
				{	
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_2_rbSpotProjector,  ATTR_CTRL_VAL, 1);
					SetCtrlAttribute(panel, TABPANEL_2_rbPhaseRetrieval, ATTR_CTRL_VAL, 0);
					SetCtrlAttribute(panel, TABPANEL_2_rbPhaseGrating,   ATTR_CTRL_VAL, 0);
					
					// bring the corresponding tabpage to the front
					SetCtrlAttribute(panel, TABPANEL_2_TAB, ATTR_CTRL_INDEX, 0);
					
					// get the spot pattern parameters and set the SLM pattern
					setSpotPatternFromParameters();
					
					break;
				}	
				case TABPANEL_2_rbPhaseRetrieval:
				{	
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_2_rbSpotProjector,  ATTR_CTRL_VAL, 0);
					SetCtrlAttribute(panel, TABPANEL_2_rbPhaseRetrieval, ATTR_CTRL_VAL, 1);
					SetCtrlAttribute(panel, TABPANEL_2_rbPhaseGrating,   ATTR_CTRL_VAL, 0);
					
					// bring the corresponding tabpage to the front
					SetCtrlAttribute(panel, TABPANEL_2_TAB, ATTR_CTRL_INDEX, 2);
					
					break;
				}	
				case TABPANEL_2_rbPhaseGrating: 
				{	
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_2_rbSpotProjector,  ATTR_CTRL_VAL, 0);
					SetCtrlAttribute(panel, TABPANEL_2_rbPhaseRetrieval, ATTR_CTRL_VAL, 0);
					SetCtrlAttribute(panel, TABPANEL_2_rbPhaseGrating,   ATTR_CTRL_VAL, 1);
					
					// bring the corresponding tabpage to the front
					SetCtrlAttribute(panel, TABPANEL_2_TAB, ATTR_CTRL_INDEX, 1);
					
					// get the phase grating parameters and set the SLM pattern
					setPhaseGratingFromParameters();
													
					break;
				}
			}
			
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


/// HIFN Callback for all controls for the phase grating pattern
int CVICALLBACK PhaseGrating_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:
			
			// get the phase grating parameters and set the SLM pattern
			setPhaseGratingFromParameters();

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


/// HIFN Callback function for toggling the phase grating mode to a pattern Gaussian peaks
int CVICALLBACK GaussianPeaks_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// check if the button is on or off
			int gptoggle;
			GetCtrlVal(panel, control, &gptoggle);
			
			// toggle the Gaussian peaks setting
			SLM_togglePhaseGratingGaussianPeaks(gptoggle);
			
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



/// HIFN callback function when the phase retrieval settings have been changed
int CVICALLBACK PhaseRetrieval_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:

			// one of the control parameters was changed, update the preview
			
			break;
			
		case EVENT_COMMIT:
		{	
			// read out all the control parameters
			int spotxoffset, spotyoffset, sigxoffset, sigyoffset;
			int spotxspacing, spotyspacing, numxspots, numyspots;
			int phaseconstraint;
			double phasestep;
			GetCtrlVal(panel, TABPANEL_5_SpotYSpacing, &spotyspacing);
			GetCtrlVal(panel, TABPANEL_5_SpotXSpacing, &spotxspacing);
			GetCtrlVal(panel, TABPANEL_5_NumYSpots, &numyspots);
			GetCtrlVal(panel, TABPANEL_5_NumXSpots, &numxspots);
			GetCtrlVal(panel, TABPANEL_5_SpotXOffset, &spotxoffset);
			GetCtrlVal(panel, TABPANEL_5_SpotYOffset, &spotyoffset);
			GetCtrlVal(panel, TABPANEL_5_SigYOffset, &sigyoffset);
			GetCtrlVal(panel, TABPANEL_5_SigXOffset, &sigxoffset);
			GetCtrlVal(panel, TABPANEL_5_PhaseConstraint, &phaseconstraint);
			GetCtrlVal(panel, TABPANEL_5_PhaseStep, &phasestep);
			
			// check if it was the GO-button that was pressed
			if (control == TABPANEL_5_PhaseRetrievalGo)
			{
				// yes, we have to start phase retrieval with the current parameters
				// check if the asked spotpattern fits the SLM
				if ((sigxoffset + spotxoffset + numxspots * (spotxspacing+1)) >= gXsize)
				{
					char msg[1024];
					sprintf(msg, "The spotpattern exceeds the SLM size in the x-dimension. Adjust X spots or X spacing.");
					InsertTextBoxLine (panel, TABPANEL_5_Output, -1, msg);
				}
				else if ((sigyoffset + spotyoffset + numyspots * (spotyspacing +1)) >= gYsize)
				{
					char msg[1024];
					sprintf(msg, "The spotpattern exceeds the SLM size in the y-dimension. Adjust Y spots or Y spacing.");
					InsertTextBoxLine (panel, TABPANEL_5_Output, -1, msg);
				}
				else
				{
					SLM_generateSpotArray(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, panel, TABPANEL_5_Output,
	                       numxspots, numyspots, spotxspacing, spotyspacing, sigxoffset, sigyoffset, spotxoffset, spotyoffset, 
						   phaseconstraint, phasestep * PI);
				
					// check if we need to update
					if (eventData1 != SLM_NO_UPDATE)
					{
						// update the SLM panel and the simulation Panel (if toggled)
						SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
					}
				}
			}
			
			
			
			break;
		}
	}
	return 0;

}
