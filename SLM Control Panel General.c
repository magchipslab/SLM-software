//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Callback functions associated with the general settings tab
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


// the current pattern on the SLM
int gCurrentPattern;


//==============================================================================
// UI callback functions for the general SLM controls


/// HIFN when the user changes the value of the HORIZONTAL slider bar
int CVICALLBACK HorizTrans_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:
		{
			// get slider value
			double sliderpos;
			GetCtrlVal (panel, control, &sliderpos);
			
			// convert from micron to meter
			sliderpos *= 1e-6;
			
			// update the horizontal translation pattern
			SLM_setHorizTrans(sliderpos);
			
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


/// HIFN when the user changes the value of the VERTICAL slider bar
int CVICALLBACK VertTrans_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:
		{
			// get slider value
			double sliderpos;
			GetCtrlVal (panel, control, &sliderpos);
			
			// convert from micron to meter
			sliderpos *= 1e-6;
			
			// update the vertical translation pattern
			SLM_setVertTrans(sliderpos);
			
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


/// HIFN when the user changes the value of the LENS slider bar
int CVICALLBACK Lens_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_VAL_CHANGED:
		{
			// get the lens settings
			int linklensphases;
			double lensxphase, lensyphase;
			GetCtrlVal (panel, TABPANEL_LensXphase, &lensxphase);
			GetCtrlVal (panel, TABPANEL_LensYphase, &lensyphase);
			GetCtrlVal (panel, TABPANEL_LinkLensPhases, &linklensphases);
			
			// convert lens coordinates to m
			lensxphase *= 1e-3;
			lensyphase *= 1e-3;
			
			// check if the lens phases should be linked (i.e. the same)
			if (linklensphases)
				lensyphase = lensxphase;
			
			// update the lens phase pattern 
			// (note that this also affects the aberration correction)
			SLM_setLensPhase(lensxphase, lensyphase);
			SLM_setLensCoordinates(0, 0);
			
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


/// HIFN Callback function for changing the focal length of the lens in the simulation
int CVICALLBACK LensFocalLength_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			double val;
			GetCtrlVal(panel, control, &val);
			
			// set the lens focal length, convert it from mm to m first
			SLM_setLensFocalLength(val * 1.0e-3);
			
			// update the displayed information about the SLM 
			UpdateSLMinfo();
			/*
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// do we also need to recalculate the pattern?
				if (gCurrentPattern == SLM_BEAMSHAPE_STD)
					BeamShape_Callback (TabPage_2_0, TABPANEL_7_BeamShapeOK, EVENT_COMMIT, NULL, 0, 0);
				
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			*/
			break;
		}
	}
	return 0;
}

/// HIFN Callback function for changing the telescope correction
int CVICALLBACK Magnification_Callback(int panel, int control, int event, 
		void *callbackData, int eventData1, int eventData2)
{
	switch(event)
	{
		case EVENT_COMMIT:
		{
			double val;
			GetCtrlVal(panel, control, &val);
			
			// set the correction, update the focal units
			SLM_setMagnification(val);
			
			// update the displayed information about the SLM
			UpdateSLMinfo();
			
			/*
			//check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				//do we need to recalculate the pattern?
				if (gCurrentPattern == SLM_BEAMSHAPE_STD)
					BeamShape_Callback (TabPage_2_0, TABPANEL_7_BeamShapeOK, EVENT_COMMIT, NULL, 0, 0);
				
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}  */
			break;
		}
	}
	return 0;
}
	



/// HIFN Callback function for changing the wavelength of the simulation
int CVICALLBACK Wavelength_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			double val;
			GetCtrlVal(panel, control, &val);
			
			// set the wavelength (convert it from nm to m first)
			SLM_setWavelength(val * 1.0e-9);
			
			// update the displayed information about the SLM 
			UpdateSLMinfo();
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// do we also need to recalculate the pattern?
				if (gCurrentPattern == SLM_BEAMSHAPE_STD)
					BeamShape_Callback (TabPage_2_0, TABPANEL_7_BeamShapeOK, EVENT_COMMIT, NULL, 0, 0);
				
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			
			break;
		}
	}
	return 0;
}


/// HIFN Callback function for when the input intensity parameters are modified
int CVICALLBACK InputIntensity_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// get all parameter values
			double xpos, ypos, sigmax, sigmay;
			GetCtrlVal(panel, TABPANEL_Xpos, &xpos);
			GetCtrlVal(panel, TABPANEL_Ypos, &ypos);
			GetCtrlVal(panel, TABPANEL_SigmaX, &sigmax);
			GetCtrlVal(panel, TABPANEL_SigmaY, &sigmay);
			
			// convert all input values to m
			xpos = xpos * 1.0e-3;
			ypos = ypos * 1.0e-3;
			sigmax = sigmax * 1.0e-6;
			sigmay = sigmay * 1.0e-6;
			
			// get the subsampling factor of the SLM
			int subsamplefactor;
			GetCtrlVal(panel, TABPANEL_SubSample, &subsamplefactor);
			
			// determine the type of input intensity desired
			int gaussian, annular, loadfromfile;
			GetCtrlVal(panel, TABPANEL_rbInputGaussian, &gaussian);
			GetCtrlVal(panel, TABPANEL_rbInputAnnular, &annular);
			GetCtrlVal(panel, TABPANEL_rbInputLoadFromFile, &loadfromfile);
			
			if (loadfromfile)
			{
				// use the intensity pattern loaded from a file (stored in gInputCameraData)
				SLM_setInputIntensityFromCameraData(xpos, ypos, sigmax, sigmay, subsamplefactor);
				
			}
			else if (annular)
			{
				// set an annular intensity profile
				SLM_setInputIntensityAnnular(xpos, ypos, sigmax, sigmay, subsamplefactor);
			}
			else 
			{
				// set a Gaussian intensity profile
				SLM_setInputIntensityGaussian(xpos, ypos, sigmax, sigmay, subsamplefactor);
			}	
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// do we also need to recalculate the pattern?
				if (gCurrentPattern == SLM_BEAMSHAPE_STD)
					BeamShape_Callback (TabPage_2_0, TABPANEL_7_BeamShapeOK, EVENT_COMMIT, NULL, 0, 0);
				
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			
			break;
		}
	}
	return 0;
}


/// HIFN Callback function for setting a uniform bias phase, a phase that is added to all pixels 
///      of the SLM
int CVICALLBACK Bias_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			int val;
			GetCtrlVal(panel, control, &val);
			
			// set the wavelength (convert it from nm to m first)
			SLM_setBias(val);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// do we also need to recalculate the pattern?
				//if (gCurrentPattern == SLM_BEAMSHAPE_STD)
				//	BeamShape_Callback (TabPage_2_0, TABPANEL_7_BeamShapeOK, EVENT_COMMIT, NULL, 0, 0);
				
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
			}

			break;
		}
	}
	return 0;
}


/// HIFN Callback function for changing the subsample factor
int CVICALLBACK SubSample_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// get the dimensions of the SLM pixel panel
			int xsize, ysize;
			GetPanelAttribute(pnlSLMpixels, ATTR_WIDTH,  &xsize);
			GetPanelAttribute(pnlSLMpixels, ATTR_HEIGHT, &ysize);
	
			// get the size of the simulation canvas
			int SimX, SimY;
			GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_WIDTH,  &SimX);
			GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_HEIGHT, &SimY);

			// get the subsampling factor
			int subsample;
			GetCtrlVal(TabPage_0, TABPANEL_SubSample, &subsample);
	
			// initialise the SLM and simulation
			SLM_initialise(xsize, ysize, SimX, SimY, subsample);
			
			// update the SLM panel and the simulation Panel (if toggled)
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
			// update the information about the SLM
			UpdateSLMinfo();
			
			break;
		}
	}
	return 0;
}


/// HIFN loads input intensity from a camera file
void InputLoadFromFile(char camfilename[])
{
	// create a file pointer and open the specified file
	MATFile *pmat = matOpen(camfilename, "r");

	// read out the camera resolution
	int camxres = (int) readMatDoubleScalar(pmat, "cam_x_size");
	int camyres = (int) readMatDoubleScalar(pmat, "cam_y_size");
	
	// read the corner coordinates of the current zoom on the canvas
	int camxmin = (int) readMatDoubleScalar(pmat, "cam_x_min");
	int camxmax = (int) readMatDoubleScalar(pmat, "cam_x_max");
	int camymin = (int) readMatDoubleScalar(pmat, "cam_y_min");
	int camymax = (int) readMatDoubleScalar(pmat, "cam_y_max");
	
	// store the size of the input data
	int camxsize = camxmax - camxmin;
	int camysize = camymax - camymin;
	
	// compute the physical size of this pattern, and store it in the 'sigma' controls
	SetCtrlVal(TabPage_0, TABPANEL_SigmaX, camxsize * 4.65);
	SetCtrlVal(TabPage_0, TABPANEL_SigmaY, camysize * 4.65);

	// allocate memory for holding the camera frame
	unsigned char* tmpcamframe = (unsigned char*) malloc(camxres * camyres * sizeof(unsigned char));

	// read out the camera data			
	mxArray* mxTemp = matGetVariable(pmat, "cam_frame");

	// copy the data to the average frame array
	memcpy(tmpcamframe, mxGetData(mxTemp), camxres * camyres * sizeof(unsigned char));

	// store the camera data
	SLM_setInputIntensityCameraData(tmpcamframe, camxres, camyres, camxmin, camymin, camxsize, camysize);
	
	// free intermediate memory
	mxDestroyArray(mxTemp);
	free(tmpcamframe);
	
	// close the mat-file
	matClose(pmat);
	
	// set the file textbox
	SetCtrlVal(TabPage_0, TABPANEL_txtInputLoadFromFile, camfilename);
}


/// HIFN Callback for the input intensity type, only really updates the status of the radio buttons, then calls
///      the input intensity callback (except the first time load from file is toggled, calls the file dialog callback)
int CVICALLBACK rbInput_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// switch off all three radio buttons
			SetCtrlVal(panel, TABPANEL_rbInputGaussian, 0);
			SetCtrlVal(panel, TABPANEL_rbInputAnnular, 0);
			SetCtrlVal(panel, TABPANEL_rbInputLoadFromFile, 0);
			
			// determine which of the radiobuttons is clicked
			switch (control)
			{
				case TABPANEL_rbInputAnnular:
				{	
					// update the radio button status
					SetCtrlVal(panel, TABPANEL_rbInputAnnular, 1);
					
					// dim the file dialog controls
					SetCtrlAttribute(panel, TABPANEL_btnInputLoadFromFile, ATTR_DIMMED, 1);
				
					break;
				}	
				case TABPANEL_rbInputLoadFromFile:
				{	
					// update the radio button status
					SetCtrlVal(panel, TABPANEL_rbInputLoadFromFile, 1);
					
					// undim the file dialog controls
					SetCtrlAttribute(panel, TABPANEL_btnInputLoadFromFile, ATTR_DIMMED, 0);
					
					// check if the file input field is filled
					char camfilename[1024];
					GetCtrlVal(panel, TABPANEL_txtInputLoadFromFile, camfilename);
					if (strlen(camfilename) == 0)
					{	
						// no, pop up the load file dialog to load a camera file
						int SelectionStatus = FileSelectPopup("", "*.mat", "*.mat", "Load input intensity pattern", 
									VAL_OK_BUTTON, 0, 1, 1, 1, camfilename);
						if (SelectionStatus)
						{
							InputLoadFromFile(camfilename);	
						}
						else
						{
							// the load failed, revert to Gaussian input intensity
							SetCtrlVal(panel, TABPANEL_rbInputLoadFromFile, 0);
							SetCtrlVal(panel, TABPANEL_rbInputGaussian, 1);
							
							// dim the file dialog controls
							SetCtrlAttribute(panel, TABPANEL_btnInputLoadFromFile, ATTR_DIMMED, 1);
						}
					}
					
					break;
				}	
				case TABPANEL_rbInputGaussian:
				default:
				{	
					// update the radio button status
					SetCtrlVal(panel, TABPANEL_rbInputGaussian, 1);
				
					// dim the file dialog controls
					SetCtrlAttribute(panel, TABPANEL_btnInputLoadFromFile, ATTR_DIMMED, 1);
					
					break;
				}
			}
			
			// recalculate the input intensity (it uses the radio buttons to determine the type)
			InputIntensity_Callback (panel, 0, EVENT_COMMIT, NULL, eventData1, eventData2);
			
			break;
		}
	}
	return 0;
}


/// HIFN Callback function for the button that loads a camera file as input intensity
int CVICALLBACK InputLoadFromFile_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// pop up the load file dialog to load a camera file
			char camfilename[1024];
			int SelectionStatus = FileSelectPopup("", "*.mat", "*.mat", "Load input intensity pattern", 
						VAL_OK_BUTTON, 0, 1, 1, 1, camfilename);
			if (SelectionStatus)
			{
				InputLoadFromFile(camfilename);	
			}
			
			break;
		}
	}
	return 0;
}

int CVICALLBACK ClearSLM_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			SLM_Clear();
			
			// update with the 'new pattern'
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

			break;
	}
	return 0;
}
