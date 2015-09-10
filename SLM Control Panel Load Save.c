//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Routines for loading/saving SLM and camera states.
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
#include "SLM Camera Panel.h"
#include "SLM_internal.h"


 
// Make the file name for the correction file
char FactoryCorrectionFilename[1024];
char SHCorrectionFilename[1024];


/// HIFN Callback function for the button for saving the current SLM phase pattern
int CVICALLBACK SaveSLMPattern_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// let the user specify a filename
			char slmfilename[1024];
		    int SelectionStatus = FileSelectPopup("", "*.slm", "*.slm", "Load SLM pattern", 
									VAL_OK_BUTTON, 0, 1, 1, 1, slmfilename);
			
			if (SelectionStatus > 0)
			{
				// display the name of the selected file in the textbox
				SetCtrlVal(panel, TABPANEL_9_LoadFilename, slmfilename);
				
				// create a Matlab file pointer for saving the data
				MATFile *pmat;

				// open the .mat file
				pmat = matOpen(slmfilename, "wz");

				// write control panel pattern settings
				WritePatternSettings(pmat);
				
				// write internal SLM settings
				// NOTE: we could have gotten these also from the control panel
				SLM_WriteSettingsToFile(pmat);

				// close the mat-file
				matClose(pmat);
				
				// did we also want to save the camera state?
				int savecam;
				GetCtrlVal(panel, TABPANEL_9_LoadSaveCameraState, &savecam);
				if (savecam)
				{
					// construct a filename from the .slm filename,
					char camfilename[1024];
					strcpy(camfilename, slmfilename);
					strcpy((camfilename + StringLength(slmfilename) - 3), "cam");	
					
					// save the camera state
					CameraSaveMatFile(camfilename, 0);
				}

			}

			break;
		}
	}
	return 0;
}


/// HIFN Callback function for button for loading an SLM phase pattern
int CVICALLBACK LoadSLMPattern_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// let the user specify a filename
			char slmfilename[1024];
		    int SelectionStatus = FileSelectPopup("", "*.slm", "*.slm", "Load SLM pattern", 
									VAL_OK_BUTTON, 0, 1, 1, 1, slmfilename);
			
			// check if a file has been succesfully selected
			if (SelectionStatus > 0)
			{
				// load the SLM settings
				LoadSLMSettings(slmfilename);
				
				// did we also want to load the camera state?
				int loadcam;
				GetCtrlVal(panel, TABPANEL_9_LoadSaveCameraState, &loadcam);
				if (loadcam)
				{
					// construct a filename from the slm data filename
					char camfilename[1024];
					strcpy(camfilename, slmfilename);
					strcpy((camfilename + StringLength(slmfilename) - 3), "cam");	
					
					// load the camera state
					if (FileExists(camfilename, 0))
						CameraLoadMatFile(camfilename);				
				}
			}
			
			break;
		}
	}
	return 0;
}


/// HIFN loads SLM settings from a file
void LoadSLMSettings(char slmfilename[])
{
	// create a file pointer and open the specified file
	MATFile *pmat = matOpen(slmfilename, "r");
	
	// load the SLM resolutions that are specified in the file
	double SLM_x_size = readMatDoubleScalar(pmat, "SLM_x_size");
	double SLM_y_size = readMatDoubleScalar(pmat, "SLM_y_size");
	
	// load the subsampling factor
	int SLM_subsampling = readMatDoubleScalar(pmat, "SLM_SubSampleFactor");
	SetCtrlVal(TabPage_0, TABPANEL_SubSample,	SLM_subsampling);
	
	// check if these resolutions match the current resolution
	if (!(((int) SLM_x_size == SLM_getXres()) && ((int) SLM_y_size == SLM_getYres())))
	{				
		// no, set the SLM size to those of the file
		
		// first resize the panel
		SetPanelSize (pnlSLMpixels, (int) SLM_y_size * SLM_subsampling, (int) SLM_x_size * SLM_subsampling);
	
		// adjust the canvas size accordingly
		SetCtrlAttribute(pnlSLMpixels, SLMpixels_SLMcanvas, ATTR_WIDTH,  (int) SLM_x_size * SLM_subsampling);
		SetCtrlAttribute(pnlSLMpixels, SLMpixels_SLMcanvas, ATTR_HEIGHT, (int) SLM_y_size * SLM_subsampling);

		// get the size of the simulation canvas, we need that to reinitialise the SLM
		int SimX, SimY;
		GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_WIDTH,  &SimX);
		GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_HEIGHT, &SimY);

		// all data structures of the SLM need to be resized to acommodate the new 
		// number of pixels, so we just reintialise the SLM in its entirety			
		SLM_initialise((int) SLM_x_size * SLM_subsampling, (int) SLM_y_size * SLM_subsampling, SimX, SimY, SLM_subsampling);
	}
	
	// load and set the slm phase pattern
	mxArray* mxTemp = matGetVariable(pmat, "SLM_phase");
	SLM_setPhase(mxGetPr(mxTemp));
	mxDestroyArray(mxTemp);
	
	// load the signal + window
	mxArray* mxTempSignal = matGetVariable(pmat, "SLM_signal");
	gSignal = (double*) realloc(gSignal, SLM_getXres() * SLM_getYres() * sizeof(double));
	memcpy(gSignal, mxGetPr(mxTempSignal), SLM_getXres() * SLM_getYres() * sizeof(double));
	mxDestroyArray(mxTempSignal);
	gWindowX = (int) readMatDoubleScalar(pmat, "SLM_window_x");
	gWindowY = (int) readMatDoubleScalar(pmat, "SLM_window_y");
	gWindowXsize = (int) readMatDoubleScalar(pmat, "SLM_window_x_size");
	gWindowYsize = (int) readMatDoubleScalar(pmat, "SLM_window_y_size");
	gFocalUnitX  = readMatDoubleScalar(pmat, "SLM_focal_x");
	gFocalUnitY  = readMatDoubleScalar(pmat, "SLM_focal_y");
	
	// check if the current general settings are to be retained
	int KeepSettings; 
	GetCtrlVal(TabPage_5, TABPANEL_9_LoadKeepSettings, &KeepSettings);
	if (!KeepSettings)
	{
		// set all the control values to those from the file
		SetCtrlVal(TabPage_0, TABPANEL_LensXphase,	readMatDoubleScalar(pmat, "SLM_LensXphase") / 1.0e-3);
		SetCtrlVal(TabPage_0, TABPANEL_LensYphase,	readMatDoubleScalar(pmat, "SLM_LensYphase") / 1.0e-3);
		SetCtrlVal(TabPage_0, TABPANEL_VertTrans,	readMatDoubleScalar(pmat, "SLM_VertTrans")  / 1.0e-6);
		SetCtrlVal(TabPage_0, TABPANEL_HorizTrans,	readMatDoubleScalar(pmat, "SLM_HorizTrans") / 1.0e-6);
		SetCtrlVal(TabPage_0, TABPANEL_SubSample,	(int) readMatDoubleScalar(pmat, "SLM_SubSampleFactor"));
		
		// the bias value also needs to be converted
		double bias = 256.0 * readMatDoubleScalar(pmat, "SLM_Bias") / (2 * PI);
		SetCtrlVal(TabPage_0, TABPANEL_Bias, (int) bias);
	}
	
	// check if the current input intensity settings are to be retained
	int KeepInputIntensity;
	GetCtrlVal(TabPage_5, TABPANEL_9_LoadKeepInpInt, &KeepInputIntensity);
	if (!KeepInputIntensity)
	{
		// set all the control values to those of the file
		SetCtrlVal(TabPage_0, TABPANEL_Xpos, 			1.0e3 * readMatDoubleScalar(pmat, "input_x_center"));
		SetCtrlVal(TabPage_0, TABPANEL_Ypos, 			1.0e3 * readMatDoubleScalar(pmat, "input_y_center"));
		SetCtrlVal(TabPage_0, TABPANEL_SigmaX, 			1.0e6 * readMatDoubleScalar(pmat, "input_x_sigma"));
		SetCtrlVal(TabPage_0, TABPANEL_SigmaY, 			1.0e6 * readMatDoubleScalar(pmat, "input_y_sigma"));
		SetCtrlVal(TabPage_0, TABPANEL_LensFocalLength, 1.0e3 * readMatDoubleScalar(pmat, "SLM_focallength"));
		SetCtrlVal(TabPage_0, TABPANEL_Wavelength, 		1.0e9 * readMatDoubleScalar(pmat, "SLM_wavelength"));
	}
	
	// update the SLM settings from the control values (also updates the SLM itself)
	setDefaultValuesFromPanel();
	
	// close the mat-file
	matClose(pmat);	
}

// Callback for loading a factory phase correction
int CVICALLBACK LoadFactoryCorrection_Callback(int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{  
			// let the user specify a filename
		    int SelectionStatus = FileSelectPopup("", "*.png", "*.jpg;*.tif;*.pcx;*.bmp;*.dib;*.rle;*.ico;*.png;*.wfm;*.emf", "Load Factory Correction", 
									VAL_OK_BUTTON, 0, 0, 1, 1, FactoryCorrectionFilename);
			
			// check the check box for loading the correction
			int FactoryCorrectionCheck;
			GetCtrlVal(panel, TABPANEL_9_FactorySwitch, &FactoryCorrectionCheck);
			
			// check if a file has been succesfully selected
			if (SelectionStatus > 0)
			{
				// load the SLM settings if the factory correction is switched on
				if (FactoryCorrectionCheck) 
				LoadFactoryCorrection(FactoryCorrectionFilename, FactoryCorrectionCheck);
				
			}						 
			// display the file name onto the load panel
			SetCtrlVal(TabPage_5, TABPANEL_9_FactoryCorrection, FactoryCorrectionFilename);
			
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


int CVICALLBACK FactorySwitch_Callback(int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			{
				int FactoryCorrectionCheck;
				GetCtrlVal(panel, TABPANEL_9_FactorySwitch, &FactoryCorrectionCheck);
				if (FactoryCorrectionCheck == FALSE)
				{
					LoadFactoryCorrection(0, 0);
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
				}
				if (FactoryCorrectionCheck == TRUE)
				{
					LoadFactoryCorrection(FactoryCorrectionFilename, FactoryCorrectionCheck);
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
				}
				break;
			}
	}
	return 0;
}

// Callback for loading a SH lens surface phase correction
int CVICALLBACK LoadSHcorrection_Callback(int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{  
			// let the user specify a filename
		    int SelectionStatus = FileSelectPopup("", "*.png", "*.jpg;*.tif;*.pcx;*.bmp;*.dib;*.rle;*.ico;*.png;*.wfm;*.emf", "Load SH Correction", 
									VAL_OK_BUTTON, 0, 0, 1, 1, SHCorrectionFilename);
			
			// check the check box for loading the correction
			int SHcorrectionCheck;
			GetCtrlVal(panel, TABPANEL_9_SHswitch, &SHcorrectionCheck);
			
			// check if a file has been succesfully selected
			if (SelectionStatus > 0)
			{
				// load the SLM settings if the factory correction is switched on
				if (SHcorrectionCheck) 
				LoadSHcorrection(SHCorrectionFilename, SHcorrectionCheck);
				
			}
			// display the file name onto the load panel
			SetCtrlVal(TabPage_5, TABPANEL_9_SHcorrection, SHCorrectionFilename);
			
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

int CVICALLBACK SHswitch_Callback(int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			{
				int SHcorrectionCheck;
				GetCtrlVal(panel, TABPANEL_9_SHswitch, &SHcorrectionCheck);
				if (SHcorrectionCheck == FALSE)
				{
					LoadSHcorrection(0, 0);
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
				}
				if (SHcorrectionCheck == TRUE)
				{
					LoadSHcorrection(SHCorrectionFilename, SHcorrectionCheck);
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
				}
				break;
			}
	}
	return 0;
}
	

int CVICALLBACK IntensityMatrix_Callback(int panel, int control, int event, void *callbackData, 
		int eventData1, int eventData2)   
{
	switch(event)
	{
		case EVENT_COMMIT:
			{
				// Get location of the box of the intensities wanted (in pixels)
				int Xmin, Xmax, Ymin, Ymax;
				GetCtrlVal(panel, TABPANEL_9_Xmax, &Xmax);
				GetCtrlVal(panel, TABPANEL_9_Xmin, &Xmin);
				GetCtrlVal(panel, TABPANEL_9_Ymax, &Ymax);
				GetCtrlVal(panel, TABPANEL_9_Ymin, &Ymin);
				double *intensitymatrix = GiveIntensities(Xmax, Xmin, Ymax, Ymin);
				// write to a matlab file
				// let the user specify a filename
				char matrixfilename[1024];
				
			    int SelectionStatus = FileSelectPopup("", "*.mat", "*.mat", "Save Intensity Matrix", 
										VAL_OK_BUTTON, 0, 1, 1, 1, matrixfilename);
			
				if (SelectionStatus > 0)
				{
					// create and open a Matlab file for saving the data
					// create a Matlab file pointer for saving the data
					MATFile *pmat;

					// open the .mat file
					pmat = matOpen(matrixfilename, "wz");
					
					//write intensity matrix to the mat file
					writeMatDoubleArray(pmat, "intensitymatrix", intensitymatrix,(Xmax - Xmin),(Ymax - Ymin));
					
					// close the mat-file
					matClose(pmat);
				}
				break;
			}
	}
	return 0;
}

/// HIFN writes the control panel pattern settings to file
void WritePatternSettings(MATFile *pmat)
{
	// write the SLM pattern settings
	writeMatDoubleScalar(pmat, "SLM_PatternMode", (double) SLM_getCurrentPattern());	         
	switch (SLM_getCurrentPattern())
	{
	
		case SLM_BEAMSHAPE_STD:
		{
			// read out the current settings from the control panel
			double width, height, sigmax, sigmay;
			int issquare, isgauss;
			GetCtrlVal(TabPage_2_0, TABPANEL_7_SignalWidth,  &width);
			GetCtrlVal(TabPage_2_0, TABPANEL_7_SignalHeight, &height);
			GetCtrlVal(TabPage_2_0, TABPANEL_7_SigmaX, &sigmax);
			GetCtrlVal(TabPage_2_0, TABPANEL_7_SigmaY, &sigmay);
			GetCtrlVal(TabPage_2_0, TABPANEL_7_BeamTypeSquare,   &issquare);
			GetCtrlVal(TabPage_2_0, TABPANEL_7_BeamTypeGaussian, &isgauss);
		
			// write them to file
			writeMatDoubleScalar(pmat, "pattern_width",  width);
			writeMatDoubleScalar(pmat, "pattern_height", height);	         
			writeMatDoubleScalar(pmat, "pattern_sigmax", sigmax);	         
			writeMatDoubleScalar(pmat, "pattern_sigmay", sigmay);	         
			writeMatDoubleScalar(pmat, "pattern_square", (double) issquare);
			writeMatDoubleScalar(pmat, "pattern_gauss",  (double) isgauss);	         
			
			break;
		}
		case SLM_BEAMSHAPE_ARB:
		{
			// read out the current settings from the control panel
			double awidth, aheight;
			int numPF, numAF;
			GetCtrlVal(TabPage_2_1, TABPANEL_8_PictureWidth,  &awidth);
			GetCtrlVal(TabPage_2_1, TABPANEL_8_PictureHeight, &aheight);
			GetCtrlVal(TabPage_2_1, TABPANEL_8_IterationsPF,  &numPF);
			GetCtrlVal(TabPage_2_1, TABPANEL_8_IterationsAF,  &numAF);
		
			// write them to file
			writeMatDoubleScalar(pmat, "pattern_width",  awidth);
			writeMatDoubleScalar(pmat, "pattern_height", aheight);	         
			writeMatDoubleScalar(pmat, "pattern_numPF", (double) numPF);	         
			writeMatDoubleScalar(pmat, "pattern_numAF", (double) numAF);
		
			break;
		}
		case SLM_BEAMSPLIT_PROJECT:
		
			break;
		
		case SLM_BEAMSPLIT_IFTA:
		{
			// read out the current settings from the control panel
			int Nx, Ny, sigxoffset, sigyoffset, spotxoffset, spotyoffset, xspacing, yspacing, phaseconstraint;
			double phasestep;
			GetCtrlVal(TabPage_1_2, TABPANEL_5_NumXSpots,    &Nx);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_NumYSpots,    &Ny);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXOffset,  &spotxoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYOffset,  &spotyoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SigXOffset,   &sigxoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SigYOffset,   &sigyoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXSpacing, &xspacing);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYSpacing, &yspacing);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_PhaseConstraint, &phaseconstraint);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_PhaseStep,    &phasestep);
		
			// write them to file
			writeMatDoubleScalar(pmat, "pattern_Nx", 		  (double) Nx);
			writeMatDoubleScalar(pmat, "pattern_Ny", 		  (double) Ny);
			writeMatDoubleScalar(pmat, "pattern_spotxoffset", (double) spotxoffset);
			writeMatDoubleScalar(pmat, "pattern_spotyoffset", (double) spotyoffset);
			writeMatDoubleScalar(pmat, "pattern_sigxoffset",  (double) sigxoffset);
			writeMatDoubleScalar(pmat, "pattern_sigyoffset",  (double) sigyoffset);
			writeMatDoubleScalar(pmat, "pattern_xspacing",    (double) xspacing);
			writeMatDoubleScalar(pmat, "pattern_yspacing",    (double) yspacing);
			writeMatDoubleScalar(pmat, "pattern_phaseconstraint", (double) phaseconstraint);
			writeMatDoubleScalar(pmat, "pattern_phasestep",   phasestep);
			
			// write the theoretical spot intensities and diffraction efficiency
			writeMatDoubleArray(pmat, "IFTASimulatedSpotIntensities", gIFTASimulatedSpotIntensities, Nx, Ny);
			writeMatDoubleScalar(pmat, "IFTASimulatedDiffractionEfficiency", gIFTASimulatedDiffractionEfficiency);
														 
			break;
		}
		case SLM_BEAMSPLIT_PHASEGRATING:
		
			break;
	}
}

