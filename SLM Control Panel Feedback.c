//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Callback functions for the Feedback panel
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
#include "SLM Camera Panel.h"

// array holding the signal window with the *desired* signal
static double* gTargetIntensity;

// array holding the camera window with the actual (measured) signal
static double* gCameraSignal;

// array holding the difference of the signal and camera windows
static double* gTargetDifference;

// array holding a new, adjusted signal based on the current signal and the 
// corresponding camera data
static double* gNewSignal;

// array holding the signal in SLM space
static double* gSLMSignal;

// array holding the current correction to the signal
static double* gCorrection;
static double* gPreliminaryCorrection;

// array holding the original phase pattern
static double* gOriginalPhase;

// the current signal quality
static double gSignalQuality = 0.0;

// the signal window dimensions and position, we store both the SLM / camera window size (gSLMWinXsize, gSLMWinYsize), 
// as well as the internally used window (gWXsize, gWYsize), which contains a resampled version of the camera signal, 
// to allow for more accurate image matching
static int gWx, gWy, gSLMWinXsize, gSLMWinYsize, gWXsize, gWYsize;

// counter for the number of iterations for a particular run
static int gIterationCounter, gStepCounter;


void transformCameraSignal(double* SigFrom, int NxFrom, int NyFrom, double* SigTo, int NxTo, int NyTo);
void constructCameraSignalTransformationMatrix(int NxFrom, int NyFrom, int NxTo, int NyTo, int OffsetX, int OffsetY);


void UpdateFeedbackData()
{
	// check if the camera has been calibrated
	if (gCamExtendedCalibration == 1)
	{
		// get the camera window
		getCameraWindow(gWx, gWy, gWXsize, gWYsize, gCameraSignal);	
		SLM_normalise(gCameraSignal, gWXsize * gWYsize, 0, 1);
		
		//update the Target Intensity (for beam splitting feedback)
		// get the desired signal intensity
		double* gTargetIntensity;
		gTargetIntensity = SLM_getUnresampledSignal();
		
		// compute the difference with the target intensity
		int border;
		GetCtrlVal(TabPage_3, TABPANE_10_EdgeBorder, &border);
		for (int l = border; l < gWYsize - border; l++)
		for (int k = border; k < gWXsize - border; k++)			
				gTargetDifference[k + l * gWXsize] = 
					((gCameraSignal[k + l * gWXsize]) - (gTargetIntensity[k + l * gWXsize]));
		
		// compute the rms deviation
		gSignalQuality = 0.0;
		for (int l = border; l < gWYsize - border; l++)
		for (int k = border; k < gWXsize - border; k++)			
			gSignalQuality += gTargetDifference[k + l * gWXsize] * gTargetDifference[k + l * gWXsize];
		gSignalQuality = sqrt(gSignalQuality / (gWXsize * gWYsize));
		
		// print the signal quality
		char val[256];
		sprintf(val, "%.2f (rms)", gSignalQuality);
		SetCtrlVal(TabPage_3, TABPANE_10_txtSignalQualityVal, val);
		
		// by what factor should we correct the signal?
		double beta;
		GetCtrlVal(TabPage_3, TABPANE_10_AdjustmentFactor, &beta);
		
		// the correction should be weighted by the inverse of the signal amplitude,
		// but this can cause trouble when the signal amplitude is small compared to the error,
		// so we ask the user by how much we should weigh the inverse signal into the correction
		// (in the terminology of the paper, we make \beta position dependent, using the original signal Us)
		double sigweight;
		GetCtrlVal(TabPage_3, TABPANE_10_SignalWeight, &sigweight);
		
		// get the maximum value of the signal
		double sigmax = 0.0;
		for (int k = 0; k < gWXsize * gWYsize; k++)
			if (gTargetIntensity[k] > sigmax)
				sigmax = gTargetIntensity[k];
			
		// should we do edge detection?
		int edgedetect;
		GetCtrlVal(TabPage_3, TABPANE_10_EdgeDetection, &edgedetect);
		
		// what kind of drop qualifies as an edge?
		double edgeheight;
		GetCtrlVal(TabPage_3, TABPANE_10_EdgeTreshold, &edgeheight);
		edgeheight *= sigmax;

		// put together a correction to the signal
		for (int l = edgedetect; l < gWYsize - edgedetect; l++)
		for (int k = edgedetect; k < gWXsize - edgedetect; k++)			
		{
			// check if we are close to an edge inside the signal
			int edgedist = edgedetect + 1;
			for (int e = edgedetect; e >= 1; e--)
			{
				double val = gTargetIntensity[k + l * gWXsize];
				if ((fabs(gTargetIntensity[(k + e) + (l + 0) * gWXsize] - val) > edgeheight) ||
					(fabs(gTargetIntensity[(k - e) + (l + 0) * gWXsize] - val) > edgeheight) ||
					(fabs(gTargetIntensity[(k + 0) + (l + e) * gWXsize] - val) > edgeheight) ||
					(fabs(gTargetIntensity[(k + 0) + (l - e) * gWXsize] - val) > edgeheight))
				{
					// yes, store the distance
					edgedist = e;
				}
				
			}
			
			// compute the correction factor for the AMPLITUDE
			double delta;
			if (edgedist < edgedetect + 1)
				delta = 0;
			else
			{
				delta = gCorrection[k + l * gWXsize] 
					- beta * (gTargetDifference[k + l * gWXsize] / (2.0 * (sigweight * sqrt(gTargetIntensity[k + l * gWXsize]) + (1.0 - sigweight))));
			}
			
			// we store the correction as a *temporary* correction, it is only made final by the RefinePhase_Callback
			gPreliminaryCorrection[k + l * gWXsize] = delta;
		}
		
		// compute the adjusted signal based on the correction
		for (int k = 0; k < gWXsize * gWYsize; k++)
		{	
			// calculate the adjusted intensity of this pixel of the signal
			double pixval = (sqrt(gTargetIntensity[k]) + gPreliminaryCorrection[k]) * (sqrt(gTargetIntensity[k]) + gPreliminaryCorrection[k]);
			
			// make sure the total signal is positive
			pixval = (pixval > 0.0 ? pixval : 0.0);

			// apply the signal mask
			gNewSignal[k] = (gTargetIntensity[k] > 0.075 * sigmax ? pixval : 0.0);
		}
		
		// after the first time an update is done, we can refine the phase
		SetCtrlAttribute(TabPage_3, TABPANE_10_RefinePhase, ATTR_DIMMED, 0);
		SetCtrlAttribute(TabPage_3, TABPANE_10_UpdateFeedbackWindow, ATTR_DIMMED, 0);
	}
}


/// HIFN updates the feedback canvas with the current data
void UpdateFeedbackCanvas(int Panel, int Canvas)
{
	// read out the radio buttons to find out what the user wants to see
	int rbSignal, rbCamera, rbDifference, rbAdjustedSignal, rbCorrection;
	GetCtrlVal(Panel, TABPANE_10_rbFeedbackSignal, 	   &rbSignal);
	GetCtrlVal(Panel, TABPANE_10_rbFeedbackCamera, 	   &rbCamera);
	GetCtrlVal(Panel, TABPANE_10_rbFeedbackDifference, &rbDifference);
	GetCtrlVal(Panel, TABPANE_10_rbFeedbackAdjusted,   &rbAdjustedSignal);
	GetCtrlVal(Panel, TABPANE_10_rbFeedbackCorrection, &rbCorrection);

	// plot the requested data on the canvas
	if (rbSignal)
	{
		// plot the signal window on the canvas
		if (gTargetIntensity != NULL)
			plotField(gTargetIntensity, gWXsize, gWYsize, Panel, Canvas);
	}
	if (rbCamera)
	{
		// plot the camera window on the canvas
		if (gCameraSignal != NULL)
			plotField(gCameraSignal, gWXsize, gWYsize, Panel, Canvas);
	}
	if (rbDifference)
	{
		// plot the difference of the signal and camera on the canvas
		if (gTargetDifference != NULL)
			plotField(gTargetDifference, gWXsize, gWYsize, Panel, Canvas);
	}
	if (rbAdjustedSignal)
	{
		// plot the adjusted signal
		if (gNewSignal != NULL)
			plotField(gNewSignal, gWXsize, gWYsize, Panel, Canvas);
	}
	if (rbCorrection)
	{
		// plot the current preliminary correction	
		if (gPreliminaryCorrection != NULL)
			plotField(gPreliminaryCorrection, gWXsize, gWYsize, Panel, Canvas);
	}
}


/// HIFN callback for updating the feedback window canvas
int CVICALLBACK UpdateFeedbackWindow_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// update the data
			UpdateFeedbackData();
			
			// display the data on the canvas
			UpdateFeedbackCanvas(panel, TABPANE_10_FeedbackCanvas);
			
			break;
	}
	return 0;
}


/// HFIN 
int CVICALLBACK rbFeedbackToggle_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		case EVENT_VAL_CHANGED:
			
			// toggle the correct radio button, switch off the others
			SetCtrlAttribute(panel, TABPANE_10_rbFeedbackSignal, 		ATTR_CTRL_VAL, 0);
			SetCtrlAttribute(panel, TABPANE_10_rbFeedbackCamera, 		ATTR_CTRL_VAL, 0);
			SetCtrlAttribute(panel, TABPANE_10_rbFeedbackDifference, 	ATTR_CTRL_VAL, 0);
			SetCtrlAttribute(panel, TABPANE_10_rbFeedbackAdjusted, 		ATTR_CTRL_VAL, 0);
			SetCtrlAttribute(panel, TABPANE_10_rbFeedbackCorrection,	ATTR_CTRL_VAL, 0);
			switch (control)
			{
				case TABPANE_10_rbFeedbackSignal:
					
					SetCtrlAttribute(panel, TABPANE_10_rbFeedbackSignal, 	   ATTR_CTRL_VAL, 1);
					break;
				
				case TABPANE_10_rbFeedbackCamera:
					
					SetCtrlAttribute(panel, TABPANE_10_rbFeedbackCamera,     ATTR_CTRL_VAL, 1);
					break;
					
				case TABPANE_10_rbFeedbackDifference:
					
					SetCtrlAttribute(panel, TABPANE_10_rbFeedbackDifference, ATTR_CTRL_VAL, 1);
					break;
					
				case TABPANE_10_rbFeedbackAdjusted:
					
					SetCtrlAttribute(panel, TABPANE_10_rbFeedbackAdjusted,   ATTR_CTRL_VAL, 1);
					break;
					
				case TABPANE_10_rbFeedbackCorrection:
					
					SetCtrlAttribute(panel, TABPANE_10_rbFeedbackCorrection,   ATTR_CTRL_VAL, 1);
					break;					
			}

			// update the canvas such that it displays the data of choice
			UpdateFeedbackCanvas(panel, TABPANE_10_FeedbackCanvas);
			
			break;
	}
	return 0;
}


/// HIFN Callback function for the refine phase button on the feedback panel
int CVICALLBACK RefinePhase_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// update the step counter (indicates the number of times the user presses 'Refine Phase')
			gStepCounter++;
			
			// should we use MRAF?
			int MRAF;
			GetCtrlVal(panel, TABPANE_10_MRAF, &MRAF);				
		
			// get the number of PF and AF iterations NOTE: check if we need only Amplitude freedom iterations
			int numPF, numAF;
			numPF = 0;
			GetCtrlVal(panel, TABPANE_10_IterationsAF, &numAF);
			
			// get the number of repeats
			int repeat;
			GetCtrlVal(panel, TABPANE_10_Repeat, &repeat);
			
			// allocate array for storing the signal quality as a function of the number of iterations
			double* sigQ = (double*) malloc(repeat * sizeof(double));
			
			// check what we should do with the phase, reset or update?
			int phasereset, phasesignal, phaseupdate;
			GetCtrlVal(panel, TABPANE_10_rbPhaseReset, &phasereset);
			GetCtrlVal(panel, TABPANE_10_rbPhaseSignal, &phasesignal);
			phaseupdate = (phasereset ? 0 : 1);
			/*
			// open a datafile
			char filename[1200];
			char filenamebase[1024];
			GetCtrlVal(panel, TABPANE_10_FeedbackFilename, filenamebase);
			sprintf(filename, "%s_%i.mat", filenamebase, gStepCounter);
			MATFile *pmat = matOpen(filename, "wz");
			
			// write the target signal to the datafile
			writeMatDoubleArray(pmat, "TargetSignal", gTargetIntensity, gWYsize, gWXsize);
			
			// write the original SLM phase to the datafile
			writeMatDoubleArray(pmat, "OriginalPhase", gOriginalPhase, gYsize, gXsize);
			
			// write the feedback settings to the file (see UpdateFeedbackData for description)
			double beta;
			GetCtrlVal(TabPage_3, TABPANE_10_AdjustmentFactor, &beta);
			double sigweight;
			GetCtrlVal(TabPage_3, TABPANE_10_SignalWeight, &sigweight);
			int edgedetect;
			GetCtrlVal(TabPage_3, TABPANE_10_EdgeDetection, &edgedetect);
			double edgeheight;
			GetCtrlVal(TabPage_3, TABPANE_10_EdgeTreshold, &edgeheight);
			int border;
			GetCtrlVal(TabPage_3, TABPANE_10_EdgeBorder, &border);
			writeMatDoubleScalar(pmat, "beta",        (double) beta);		
			writeMatDoubleScalar(pmat, "sigweight",   (double) sigweight);		
			writeMatDoubleScalar(pmat, "edgedetect",  (double) edgedetect);		
			writeMatDoubleScalar(pmat, "edgeheight",  (double) edgeheight);		
			writeMatDoubleScalar(pmat, "edgeborder",  (double) border);		
			writeMatDoubleScalar(pmat, "MRAF",  	  (double) MRAF);		
			writeMatDoubleScalar(pmat, "numPF", 	  (double) numPF);		
			writeMatDoubleScalar(pmat, "numAF", 	  (double) numAF);		
			writeMatDoubleScalar(pmat, "repeat", 	  (double) repeat);		
			writeMatDoubleScalar(pmat, "phasereset",  (double) phasereset);		
			writeMatDoubleScalar(pmat, "phasesignal", (double) phasesignal);		
			writeMatDoubleScalar(pmat, "phaseupdate", (double) phaseupdate);		
			
			// close the file
			matClose(pmat);
			*/
			// perform the feedback iteratively
			for (int n = 0; n < repeat; n++)
			{
				// update the iteration counter
				gIterationCounter++;
				
				// update the iteration counter display				
				char itup[256];
				sprintf(itup, "(%i / %i / %i)", n + 1, gIterationCounter, repeat);
				SetCtrlVal(panel, TABPANE_10_txtIterationUpdate, itup);
				
				// open a datafile
				//sprintf(filename, "%s_%i_%i.mat", filenamebase, gStepCounter, gIterationCounter);
				//MATFile *pmatdata = matOpen(filename, "wz");
				
				//if (n == repeat - 1)
				//{
					// write the modified signal to the datafile
				//	writeMatDoubleArray(pmatdata, "Correction", gCorrection, gWYsize, gWXsize);
				//}
				
				// take a measurement and compute the new signal
				UpdateFeedbackData();
				
				// write the measured signal to the datafile
				//writeMatDoubleArray(pmatdata, "TargetDifference", gTargetDifference, gWYsize, gWXsize);
				
				// write the (estimated) signal quality
				//writeMatDoubleScalar(pmatdata, "SignalQuality", gSignalQuality);
				
				// close the file
				//matClose(pmatdata);
				
				// store the current signal quality, and plot in a debuggraph
				sigQ[n] = gSignalQuality;
				DeleteGraphPlot(gDebugPanel, gDebugGraph2, -1, VAL_DELAYED_DRAW);
				PlotY(gDebugPanel, gDebugGraph2, sigQ, n + 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0xFF);			
				RefreshGraph(gDebugPanel, gDebugGraph2);
				
				// plot the latest camera signal etc. in the debug panel
				plotField(gTargetIntensity, gWXsize, gWYsize, gDebugPanel, gDebugCanvas1);
				plotField(gCameraSignal, gWXsize, gWYsize, gDebugPanel, gDebugCanvas2);
				plotField(gPreliminaryCorrection, gWXsize, gWYsize, gDebugPanel, gDebugCanvas3);
				plotField(gNewSignal, gWXsize, gWYsize, gDebugPanel, gDebugCanvas4);
				
				// make the preliminary correction calculated by UpdateFeedbackData final
				memcpy(gCorrection, gPreliminaryCorrection, gWXsize * gWYsize * sizeof(double));
				
				// should we reset the phase to the original phase for the unaltered signal?
				if (phasesignal)
					SLM_setPhase(gOriginalPhase);
				
				// transform the signal from camera space to SLM space	
				transformCameraSignal(gNewSignal, gWXsize, gWYsize, gSLMSignal, gSLMWinXsize, gSLMWinYsize);
				
				// call the phase generation code with the transformed, corrected signal
				SLM_generatePhase(gSLMSignal, gSLMWinXsize, gSLMWinYsize, 0, 0, gSLMWinXsize, gSLMWinYsize, numPF, numAF, 1, phaseupdate, 1, MRAF, 0);

				//We want to do beamsplitting hence:
				//SLM_UpdateSpotarray(gSLMSignal, pnlSLMpixels, pnlSimPanel, SimPanel_CANVAS, panel, OutputTextBox, 
				
				// update the SLM
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
				
				// wait a short while to make sure the image is adjusted to the new situation
				Delay(1.0);
				
				// call the camera relevant part of the camera timer callback function manually, since the camera 
				// timer events are not processed when the feedback routine is running
				for (int c = 0; c < getCameraNumAveragingFrames(); c++)
				{
					// first wait one frame
					Delay(1.0 / getCameraFrameRate());
					
					// grab a frame
					grabCameraFrame();
				}
				
				// average the camera frames
				performCameraAveraging();
			}
			
			break;
		}
	}
	return 0;
}


int CVICALLBACK FixSignal_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// reset the iteration and step counters
			gIterationCounter = 0;
			//gStepCounter = 0;
			
			// for the moment, feedback only works for beam shaping
			if (SLM_getCurrentPattern() == SLM_BEAMSHAPE_ARB) 
			{	
				/*// debug create artificial mapping
				double xf[9], yf[9];
				xf[0] = -1; xf[3] = -1; xf[6] = -1;
				xf[1] =  0; xf[4] =  0; xf[7] =  0;
				xf[2] =  1; xf[5] =  1; xf[8] =  1;
				
				yf[0] = -1; yf[1] = -1; yf[2] = -1;
				yf[3] =  0; yf[4] =  0; yf[5] =  0;
				yf[6] =  1; yf[7] =  1; yf[8] =  1;
				
				double xto[9], yto[9];
				for (int k = 0; k < 9; k++)
				{
					xf[k] /= 2.0;
					yf[k] /= 2.0;
					xto[k] = 3 * xf[k];
					yto[k] = -2 * yf[k];
				}
				
				double xt, yt;
				transformCoordinate(xf, yf, xto, yto, 3, 3, -0.1, 0.501, &xt, &yt);
				
				*/
				
			
				// get the desired signal intensity
				double* tmpTargetIntensity;
				int ttix, ttiy;
				tmpTargetIntensity = SLM_getUnresampledSignal();
				ttix = gUnresampledSignalXsize;
				ttiy = gUnresampledSignalYsize;
			
				// get the desired signal dimensions from the control panel
				double width, height;
				GetCtrlVal(TabPage_2_1, TABPANEL_8_PictureWidth,  &width);
				GetCtrlVal(TabPage_2_1, TABPANEL_8_PictureHeight, &height);
				
				// convert to meters
				width *= 1.0e-6;
				height *= 1.0e-6;
			
				// get the number of camera pixels the signal comprises
				gWXsize =  width / gCamPixelSize;
				gWYsize = height / gCamPixelSize;
				
				// realloc and resample the target signal to the camera pixel window dimensions
				gTargetIntensity = (double*) realloc(gTargetIntensity, gWXsize * gWYsize * sizeof(double));
				SLM_resampleBitmapInPlace(tmpTargetIntensity, gTargetIntensity, ttix, ttiy, gWXsize, gWYsize);
				SLM_normalise(gTargetIntensity, gWXsize * gWYsize, 0, 1);
				
				// get the number of focal units that the window (approximately) comprises in SLM space
				// NOTE: we cheat a little by extending the size as long as the 'border' option has bugs
				gSLMWinXsize = 1.1 * width / SLM_getFocalUnitX();				
				gSLMWinYsize = 1.1 * height / SLM_getFocalUnitY();
				
				// do we want to add a border in SLM space?
				double border;
				GetCtrlVal(panel, TABPANE_10_WinBorder, &border);
				border *= 1.0e-6;
				int BorderX = border / SLM_getFocalUnitX();
				int BorderY = border / SLM_getFocalUnitY();
				gSLMWinXsize += 2 * BorderX;
				gSLMWinYsize += 2 * BorderY;
				
				// get the camera window origin
				// (this is the mapped image of the origin of the SLM coordinate system)
				getCameraWindowOriginPixel(&gWx, &gWy);
				
				// debug plot the target intensity
				plotField(gTargetIntensity, gWXsize, gWYsize, gDebugPanel, gDebugCanvas1);
			
				// create a suitable window area array in SLM space
				gSLMSignal = (double*) realloc(gSLMSignal, gSLMWinXsize * gSLMWinYsize * sizeof(double));
				
				/*// debug test the transformation camera -> SLM
				memset(gCameraSignal, 0, gWXsize * gWYsize * sizeof(double));
				#define NUMP 20
				double xwidth = 0.6 * gSLMWinXsize * SLM_getFocalUnitX();
				double ywidth = 0.6 * gSLMWinYsize * SLM_getFocalUnitY();
				double x0 = 0.2 * gSLMWinXsize * SLM_getFocalUnitX();
				double y0 = 0.2 * gSLMWinYsize * SLM_getFocalUnitY();
				for (int nx = 0; nx < NUMP; nx++)
				for (int ny = 0; ny < NUMP; ny++)
				{
					double x = nx * xwidth / (NUMP - 1) + x0;
					double y = ny * ywidth / (NUMP - 1) + y0;
					
					double xc, yc;
					transformCoordinateSLMToCamera(x, y, &xc, &yc);
					int ix = (int) (xc / gCamPixelSize) - gWx;
					int iy = (int) (yc / gCamPixelSize) - gWy;
					if ((ix >= 0) && (ix < gWXsize) && (iy >= 0) && (iy < gWYsize))
						gCameraSignal[ix + iy * gWXsize] = 1.0;
					
				}
				plotField(gCameraSignal, gWXsize, gWYsize, gDebugPanel, gDebugCanvas3); */
				
				// construct the mapping matrix between the SLM space and the camera space
				constructCameraSignalTransformationMatrix(gWXsize, gWYsize, gSLMWinXsize, gSLMWinYsize, BorderX, BorderY);
			
				// map the signal from camera space to SLM space
				transformCameraSignal(gTargetIntensity, gWXsize, gWYsize, gSLMSignal, gSLMWinXsize, gSLMWinYsize);
				
				// DEBUG save to file
				MATFile *pmat;
				pmat = matOpen("TransformedSignal.mat", "wz");
				writeMatDoubleArray(pmat, "gSLMSignal", gSLMSignal, gSLMWinXsize, gSLMWinYsize);
				writeMatDoubleArray(pmat, "gTargetIntensity", gTargetIntensity, gWXsize, gWYsize);
				matClose(pmat);				

				
				
				// debug plot the transformed intensity
				plotField(gSLMSignal, gSLMWinXsize, gSLMWinYsize, gDebugPanel, gDebugCanvas2);
			
				// get the beam shaping parameters from the panel
				int numPF, numAF;
				GetCtrlVal(TabPage_2_1, TABPANEL_8_IterationsPF, &numPF);
				GetCtrlVal(TabPage_2_1, TABPANEL_8_IterationsAF, &numAF);
				
				// should we use MRAF?
				int MRAF;
				GetCtrlVal(panel, TABPANE_10_MRAF, &MRAF);				

				// redo the beam shaping with the transformed signal
				SLM_generatePhase(gSLMSignal, gSLMWinXsize, gSLMWinYsize, BorderX, BorderY, gSLMWinXsize, gSLMWinYsize, numPF, numAF, 1, 0, 1, MRAF, 0);				
				
				// update the SLM
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
				
				// allocate memory for the new signal composed at each iteration, the camera signal, the difference 
				// between the target and camera intensity, and two arrays for holding the corrections to the signal amplitude (not the intensity)
				// and finally the (camera) signal mapped onto SLM space
				gNewSignal        	   = (double*) realloc(gNewSignal,        	   gWXsize * gWYsize * sizeof(double));
				gCameraSignal     	   = (double*) realloc(gCameraSignal,     	   gWXsize * gWYsize * sizeof(double));
				gTargetDifference 	   = (double*) realloc(gTargetDifference, 	   gWXsize * gWYsize * sizeof(double));
				gCorrection		  	   = (double*) realloc(gCorrection,       	   gWXsize * gWYsize * sizeof(double));
				gPreliminaryCorrection = (double*) realloc(gPreliminaryCorrection, gWXsize * gWYsize * sizeof(double));
				memset(gCorrection, 		   0, gWXsize * gWYsize * sizeof(double));
				memset(gPreliminaryCorrection, 0, gWXsize * gWYsize * sizeof(double));
				memset(gNewSignal, 			   0, gWXsize * gWYsize * sizeof(double));
				memset(gTargetDifference,      0, gWXsize * gWYsize * sizeof(double));

				// store the current SLM phase pattern, in case we want each new adjusted phase pattern calculation to be started from this phase pattern
				gOriginalPhase = (double*) realloc(gOriginalPhase, gXsize * gYsize * sizeof(double));
				memcpy(gOriginalPhase, SLM_getPhase(), gXsize * gYsize * sizeof(double));
				
				UpdateFeedbackData();
			}			
			   /*
				// Feedback for Beamsplitting
			if (SLM_getCurrentPattern() == SLM_BEAMSPLIT_IFTA) 
			{	
			
				// get the desired signal intensity
				double* tmpTargetIntensity;
				int ttix, ttiy;
				tmpTargetIntensity = SLM_getUnresampledSignal();
				ttix = gUnresampledSignalXsize;
				ttiy = gUnresampledSignalYsize;
				
				// get some parameters from beam splitting
				int sigxoffset, sigyoffset, numxspots, numyspots, spotxspacing, spotyspacing;
				GetCtrlVal(TabPage_1_2, TABPANEL_5_SigYOffset, &sigyoffset);
				GetCtrlVal(TabPage_1_2, TABPANEL_5_SigXOffset, &sigxoffset);
				GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYSpacing, &spotyspacing);
				GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXSpacing, &spotxspacing);
				GetCtrlVal(TabPage_1_2, TABPANEL_5_NumYSpots, &numyspots);
				GetCtrlVal(TabPage_1_2, TABPANEL_5_NumXSpots, &numxspots);
				
				// get the desired signal dimensions from the control panel
				double width, height;
				width = 2 * sigxoffset + numxspots + (numxspots - 1) * spotxspacing;
				height= 2 * sigyoffset + numyspots + (numyspots - 1) * spotyspacing;
				
				// convert to meters
				width *= 1.0e-6;
				height *= 1.0e-6;
			
				// get the number of camera pixels the signal comprises
				gWXsize =  width / gCamPixelSize;
				gWYsize = height / gCamPixelSize;
				
				// realloc and resample the target signal to the camera pixel window dimensions
				gTargetIntensity = (double*) realloc(gTargetIntensity, gWXsize * gWYsize * sizeof(double));
				SLM_resampleBitmapInPlace(tmpTargetIntensity, gTargetIntensity, ttix, ttiy, gWXsize, gWYsize);
				SLM_normalise(gTargetIntensity, gWXsize * gWYsize, 0, 1);
				
				// get the number of focal units that the window (approximately) comprises in SLM space
				// NOTE: we cheat a little by extending the size as long as the 'border' option has bugs
				gSLMWinXsize = 1.1 * width / SLM_getFocalUnitX();				
				gSLMWinYsize = 1.1 * height / SLM_getFocalUnitY();
				
				// do we want to add a border in SLM space?
				double border;
				GetCtrlVal(panel, TABPANE_10_WinBorder, &border);
				border *= 1.0e-6;
				int BorderX = border / SLM_getFocalUnitX();
				int BorderY = border / SLM_getFocalUnitY();
				gSLMWinXsize += 2 * BorderX;
				gSLMWinYsize += 2 * BorderY;
				
				// get the camera window origin
				// (this is the mapped image of the origin of the SLM coordinate system)
				getCameraWindowOriginPixel(&gWx, &gWy);
				
				// debug plot the target intensity
				plotField(gTargetIntensity, gWXsize, gWYsize, gDebugPanel, gDebugCanvas1);
			
				// create a suitable window area array in SLM space
				gSLMSignal = (double*) realloc(gSLMSignal, gSLMWinXsize * gSLMWinYsize * sizeof(double));
				
				/ debug test the transformation camera -> SLM
				memset(gCameraSignal, 0, gWXsize * gWYsize * sizeof(double));
				#define NUMP 20
				double xwidth = 0.6 * gSLMWinXsize * SLM_getFocalUnitX();
				double ywidth = 0.6 * gSLMWinYsize * SLM_getFocalUnitY();
				double x0 = 0.2 * gSLMWinXsize * SLM_getFocalUnitX();
				double y0 = 0.2 * gSLMWinYsize * SLM_getFocalUnitY();
				for (int nx = 0; nx < NUMP; nx++)
				for (int ny = 0; ny < NUMP; ny++)
				{
					double x = nx * xwidth / (NUMP - 1) + x0;
					double y = ny * ywidth / (NUMP - 1) + y0;
					
					double xc, yc;
					transformCoordinateSLMToCamera(x, y, &xc, &yc);
					int ix = (int) (xc / gCamPixelSize) - gWx;
					int iy = (int) (yc / gCamPixelSize) - gWy;
					if ((ix >= 0) && (ix < gWXsize) && (iy >= 0) && (iy < gWYsize))
						gCameraSignal[ix + iy * gWXsize] = 1.0;
					
				}
				plotField(gCameraSignal, gWXsize, gWYsize, gDebugPanel, gDebugCanvas3); 
				
				// construct the mapping matrix between the SLM space and the camera space
				constructCameraSignalTransformationMatrix(gWXsize, gWYsize, gSLMWinXsize, gSLMWinYsize, BorderX, BorderY);
			
				// map the signal from camera space to SLM space
				transformCameraSignal(gTargetIntensity, gWXsize, gWYsize, gSLMSignal, gSLMWinXsize, gSLMWinYsize);
				
				// DEBUG save to file
				MATFile *pmat;
				pmat = matOpen("TransformedSignal.mat", "wz");
				writeMatDoubleArray(pmat, "gSLMSignal", gSLMSignal, gSLMWinXsize, gSLMWinYsize);
				writeMatDoubleArray(pmat, "gTargetIntensity", gTargetIntensity, gWXsize, gWYsize);
				matClose(pmat);				

				
				
				// debug plot the transformed intensity
				plotField(gSLMSignal, gSLMWinXsize, gSLMWinYsize, gDebugPanel, gDebugCanvas2);
			
				// get the beam shaping parameters from the panel
				int numPF, numAF;
				GetCtrlVal(TabPage_2_1, TABPANEL_8_IterationsPF, &numPF);
				GetCtrlVal(TabPage_2_1, TABPANEL_8_IterationsAF, &numAF);
				
				// should we use MRAF?
				int MRAF;
				GetCtrlVal(panel, TABPANE_10_MRAF, &MRAF);				

				// redo the beam shaping with the transformed signal
				SLM_generatePhase(gSLMSignal, gSLMWinXsize, gSLMWinYsize, BorderX, BorderY, gSLMWinXsize, gSLMWinYsize, numPF, numAF, 1, 0, 1, MRAF, 0);				
				
				// update the SLM
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
				
				// allocate memory for the new signal composed at each iteration, the camera signal, the difference 
				// between the target and camera intensity, and two arrays for holding the corrections to the signal amplitude (not the intensity)
				// and finally the (camera) signal mapped onto SLM space
				gNewSignal        	   = (double*) realloc(gNewSignal,        	   gWXsize * gWYsize * sizeof(double));
				gCameraSignal     	   = (double*) realloc(gCameraSignal,     	   gWXsize * gWYsize * sizeof(double));
				gTargetDifference 	   = (double*) realloc(gTargetDifference, 	   gWXsize * gWYsize * sizeof(double));
				gCorrection		  	   = (double*) realloc(gCorrection,       	   gWXsize * gWYsize * sizeof(double));
				gPreliminaryCorrection = (double*) realloc(gPreliminaryCorrection, gWXsize * gWYsize * sizeof(double));
				memset(gCorrection, 		   0, gWXsize * gWYsize * sizeof(double));
				memset(gPreliminaryCorrection, 0, gWXsize * gWYsize * sizeof(double));
				memset(gNewSignal, 			   0, gWXsize * gWYsize * sizeof(double));
				memset(gTargetDifference,      0, gWXsize * gWYsize * sizeof(double));

				// store the current SLM phase pattern, in case we want each new adjusted phase pattern calculation to be started from this phase pattern
				gOriginalPhase = (double*) realloc(gOriginalPhase, gXsize * gYsize * sizeof(double));
				memcpy(gOriginalPhase, SLM_getPhase(), gXsize * gYsize * sizeof(double));
				
				UpdateFeedbackData();
			}
			*/ 
			
			break;
	}
	return 0;
}



int CVICALLBACK PhaseOption_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		case EVENT_VAL_CHANGED:

			// toggle the correct radio button, switch off the others
			SetCtrlAttribute(panel, TABPANE_10_rbPhaseReset,  ATTR_CTRL_VAL, 0);
			SetCtrlAttribute(panel, TABPANE_10_rbPhaseSignal, ATTR_CTRL_VAL, 0);
			SetCtrlAttribute(panel, TABPANE_10_rbPhaseUpdate, ATTR_CTRL_VAL, 0);
			
			switch (control)
			{
				case TABPANE_10_rbPhaseReset:
					
					SetCtrlAttribute(panel, TABPANE_10_rbPhaseReset,  ATTR_CTRL_VAL, 1);
					break;
				
				case TABPANE_10_rbPhaseSignal:
					
					SetCtrlAttribute(panel, TABPANE_10_rbPhaseSignal, ATTR_CTRL_VAL, 1);
					break;
					
				case TABPANE_10_rbPhaseUpdate:
					
					SetCtrlAttribute(panel, TABPANE_10_rbPhaseUpdate, ATTR_CTRL_VAL, 1);
					break;
					
			}
			break;
	}
	return 0;
}


// clear the current correction
int CVICALLBACK ClearCorr_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// reset the iteration and step counters
			gIterationCounter = 0;
			//gStepCounter = 0;
			
			// clear the correction
			memset(gCorrection, 		   0, gWXsize * gWYsize * sizeof(double));
			memset(gPreliminaryCorrection, 0, gWXsize * gWYsize * sizeof(double));
			
			// set the original signal phase back
			SLM_setPhase(gOriginalPhase);
			
			// update the SLM
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

			break;
	}
	return 0;
}


// arrays describing a sparse mapping matrix transforming camera space to SLM space
static int* gMapMatrixColIndices;
static int* gMapMatrixNumIndicesPerRow;
static int* gMapMatrixRowStartIndices;
double* gMapMatrixValues;

double minFour(double* a)
{
	double min = a[0];
	for (int k = 1; k < 4; k++)
		if (a[k] < min)
			min = a[k];
	
	return min;
}


double maxFour(double* a)
{
	double max = a[0];
	for (int k = 1; k < 4; k++)
		if (a[k] > max)
			max = a[k];
	
	return max;
}

void constructCameraSignalTransformationMatrix(int NxFrom, int NyFrom, int NxTo, int NyTo, int OffsetX, int OffsetY)
{
	
	//TODO: ADD PROGRESS INDICATOR
	
	// construct arrays with the row and column indices and corresponding values
	int blocksize = 1024 * 1024;
	int numblocks = 1;
	int* colindices = (int*) malloc(blocksize * sizeof(int));
	int* rowindices = (int*) malloc(blocksize * sizeof(int));
	double* values  = (double*) malloc(blocksize * sizeof(double));
	
	// determine the matrix dimensions
	int numrows = NxTo * NyTo;
	int numcols = NxFrom * NyFrom;
	int numvalues = 0;
	
	// SLM space focal unit size
	double SX = SLM_getFocalUnitX();
	double SY = SLM_getFocalUnitY();
	
	// offset of SLM pixels
	// (we take half a pixel, such that the origin pixel is centered around (0, 0) in SLM space)
	double pixeloffsetx = -0.5;
	double pixeloffsety = -0.5;
	
	// DEBUG store transformed pixels
	/*double* cx0 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cx1 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cx2 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cx3 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cy0 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cy1 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cy2 = (double*) malloc(NxTo * NyTo * sizeof(double));
	double* cy3 = (double*) malloc(NxTo * NyTo * sizeof(double));	   */
	
	double cx[4], cy[4]; 
	
	// loop over all the pixels in the destination array (these are SLM space pixels)
	for (int l = 0; l < NyTo; l++)
	for (int k = 0; k < NxTo; k++)		
	{
		// get the four corner coordinates of this pixel in SLM space (starting upper right, going counterclockwise)
		double sx[4], sy[4];
		sx[0] = (((double) k) + 1.0 + pixeloffsetx - OffsetX) * SX;
		sx[1] = (((double) k) + 0.0 + pixeloffsetx - OffsetX) * SX;
		sx[2] = sx[1];
		sx[3] = sx[0];
		sy[0] = (((double) l) + 0.0 + pixeloffsety - OffsetY) * SY;
		sy[1] = sy[0];
		sy[2] = (((double) l) + 1.0 + pixeloffsety - OffsetY) * SY;
		sy[3] = sy[2];
	
		// get the corresponding coordinates in camera space
		// re use half the coordinates from the previous iteration!
		if (k > 0)
		{
			cx[1] = cx[0];
			cy[1] = cy[0];
			cx[2] = cx[3];
			cy[2] = cy[3];
		}
		else
		{
			transformCoordinateSLMToCamera(sx[1], sy[1], &(cx[1]), &(cy[1]));	
			transformCoordinateSLMToCamera(sx[2], sy[2], &(cx[2]), &(cy[2]));    
		}
		transformCoordinateSLMToCamera(sx[0], sy[0], &(cx[0]), &(cy[0]));
		transformCoordinateSLMToCamera(sx[3], sy[3], &(cx[3]), &(cy[3]));
		
		
		// DEBUG variables
		/*cx0[k + l * NxTo] = cx[0];
		cx1[k + l * NxTo] = cx[1];
		cx2[k + l * NxTo] = cx[2];
		cx3[k + l * NxTo] = cx[3];
		cy0[k + l * NxTo] = cy[0];
		cy1[k + l * NxTo] = cy[1];
		cy2[k + l * NxTo] = cy[2];
		cy3[k + l * NxTo] = cy[3];  */
		
		
		// create a bounding box for the mapped pixel
		double ulx = minFour(cx);
		double uly = minFour(cy);
		double lrx = maxFour(cx);
		double lry = maxFour(cy);
		
		// get the coordinates of the pixel at the origin of the camera window
		int cwox, cwoy;
		getCameraWindowOriginPixel(&cwox, &cwoy);
		
		
		// TODO: From the bounding box we can easily find the overlapping pixels immediately,
		//       without having to loop over all of them!!
		
		int ixmin = (ulx / gCamPixelSize) - 2 - cwox;
		int ixmax = (lrx / gCamPixelSize) + 1 - cwox;
		ixmin = (ixmin >= 0 ? ixmin : 0);
		ixmax = (ixmax <= NxFrom ? ixmax : NxFrom);
		
		int jymax = (lry / gCamPixelSize) + 1 - cwoy;
		int jymin = (uly / gCamPixelSize) - 2 - cwoy;
		jymin = (jymin >= 0 ? jymin : 0);
		jymax = (jymax <= NyFrom ? jymax : NyFrom);
		
		// loop over all camera space pixels
		for (int i = ixmin; i < ixmax; i++)
		for (int j = jymin; j < jymax; j++)
		{
			// get the current camera pixel coordinates
			double px[4], py[4];
			px[0] = (i + 1 + cwox) * gCamPixelSize; 
			px[1] = (i + 0 + cwox) * gCamPixelSize; 
			px[2] = (i + 0 + cwox) * gCamPixelSize;
			px[3] = (i + 1 + cwox) * gCamPixelSize;
			py[0] = (j + 0 + cwoy) * gCamPixelSize;
			py[1] = (j + 0 + cwoy) * gCamPixelSize;
			py[2] = (j + 1 + cwoy) * gCamPixelSize;
			py[3] = (j + 1 + cwoy) * gCamPixelSize;
			
			
			
			// DEBUG
			
				// create a Matlab file pointer for saving the data
				//MATFile *pmat;

				// open the .mat file
				//pmat = matOpen("debug.mat", "wz");
				
				//writeMatDoubleArray(pmat, "px", px, 1, 4);
				//writeMatDoubleArray(pmat, "py", py, 1, 4);
			

			
			// check if the current camera pixel has any chance of overlapping with the mapped pixel, based on the bounding box
			if (   !(px[0] < ulx)     
				&& !(px[1] > lrx) 
				&& !(py[3] < uly) 
				&& !(py[0] > lry) )
			{		
				// determine its overlap with the mapped pixel
				double vx[16];
				double vy[16];
				double tmpvx[16];
				double tmpvy[16];
				
				// copy the four original points to the list of vertices, and add the first point also to the end (!)
				int nv = 5, nvt;
				memcpy(vx, cx, 4 * sizeof(double));
				memcpy(vy, cy, 4 * sizeof(double));
				vx[4] = vx[0];
				vy[4] = vy[0];
				
				//writeMatDoubleArray(pmat, "vx0", vx, 1, nv);
				//writeMatDoubleArray(pmat, "vy0", vy, 1, nv);
				
				// variables indicating whether a point is inside or outside the half space defined by the clipping edge
				int previn, currin, firstin;
		
				// note that the surrounding pixel (in camera space) is a square, so the Sutherland-Hodgman algorithm
				// greatly simplifies

				// we will continue by treating each edge of the camera pixel in turn
				
				// RIGHT VERTICAL
				{
					// initialise: if the first point lies inside, add it
					previn = 0;
					firstin = 0;
					nvt = 0;
					if (vx[0] < px[0])
					{
						// this point lies inside compared to the edge, move to the list of 'in'
						tmpvx[nvt] = vx[0];
						tmpvy[nvt] = vy[0];
						nvt++;
						previn = 1;
						firstin = 1;
					}
					
					// loop over the remaining points in vx
					for (int v = 1; v < nv; v++)
					{				
						// is this point inside?
						currin = 0;
						if (vx[v] < px[0])
							currin = 1;
						
						// compute the intersection (if necessary)
						double ix, iy;
						if (previn != currin)
						{
							double slope = (vy[v] - vy[v - 1]) / (vx[v] - vx[v - 1]);
							iy = (px[0] - vx[v - 1]) * slope + vy[v - 1];
							ix = px[0];
						}
						
						// check the 4 possible cases
						if ((previn == 1) && (currin == 1))
						{
							// both in, just add the current point
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 0) && (currin == 1))
						{
							// moving from out to in, add the intersection AND the point (in that order)
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 1) && (currin == 0))
						{
							// moving from in to out, add the intersection
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							
							previn = 0;
						}
						else if ((previn == 0) && (currin == 0))
						{
							// both out, do nothing
							previn = 0;
						}
						
					}
					
					// check if the first point was added
					// if not, a new (intersection) point was added as the first point, and we have
					// to close the loop with that point instead
					if (!firstin)
					{
						tmpvx[nvt] = tmpvx[0];
						tmpvy[nvt] = tmpvy[0];
						nvt++;
					}
					
					// copy tmpvx to vx
					memcpy(vx, tmpvx, nvt * sizeof(double));
					memcpy(vy, tmpvy, nvt * sizeof(double));
					nv = nvt;
					
				}   // (done with RIGHT VERTICAL edge)
				
				
				//writeMatDoubleArray(pmat, "vx1", vx, 1, nv);
				//writeMatDoubleArray(pmat, "vy1", vy, 1, nv);
				
				
				
				
				
				
				
				// UP HORIZONTAL
				if (nv > 2)
				{
					// initialise: if the first point lies inside, add it
					previn = 0;
					firstin = 0;
					nvt = 0;
					if (vy[0] > py[0])
					{
						// this point lies inside compared to the edge, move to the list of 'in'
						tmpvx[nvt] = vx[0];
						tmpvy[nvt] = vy[0];
						nvt++;
						previn = 1;
						firstin = 1;
					}
					
					// loop over the remaining points in vx
					for (int v = 1; v < nv; v++)
					{				
						// is this point inside?
						currin = 0;
						if (vy[v] > py[0])
							currin = 1;
						
						// compute the intersection (if necessary)
						double ix, iy;
						if (previn != currin)
						{
							double slope = (vx[v] - vx[v - 1]) / (vy[v] - vy[v - 1]);
							ix = (py[0] - vy[v - 1]) * slope + vx[v - 1];
							iy =  py[0];
						}
						
						// check the 4 possible cases
						if ((previn == 1) && (currin == 1))
						{
							// both in, just add the current point
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 0) && (currin == 1))
						{
							// moving from out to in, add the intersection AND the point (in that order)
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 1) && (currin == 0))
						{
							// moving from in to out, add the intersection
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							
							previn = 0;
						}
						else if ((previn == 0) && (currin == 0))
						{
							// both out, do nothing
							previn = 0;
						}
						
					}
					
					// check if the first point was added
					// if not, a new (intersection) point was added as the first point, and we have
					// to close the loop with that point instead
					if (!firstin)
					{
						tmpvx[nvt] = tmpvx[0];
						tmpvy[nvt] = tmpvy[0];
						nvt++;
					}
					
					// copy tmpvx to vx
					memcpy(vx, tmpvx, nvt * sizeof(double));
					memcpy(vy, tmpvy, nvt * sizeof(double));
					nv = nvt;
					
				}   // (done with UP HORIZONTAL edge)
				
				
					//writeMatDoubleArray(pmat, "vx2", vx, 1, nv);
				//writeMatDoubleArray(pmat, "vy2", vy, 1, nv);
				
				
				
				
				
				// LEFT VERTICAL
				if (nv > 2)
				{
					// initialise: if the first point lies inside, add it
					previn = 0;
					firstin = 0;
					nvt = 0;
					if (vx[0] > px[1])
					{
						// this point lies inside compared to the edge, move to the list of 'in'
						tmpvx[nvt] = vx[0];
						tmpvy[nvt] = vy[0];
						nvt++;
						previn = 1;
						firstin = 1;
					}
					
					// loop over the remaining points in vx
					for (int v = 1; v < nv; v++)
					{				
						// is this point inside?
						currin = 0;
						if (vx[v] > px[1])
							currin = 1;
						
						// compute the intersection (if necessary)
						double ix, iy;
						if (previn != currin)
						{
							double slope = (vy[v] - vy[v - 1]) / (vx[v] - vx[v - 1]);
							iy = (px[1] - vx[v - 1]) * slope + vy[v - 1];
							ix = px[1];
						}
						
						// check the 4 possible cases
						if ((previn == 1) && (currin == 1))
						{
							// both in, just add the current point
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 0) && (currin == 1))
						{
							// moving from out to in, add the intersection AND the point (in that order)
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 1) && (currin == 0))
						{
							// moving from in to out, add the intersection
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							
							previn = 0;
						}
						else if ((previn == 0) && (currin == 0))
						{
							// both out, do nothing
							previn = 0;
						}
						
					}
					
					// check if the first point was added
					// if not, a new (intersection) point was added as the first point, and we have
					// to close the loop with that point instead
					if (!firstin)
					{
						tmpvx[nvt] = tmpvx[0];
						tmpvy[nvt] = tmpvy[0];
						nvt++;
					}
					
					// copy tmpvx to vx
					memcpy(vx, tmpvx, nvt * sizeof(double));
					memcpy(vy, tmpvy, nvt * sizeof(double));
					nv = nvt;
					
				}   // (done with LEFT VERTICAL edge)
				
				
				
									  //writeMatDoubleArray(pmat, "vx3", vx, 1, nv);
				//writeMatDoubleArray(pmat, "vy3", vy, 1, nv);
				
				
				
				// DOWN HORIZONTAL
				if (nv > 2)
				{
					// initialise: if the first point lies inside, add it
					previn = 0;
					firstin = 0;
					nvt = 0;
					if (vy[0] < py[2])
					{
						// this point lies inside compared to the edge, move to the list of 'in'
						tmpvx[nvt] = vx[0];
						tmpvy[nvt] = vy[0];
						nvt++;
						previn = 1;
						firstin = 1;
					}
					
					// loop over the remaining points in vx
					for (int v = 1; v < nv; v++)
					{				
						// is this point inside?
						currin = 0;
						if (vy[v] < py[2])
							currin = 1;
						
						// compute the intersection (if necessary)
						double ix, iy;
						if (previn != currin)
						{
							double slope = (vx[v] - vx[v - 1]) / (vy[v] - vy[v - 1]);
							ix = (py[2] - vy[v - 1]) * slope + vx[v - 1];
							iy =  py[2];
						}
						
						// check the 4 possible cases
						if ((previn == 1) && (currin == 1))
						{
							// both in, just add the current point
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 0) && (currin == 1))
						{
							// moving from out to in, add the intersection AND the point (in that order)
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							tmpvx[nvt] = vx[v];
							tmpvy[nvt] = vy[v];
							nvt++;
							
							previn = 1;
						}
						else if ((previn == 1) && (currin == 0))
						{
							// moving from in to out, add the intersection
							tmpvx[nvt] = ix;
							tmpvy[nvt] = iy;
							nvt++;
							
							previn = 0;
						}
						else if ((previn == 0) && (currin == 0))
						{
							// both out, do nothing
							previn = 0;
						}
						
					}
					
					// check if the first point was added
					// if not, a new (intersection) point was added as the first point, and we have
					// to close the loop with that point instead
					if (!firstin)
					{
						tmpvx[nvt] = tmpvx[0];
						tmpvy[nvt] = tmpvy[0];
						nvt++;
					}
					
					// copy tmpvx to vx
					memcpy(vx, tmpvx, nvt * sizeof(double));
					memcpy(vy, tmpvy, nvt * sizeof(double));
					nv = nvt;
					
				}   // (done with DOWN HORIZONTAL edge)
				
				
				//writeMatDoubleArray(pmat, "vx4", vx, 1, nv);
				//writeMatDoubleArray(pmat, "vy4", vy, 1, nv);
				
				
				// the list nv should now contain all the vertices of the polygon that
				// resulted from clipping the mapped pixel with the current (square) pixel
				
				
				// compute the area of the resulting polygon: https://en.wikipedia.org/wiki/Polygon#Area_and_centroid
				double A = 0.0;
				for (int n = 0; n < nv - 1; n++)
					A -= 0.5 * (vx[n] * vy[n + 1] - vx[n + 1] * vy[n]);	
		
				// was there any overlap in the end?
				if (A > 0.0)
				{
					// yes, add the pixel indices and overlap to the mapping matrix arrays
					rowindices[numvalues] = k + l * NxTo;
					colindices[numvalues] = i + j * NxFrom;
					values[numvalues] = A / (gCamPixelSize * gCamPixelSize);
					numvalues++;
		
					// increase the array sizes if necessary
					if (numvalues > numblocks * blocksize)
					{
						numblocks++;
						rowindices = (int*)    realloc(rowindices, numblocks * blocksize * sizeof(int));
						colindices = (int*)    realloc(colindices, numblocks * blocksize * sizeof(int));
						values     = (double*) realloc(values,     numblocks * blocksize * sizeof(double));
					}
				}
			}
			
			//matClose(pmat);
		}
	}
	
	
	
	// we have now created arrays with all the row and column indices and the corresponding values
	// let's turn that into a sparse matrix (borrowed tested code from Spin1 project)
	
	
	
	// allocate memory for the indices and values of the sparse matrix
    gMapMatrixColIndices       = (int*) realloc(gMapMatrixColIndices,       numvalues * sizeof(int));
    gMapMatrixNumIndicesPerRow = (int*) realloc(gMapMatrixNumIndicesPerRow, numrows * sizeof(int));
    gMapMatrixRowStartIndices  = (int*) realloc(gMapMatrixRowStartIndices,  numrows * sizeof(int));
	gMapMatrixValues           = (double*) realloc(gMapMatrixValues, numvalues * sizeof(double));

    // copy the column indices
    memcpy(gMapMatrixColIndices, colindices, numvalues * sizeof(int));
	
	// copy the values
	memcpy(gMapMatrixValues, values, numvalues * sizeof(double));

    // copy and arrange the row indices, looping over all rows of the matrix
    int currrowindex = 0;
    for (int k = 0; k < numrows; k++)
    {
        // write the start of the row index
        gMapMatrixRowStartIndices[k] = currrowindex;

        // count the number of entries for this row
        int count = 0;
        if (numvalues > 0)
            while (rowindices[currrowindex + count] == k)
                count++;

        // store the count
        gMapMatrixNumIndicesPerRow[k] = count;

        // update the current row index for the next iteration
        currrowindex += count;
    }
	
	
	/// DEBUG SAVE TO FILE
	
	// create a Matlab file pointer for saving the data
	MATFile *pmat;

	// open the .mat file
	pmat = matOpen("MapMatrix2.mat", "wz");
	
	// write the matrix initialisation list
	writeMatIntArray(pmat, "rowindices", rowindices, numvalues);
	writeMatIntArray(pmat, "colindices", colindices, numvalues);
	writeMatDoubleArray(pmat, "values",     values,  numvalues, 1);
	/*writeMatDoubleArray(pmat, "cx0", cx0, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cx1", cx1, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cx2", cx2, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cx3", cx3, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cy0", cy0, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cy1", cy1, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cy2", cy2, NxTo, NyTo);
	writeMatDoubleArray(pmat, "cy3", cy3, NxTo, NyTo);  */

	// write the input arguments to the function
	writeMatDoubleScalar(pmat, "NxFrom", NxFrom);
	writeMatDoubleScalar(pmat, "NyFrom", NyFrom);
	writeMatDoubleScalar(pmat, "NxTo", NxTo);
	writeMatDoubleScalar(pmat, "NyTo", NyTo);
	
	matClose(pmat);
	
	
				   /*
	
	// create mxArray for holding the values of the matrix
    mxArray* mxTemp = mxCreateSparse((mwSize) numrows, (mwSize) numcols, (mwSize) numvalues, mxREAL);

    // store the number of nonzero elements
    mxSetNzmax(mxTemp, (mwSize) numvalues);

    // get pointers to the row indices (ir), and column indices (jc) arrays
    mwIndex* mxJc = mxGetJc(mxTemp);
    mwIndex* mxIr = mxGetIr(mxTemp);

    // get pointer to the values (pr) array
    double* mxPr = mxGetPr(mxTemp);

    // fill the mxArray's ROW indices (ir) with the Matrix' COLUMN indices, essentially storing the matrix
    // in a transposed fashion; we do this because Matlab uses column-major storage, whereas we use row-major
    // FIXME: store matrix untransposed; fttb we only expect to use symmetric matrices so it doesn't matter
    int count = 0;
    for (int k = 0; k < numrows; k++)
    for (int l = 0; l < gMapMatrixNumIndicesPerRow[k]; l++)
        mxIr[count++] = (mwIndex) gMapMatrixColIndices[l + gMapMatrixRowStartIndices[k]];

    // also fill the column indices (jc) array
    for (int k = 0; k < numrows; k++)
        mxJc[k] = (mwIndex) gMapMatrixRowStartIndices[k];
    mxJc[numrows] = numvalues;

    // finally, copy the values
	memcpy(mxPr, gMapMatrixValues, numvalues * sizeof(double));

    // write mxArray to file
    matPutVariable(pmat, "MapMatrix", mxTemp);

    // clean up the mxArray
    //mxDestroyArray(mxTemp);
	
	matClose(pmat);	
	
					*/
	
	// clean up the temporary arrays
	free(rowindices);
	free(colindices);
	free(values);
}


// transforms the signal from Camera space to SLM space
void transformCameraSignal(double* SigFrom, int NxFrom, int NyFrom, double* SigTo, int NxTo, int NyTo)
{
	// the transformation is done using a multiplication by a sparse matrix, with each row corresponding
	// to a pixel in the destination array, and each column corresponding with a pixel in the source
	// array, and the value at (row, col) indicating how much of pixel 'col' contributed to pixel 'row'
	
	// loop over all the pixels in the destination array, i.e. the rows of the matrix
    for (int k = 0; k < NxTo; k++)
	for (int l = 0; l < NyTo; l++)
    {
		// compute the row index
		int row = k + l * NxTo;
		
        // loop over all columns
        double destval = 0.0;
        int rowstart = gMapMatrixRowStartIndices[row];
        for (int c = 0; c < gMapMatrixNumIndicesPerRow[row]; c++)
        {
            // read data elements
            int colindex = gMapMatrixColIndices[c + rowstart];
            double src   = SigFrom[colindex];
            double mul   = gMapMatrixValues[rowstart + c];

            // perform multiplication and store
            destval += mul * src;
        }

        // store the final result in the destination vector
        SigTo[row] = destval;
    }
}

int CVICALLBACK RunBatchFeedback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// get the filename
			char filenamebase[200];
			GetCtrlVal(panel, TABPANE_10_FeedbackFilename, filenamebase);
			
			// get the min and max adj and the steps
			double adjmin, adjmax;
			int nadj;
			GetCtrlVal(panel, TABPANE_10_adjmin, &adjmin);
			GetCtrlVal(panel, TABPANE_10_adjmax, &adjmax);
			GetCtrlVal(panel, TABPANE_10_adjmax, &nadj);
			
			// the number of iterations to be taken for adjmax
			int itadjmax;
			GetCtrlVal(panel, TABPANE_10_itadjmax, &itadjmax);
			
			for (int k = 0; k < nadj; k++)
			{
				// compute the adjustment factor
				double adj = adjmin + k * ((adjmax - adjmin) / (((double) nadj) - 1.0));
				
				// compute the number of iterations
				int it = (int) ((adjmax / adj) * ((double) itadjmax));
				
				// clear the current correction
				ClearCorr_Callback (panel, control, EVENT_COMMIT, NULL, eventData1, eventData2);
				
				// update the filename
				char filename[220];
				sprintf(filename, "%s_adj%.2f.mat", filenamebase, adj);
				SetCtrlVal(panel, TABPANE_10_FeedbackFilename, filename); 
				
				// change the control settings for the adjustment factor and the iterations
				SetCtrlVal(panel, TABPANE_10_AdjustmentFactor, adj);
				SetCtrlVal(panel, TABPANE_10_Repeat, it);
				
				// delay a little bit to let the camera adapt to the cleared correction
				Delay(1.0);
				
				// call the camera relevant part of the camera timer callback function manually, since the camera 
				// timer events are not processed when the feedback routine is running
				for (int c = 0; c < getCameraNumAveragingFrames(); c++)
				{
					// first wait one frame
					Delay(1.0 / getCameraFrameRate());
					
					// grab a frame
					grabCameraFrame();
				}
				
				// average the camera frames
				performCameraAveraging();
				
				// call the 'refine' callback	
				RefinePhase_Callback (panel, control, EVENT_COMMIT, NULL, eventData1, eventData2);
			}
			

			break;
	}
	return 0;
}

int CVICALLBACK itcalc (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// get the min and max adj and the steps
			double adjmin, adjmax;
			int nadj;
			GetCtrlVal(panel, TABPANE_10_adjmin, &adjmin);
			GetCtrlVal(panel, TABPANE_10_adjmax, &adjmax);
			GetCtrlVal(panel, TABPANE_10_nadj, &nadj);
			
			// the number of iterations to be taken for adjmax
			int itadjmax;
			GetCtrlVal(panel, TABPANE_10_itadjmax, &itadjmax);
			
			// count the total number of iterations ('repeats')
			int it = 0;
			for (int k = 0; k < nadj; k++)
			{
				// compute the adjustment factor
				double adj = adjmin + k * ((adjmax - adjmin) / (((double) nadj) - 1.0));
				
				// compute the number of iterations
				it += (int) ((adjmax / adj) * ((double) itadjmax));
			}
			int itadjmin = (int) ((adjmax / adjmin) * ((double) itadjmax));
			
			// update the label
			char itup[256];
			sprintf(itup, "total: %i (%i)", it, itadjmin);
			SetCtrlVal(panel, TABPANE_10_totalit, itup);

			break;
	}
	return 0;
}
