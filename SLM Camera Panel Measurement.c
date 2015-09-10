//==============================================================================
//
// Title:       Camera Controller
// Purpose:     Library of functions to perform some specific camera measurements
//
// Created on:  12-12-2011 at 12:40:12 by Rick.
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
#include "toolbox.h"
#include "SLM.h"
#include "Flycapture2/FlyCapture2_C.h"
#include "SLM Camera Panel.h"
#include "libmat/mat.h" 
#include "libmat/matrix.h"
#include "SLM_internal.h"
#include "SLM Control Panel Internal.h"  




// filename prefix for storing the amplitude modulation measurements
static char gAmplModFilePrefix[1024];

// indicator variable to indicate whether the measurent is running
static int gRunAmplMod = 0;

// variables for the amplitude measurement
static int gCurrBias = 0, gFrameBiasSet, *gBiases;

// parameters for the spot sequence measurements: all are arrays of values, through which the measurement then cycles
static int *gMnumxspots, *gMnumyspots, *gMspotxspacing, *gMspotyspacing, *gMsigxoffset, *gMsigyoffset, *gMspotxoffset, *gMspotyoffset, gMphaseconstraint = 0; 
static double *gMaberration, *gMlensxphase, *gMlensyphase;
static double *gMxaberration;
static double *gMyaberration;
static double *gMaberranisotr;
static double gMphasestep = 0.0;

// array of randomised measurement indices
static int *gMindices;

// the spot sequence counters and the total number of measurements in the sequence
static int gSpotSeqCounter, gFrameSpotSeqSet, gNumMeasurements;

// filename prefix for storing the amplitude modulation measurements
static char gSpotSeqFilePrefix[1024];

// indicator variable to indicate whether the spot sequence measurements are running
static int gRunSpotSeq = 0;



/// HIFN generates an array of randomly permuted indices from 0..N
int* RandomIndices(int N)
{
	// allocate an array with the indices that we have to permute randomly
	int* Indices = (int*) malloc(N * sizeof(int));
	for (int k = 0; k < N; k++)
		Indices[k] = k;
	
	// now we loop over all the indices, and swap each with a random other
	for (int k = 0; k < N; k++)
	{
		// generate the index to swap this value with
		int swapindex = (rand() % N);		
		
		// swap the values
		int tmp = Indices[k];
		Indices[k] = Indices[swapindex];
		Indices[swapindex] = tmp;
	}
	
	return Indices;
}


/// HIFN save the camera frame as a Matlab .mat file
void CameraSaveMatFile(char matfilename[], int saveslmdata)
{
	// create a Matlab file pointer for saving the raw data
	MATFile *pmat;

	// open the .mat file
	pmat = matOpen(matfilename, "wz");

	// create an array for the camera data
	mxArray *paCamData;
	paCamData = mxCreateNumericMatrix(gCamX, gCamY, mxUINT8_CLASS, mxREAL);

	// copy the camera data to the array
	memcpy((void*) mxGetData(paCamData), (void*) gAvgFrame, gCamX * gCamY * sizeof(unsigned char));			

	// write the data to the mat file
	matPutVariable(pmat, "cam_frame", paCamData);

	// destroy the mxArray
	mxDestroyArray(paCamData);

	// write the pixel numbers of the camera feed
	writeMatDoubleScalar(pmat, "cam_x_size", (double) gCamX);			
	writeMatDoubleScalar(pmat, "cam_y_size", (double) gCamY);			

	// write the corner coordinates of the current zoom on the canvas
	writeMatDoubleScalar(pmat, "cam_x_min", (double) gCamXmin);
	writeMatDoubleScalar(pmat, "cam_x_max", (double) gCamXmax);
	writeMatDoubleScalar(pmat, "cam_y_min", (double) gCamYmin);
	writeMatDoubleScalar(pmat, "cam_y_max", (double) gCamYmax);

	// write the center positions
	writeMatDoubleScalar(pmat, "cam_x_center", (double) gCamCenterX);			
	writeMatDoubleScalar(pmat, "cam_y_center", (double) gCamCenterY);

	// write the number of frames that is being averaged over
	writeMatDoubleScalar(pmat, "cam_numframes", (double) gNumFrames);
	
	// write the zoom factor
	writeMatDoubleScalar(pmat, "cam_zoomfactor", (double) gCamZoomFactor);

	// write the center positions
	writeMatDoubleScalar(pmat, "cam_x_center", (double) gCamCenterX);			
	writeMatDoubleScalar(pmat, "cam_y_center", (double) gCamCenterY);

	// write the total frame number
	writeMatDoubleScalar(pmat, "cam_framenumber", (double) gFrameNumber);
		
	// write the calibration data
	writeMatDoubleScalar(pmat, "cam_calibrated", (double) gCamCalibrated);
	if (gCamCalibrated)
	{
		writeMatDoubleScalar(pmat, "cam_calibration_Ox", (double) gOx);
		writeMatDoubleScalar(pmat, "cam_calibration_Oy", (double) gOy);
		writeMatDoubleScalar(pmat, "cam_calibration_Mx", gCamMx);
		writeMatDoubleScalar(pmat, "cam_calibration_My", gCamMy);
		writeMatDoubleScalar(pmat, "cam_calibration_MPx", gCamMPx);
		writeMatDoubleScalar(pmat, "cam_calibration_MPy", gCamMPy);
	}

	// do we also want to save the SLM data?
	if (saveslmdata)
	{
		// SLM internal settings
		SLM_WriteSettingsToFile(pmat);
		
		// SLM settings from the control panel
		WritePatternSettings(pmat);
	}

	// close the mat-file
	matClose(pmat);
}


/*/// HIFN starts a measurement sequence for the amplitude modulation
int CVICALLBACK MeasureAmplitude_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			if (gRunAmplMod)
			{
				// the measurement was already running, stop it	
				SetCtrlAttribute(pnlCamera, PANEL_AmplModTimer, ATTR_ENABLED, 0);
				gRunAmplMod = 0;
			}
			else
			{
				gRunAmplMod = 1;
				gCurrBias = 0;
				
				// initialise an array of randomly permuted bias values
				gBiases = RandomIndices(256);
				
				// get the filename prefix
				GetCtrlVal(pnlCamera, PANEL_AmplModFilePrefix, gAmplModFilePrefix);

				// start the timer
				SetCtrlAttribute(pnlCamera, PANEL_AmplModTimer, ATTR_ENABLED, 1);
			}
			
			break;
	}
	return 0;
}		*/


	/*
/// HIFN the callback for the amplitude modulation timer, each time it fires
///      it measures the amplitude modulation and stores a datafile
int CVICALLBACK AmplModTimer_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_TIMER_TICK:
			
			// check if the current frame is at least gNumFrames frames
			// after the bias was set
			if (gFrameNumber > gFrameBiasSet + 3 * gNumFrames)
			{
				// yes, we can assume that the current camera frame has adjusted to the state of the SLM
				
				// generate a new filename
				char filename[1024];
				sprintf(filename, "%s%i.mat", gAmplModFilePrefix, gCurrBias);
				
				// store the frame
				CameraSaveMatFile(filename, 1);
				
				// update the bias 
				gCurrBias++;
				SLM_setBias(gBiases[gCurrBias % 256]); 
				
				// update the SLM with the new bias
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
				gFrameBiasSet = gFrameNumber;       
				
				// update the camera panel bias text label
				char tmp[100];
				sprintf(tmp, "%i", gCurrBias);
				SetCtrlVal(pnlCamera, PANEL_CurrAmplModBias, tmp);
			}

			break;
	}
	return 0;
}			  */


/// HIFN Callback for starting / stopping the spot sequence measurement
int CVICALLBACK SpotSequence_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			if (gRunSpotSeq)
			{
				// the measurement was already running, stop it	
				SetCtrlAttribute(pnlCamera, PANEL_SpotSeqTimer, ATTR_ENABLED, 0);
				gRunSpotSeq = 0;
			}
			else
			{
				// start the spot measurement sequence
				gRunSpotSeq = 1;
				
				// the array dimensions
				int numnumspots = 1;
				int numspots[1] = {6};//{2, 4, 8, 16, 32};
				
				// the number of steps for the offsets
				int numoffsets = 1;
				int offsets[1] = {0};//{2, 4, 8, 16, 32, 64};
				
				// the number of spot spacings (NOTE: this is the y-spacing, the x-spacing is 2 * y-spacing)
				int numspacings = 1;
				int spacings[1] = {20};//{2, 4, 8, 16};
				
				// the number of steps for the lens focal length corrections (2 of them!)
				int nlx = 1;
				int nly = 100;
				
				// the ranges and values of the lens parameters
				double mm = 1e-3;
				double lrange = 0.2 * mm;
				double lxcenter = 0 * mm;
				double lycenter = 0 * mm;
				double lxmin = lxcenter;// - lrange;
				double lxmax = lxcenter ;//+ lrange;
				double lymin = lycenter - lrange;
				double lymax = lycenter + lrange;
				
				// the number of times a measurement should be repeated
				int numrepeats = 1;
				
				// calculate the number of measurements
				gNumMeasurements = numnumspots * numoffsets * numspacings * nlx * nly * numrepeats;
				
				// check if there is a sequence to resume?
				if ((gSpotSeqCounter > 0) && (gSpotSeqCounter < gNumMeasurements))
				{
					// yes ....	
				}
				else
				{
					// no sequence to resume, we want to start a new one
				
					// reset the spot sequence counter
					gSpotSeqCounter = 0;
				
					// allocate the arrays with parameter values
					gMnumxspots    = (int*) realloc(gMnumxspots,    gNumMeasurements * sizeof(int));
					gMnumyspots    = (int*) realloc(gMnumyspots,    gNumMeasurements * sizeof(int));
					gMspotxspacing = (int*) realloc(gMspotxspacing, gNumMeasurements * sizeof(int));
					gMspotyspacing = (int*) realloc(gMspotyspacing, gNumMeasurements * sizeof(int));
					gMspotxoffset  = (int*) realloc(gMspotxoffset,  gNumMeasurements * sizeof(int));
					gMspotyoffset  = (int*) realloc(gMspotyoffset,  gNumMeasurements * sizeof(int));
					gMsigxoffset   = (int*) realloc(gMsigxoffset,  gNumMeasurements * sizeof(int));
					gMsigyoffset   = (int*) realloc(gMsigyoffset,  gNumMeasurements * sizeof(int));
					gMlensxphase   = (double*) realloc(gMlensxphase,  gNumMeasurements * sizeof(double));
					gMlensyphase   = (double*) realloc(gMlensyphase,  gNumMeasurements * sizeof(double));
				
					// fill the arrays with parameter values
					for (int n = 0; n < numnumspots; n++)
					for (int o = 0; o < numoffsets; o++)
					for (int s = 0; s < numspacings; s++)
					for (int lx = 0; lx < nlx; lx++)
					for (int ly = 0; ly < nly; ly++)
					for (int r = 0; r < numrepeats; r++)
					{
						// calculate the current index into the array of measurements
						int index = n  +
									o  * numnumspots +
									s  * numnumspots * numoffsets +
									lx * numnumspots * numoffsets * numspacings +
									ly * numnumspots * numoffsets * numspacings * nlx + 
									r  * numnumspots * numoffsets * numspacings * nlx * nly;
					
						// the number of spots
						gMnumxspots[index] = numspots[n];
						gMnumyspots[index] = numspots[n];
						
						// the spot offsets within the signal area
						gMspotxoffset[index] = 2 * offsets[o];
						gMspotyoffset[index] = 1 * offsets[o];
					
						// the spot spacing
						gMspotxspacing[index] = 2 * spacings[s];
						gMspotyspacing[index] = 1 * spacings[s];
					
						// the signal area offsets
						gMsigxoffset[index] = 10;
						gMsigyoffset[index] = 10;
						
						// the lens phase
						gMlensxphase[index] = lxmin + lx * (lxmax - lxmin) / (nlx > 1 ? ((double) (nlx - 1)) : 1.0);
						gMlensyphase[index] = lymin + ly * (lymax - lymin) / (nly > 1 ? ((double) (nly - 1)) : 1.0);
					}
				
					// generate an array of randomised measurement indices
					gMindices = RandomIndices(gNumMeasurements);
				
					// get the filename prefix
					GetCtrlVal(pnlCamera, PANEL_SpotSeqFilePrefix, gSpotSeqFilePrefix);
				
					// write a separate file with the measurement sequence, linked with the datafile-names
					char matfilename[1024];
					sprintf(matfilename, "%s_sequence.mat", gSpotSeqFilePrefix);
					MATFile *pmat;
					pmat = matOpen(matfilename, "wz");
					writeMatIntArray(pmat,  "gMnumxspots",    gMnumxspots,    gNumMeasurements);
					writeMatIntArray(pmat,  "gMnumyspots",    gMnumyspots,    gNumMeasurements);
					writeMatIntArray(pmat,  "gMspotxspacing", gMspotxspacing, gNumMeasurements);
					writeMatIntArray(pmat,  "gMspotyspacing", gMspotyspacing, gNumMeasurements);
					writeMatIntArray(pmat,  "gMsigxoffset",   gMsigxoffset,   gNumMeasurements);
					writeMatIntArray(pmat,  "gMsigyoffset",   gMsigyoffset,   gNumMeasurements);
					writeMatIntArray(pmat,  "gMspotxoffset",  gMspotxoffset,  gNumMeasurements);
					writeMatIntArray(pmat,  "gMspotyoffset",  gMspotyoffset,  gNumMeasurements);
					writeMatDoubleArray(pmat, "gMlensxphase", gMlensxphase,   1, gNumMeasurements);
					writeMatDoubleArray(pmat, "gMlensyphase", gMlensyphase,   1, gNumMeasurements);
					writeMatIntArray(pmat,  "gMindices",      gMindices,      gNumMeasurements);
				
					// close the sequence file
					matClose(pmat);
				}
				
				// set the first measurement of the sequence
				int i = gMindices[gSpotSeqCounter];
				SLM_setLensPhase(gMlensxphase[i], gMlensyphase[i]);
				SLM_generateSpotArray(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, -1, -1,
	                       gMnumxspots[i], gMnumyspots[i], gMspotxspacing[i], gMspotyspacing[i], gMsigxoffset[i], 
						   gMsigyoffset[i], gMspotxoffset[i], gMspotyoffset[i], gMphaseconstraint, gMphasestep);
				
				// also update the control panel values
				SetCtrlVal(TabPage_1_2, TABPANEL_5_NumXSpots,    gMnumxspots[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_NumYSpots,    gMnumyspots[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXOffset,  gMspotxoffset[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYOffset,  gMspotyoffset[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_SigXOffset,   gMsigxoffset[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_SigYOffset,   gMsigyoffset[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXSpacing, gMspotxspacing[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYSpacing, gMspotyspacing[i]);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_PhaseConstraint, gMphaseconstraint);
				SetCtrlVal(TabPage_1_2, TABPANEL_5_PhaseStep,    gMphasestep);
			
				// update the SLM with the new pattern
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
				// record the framenumber at which we have reset the SLM
				gFrameSpotSeqSet = gFrameNumber;       
			
				// update the camera panel spot seq. text label
				char tmp[100];
				sprintf(tmp, "%i / %i", gSpotSeqCounter, gNumMeasurements);
				SetCtrlVal(pnlCamera, PANEL_CurrSpotSeq, tmp);

				// start the timer
				SetCtrlAttribute(pnlCamera, PANEL_SpotSeqTimer, ATTR_ENABLED, 1);
			}

			break;
	}
	return 0;
}

int CVICALLBACK SpotSeqTimer_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_TIMER_TICK:
			
			// check if the current frame is at least gNumFrames frames
			// after the SLM pattern was set
			if (gFrameNumber > gFrameSpotSeqSet + 2 * gNumFrames)
			{
				// yes, we can assume that the current camera frame has adjusted to the state of the SLM
				
				// generate a new filename
				char filename[1024];
				sprintf(filename, "%s%i.mat", gSpotSeqFilePrefix, gMindices[gSpotSeqCounter]);
				
				// store the frame
				CameraSaveMatFile(filename, 1);
				
				// advance to the next measurement
				gSpotSeqCounter++;
				
				// check if there are still measurements left in the sequence
				if (gSpotSeqCounter < gNumMeasurements)
				{
					// set the new SLM pattern
					int i = gMindices[gSpotSeqCounter];
					SLM_setLensPhase(gMlensxphase[i], gMlensyphase[i]);
					//SLM_generateSpotArray(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, -1, -1,
	                //       gMnumxspots[i], gMnumyspots[i], gMspotxspacing[i], gMspotyspacing[i], gMsigxoffset[i], 
					//	   gMsigyoffset[i], gMspotxoffset[i], gMspotyoffset[i], gMphaseconstraint, gMphasestep);
				
					// update the SLM with the new pattern
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
				
					// record the framenumber at which we reset the SLM
					gFrameSpotSeqSet = gFrameNumber;       
				
					// update the camera panel spot seq. text label
					char tmp[100];
					sprintf(tmp, "%i / %i", gSpotSeqCounter, gNumMeasurements);
					SetCtrlVal(pnlCamera, PANEL_CurrSpotSeq, tmp);
				}
				else
				{
					// no measurements left, finish the sequence
					gRunSpotSeq = 0;
					SetCtrlAttribute(pnlCamera, PANEL_SpotSeqTimer, ATTR_ENABLED, 0);
				}
					
			}

			break;
	}
	
	return 0;
}

/*
int CVICALLBACK TestAmplModCorr_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// set the number of measurements
			int NumMeasurements;
			GetCtrlVal(panel, PANEL_NumAMCMeasurements, &NumMeasurements);
			
			// get the filename prefix
			GetCtrlVal(pnlCamera, PANEL_AmplModFilePrefix, gAmplModFilePrefix);
			
			// perform the measurements in sequence
			for (int k = 0; k < NumMeasurements; k++)
			{
				// update the camera panel bias text label
				char tmp[100];
				sprintf(tmp, "%i / %i", k, NumMeasurements);
				SetCtrlVal(pnlCamera, PANEL_CurrAmplModBias, tmp);
				
				// set the current amplitude modulation
				double am = 2.0 * ((double) k) / ((double) (NumMeasurements - 1)) - 0.5;
				SLM_setAmplitudeModulation(am);
				
				// recompute the current beam shaping pattern using the current control panel settings
				PicPattern_Callback (pnlControl, TABPANEL_8_PicPatternGo, EVENT_COMMIT,
					NULL, SLM_NO_UPDATE, 0);
				
				// update the SLM with the new pattern
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
				// wait a little bit for the SLM pixels to adjust to their new setting
				Delay(1.0);
				
				// call the relevant part of the camera timer callback function manually, since the camera 
				// timer events are not processed when this routine is running
				for (int c = 0; c < getCameraNumAveragingFrames(); c++)
				{
					// first wait one frame
					Delay(1.0 / getCameraFrameRate());
					
					// grab a frame
					grabCameraFrame();
				}
				
				// average the camera frames
				performCameraAveraging();
				
				// update the canvas too (we need this to save a .png screenshot)
				CameraTimer_Callback (panel, PANEL_CameraTimer, EVENT_TIMER_TICK, NULL, 0, 0);
			
				// generate a new filename
				char filename[1024];
				sprintf(filename, "%s%.3f.mat", gAmplModFilePrefix, am);
				
				// store the frame
				CameraSaveMatFile(filename, 1);
				
				// generate a bitmap from the canvas			
				int TmpBitmap = -1;
				GetCtrlDisplayBitmap (panel, PANEL_CameraCanvas, 0, &TmpBitmap);
				
				// save it to file
				sprintf(filename, "%s%.3f.png", gAmplModFilePrefix, am);
				SaveBitmapToPNGFile (TmpBitmap, filename);
				
				// discard the bitmap 
				DiscardBitmap(TmpBitmap);
			}

			break;
	}
	return 0;
}			*/
