//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Routines for determining the phase aberration with a modified IFTA.
//
// Created on:  02-04-2012 by Rick.
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
#include "SLM_Zernike.h"


// the maximum number of Zernike polynomials to include
static int NPmax = 28;

// the current number of Zernike polynomials
static int gNZP = 28;

// the Zernike coefficients
static double* gCj;

// should we preallocate the Zernike polynomials?
static int gPreallocZernike;


void SaveBatchToFile(char filename[1024]);

/// HIFN Initialise the Zernike library: grab a camera frame, copy the signal, 
///      and precalculate the Zernike polynomials (if required)
int CVICALLBACK InitPhaseError_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// get the values of the supersampling ratio
			int ssrnumerator, ssrdenominator;
			GetCtrlVal(panel, TABPANE_11_SSRDenominator, &ssrdenominator);
			GetCtrlVal(panel, TABPANE_11_SSRNumerator,   &ssrnumerator);
			
			// compute the number of sampling points
			int NumSamplingPoints = SLM_getXres() * ssrnumerator / ssrdenominator;
			
			// compute the resolutions
			// NOTE: the camera pixels are square, so xres == yres
			int xres = SLM_getXres() * ssrnumerator / ssrdenominator;
			int yres = SLM_getXres() * ssrnumerator / ssrdenominator;
			
			// should we preallocate the Zernike polynomials?
			GetCtrlVal(panel, TABPANE_11_PreallocateZernike, &gPreallocZernike);
			
			// read out the cutoff level for noise			  
			double NoiseCutoff;
			GetCtrlVal(panel, TABPANE_11_PhaseErrorNoiseCutoff, &NoiseCutoff);
			
			// read out the number of Zernike polynomials to include
			int NZP;
			GetCtrlVal(panel, TABPANE_11_PhaseErrorNumZernike, &NZP);
				
			// initialise the phase error characterisation part of the SLM library
			SLM_initPhaseErrorCharacterisation(xres, yres, NZP, NoiseCutoff, gPreallocZernike, 
				pnlDebug, DebugPanel_CANVAS,   DebugPanel_CANVAS_2, DebugPanel_CANVAS_3, 
						  DebugPanel_CANVAS_4, DebugPanel_CANVAS_5, DebugPanel_CANVAS_6);
			
			// allocate and clear the array of Zernike polynomials, if it wasn't allocated already
			if (gCj == NULL)
				gCj = (double*) calloc(NPmax, sizeof(double));

			break;
		}
	}
	return 0;
}


/// HIFN start the phase error calculation with the current settings
int CVICALLBACK PhaseErrorGoCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// read out the number of iterations
			int Nit;
			GetCtrlVal(panel, TABPANE_11_PhaseErrorIterations, &Nit);
			
			// read out the number of Zernike polynomials to include
			int NZP;
			GetCtrlVal(panel, TABPANE_11_PhaseErrorNumZernike, &NZP);
			
			// read out the Zernike projection interval
			int prji;
			GetCtrlVal(panel, TABPANE_11_PhaseErrorPrjInterval, &prji);
			
			// read out the adjustment factor beta
			double beta;				
			GetCtrlVal(panel, TABPANE_11_PhaseErrorBeta, &beta);
			
			// start the IFTA iterations
			double refampmatch;
			double* tmpcj = SLM_IFTA_PE(gCj, Nit, 1, NZP, prji, beta, &refampmatch, pnlDebug,  DebugPanel_CANVAS,   DebugPanel_CANVAS_2, DebugPanel_CANVAS_3, 
											    	DebugPanel_CANVAS_4, DebugPanel_CANVAS_5, DebugPanel_CANVAS_6, DebugPanel_DBGGRAPH, DebugPanel_DBGGRAPH2, DebugPanel_DBGGRAPH3);
			
			// copy the values of the cj
			for (int j = 1; j <= NZP; j++)
				gCj[j - 1] = tmpcj[j - 1];
			
			// clean up
			free(tmpcj);

			break;
		}
	}
	return 0;
}


/// HIFN update the memory requirement and pixel dimensions for the supersampling settings
int CVICALLBACK SSRUpdate (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// get the values of the supersampling ratio
			int numerator, denominator;
			GetCtrlVal(panel, TABPANE_11_SSRDenominator, &denominator);
			GetCtrlVal(panel, TABPANE_11_SSRNumerator,   &numerator);
			
			// update the text label displaying the number of sampling points 
			// and the memory requirements for a single FFT structure
			// NOTE: the sampling points form a SQUARE area, so xres == yres
			char txtSP[1024];
			int xres = SLM_getXres() * numerator / denominator;
			int yres = SLM_getXres() * numerator / denominator;
			sprintf(txtSP, "%i x %i (%i mb)", xres, yres, xres * yres * 8 * 2 / (1024 * 1024));
			SetCtrlVal(panel, TABPANE_11_SamplingPointsValue, txtSP);
			
			// also display the number of sampling points needed to exactly sample the camera pixels
			char txtCP[1024];
			int cxres = SLM_getLensFocalLength() * SLM_getWavelength() / (4.65e-6 * SLM_getLxSLM()) * SLM_getXres();
			int cyres = SLM_getLensFocalLength() * SLM_getWavelength() / (4.65e-6 * SLM_getLySLM()) * SLM_getYres();
			sprintf(txtCP, "%i x %i (%i mb)", cxres, cyres, cxres * cyres * 8 * 2 / (1024 * 1024));
			SetCtrlVal(panel, TABPANE_11_SamplingPointsValCam, txtCP);
			
			break;
		}
	}
	return 0;
}


/// HIFN reset the current Zernike coefficients to some random initial guess
int CVICALLBACK RandomizeCJ_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			for (int k = 0; k < NPmax; k++)
				gCj[k] = (Random(0.0, 2.0) - 1.0);

			break;
	}
	return 0;
}


/// HIFN reset the current Zernike coefficients to zero
int CVICALLBACK ClearCJ_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			for (int k = 0; k < NPmax; k++)
				gCj[k] = 0.0;

			break;
	}
	return 0;
}


/// HIFN Apply the current phase error correction pattern that was found through iteration
int CVICALLBACK PhaseErrorApplyCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// read out the number of Zernike polynomials to include
			int NZP;
			GetCtrlVal(panel, TABPANE_11_PhaseErrorNumZernike, &NZP);
			
			// check if we should mirror the aberration correction
			int mx, my;
			GetCtrlVal(panel, TABPANE_11_AbbCorrMirrorX, &mx);
			GetCtrlVal(panel, TABPANE_11_AbbCorrMirrorY, &my);
			
			// by how much should the Correction be multiplied?
			double zcf;
			GetCtrlVal(panel, TABPANE_11_ZerCorrFactor, &zcf);
			
			// create a copy of the zernike coefficients, multiplied by the specified factor
			double* tmp = (double*) malloc(NZP * sizeof(double));
			for (int k = 0; k < NZP; k++)
				tmp[k] = zcf * gCj[k];
			
			// set the tilt and piston to zero
			tmp[0] = 0.0;
			tmp[1] = 0.0;
			tmp[2] = 0.0;
			
			// set the aberration correction
			SLM_setAberrationCorrection(tmp, NZP, mx, my);
			
			// clean up
			free(tmp);
			
			// update the SLM
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
			break;
		}
	}
	return 0;
}


/// HIFN clear the current phase error correction pattern
int CVICALLBACK PhaseErrorClearCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// set the aberration correction
			SLM_setAberrationCorrection(NULL, 0, 0, 0);
			
			// update the SLM
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

			break;
		}
	}
	return 0;
}


/// HIFN Load a single set of Zernike coefficients from file
int CVICALLBACK LoadCJ_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{												   
	switch (event)
	{
		case EVENT_COMMIT:
			
			// let the user specify a filename for the batch
			char filename[1024];
		    int SelectionStatus = FileSelectPopup("", "*.mat", "*.mat", "Load Batch", 
									VAL_OK_BUTTON, 0, 1, 1, 1, filename);
			
			if (SelectionStatus > 0)
			{
				// create a Matlab file pointer
				MATFile *pmat;

				// open the .mat file
				pmat = matOpen(filename, "r");
				
				// the array dimensions
				int N, M;
				
				// read the coefficients from the file
				double* tmp = readMatDoubleArray(pmat, "Cj", &M, &N);
				
				// set the current number of Zernike polynomials
				int NZP = (M > N ? M : N);
				SetCtrlVal(panel, TABPANE_11_PhaseErrorNumZernike, NZP);
				
				// make sure the Zernike library is initialised to the correct number
				InitPhaseError_Callback (panel, 0, EVENT_COMMIT, NULL, 0, 0);
				
				// copy the Zernike coefficients from the temporary array
				memcpy(gCj, tmp, NZP * sizeof(double));
				
				// clean up the temporary array
				free(tmp);
				
				// close the mat-file
				matClose(pmat);
			}

			break;
	}
	return 0;
}



// indicator variable whether a batch is currently running
static int gBatchRunning = 0;

// current calculation number within batch
static int gBatchCurrCalculation = 0;

// the total number of calculations in the current batch
static int gBatchNumCalculations = 0;

// the coefficients calculated in a batch (1 row == 1 calculation of Zernike coeffs)
static double* gBatchZernikeCoefficients;

// the evaluation score of the coefficients, each row contains 3 evaluation criteria, 1 row per calculation
static double* gBatchEvaluationScore;

// the number of Zernike coefficients to be calculated in each calculation of the current batch
static int gBatchNumZernikeCoeffs;

// the number of iterations per set of Zernike coefficients in the current batch
static int gBatchNumIterations;

// the adjustment factor beta of the current batch
static double gBatchBeta;

// the noise cutoff that was used for the reference amplitude
static double gBatchNoiseCutoff;

static char gBatchFilename[1024];

/// HIFN start / stop a batch calculation
int CVICALLBACK PhaseBatchGoCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			if (!gBatchRunning)
			{
				// let the user specify a filename for the batch
			    int SelectionStatus = FileSelectPopup("", "*.zcj", "*.zcj", "Save Batch", 
										VAL_OK_BUTTON, 0, 1, 1, 1, gBatchFilename);
			
				if (SelectionStatus > 0)
				{
					// start a new batch
			
					// get the number of Zernike coefficients to calculate per batch
					GetCtrlVal(panel, TABPANE_11_PhaseErrorNumZernike, &gBatchNumZernikeCoeffs);
				
					// get the number of calculations to be performed
					GetCtrlVal(panel, TABPANE_11_PhaseErrorBatchNumIt, &gBatchNumCalculations);
				
					// get the number of iterations to be performed per set of Zernike coeffs
					GetCtrlVal(panel, TABPANE_11_PhaseErrorIterations, &gBatchNumIterations);
				
					// get the beta adjustment coefficient
					GetCtrlVal(panel, TABPANE_11_PhaseErrorBeta, &gBatchBeta);
				
					// get the noise cutoff
					GetCtrlVal(panel, TABPANE_11_PhaseErrorNoiseCutoff, &gBatchNoiseCutoff);
				
					// was there already memory allocated?
					if (gBatchZernikeCoefficients != NULL)
					{
						// yes, free it first
						free(gBatchZernikeCoefficients);
						free(gBatchEvaluationScore);
					}
				
					// allocate new memory
					gBatchZernikeCoefficients = (double*) calloc(gBatchNumCalculations * gBatchNumZernikeCoeffs, sizeof(double));
					gBatchEvaluationScore     = (double*) calloc(gBatchNumCalculations * 3, sizeof(double));
				
					// reset counting variables
					gBatchCurrCalculation = 0;
			
					// start timer
					SetCtrlAttribute(panel, TABPANE_11_PhaseErrorBatchTimer, ATTR_ENABLED, 1);
				
					// indicate that the batch is running
					gBatchRunning = 1;
				}
			}
				
			break;
	}
	return 0;
}


/// HIFN perform one calculation of a batch
int CVICALLBACK PhaseErrorBatchTimer_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_TIMER_TICK:
		{	
			// update progress counter label
			char txtBP[1024];
			sprintf(txtBP, "Current Calculation: %i / %i", gBatchCurrCalculation + 1, gBatchNumCalculations);
			SetCtrlVal(panel, TABPANE_11_BatchProgress, txtBP);
			
			// randomise coefficients
			for (int k = 0; k < NPmax; k++)
				gCj[k] = 5 * 4 * (Random(0.0, 2.0) - 1.0) / (k + 4);
			
			// start iteration
			double* tmpcj = SLM_IFTA_PE(gCj, gBatchNumIterations, 1, gBatchNumZernikeCoeffs, 1, gBatchBeta, &gBatchEvaluationScore[gBatchCurrCalculation * 3 + 0], 
								pnlDebug,  DebugPanel_CANVAS,   DebugPanel_CANVAS_2, DebugPanel_CANVAS_3, DebugPanel_CANVAS_4, 
								DebugPanel_CANVAS_5, DebugPanel_CANVAS_6, DebugPanel_DBGGRAPH, DebugPanel_DBGGRAPH2, DebugPanel_DBGGRAPH3);
			
			// copy the values of the coefficients to the batch structure
			for (int j = 1; j <= gBatchNumZernikeCoeffs; j++)
				gBatchZernikeCoefficients[gBatchCurrCalculation * gBatchNumZernikeCoeffs + j - 1] = tmpcj[j - 1];
			
			// clean up
			free(tmpcj);
			
			// DEBUG plot the current batch structure
			plotField(gBatchZernikeCoefficients, gBatchNumZernikeCoeffs, gBatchNumCalculations, pnlDebug, DebugPanel_CANVAS_6);
			
			// save the latest result to file
			char tmpfilename[1024];
			sprintf(tmpfilename, "%s%i", gBatchFilename, gBatchCurrCalculation);
			SaveBatchToFile(tmpfilename);
			
			// update counter
			gBatchCurrCalculation++;
			
			// are we finished?
			if (gBatchCurrCalculation == gBatchNumCalculations)
			{
				// yes, we can stop the timer
				SetCtrlAttribute(panel, TABPANE_11_PhaseErrorBatchTimer, ATTR_ENABLED, 0);
			
				// indicate that no batch is running
				gBatchRunning = 0;
				
				// change the label of the button
				SetCtrlAttribute(panel, TABPANE_11_PhaseErrorBatchGo, ATTR_LABEL_TEXT, "__Start Batch");
			}

			break;
		}
	}
	return 0;
}


/// HIFN save the results of a batch calculation
int CVICALLBACK PhaseErrorBatchSaveCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// let the user specify a filename for the bitmap
			char filename[1024];
		    int SelectionStatus = FileSelectPopup("", "*.zcj", "*.zcj", "Save Batch", 
									VAL_OK_BUTTON, 0, 1, 1, 1, filename);
			
			if (SelectionStatus > 0)
			{
				SaveBatchToFile(filename);
			}

			break;
		}
	}
	return 0;
}


void SaveBatchToFile(char filename[1024])
{
	// create a Matlab file pointer
	MATFile *pmat;

	// open the .mat file
	pmat = matOpen(filename, "wz");
	
	// write the batch settings
	writeMatDoubleScalar(pmat, "batch_CurrCalc", (double) gBatchCurrCalculation);			
	writeMatDoubleScalar(pmat, "batch_NumCalc", (double) gBatchNumCalculations);
	writeMatDoubleScalar(pmat, "batch_NumZernikeCoeffs", (double) gBatchNumZernikeCoeffs);
	writeMatDoubleScalar(pmat, "batch_NumIterations", (double) gBatchNumIterations);
	writeMatDoubleScalar(pmat, "batch_Beta", (double) gBatchBeta);
	writeMatDoubleScalar(pmat, "batch_NoiseCutoff", (double) gBatchNoiseCutoff);
	
	// retrieve the reference amplitude
	double* refampl;
	int wx, wy, wxpos, wypos;
	SLM_getCharacterisationReferenceAmplitude(&refampl, &wx, &wy, &wxpos, &wypos);
	
	// write the reference amplitude parameters
	writeMatDoubleScalar(pmat, "batch_RefAmpWx", (double) wx);
	writeMatDoubleScalar(pmat, "batch_RefAmpWy", (double) wy);
	writeMatDoubleScalar(pmat, "batch_RefAmpWxPos", (double) wxpos);
	writeMatDoubleScalar(pmat, "batch_RefAmpWxPos", (double) wypos);
	
	// create an array for the reference amplitude
	mxArray *mxTmpRA;
	mxTmpRA = mxCreateNumericMatrix(wx, wy, mxDOUBLE_CLASS, mxREAL);

	// copy the coefficients to the array
	memcpy((void*) mxGetData(mxTmpRA), (void*) refampl, wx * wy * sizeof(double));			

	// write the data to the mat file
	matPutVariable(pmat, "batch_ReferenceAmplitude", mxTmpRA);

	// destroy the mxArray
	mxDestroyArray(mxTmpRA);
	
	// create an array for the coefficients
	mxArray *mxTmp;
	mxTmp = mxCreateNumericMatrix(gBatchNumZernikeCoeffs, gBatchNumCalculations, mxDOUBLE_CLASS, mxREAL);

	// copy the coefficients to the array
	memcpy((void*) mxGetData(mxTmp), (void*) gBatchZernikeCoefficients, gBatchNumCalculations * gBatchNumZernikeCoeffs * sizeof(double));			

	// write the data to the mat file
	matPutVariable(pmat, "batch_ZernikeCoefficients", mxTmp);

	// destroy the mxArray
	mxDestroyArray(mxTmp);
	
	// create an array for the evaluation scores
	mxArray *mxTmpES;
	mxTmpES = mxCreateNumericMatrix(3, gBatchNumCalculations, mxDOUBLE_CLASS, mxREAL);

	// copy the coefficients to the array
	memcpy((void*) mxGetData(mxTmpES), (void*) gBatchEvaluationScore, 3 * gBatchNumCalculations * sizeof(double));			

	// write the data to the mat file
	matPutVariable(pmat, "batch_EvaluationScore", mxTmpES);

	// destroy the mxArray
	mxDestroyArray(mxTmpES);
	
	// close the mat-file
	matClose(pmat);
}


/// HIFN load the results of a batch calculation
int CVICALLBACK PhaseErrorBatchLoadCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// let the user specify a filename for the batch
			char filename[1024];
		    int SelectionStatus = FileSelectPopup("", "*.zcj", "*.zcj", "Load Batch", 
									VAL_OK_BUTTON, 0, 1, 1, 1, filename);
			
			if (SelectionStatus > 0)
			{
				// was there already memory allocated for the Zernike coefficients?
				if (gBatchZernikeCoefficients != NULL)
				{
					// yes, free it first
					free(gBatchZernikeCoefficients);
				}
				
				// create a Matlab file pointer
				MATFile *pmat;

				// open the .mat file
				pmat = matOpen(filename, "r");
				
				// read the batch settings
				gBatchCurrCalculation  = readMatDoubleScalar(pmat, "batch_CurrCalc");			
				gBatchNumCalculations  = readMatDoubleScalar(pmat, "batch_NumCalc");
				gBatchNumZernikeCoeffs = readMatDoubleScalar(pmat, "batch_NumZernikeCoeffs");
				gBatchNumIterations    = readMatDoubleScalar(pmat, "batch_NumIterations");
				gBatchBeta			   = readMatDoubleScalar(pmat, "batch_Beta");
				
				// TODO: resolutions etc.
				//readMatDoubleScalar(pmat, "batch_");
				
				
				// allocate new memory for the Zernike coefficients
				gBatchZernikeCoefficients = (double*) malloc(gBatchNumCalculations * gBatchNumZernikeCoeffs * sizeof(double));
				
				// create an array for reading out the coefficients
				mxArray* mxTmp = matGetVariable(pmat, "batch_ZernikeCoefficients");
				
				// copy the data to the array of Zernike coefficients
				memcpy(gBatchZernikeCoefficients, mxGetPr(mxTmp), gBatchNumZernikeCoeffs * gBatchNumCalculations * sizeof(double));
					
				// clean up temporary array
				mxDestroyArray(mxTmp);
				
				// create an array for reading out the evaluation score
				mxArray* mxTmpES = matGetVariable(pmat, "batch_EvaluationScore");
				
				// allocate new memory for the evaluation score
				gBatchEvaluationScore = (double*) realloc(gBatchEvaluationScore, 3 * gBatchNumCalculations * sizeof(double));
				
				// copy the data to the array of evaluation scores
				memcpy(gBatchEvaluationScore, mxGetPr(mxTmpES), 3 * gBatchNumCalculations * sizeof(double));
					
				// clean up temporary array
				mxDestroyArray(mxTmpES);
				
				// close the mat-file
				matClose(pmat);
				
				// DEBUG plot the current batch structure
				plotField(gBatchZernikeCoefficients, gBatchNumZernikeCoeffs, gBatchNumCalculations, pnlDebug, DebugPanel_CANVAS_6);
			}

			break;
		}
	}
	return 0;
}

int CVICALLBACK PhaseErrorBatchEvalCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			break;
	}
	return 0;
}


/// HIFN selects a calculation from the batch and loads the Zernike coefficients
int CVICALLBACK PhaseErrorBatchSelect_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// which calculation is selected?
			int selected;
			GetCtrlVal(panel, control, &selected);
			
			// display the 'score' 
			char txt[1024];
			sprintf(txt, "Score: %f", gBatchEvaluationScore[(selected - 1) * 3 + 0]);
			SetCtrlVal(panel, TABPANE_11_BatchScore, txt);
			
			// clear the current coefficients
			memset(gCj, 0, NPmax * sizeof(double));
			
			// copy the Zernike coefficients
			for (int j = 1; j <= gBatchNumZernikeCoeffs; j++)
				gCj[j - 1] = gBatchZernikeCoefficients[gBatchNumZernikeCoeffs * (selected - 1) + (j - 1)];

			break;
		}
	}
	return 0;
}


/// HIFN apply the Zernike coefficients of the currently selected batch calculation
int CVICALLBACK PhaseErrorBatchApplyCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// check if we should mirror the aberration correction
			int mx, my;
			GetCtrlVal(panel, TABPANE_11_AbbCorrMirrorX, &mx);
			GetCtrlVal(panel, TABPANE_11_AbbCorrMirrorY, &my);
			
			// set the aberration correction
			SLM_setAberrationCorrection(gCj, gBatchNumZernikeCoeffs, mx, my);
			
			// update the SLM
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

			break;
	}
	return 0;
}

int CVICALLBACK PhaseBatchStopCallback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// stop the timer
			SetCtrlAttribute(panel, TABPANE_11_PhaseErrorBatchTimer, ATTR_ENABLED, 0);
		
			// indicate that no batch is running
			gBatchRunning = 0;
			
			break;
	}
	return 0;
}

int CVICALLBACK ManualSlide_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// get the j-index
			int j;
			GetCtrlVal(panel, TABPANE_11_ManualCj, &j);
			
			// get the value of cj
			double cj;
			GetCtrlVal(panel, control, &cj);
			
			// set the value of cj in the array of Zernike coefficients
			gCj[j - 1] = cj;
			
			// check if we should mirror the aberration correction
			int mx, my;
			GetCtrlVal(panel, TABPANE_11_AbbCorrMirrorX, &mx);
			GetCtrlVal(panel, TABPANE_11_AbbCorrMirrorY, &my);
			
			// apply the coefficients
			SLM_setAberrationCorrection(gCj, gNZP, mx, my);
			
			// update the SLM
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

			break;
		}
	}
	return 0;
}

int CVICALLBACK ManualCj_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{	
			// get the current j-index
			int j;
			GetCtrlVal(panel, control, &j);
			
			// set the sliderbar to the current value of the j-th coefficient
			SetCtrlVal(panel, TABPANE_11_ManualSlide, gCj[j - 1]);

			break;
		}
	}
	return 0;
}


