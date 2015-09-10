//==============================================================================
//
// Title:       SLM_characterisation.c
// Purpose:     Library of functions that are used for characterising the input
//              intensity incident onto the SLM, and the total phase error 
//              induced by the optical system, via a modified IFTA.
//
// Created on:  23-3-2012 at 13:31:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//============================================================================== 

#include "SLM_internal.h"
#include "SLM_Zernike.h"
#include "SLM Camera Panel.h"

// arrays for describing the light field in the SLM and Focal plane, for phase error characterisation
fftw_complex *gFFTpeSLM, *gFFTpeFocal;

// FFT plans for propagation between the two planes (SLM and Focal)
fftw_plan gFFTpeSF, gFFTpeFS;

// the number of sampling points in the SLM and focal plane for phase error characterisation
static int gNxPE, gNyPE;

// the window size and position in the focal plane
static int gWxPE, gWyPE, gWXsizePE, gWYsizePE;

// the window size in camera units
static int gWXsizeCam, gWYsizeCam;

// the SLM coordinates in the SLM plane
static int gSLMxPE, gSLMyPE;

// the maximum number of Zernike polynomials to include
static int gNPmax;

// the SLM base phase pattern, without phase errors
static double* gSLMBasePhase;

// the measured amplitude, to be used as constraint in the IFTA
static double* gReferenceAmplitude;

// the total intensity contained in the SLM plane (and thus the total
// intensity in the focal plane too)
static double gTotalIntensity;

// cutoff factor for noise: any amplitude below this value is regarded as noise
static double gNoiseCutoff;


/// HIFN performs a Aitken transformation of a sequence, accelerating its convergence
double* AitkenTransform(double *seq, int N);

/// HIFN Initialise the data structures for the phase error characterisation
void SLM_initPhaseErrorCharacterisation(int NxPE, int NyPE, int NZP, double NoiseCutoff, int prealloc, 
	int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6)
{
	// should we precalculate the Zernike polynomials?
	if (prealloc)
		PreCalculateZernikePolynomials(gXsize, gYsize, NZP);
	
	// store the noise cutoff level
	gNoiseCutoff = NoiseCutoff;
	
	// copy the current base SLM phase pattern
	gSLMBasePhase = SLM_getPhase();
	
	// unwrap the base phase pattern
	SLM_unwrapPhase(gSLMBasePhase, gXsize, gYsize);
	
	// get the window size and position in focal units
	int wxf, wyf, wxsizef, wysizef;
	SLM_getSignalWindowParameters(&wxf, &wyf, &wxsizef, &wysizef);
	
	// get the window size in camera pixels
	int ulxwc = TransformFocalPlaneToBitmapX(wxf);
	int ulywc = TransformFocalPlaneToBitmapY(wyf);
	int lrxwc = TransformFocalPlaneToBitmapX(wxf + wxsizef);
	int lrywc = TransformFocalPlaneToBitmapY(wyf + wysizef);
	gWXsizeCam = lrxwc - ulxwc;
	gWYsizeCam = lrywc - ulywc;
	
	// set the number of sampling points
	gNxPE = NxPE;
	gNyPE = NyPE;
	
	// calculate the supersampling factor, i.e. how many sampling points fall within one focal unit?
	double SSFx = ((double) gNxPE) / ((double) SLM_getXres());
	double SSFy = ((double) gNyPE) / ((double) SLM_getYres());
	
	// set the position and number of sampling points of the window within the focal plane
	gWxPE = SSFx * wxf;
	gWyPE = SSFy * wyf;
	gWXsizePE = SSFx * wxsizef;
	gWYsizePE = SSFy * wysizef;
	
	// check if the FFT arrays were already allocated, if so, we have to free them first
	if (gFFTpeSLM != NULL)
	{
		fftw_free(gFFTpeSLM);
		fftw_free(gFFTpeFocal);
		fftw_destroy_plan(gFFTpeSF);
		fftw_destroy_plan(gFFTpeFS);
	}
	
	// initialise the FFT arrays
	gFFTpeSLM   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gNxPE * gNyPE);
	gFFTpeFocal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gNxPE * gNyPE);
	
	// create the FFT plans
	gFFTpeSF = fftw_plan_dft_2d(gNyPE, gNxPE, gFFTpeSLM, gFFTpeFocal, FFTW_FORWARD,  FFTW_ESTIMATE);
	gFFTpeFS = fftw_plan_dft_2d(gNyPE, gNxPE, gFFTpeFocal, gFFTpeSLM, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	// determine the coordinates of the SLM in the SLM plane (the FFT array describing the SLM plane)
	gSLMxPE = 0;
	gSLMyPE = 0;
	
	// grab the window 1:1 from the camera
	double* tmpwindow = getSignalWindowFromCameraResampled(NULL, wxf, wyf, wxsizef, wysizef, 0, 0, gWXsizeCam, gWYsizeCam);
	
	plotField(tmpwindow, gWXsizeCam, gWYsizeCam, debugpanel, dbc1);
	
	double tmpwinmax = 0.0;
	for (int k = 0; k < gWXsizeCam * gWYsizeCam; k++)
		tmpwinmax = (tmpwindow[k] > tmpwinmax ? tmpwindow[k] : tmpwinmax);
	
	// convert the intensity to an amplitude
	for (int k = 0; k < gWXsizeCam * gWYsizeCam; k++)
		tmpwindow[k] = (sqrt(tmpwindow[k]) > gNoiseCutoff * sqrt(tmpwinmax) ? sqrt(tmpwindow[k]) : 0.0);
	
	// resample the camera window to PE units
	gReferenceAmplitude = SLM_resampleBitmap(tmpwindow, gWXsizeCam, gWYsizeCam, gWXsizePE, gWYsizePE);
	free(tmpwindow);
	
	// compute the total intensity contained in the signal
	gTotalIntensity = 0.0;
	for (int k = 0; k < gWXsizePE * gWYsizePE; k++)
		gTotalIntensity += gReferenceAmplitude[k] * gReferenceAmplitude[k];
	
	plotField(gReferenceAmplitude, gWXsizePE, gWYsizePE, debugpanel, dbc2);
	
	// normalise the amplitude to that contained in the SLM window (integrate over the entire (simulated) focal plane??)
	// ################# TODO!
}


/// HIFN Set the value (amplitude + phase) inside a specified window
void SetValuesInWindow(fftw_complex* U, int NxU, int NyU, double* amplitude, double* phase,
	                       int xoffset, int yoffset, int Wx, int Wy, int ClearComplement, double amplc)
{
	// check if we want to clear the area outside the window
	if (ClearComplement)
	{
		for (int k = 0; k < NxU; k++)
		for (int l = 0; l < NyU; l++)
		{
			// check if we are in the window area 
			if (   (k >= xoffset) && (k < xoffset + Wx)
				&& (l >= yoffset) && (l < yoffset + Wy) )
			{
				// calculate the index into the window
				int winindex = (k - xoffset) + (l - yoffset) * Wx;
				
				// yes. set the value of the light field to the desired amplitude and phase
				U[k + l * NxU][0] = amplc * amplitude[winindex] * cos(phase[winindex]);
				U[k + l * NxU][1] = amplc * amplitude[winindex] * sin(phase[winindex]);
			}
			else
			{
				// no, set to zero: there is no light here
				gFFTpeSLM[k + l * NxU][0] = 0.0;
				gFFTpeSLM[k + l * NxU][1] = 0.0;
			}
		}
	}
	else
	{
		// no, we can iterate only over the window points
		for (int k = 0; k < Wx; k++)
		for (int l = 0; l < Wy; l++)
		{
			// set amplitude and phase
			U[k + xoffset + (l + yoffset) * NxU][0] = amplc * amplitude[k + l * Wx] * cos(phase[k + l * Wx]);
			U[k + xoffset + (l + yoffset) * NxU][1] = amplc * amplitude[k + l * Wx] * sin(phase[k + l * Wx]);
		}
	}
	
}


/// HIFN Obtain the phase of the complex valued field in a specified window,
///      the output array 'phase' should be pre-allocated by the user.
void GetPhaseFromWindow(fftw_complex* U, int NxU, int NyU, double* phase, int xoffset, int yoffset, int Wx, int Wy)
{
	for (int k = 0; k < Wx; k++)
	for (int l = 0; l < Wy; l++)
	{
		phase[k + l * Wx] = atan2(U[k + xoffset + (l + yoffset) * NxU][1], 
			                      U[k + xoffset + (l + yoffset) * NxU][0]);
		if (phase[k + l * Wx] < 0.0)
			phase[k + l * Wx] += 2 * PI;
	}
}


/// HIFN Obtain the amplitude of the complex valued field in a specified window,
///      the output array 'ampl' should be pre-allocated by the user.
void GetAmplitudeFromWindow(fftw_complex* U, int NxU, int NyU, double* ampl, int xoffset, int yoffset, int Wx, int Wy)
{
	for (int k = 0; k < Wx; k++)
	for (int l = 0; l < Wy; l++)
	{
		ampl[k + l * Wx] = sqrt(U[k + xoffset + (l + yoffset) * NxU][0] * U[k + xoffset + (l + yoffset) * NxU][0] + 
							U[k + xoffset + (l + yoffset) * NxU][1] * U[k + xoffset + (l + yoffset) * NxU][1]);
	}
}


/// HIFN Performs 'N' IFTA iterations, where the projection step in the SLM plane consists of 
///      projecting onto a subspace of a fixed SLM pattern + a linear combination of the 
///      lowest few Zernike rectangle polynomials.
///      Function returns the coefficients of the Zernike polynomials.
double* SLM_IFTA_PE(double* cj0, int N, int UseSoftOp, int NPmax, int prjinterval, double beta, double* RefAmpMatch,
	int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6, int debuggraph, int debuggraph2, int debuggraph3)
{
	// store the number of Zernike polynomials
	gNPmax = NPmax;
	
	// arrays to keep track of convergence
	double* convergence = (double*) malloc((N + 1) * sizeof(double));
	double* zdists      = (double*) malloc((N + 0) * sizeof(double));
	
	// initialise array of array of Zernike coefficients
	double** cj = (double**) malloc((N + 1) * sizeof(double*));
	for (int n = 0; n <= N; n++)
		cj[n] = (double*) malloc(NPmax * sizeof(double));

	// are there initial values provided for the Zernike coefficients?
	if (cj0 != NULL)
		memcpy(cj[0], cj0, gNPmax * sizeof(double));
	else
	{
		// no, initialise with random values
		double ampl = 0.5;
		for (int j = 1; j <= gNPmax; j++)
			cj[0][j - 1] = ampl * (Random(0.0, 2.0) - 1.0);
		
		// we manually set the coefficients associated with translation to zero
		cj[0][1] = 0.0;
		cj[0][2] = 0.0;
	}
	
	// create an array for holding the projected SLM phase
	double* prjphase = (double*) malloc(gXsize * gYsize * sizeof(double));
	
	// create an array for temporarily holding the phase in the SLM plane while operating on it
	double* slmphase = (double*) malloc(gXsize * gYsize * sizeof(double));
	
	// create an array for temporarily holding the phase in the focal plane window (for applying the amplitude constraints)
	double* focalphase = (double*) malloc(gWXsizePE * gWYsizePE * sizeof(double));
	
	// create an array for temporarily holding the amplitude in the focal plane window (debug purposes!)
	double* focalamp   = (double*) malloc(gWXsizePE * gWYsizePE * sizeof(double));
	
	// initialise the phase array with the base phase + the random vector in Zernike space
	memcpy(prjphase, gSLMBasePhase, gXsize * gYsize * sizeof(double));
	//plotField(prjphase, gXsize, gYsize, debugpanel, dbc2);
	for (int j = 1; j <= gNPmax; j++)
	{
		// sample the j-th Zernike polynomial
		double* tmp = SampleZernikeJ(j, gXsize, gYsize);

		// add the j-th Zernike polynomial to the phase pattern
		for (int k = 0; k < gXsize * gYsize; k++)
			prjphase[k] += cj[0][j - 1] * tmp[k];	
		
		// free the temporary array
		free(tmp);
	}
	plotField(prjphase, gXsize, gYsize, debugpanel, dbc3);
	
	// initialise the FFT array of the SLM plane
	SetValuesInWindow(gFFTpeSLM, gNxPE, gNyPE, gInputAmplitude, prjphase, gSLMxPE, gSLMyPE, gXsize, gYsize, 1, 1.0);
	
	plotFFTFieldPhase(gFFTpeSLM, gNxPE, gNyPE, debugpanel, dbc5);
	DeleteGraphPlot (debugpanel, debuggraph, -1, VAL_IMMEDIATE_DRAW); 
	
	// perform the IFTA iteration
	for (int n = 1; n <= N; n++)
	{
		// propagate from the SLM plane to the focal plane
		fftw_execute(gFFTpeSF);
		
		// compute amplitude correction factor
		double amplc = 0.0;
		for (int k = 0; k < gWXsizePE; k++)
		for (int l = 0; l < gWYsizePE; l++)
		{
			int index = k + gWxPE + (l + gWyPE) * gNxPE;
			amplc += sqrt(gFFTpeFocal[index][0] * gFFTpeFocal[index][0] + gFFTpeFocal[index][1] * gFFTpeFocal[index][1]) * gReferenceAmplitude[k + l * gWXsizePE];
		}
	
		// divide by the total amplitude contained in the signal f
		amplc = amplc / gTotalIntensity;
		
		// get the phase inside the window in the focal plane
		GetPhaseFromWindow(gFFTpeFocal, gNxPE, gNyPE, focalphase, gWxPE, gWyPE, gWXsizePE, gWYsizePE);
		
		// DEBUG: get the amplitude inside the window in the focal plane
		GetAmplitudeFromWindow(gFFTpeFocal, gNxPE, gNyPE, focalamp, gWxPE, gWyPE, gWXsizePE, gWYsizePE);
		plotField(focalamp, gWXsizePE, gWYsizePE, debugpanel, dbc4);
		//plotFFTField(gFFTpeFocal, gNxPE, gNyPE, debugpanel, dbc6);
		
		// calculate and plot convergence
		double famax = 0.0;
		for (int k = 0; k < gWXsizePE * gWYsizePE; k++)
		{
			focalamp[k] /= amplc;
			if (focalamp[k] > famax)
				famax = focalamp[k];
		}
		for (int k = 0; k < gWXsizePE * gWYsizePE; k++)
			focalamp[k] = (focalamp[k] > gNoiseCutoff * famax ? focalamp[k] : 0.0);
		convergence[n] = calcSquaredDifference(gReferenceAmplitude, focalamp, gWXsizePE * gWYsizePE) / (gWXsizePE * gWYsizePE);
		if (n > 1)
		{
			DeleteGraphPlot(debugpanel, debuggraph2, -1, VAL_DELAYED_DRAW);
			PlotY(debugpanel, debuggraph2, &(convergence[1]), n - 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0xFF);			
			RefreshGraph(debugpanel, debuggraph2);
		}
		
		// apply amplitude constraints in the focal plane
		SetValuesInWindow(gFFTpeFocal, gNxPE, gNyPE, gReferenceAmplitude, focalphase, gWxPE, gWyPE, gWXsizePE, gWYsizePE, 0, amplc);
		
		GetAmplitudeFromWindow(gFFTpeFocal, gNxPE, gNyPE, focalamp, gWxPE, gWyPE, gWXsizePE, gWYsizePE);
		//plotField(focalphase, gWXsizePE, gWYsizePE, debugpanel, dbc2);
		plotField(focalamp,   gWXsizePE, gWYsizePE, debugpanel, dbc5);
		//plotFFTField(gFFTpeFocal, gNxPE, gNyPE, debugpanel, dbc6);
		
		// DEBUG: get the amplitude inside the window in the focal plane
		//GetAmplitudeFromWindow(gFFTpeFocal, gNxPE, gNyPE, focalamp, gWxPE, gWyPE, gWXsizePE, gWYsizePE);
		//plotField(focalamp, gWXsizePE, gWYsizePE, debugpanel, dbc4);
		
		// propagate backwards from the focal plane to the SLM plane
		fftw_execute(gFFTpeFS);
		
		// get the phase of the window in the SLM plane
		GetPhaseFromWindow(gFFTpeSLM, gNxPE, gNyPE, slmphase, gSLMxPE, gSLMyPE, gXsize, gYsize);
		
		// get amount of new phase we want to mix in with the old phase
		//double beta = 0.05;//((double) n + 1.0) / ((double) N);

		// unwrap phase
		//plotField(slmphase,   gXsize, gYsize, debugpanel, dbc2);
		//SLM_unwrapPhase(slmphase, gXsize, gYsize);
		//SLM_unwrapPhase(prjphase, gXsize, gYsize);
		//plotField(slmphase,   gXsize, gYsize, debugpanel, dbc3);
		//SLM_unwrapPhase(slmphase, gXsize, gYsize);
		plotField(slmphase,   gXsize, gYsize, debugpanel, dbc1);
		
		
		// perfrom a weighted average of the new phase with the old one, 
		// and subtract the SLM base phase pattern
		for (int k = 0; k < gXsize * gYsize; k++)
		{
			// perform the averaging
			// NOTE: after averaging, the result is guaranteed to be between -PI and PI
			slmphase[k] = SLM_averagePhase(prjphase[k], slmphase[k], beta) - gSLMBasePhase[k];
		}
		SLM_unwrapPhase(slmphase, gXsize, gYsize);
		plotField(slmphase,   gXsize, gYsize, debugpanel, dbc2);
		
		// project onto Zernike space
		if (n % prjinterval == 0)
			ProjectPhaseOntoZernikeSpace(slmphase, cj[n], gXsize, gYsize, gNPmax);
		
		// DEBUG plot
		if (n >= 1)
		{
			zdists[n - 1] = calcSquaredDifference(&(cj[n - 1][1]), &(cj[n][1]), NPmax - 1);
			DeleteGraphPlot(debugpanel, debuggraph3, -1, VAL_DELAYED_DRAW);
			if (n > 1)
				PlotY(debugpanel, debuggraph3, zdists, n - 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, 0xFF0000);			
			RefreshGraph(debugpanel, debuggraph3);
		}
		
		// DEBUG plot of the unwrapped phase
		SLM_unwrapPhase(slmphase, gXsize, gYsize);
		plotField(slmphase,   gXsize, gYsize, debugpanel, dbc3);
		
		// add the base phase back
		for (int k = 0; k < gXsize * gYsize; k++)
			prjphase[k] = slmphase[k] + gSLMBasePhase[k];
		
		// DEBUG plot
		int ncol = (255 / (N - 1 > 0 ? N - 1 : 1));
		unsigned int color = ((n - 1) * ncol) + 0xAA00;
		PlotY(debugpanel, debuggraph, cj[n], gNPmax, VAL_DOUBLE, VAL_THIN_LINE, VAL_SMALL_EMPTY_SQUARE, VAL_SOLID, 1, color);			
		
		// set the amplitude and phase in the SLM window, set the remainder of the SLM plane to zero
		SetValuesInWindow(gFFTpeSLM, gNxPE, gNyPE, gInputAmplitude, prjphase, gSLMxPE, gSLMyPE, gXsize, gYsize, 1, 1.0);
	}
	
	// clean up temporary arrays
	free(prjphase);
	free(slmphase);
	free(focalphase);
	free(focalamp);

	// copy the last row of Zernike coefficients, and clear the cj array
	double* cjtmp = (double*) malloc(gNPmax * sizeof(double));
	for (int j = 1; j <= gNPmax; j++)
		cjtmp[j - 1] = cj[N][j - 1];
	for (int n = 0; n <= N; n++)	
		free(cj[n]);
	free(cj);
	
	// return the squared difference of the final signal with the reference amplitude
	*RefAmpMatch = convergence[N];
	
	// clean up the convergence arrays
	free(convergence);
	free(zdists);
	
	// return the Zernike coefficients
	return cjtmp;
}


/// HIFN performs a Aitken transformation of a sequence, accelerating its convergence
double* AitkenTransform(double *seq, int N)
{
	// the return var
	double* aseq;
	
	// check if we have at least 3 elements in the array	
	if (N > 2)
	{
		// allocate memory for storing the accelerated sequence
		aseq = (double*) malloc((N - 2) * sizeof(double));
		
		// compute the elements of the accelerated sequence
		for (int k = 2; k < N; k++)
		{
			double dx  = seq[k - 1] -       seq[k - 2];
			double dx2 = seq[k - 2] - 2.0 * seq[k - 1] + seq[k];
			aseq[k - 2] = seq[k - 2] - (dx * dx / dx2);
		}
		
		return aseq;
	}
	else 
	{
		// not enough elements to accelerate the sequence, we simply copy the input
		aseq = (double*) malloc(N * sizeof(double));
		memcpy(aseq, seq, N * sizeof(double));
		return seq;
	}
}


/// HIFN evaluates a set of Zernike coefficients, i.e. the aberration correction is loaded to the SLM and the resulting camera image is rated
void SLM_evaluateAberrationCorrection(int SLMPixelsPanel, int SLMPixelsCanvas, int SimPanel, int SimCanvas, 
		double* ZernikeCoeffs, int NZ, double** correlation, double** sharpness)
{
	// make sure no aberration correction is applied to the SLM
	SLM_setAberrationCorrection(NULL, 0, 0, 0);
	
	// simulate the intensity in the focal plane
	
	// store the signal window of the simulation, this is our theoretical ideal situation
	
	// apply the specified aberration correction to the SLM
	//SLM_setAberrationCorrection(gCj, NZP);
			
	// update the SLM
	SLM_update(SLMPixelsPanel, SLMPixelsCanvas, SimPanel, SimCanvas, 1);
	
	// wait a short while for the signal to adapt
	
	// grab a camera frame (or more?)
	
	// retrieve the signal window
	
	// compute cross correlation with the simulated intensity
	
	// find the maximum value of the cross correlation, this is a measure for how good the camera signal matches the ideal situation
	
	// compute the Fourier transform of the amplitude
	
	// find the maximum value of the frequency components related to the spot spacing, this is a measure for the sharpness
	
}	
	


/// HIFN get the SLM base phase pattern used in the phase error characterisation
double* SLM_getCharacterisationBasePhase()
{
	return gSLMBasePhase;
}


/// HIFN get the reference amplitude and its parameters used in the phase error characterisation
void SLM_getCharacterisationReferenceAmplitude(double** refampl, int* wx, int* wy, int* wxpos, int* wypos)
{
	*wx = gWXsizePE;
	*wy = gWYsizePE;
	*wxpos = gWxPE;
	*wypos = gWyPE;
	*refampl = gReferenceAmplitude;
}
