//==============================================================================
//
// Title:       SLM.c
// Purpose:     Library of general functions for controlling the pixels of the
//              SLM, and simulating the results. Includes initialisation and
//              memory allocating etc.
//
// Created on:  1-9-2011 at 19:46:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#include "SLM_internal.h"
#include "SLM_Zernike.h"
#include <stdio.h>
#include <stdlib.h>



// Number of pixels in each dimension on the SLM pixel canvas
int gXsize;
int gYsize;

// simulate intensity (1) or not (0)
unsigned int gSimulate;

// simulate SLMpixels (1) or not (0)
unsigned int gSimulateSLMpixels;

// show the phase instead of intensity in the simulation
unsigned int gSimPhase;

// the simulation saturation treshold, relative to the maximum intensity
double gSimSaturation;

// the simulation lens focal length (in m)
double gFocalLength;

// the magnification of focal length
double gMagnification;

// the simulation wavelength (in m)
double gWavelength;

// the (Gaussian) input intensity pattern
double* gInputAmplitude;

// the separable version of the input intensity pattern
double* gInputAmplitudeX;
double* gInputAmplitudeY;

// parameters of the Gaussian input intensity
double gIxcenter, gIycenter, gIxsigma, gIysigma;

// the physical dimensions of the SLM (in m)
const double LxSLM = 1.5 * 0.01;
const double LySLM = 0.8 * 0.01;

// the zoom factor of the pattern
double gZoom = 1.0;

// the focal unit, i.e. the elementary pixel size in the focal plane
double gFocalUnitX = 1.0, gFocalUnitY = 1.0;

// the signal intensity, set by most pattern generators, and used for feedback purposes
double* gSignal;

// variables indicating the window position and size within the gSignal array
int gWindowX, gWindowY, gWindowXsize, gWindowYsize;

// the SLM base phase pattern
double* gSLMphase;

// the total SLM phase pattern, without having taken the modulus
double* gSLMtotal;

// the total SLM phase pattern, the same size as the canvas (for use when subsampling the SLM)
double* gSLMFPtotal;

// the SLM horizontal and vertical translation phase patterns
double* gSLMhtrans;
double* gSLMvtrans;

// the SH correction for the SLM
double* gSLMSHpattern;

// the factory correction and the SH correction to the SLM phase pattern
double* gSLMFactoryCorrection;
double* gSLMSHcorrection;

// the SLM lens correction phase pattern
unsigned char* gSLMLensX;
unsigned char* gSLMLensY;
double* gSLMXsquared;
double* gSLMYsquared;

// the SLM aberration correction phase pattern, that will be SUBTRACTED from the global SLM phase
double* gSLMaberration;

// coefficients of Zernike polynomials, describing the aberration correction
double* gZernikeCoefficients;

// number of Zernike polynomials used to describe the aberration correction
int gNumZernikeCoefficients;

// the current settings for translation, aberration and lens phase (focal length in m) and lens position
double gLensXphase = 0.0, gLensYphase = 0.0, gHorizTrans = 0.0, gVertTrans = 0.0, gAberration = 0.0, gLensX = 0.0, gLensY = 0.0;

// SLM pixel data in which we can directly write
double* gSLMPatternPixels;

// SLM pixel data to which we have added the translation phases, lens phases, etc.
unsigned char* gSLMPixels;

// variables for the FFTs used in the IFTAs
fftw_complex *gFFTin, *gFFTout;
fftw_plan gFFTplan_Gg, gFFTplan_gG;

// variables for the FFT used for the simulation
fftw_complex *gSimFFTin, *gSimFFTout;
fftw_plan     gSimFFTplan;

// the color table for the SLM bitmap
int gSLMColorMap[256];

// a bitmap structure to transfer SLM and simulation data to other
// information carriers (e.g. a Canvas to draw on, or a file to save in)
int gBitmap = -1;


// the bias phase of the SLM (offset for all pixels)
unsigned char gBias = 0;

// variable that indicates whether the SLM should be subsampled,
// and how many pixels should form a new unit pixel (1 = no subsampling)
int gSubSampleFactor = 1;

// the size of the canvas displaying the SLM pixels
int gCanvasXsize, gCanvasYsize;

// raw camera data for the input intensity (when loaded from file)
static double* gInputCameraData;
static int gInputCamNx, gInputCamNy;

// input intensity mode indicator
#define SLM_INPUT_GAUSSIAN 		0
#define SLM_INPUT_ANNULAR  		1
#define SLM_INPUT_LOADFROMFILE 	2
static int gInputIntensityMode = 0;

// helper functions that set the lens and aberration patterns
void SLM_updateLensPattern(void);

// debug plotting controls
int gDebugPanel, gDebugCanvas1, gDebugCanvas2, gDebugCanvas3, gDebugCanvas4,  gDebugCanvas5, gDebugCanvas6, gDebugGraph1, gDebugGraph2, gDebugGraph3;


//==============================================================================
// SLM general control functions


/// HIFN (re)initialises all arrays and the canvas and bitmap structures of the SLM panel
///      used not only at start of program, but also when the SLM panel is resized
void SLM_initialise(int Xsize, int Ysize, int SimXsize, int SimYsize, int SubSampleFactor)
{
	// store the subsampling factor
	gSubSampleFactor = SubSampleFactor;

	// store the canvas sizes
	gCanvasXsize = Xsize;
	gCanvasYsize = Ysize;

	// adjust the dimensions for the subsampling
	Xsize /= SubSampleFactor;
	Ysize /= SubSampleFactor;

	// if already initialised, try to maintain the current SLM pattern
	double* tmpPhase;
	int useOldPhase = 0;
	if (gSLMphase != NULL)
	{
		// indicate that we should use the old phase
		useOldPhase = 1;

		// resample the phase
		tmpPhase = SLM_resamplePhase(gSLMphase, gXsize, gYsize, Xsize, Ysize);
	}

	// store the (new) sizes of the SLM
	gXsize = Xsize;
	gYsize = Ysize;

	// discard previous bitmap (if present)
	if (gBitmap != -1)
	{
		DiscardBitmap(gBitmap);
		gBitmap = -1;
	}

	// initialize the SLM phase pattern and the corrections (translations and lens phases)
	gSLMphase      = (double *) realloc(gSLMphase,      gXsize * gYsize * sizeof(double));
	gSLMtotal      = (double *) realloc(gSLMtotal,      gXsize * gYsize * sizeof(double));
	gSLMhtrans     = (double *) realloc(gSLMhtrans,             gCanvasXsize * sizeof(double));
	gSLMvtrans     = (double *) realloc(gSLMvtrans,             gCanvasYsize * sizeof(double));
	gSLMLensX       = (unsigned char *) realloc(gSLMLensX, gCanvasXsize * sizeof(unsigned char));
	gSLMLensY       = (unsigned char *) realloc(gSLMLensY, gCanvasYsize * sizeof(unsigned char));
	gSLMXsquared    = (double *) realloc(gSLMXsquared,  gCanvasXsize * sizeof(double));
	gSLMYsquared    = (double *) realloc(gSLMYsquared,  gCanvasYsize * sizeof(double));
	gSLMaberration = (double *) realloc(gSLMaberration, gXsize * gYsize * sizeof(double));
	gSignal        = (double *) realloc(gSignal,        gXsize * gYsize * sizeof(double));
	gSLMSHpattern =  (double *) realloc(gSLMSHpattern, gXsize * gYsize * sizeof(double));
	
	// (re)allocate memory for the input intensity arrays
	gInputAmplitude  = (double*) realloc(gInputAmplitude, gXsize * gYsize * sizeof(double));
	gInputAmplitudeX = (double*) realloc(gInputAmplitudeX, gXsize * sizeof(double));
	gInputAmplitudeY = (double*) realloc(gInputAmplitudeY, gYsize * sizeof(double));
	


	// recalculate the input intensity
	switch (gInputIntensityMode)
	{
		case SLM_INPUT_ANNULAR:

			SLM_setInputIntensityFromCameraData(gIxcenter, gIycenter, gIxsigma, gIysigma, SubSampleFactor);
			break;

		case SLM_INPUT_LOADFROMFILE:

			SLM_setInputIntensityAnnular(gIxcenter, gIycenter, gIxsigma, gIysigma, SubSampleFactor);
			break;

		case SLM_INPUT_GAUSSIAN:
		default:

			SLM_setInputIntensityGaussian(gIxcenter, gIycenter, gIxsigma, gIysigma, SubSampleFactor);
			break;
	}

		
	// compute the x-squared and y-squared arrays for the lens phases
	double x, y;
	for (int k = 0; k < gCanvasXsize; k++)
	{
		// get the current x coordinate in m
		x = (((double) k) / ((double) gCanvasXsize) - 0.5) * LxSLM;
		gSLMXsquared[k] = x * x;
	}
	for (int l = 0; l < gCanvasYsize; l++)
	{
		// get the current y coordinate in m
		y = (((double) l) / ((double) gCanvasYsize) - 0.5) * LySLM;
		gSLMYsquared[l] = y * y;
	}

	// set the lens phase, aberration correction, and translation patterns and the correction
	SLM_setLensPhase(gLensXphase, gLensYphase);
	SLM_setAberrationCorrection(NULL, 0, 0, 0);
	SLM_setHorizTrans(gHorizTrans);
	SLM_setVertTrans(gVertTrans);
	SLM_SetShackHartmannPattern(0, 0, 1000, gHorizTrans, gVertTrans);
	
	
	// Load the factory correction
	LoadFactoryCorrection("C:/Program Files (x86)/SLM Controller UvA/780_correction.png" , 1);
	
	// fill the colormap, following the pattern AARRGGBB
	for (int k = 0; k < 256; k++)
		gSLMColorMap[k] = k * 0x00010101;
	gSLMColorMap[255] = 0x0000FF00;

	// if subsampling is enabled, allocate memory for an intermediate resampled floating point version of the phase
	// that is the same size as the canvas
	if (gSubSampleFactor > 1)
		gSLMFPtotal = (double *) realloc(gSLMFPtotal, gCanvasXsize * gCanvasYsize * sizeof(double));

	// allocate memory for the final SLM pixels
	gSLMPixels = (double *) realloc(gSLMPixels, gCanvasXsize * gCanvasYsize * sizeof(double));

	// initialize the bitmap with the SLM pattern pixel data
	gSLMPatternPixels = (double *) realloc(gSLMPatternPixels, gCanvasXsize * gCanvasYsize * sizeof(double));
	NewBitmap(-1, 8, gCanvasXsize, gCanvasYsize, gSLMColorMap, gSLMPatternPixels, NULL, &gBitmap);

	// check if the FFT variables were already initialised
	if (gFFTin != NULL)
	{
		// yes, free all that memory first
		fftw_free(gFFTin);
		fftw_free(gFFTout);
		fftw_destroy_plan(gFFTplan_Gg);
		fftw_destroy_plan(gFFTplan_gG);
	}

	// initialise the Fourier transform structures for the IFTAs
	gFFTin  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gXsize * gYsize);
	gFFTout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gXsize * gYsize);
	gFFTplan_Gg = fftw_plan_dft_2d(gYsize, gXsize, gFFTin, gFFTout, FFTW_FORWARD,  FFTW_ESTIMATE);
	gFFTplan_gG = fftw_plan_dft_2d(gYsize, gXsize, gFFTout, gFFTin, FFTW_BACKWARD, FFTW_ESTIMATE);

	// see if we already had a phase pattern that can be copied in resized form
	if (useOldPhase == 1)
	{
		// copy the old phase pattern (resampled at the beginning of this function)
		memcpy(gSLMphase, tmpPhase, gXsize * gYsize * sizeof(double));

		// clean up
		free(tmpPhase);
	}

	// initialise the simulation
	SLM_initialiseSimulation(SimXsize, SimYsize);
}


/// HIFN Render the pixels of the SLM onto a canvas
void SLM_update(int Panel, int Canvas, int SimPanel, int SimCanvas, int pattern_updated)
{
	if (pattern_updated > 0)
	{
		// add aberration corrections, and compute the pixel color (unsigned char)
		for (int k = 0; k < gXsize; k++)
		for (int l = 0; l < gYsize; l++)
		{
			// variables for determining the color of the current pixel
			double pixelval;
			
			// get the phase for this pixel
			pixelval = gSLMphase[l * gXsize + k];

			// set the current pixel, subtracting the aberration correction
			gSLMtotal [l * gXsize + k] = pixelval - gSLMaberration[l * gXsize + k];
		}
	
		// conversion factor from 2 * pi to 256
		double convf = 256.0 / (2 * PI);

		// if subsampling is enabled, resample the total phase
		if (gSubSampleFactor > 1)
		{
			SLM_resamplePhaseInPlace(gSLMtotal, gSLMFPtotal, gXsize, gYsize, gCanvasXsize, gCanvasYsize);

			// convert and copy the resampled phase to the SLM pixels
			for (int k = 0; k < gCanvasXsize * gCanvasYsize; k++)
				gSLMPatternPixels[k] = (double) (convf * gSLMFPtotal [k]);
		}
		else
		{
			// convert the floating point phase to unsigned char got rid of the convf
			for (int k = 0; k < gCanvasXsize * gCanvasYsize; k++)
				gSLMPatternPixels[k] = (double) (convf * gSLMtotal [k]);
		}
		
	}

	// add translation phases and lens correction to the SLM pixels, ugly solution b/c General protection fault without source line info.  
	for (int k = 0; k < gCanvasXsize; k++)
	for (int l = 0; l < gCanvasYsize; l++)
	{
		
		if (gSLMFactoryCorrection == NULL && gSLMSHcorrection == NULL)
		{
			gSLMPixels[k + l * gCanvasXsize] = (gSLMPatternPixels[k + l * gCanvasXsize] + gSLMhtrans[k]
										 + gSLMvtrans[l] + gSLMLensX[k] + gSLMLensY[l] + gBias + gSLMSHpattern[k + l * gCanvasXsize]);
		}
		else if (gSLMFactoryCorrection != NULL && gSLMSHcorrection != NULL)
		{
			gSLMPixels[k + l * gCanvasXsize] = (gSLMPatternPixels[k + l * gCanvasXsize] + gSLMhtrans[k]
										 + gSLMvtrans[l] + gSLMLensX[k] + gSLMLensY[l] + gBias + gSLMFactoryCorrection[k + l * gCanvasXsize] + gSLMSHcorrection[k + l * gCanvasXsize]); //(k > gCanvasXsize / 2 ? gBias : 0);
		}
		else if (gSLMFactoryCorrection != NULL && gSLMSHcorrection == NULL)
		{
			gSLMPixels[k + l * gCanvasXsize] = (gSLMPatternPixels[k + l * gCanvasXsize] + gSLMhtrans[k]
										 + gSLMvtrans[l] + gSLMLensX[k] + gSLMLensY[l] + gBias + gSLMFactoryCorrection[k + l * gCanvasXsize] + gSLMSHpattern[k + l * gCanvasXsize]);
		}
		else if (gSLMFactoryCorrection == NULL && gSLMSHcorrection != NULL)
		{
			gSLMPixels[k + l * gCanvasXsize] = (gSLMPatternPixels[k + l * gCanvasXsize] + gSLMhtrans[k]
										 + gSLMvtrans[l] + gSLMLensX[k] + gSLMLensY[l] + gBias + gSLMSHcorrection[k + l * gCanvasXsize]);
		}
	}
	
	/*
	// Convert the right side of the SLM pattern to the blue laser
	for (int k = 0.5 * gCanvasXsize; k < gCanvasXsize; k++)
	{
		for (int l = 0; l < gCanvasYsize; l++)
		{
			gSLMPixels[k + l * gCanvasXsize] = 480/780 * gSLMPixels[k + l * gCanvasXsize];
		}
	}
	*/
	
	
	// update the bitmap
	SetBitmapData(gBitmap, -1, 8, gSLMColorMap, gSLMPixels, NULL);

	// and copy the bitmap to the canvas
	CanvasDrawBitmap(Panel, Canvas, gBitmap, VAL_ENTIRE_OBJECT , VAL_ENTIRE_OBJECT);

	// update the canvas immediately
	CanvasUpdate(Panel, Canvas, VAL_ENTIRE_OBJECT);

	// check if we also have to update the simulated intensity
	if (gSimulate == 1)
		SLM_simulateIntensity(SimPanel, SimCanvas);
	// check if we want to show the SLM pixels
	if (gSimulateSLMpixels == 1)
		CanvasDrawBitmap(SimPanel, SimCanvas, gBitmap, VAL_ENTIRE_OBJECT, VAL_ENTIRE_OBJECT);
}


/// HIFN set the Gaussian input intensity (parameters are all lengths, in m)
void SLM_setInputIntensityGaussian(double x, double y, double sigmax, double sigmay, int subsamplefactor)
{
	// indicate that the current input intensity is a Gaussian
	gInputIntensityMode = SLM_INPUT_GAUSSIAN;

	// store the new parameters
	gIxcenter = x;
	gIycenter = y;
	gIxsigma = sigmax;
	gIysigma = sigmay;

	// set the x intensity to a Gaussian with the specified parameters
	for (int k = 0; k < gXsize; k++)
	{
		// get the current coordinate (in m)
		double x = LxSLM * ((double) k) / ((double) gXsize - 1.0) - 0.5 * LxSLM - gIxcenter;

		// set the x amplitude
		gInputAmplitudeX[k] = exp(-0.5 * (x * x / (sigmax * sigmax)));
	}

	// set the y intensity to a Gaussian with the specified parameters
	for (int l = 0; l < gYsize; l++)
	{
		// get the current coordinate (in m)
		double y = LySLM * ((double) l) / ((double) gYsize - 1.0) - 0.5 * LySLM - gIycenter;

		// set the y amplitude
		gInputAmplitudeY[l] = exp(-0.5 * (y * y / (sigmay * sigmay)));
	}

	// normalise the intensities (not the amplitudes!), such that the total *intensity* integrates to 1
	SLM_normalise_sq(gInputAmplitudeX, gXsize, 0, LxSLM);
	SLM_normalise_sq(gInputAmplitudeY, gYsize, 0, LySLM);

	// set the 2D amplitude as the product of 1D intensities (2D is now also normalised)
	for (int k = 0; k < gXsize; k++)
	for (int l = 0; l < gYsize; l++)
		gInputAmplitude[l * gXsize + k] = gInputAmplitudeX[k] * gInputAmplitudeY[l];
}



/// HIFN set an annular Gaussian input intensity (parameters are all lengths, in m)
void SLM_setInputIntensityAnnular(double x, double y, double sigmax, double sigmay, int subsamplefactor)
{
	// indicate that the current input intensity is an annular Gaussian
	gInputIntensityMode = SLM_INPUT_ANNULAR;
}



/// HIFN calculate input intensity from camera data (loaded with SLM_setInputIntensityCameraData)
void SLM_setInputIntensityFromCameraData(double x, double y, double sigmax, double sigmay, int subsamplefactor)
{
	// indicate that the current input intensity is loaded from camera data
	gInputIntensityMode = SLM_INPUT_LOADFROMFILE;

	// compute the number of SLM pixels that the resized camera data should cover (incl. subsampling correction)
	int camXsize = gXsize * (sigmax / LxSLM) ;
	int camYsize = gYsize * (sigmay / LySLM);

	// resample the camera data in a temporary array
	double* tmpCamData = SLM_resampleBitmap(gInputCameraData, gInputCamNx, gInputCamNy, camXsize, camYsize);

	// compute the pixel coordinates in the SLM plane of the center of the input intensity pattern
	// (note that x and y should be specified relative to the center of the SLM)
	int cx = gXsize * ((x + LxSLM * 0.5) / LxSLM);
	int cy = gYsize * ((y + LySLM * 0.5) / LySLM);

	// loop over all the pixels of the SLM, and set the correct light field amplitude in the input amplitude array
	for (int k = 0; k < gXsize; k++)
	for (int l = 0; l < gYsize; l++)
	{
		// compute the coordinates inside the resampled camera data array
		int ix = (k - cx) + camXsize / 2;
		int iy = (l - cy) + camYsize / 2;

		// see if the coordinate falls inside the array
		if ( (ix >= 0) && (iy >= 0) && (ix < camXsize) && (iy < camYsize) )
		{
			// yes, write the camera data to the input amplitude array
			gInputAmplitude[k + l * gXsize] = tmpCamData[ix + iy * camXsize];
		}
		else
		{
			// no, we just write a zero
			gInputAmplitude[k + l * gXsize] = 0;
		}
	}

	plotField(gInputAmplitude, gXsize, gYsize, 19, 2);

	// create a separable approximation of the input intensity
	memset(gInputAmplitudeX, 0, gXsize * sizeof(double));
	memset(gInputAmplitudeY, 0, gYsize * sizeof(double));
	SLM_SepOp(gInputAmplitude, gXsize, gYsize, &gInputAmplitudeX, &gInputAmplitudeY);

	// clean up the temporary array
	free(tmpCamData);
}




/// HIFN adjusts the horizontal translation of the intensity pattern
void SLM_setHorizTrans(double htrans)
{
	// store the new value
	gHorizTrans = htrans;

	// update the horizontal translation pattern
	for (int k = 0; k < gCanvasXsize; k++)
	{
		gSLMhtrans[k] = (double) (256 * (((double) k) / ((double) gCanvasXsize)) * (htrans / gFocalUnitX));
		int a = gSLMhtrans[k];
		a = a%256;
		gSLMhtrans[k] = a;
	}
					
}


/// HIFN adjusts the vertical translation of the intensity pattern, the translation is specified in meters
void SLM_setVertTrans(double vtrans)
{
	// store the new value
	gVertTrans = vtrans;

	// update the vertical translation pattern
	for (int l = 0; l < gCanvasYsize; l++)
		gSLMvtrans[l] = (double) (256.0 * (((double) l) / ((double) gCanvasYsize)) * (vtrans / gFocalUnitY));
}


/// HIFN sets the SLM lens phase pattern (characterised by two lens phases, to allow for corre
void SLM_setLensPhase(double xstrength, double ystrength)
{
	// store the new values
	gLensXphase = xstrength;
	gLensYphase = ystrength;

	// update the lens phase pattern
	SLM_updateLensPattern();
}

/// HIFN sets a pattern used for Shack-Hartmann aberration correction
void SLM_SetShackHartmannPattern(double spotx, double spoty, double spotdiameter, double strayx, double strayy)
{
	// loop over all the pixels of the SLM
	for (int k = 0; k < gXsize; k++)
	for (int l = 0; l < gYsize; l++)
	{
		// compute the current (x, y) coordinate on the SLM	in meter
		double xum = (-0.5 + (k) / ((double) gXsize - 1)) * LxSLM;
		double yum = (-0.5 + (l) / ((double) gYsize - 1)) * LySLM;
		
		// is this pixel part of the spot?
		double dx = xum - spotx;
		double dy = yum - spoty;
		if ((dx * dx + dy * dy) < (spotdiameter * spotdiameter) / 4.0)
		{
			// yes, part of the spot
			gSLMSHpattern[k + l * gXsize] = (double) (256 * ((double)k / ((double)gXsize - 1)) * ((strayx) / gFocalUnitX)) 
									  + (double) (256 * ((double)l / ((double)gYsize - 1)) * ((strayy) / gFocalUnitY));	 
		}
		else
		{
			// no, not part of the spot, but part of the stray light Switched 2 Pi for 256
			//gSLMphase[k + l * gXsize] = (double) (2 * PI * (((double) k) / ((double) gXsize - 1)) * (strayx / gFocalUnitX)) 
			//+ (2 * PI * (((double) l) / ((double) gYsize - 1)) * (strayy / gFocalUnitY));
			gSLMSHpattern[k +l * gXsize] = 0;
		}
		
	}
	
	// set the current pattern indicator
	//gCurrentPattern = SLM_SHACKHARTMANN;
}



/// HIFN calculates the lens phase pattern from the current settings
void SLM_updateLensPattern()
{
	// calculate the lens phase amplitude: gLensPhase holds a focal length,
	// which is a correction to the focal length of the true lens,
	// and from this, together with the wavelength, we can compute a phase
	// we do this for both x and y directions
	double qamplx = PI / (gWavelength * (gLensXphase + gFocalLength)) - PI / (gWavelength * (gFocalLength));
	double qamply = PI / (gWavelength * (gLensYphase + gFocalLength)) - PI / (gWavelength * (gFocalLength));


	// conversion factor from 2 * pi to 256
	double convf = 256.0 / (2 * PI);

	// loop over all pixels of the SLM, and calculate the lens phase there
	for (int k = 0; k < gCanvasXsize; k++)
		gSLMLensX[k] = (unsigned char) (convf * qamplx * gSLMXsquared[k]);
	for (int l = 0; l < gCanvasYsize; l++)
		gSLMLensY[l] = (unsigned char) (convf * qamply * gSLMYsquared[l]);
}

// HIFN loads a factory correction to the SLM pattern
void LoadFactoryCorrection(char filename[], int Check)
{
	if (filename == 0 || Check == 0)
	{
		gSLMFactoryCorrection = NULL;
	}
	else
	{
							
		int FactoryCorrectionBitmap = -1;
		GetBitmapFromFile(filename, &FactoryCorrectionBitmap);
	
		// allocate variables to hold the bitmap data
		int BytesPerRow, PixelDepth, Width, Height;
		unsigned char* Correction;
		AllocBitmapData (FactoryCorrectionBitmap, NULL, &Correction, NULL);
		// get bitmap pixel data
		GetBitmapData (FactoryCorrectionBitmap, &BytesPerRow, &PixelDepth, &Width, &Height, NULL, Correction, NULL);
	
		// allocate array to hold the signal intensity
		double* FactoryCorrection = (double*) malloc(Width * Height * sizeof(double));
		
		// extract the GREEN channel pixel data
		for (int k = 0; k < Width * Height; k++)
			FactoryCorrection[k] = ((double) Correction[k * (PixelDepth >> 3)+1]);
	
		// Resample to the SLM size
		gSLMFactoryCorrection = FactoryCorrection; //SLM_resamplePhase(FactoryCorrection, Width, Height, gCanvasXsize, gCanvasYsize);
	}
	
}


// HIFN loads a SH surface correction to the SLM pattern
void LoadSHcorrection(char filename[], int check)
{
	if (filename == 0 || check == 0)
	{
		gSLMSHcorrection = NULL;
	}
	else
	{
		int SHcorrectionBitmap = -1;
		GetBitmapFromFile(filename, &SHcorrectionBitmap);
		
		//allocate variable to hold the bitmap data
		int BytesPerRow, PixelDepth, Width, Height;
		unsigned char* Correction;
		AllocBitmapData(SHcorrectionBitmap, NULL, &Correction, NULL);
		//get bitmap pixel data
		GetBitmapData (SHcorrectionBitmap, &BytesPerRow, &PixelDepth, &Width, &Height, NULL, Correction, NULL);
		
		// allocate array to hold the signal intensity
		double* SHcorrection = (double*) malloc(Width * Height * sizeof(double));
		
		//fill gSLMSHcorrection
		for (int k = 0; k < Width * Height; k++)
			SHcorrection[k] = ((double) Correction[k * (PixelDepth >> 3)+1]);
		gSLMSHcorrection = SHcorrection;
	}
}


/// HIFN set the simulation lens focal length (in mm)
void SLM_setLensFocalLength(double f)
{
	gFocalLength = f/gMagnification;

	// compute the new focal units
	SLM_updateFocalUnit();
}

/// HIFN set the telescope correction
void SLM_setMagnification(double h)
{
	gMagnification = h;
	
	// compute the new focal length
	SLM_setLensFocalLength(gFocalLength);
}

/// HIFN computes the focal unit, and updates gFocalUnitX and gFocalUnitY
void SLM_updateFocalUnit(void)
{
	// compute the focal units from the lens focal length, wavelength, and the SLM dimensions
	gFocalUnitX = gFocalLength * gWavelength / LxSLM;
	gFocalUnitY = gFocalLength * gWavelength / LySLM;
}


/// HIFN set the simulation wavelength (in nm)
void SLM_setWavelength(double lambda)
{
	gWavelength = lambda;

	// the lens correction factor also depends on the wavelength, so
	// we also update that
	SLM_updateLensPattern();

	// compute the new focal units
	SLM_updateFocalUnit();
}

/// HIFN sets the SLM pattern zoom factor
void SLM_setZoom(double zoom)
{
	gZoom = zoom;
}


/// HIFN get the SLM x resolution
int SLM_getXres(void)
{
	return gXsize;
}


/// HIFN get the SLM y resolution
int SLM_getYres(void)
{
	return gYsize;
}


/// HIFN get the SLM canvas x resolution
int SLM_getCanvasXres(void)
{
	return gCanvasXsize;
}


/// HIFN get the SLM canvas y resolution
int SLM_getCanvasYres(void)
{
	return gCanvasYsize;
}


/// HIFN get the lens focal length (in m)
double SLM_getLensFocalLength(void)
{
	return gFocalLength;
}


/// HIFN get the wavelength (in m)
double SLM_getWavelength(void)
{
	return gWavelength;
}


/// HIFN get the focal unit in the X direction (in meters)
double SLM_getFocalUnitX(void)
{
	return gFocalUnitX;
}


/// HIFN get the focal unit in the Y direction (in meters)
double SLM_getFocalUnitY(void)
{
	return gFocalUnitY;
}


/// HIFN the amplitude modulation of the SLM as a function of the phase
double SLM_AmplitudeModulation(double phase)
{
	// reduce the phase such that it is in the interval [0 .. 2 pi]
	double rphase = fmod(phase, 2 * PI);
	if (rphase < 0.0)
		rphase = rphase + 2 * PI;

	// return the amplitude modulation as measured
	return sqrt(1.0 + gAmplitudeModulation * sin(0.6 * rphase + 0.4 * rphase * rphase / (2.0 * PI)));
}


/// HIFN sets the amplitude modulation
void SLM_setAmplitudeModulation(double amplmod)
{
	gAmplitudeModulation = amplmod;
}


/// HIFN sets the offset (bias) phase for all pixels (input: 0..255)
void SLM_setBias(int bias)
{
	gBias = (unsigned char) (bias & 255);
}


/// HIFN sets the coordinates of the focusing lens
void SLM_setLensCoordinates(double lensX, double lensY)
{
	gLensX = lensX;
	gLensY = lensY;

	// update the lens phase pattern
	SLM_updateLensPattern();
}


/// HIFN get the current base phase pattern of the SLM
double* SLM_getPhase(void)
{
	return gSLMphase;
}


/// HIFN copies the SLM phase from a specified memory location
void SLM_setPhase(double* src)
{
	memcpy((void*) gSLMphase, (void*) src, gXsize * gYsize * sizeof(double));
}

/// HIFN clear the SLM base phase pattern (i.e. not the lens, translation, corrections.. etc.)
void SLM_Clear()
{
	memset((void*) gSLMphase, 0, gXsize * gYsize * sizeof(double));
}


/// HIFN get the current zoom of the SLM pattern
double SLM_getZoom(void)
{
	return gZoom;
}


/// HIFN get the current amplitude modulation
double SLM_getAmplitudeModulation(void)
{
	return gAmplitudeModulation;
}


/// HIFN get the current pattern mode
int SLM_getCurrentPattern(void)
{
	return gCurrentPattern;
}


/// HIFN get the signal window area position and dimensions
void SLM_getSignalWindowParameters(int *wx, int *wy, int *wxsize, int* wysize)
{
	*wx = gWindowX;
	*wy = gWindowY;
	*wxsize = gWindowXsize;
	*wysize = gWindowYsize;
}


/// HIFN get the signal window contents
double* SLM_getSignalWindow()
{
	double* tmp = (double*) malloc(gWindowXsize * gWindowYsize * sizeof(double));

	for (int k = 0; k < gWindowXsize; k++)
	for (int l = 0; l < gWindowYsize; l++)
		tmp[k + l * gWindowXsize] = gSignal[((k + gWindowX + 2 * gXsize) % gXsize) + ((l + gWindowY + 2 * gYsize) % gYsize) * gXsize];

	return tmp;
}

/// HIFN get the physical size of the SLM in the x-direction (in m)
double SLM_getLxSLM()
{
	return LxSLM;
}


/// HIFN get the physical size of the SLM in the y-direction (in m)
double SLM_getLySLM()
{
	return LySLM;
}


/// HIFN set a phase pattern, specified by an array of coefficients of Zernike polynomials,
///		 to correct any phase aberrations of the beamline
void SLM_setAberrationCorrection(double* ZernikeCoefficients, int Nz, int mx, int my)
{
	// clear the aberration phase pattern
	memset(gSLMaberration, 0, gXsize * gYsize * sizeof(double));

	// store the number of Zernike coefficients
	gNumZernikeCoefficients = Nz;

	// check if there is a nonzero number of Zernike coefficients supplied
	if (Nz > 0)
	{
		// reallocate the array of Zernike coefficients and copy
		gZernikeCoefficients = (double*) realloc(gZernikeCoefficients, Nz * sizeof(double));
		memcpy(gZernikeCoefficients, ZernikeCoefficients, Nz * sizeof(double));

		// allocate a temporary array for holding the sampled Zernike polynomials
		double* Ztmp = (double*) malloc(gXsize * gYsize * sizeof(double));

		// add the contribution of each Zernike polynomial to the aberration phase pattern
		for (int j = 1; j <= Nz; j++)
		{
			// sample the j-th Zernike rectangle polynomial
			FillArrayZernikeJ(Ztmp, j, gXsize, gYsize);

			// add the polynomial to the global aberration correction
			for (int k = 0; k < gXsize; k++)
			for (int l = 0; l < gYsize; l++)
			{
				int ik = (mx == 1 ? gXsize - 1 - k : k);
				int il = (my == 1 ? gYsize - 1 - l : l);
				gSLMaberration[k + l * gXsize] += gZernikeCoefficients[j - 1] * Ztmp[ik + il * gXsize];
			}
		}

		// free the temporary array
		free(Ztmp);
	}
	else
	{
		if (gZernikeCoefficients != NULL)
			free(gZernikeCoefficients);
	}
}


/// HIFN get the coefficients of rectangular Zernike polynomials describing the aberration correction of the SLM
double* SLM_getZernikeAberrationCorrection()
{
	return gZernikeCoefficients;
}


/// HFIN store the raw data for the input intensity, when a measurement from a camera is used
void SLM_setInputIntensityCameraData(unsigned char* camdata, int camxres, int camyres,
	int camxmin, int camymin, int camNx, int camNy)
{
	// store the pixel dimensions of the camera data
	gInputCamNx = camNx;
	gInputCamNy = camNy;

	// (re)allocate memory for holding the raw camera input intensity data
	gInputCameraData = (double*) realloc(gInputCameraData, gInputCamNx * gInputCamNy * sizeof(double));

	// copy the zoomed area from the camera frame (and take the square root, to go from intensity to amplitude)
	for (int k = 0; k < gInputCamNx; k++)
	for (int l = 0; l < gInputCamNy; l++)
		gInputCameraData[k + l * gInputCamNx] = sqrt((double) camdata[(k + camxmin) + (l + camymin) * camxres]);
}


/// HIFN assign the debug plotting controls
void SLM_setDebugPlottingControls(int dbgpanel, int dbgc1, int dbgc2, int dbgc3, int dbgc4, int dbgc5, int dbgc6, int dbgg1, int dbgg2, int dbgg3)
{
	gDebugPanel = dbgpanel;
	gDebugCanvas1 = dbgc1;
	gDebugCanvas2 = dbgc2;
	gDebugCanvas3 = dbgc3;
	gDebugCanvas4 = dbgc4;
	gDebugCanvas5 = dbgc5;
	gDebugCanvas6 = dbgc6;
	gDebugGraph1 = dbgg1;
	gDebugGraph2 = dbgg2;
	gDebugGraph3 = dbgg3;
}
