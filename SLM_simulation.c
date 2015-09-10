//==============================================================================
//
// Title:       SLM Simulation
// Purpose:     Library of general functions for simulating the light field
// 				in the focal plane resulting from the current SLM phase pattern.
//
// Created on:  05-01-2012 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#include "SLM_internal.h"

// the simulation dimensions
int gSimXsize, gSimYsize;

// the dimensions of the simulation CANVAS
int gSimCanvasXsize, gSimCanvasYsize;

// simulate intensity (1) or not (0)
unsigned int gSimulate;

// simulate SLMpixels (1) or not (0)
unsigned int gSimulateSLMpixels;

// show the phase instead of intensity in the simulation
unsigned int gSimPhase;

// display grid lines?
int gSimShowGrid;

// grid spacing
double gSimGridSpacing;

// the simulation saturation treshold, relative to the maximum intensity
double gSimSaturation;

// the simulation zoom factor
double gSimZoom = 1.0;

// the minimum zoom factor such that the whole SLM is sampled
double gMinSimZoom;

// pixel data for the simulation
unsigned char* gSimPixels;

// (resampled) field at the SLM
double *ReducedSLM, *ReducedIntensity;

// Helmholtz propagation phases
double* gHelmholtzPhase1;
double* gHelmholtzPhase2;
double* gHelmholtzLens;

// the Helmholtz propagation lengths
double gSLMz = 0.09;
double gSLMzPrev = 0.0;
double gLensFocalLengthPrev = 0.0;

// SLM pixel values converted to doubles and rescaled to 0..2PI
double* gSimSLMsampling;

// convenience variables to keep track of changes in the number of sampling points on the SLM, and the Helmholtz simulation arrays
int NxSLMprev, NySLMprev, Nxprev, Nyprev;

// the color table for the simulation
int* gSimColorMap;

// the amplitude modulation
double gAmplitudeModulation = 0.0;

// variables for the FFT used for the simulation
fftw_complex *gSimFFT_SLM, *gSimFFT_focal;
fftw_plan     gSimFFT_plan;
fftw_plan	  gSimFFTinverse_plan;

// should the simulation supersample or not?
int gSuperSampling = 1;

// should we use Helmholtz propagation?
int gHelmholtzPropagation = 0;

// a bitmap structure to transfer simulation data to other 
// information carriers (e.g. a Canvas to draw on, or a file to save in)
int gSimBitmap = -1;


/// HIFN Initialises the Simulation FFT and bitmap 
void SLM_initialiseSimulation(int Xsize, int Ysize)
{
	if (gSuperSampling)
	{
		// set the simulation sizes
		gSimXsize = Xsize;
		gSimYsize = Ysize;
	}
	else
	{
		gSimXsize = SLM_getXres();
		gSimYsize = SLM_getYres();
	}
	
	// compute the minimum simulation zoom, such that (at least) 
	// the whole SLM is sampled in Fourier space 
	double Zx = ((double) gXsize) / ((double) gSimXsize);
	double Zy = ((double) gYsize) / ((double) gSimYsize);
	gMinSimZoom = (Zx > Zy ? (Zx > 1.0 ? Zx : 1.0) : (Zy > 1.0 ? Zy : 1.0));
	
	// make sure the actual zoom does not exceed that 
	gSimZoom = gMinSimZoom;
	
	// initialise the colormap
	if (gSimColorMap != NULL)
		free(gSimColorMap);
	gSimColorMap = SLM_CreateColorMap(255);
	
	// initialize the bitmap with the Simulation pixel data
	gSimPixels = (unsigned char*) realloc(gSimPixels, gSimXsize * gSimYsize * sizeof(unsigned char));
	NewBitmap(-1, 8, gSimXsize, gSimYsize, gSimColorMap, gSimPixels, NULL, &gSimBitmap);
	
	// allocate the SLM pixel sampling data
	gSimSLMsampling = (double*) realloc(gSimSLMsampling, SLM_getCanvasXres() * SLM_getCanvasYres() * sizeof(double));
	
	// check if the FFT variables were already initialised
	if (gSimFFT_SLM != NULL)
	{
		// yes, free all that memory first
		fftw_free(gSimFFT_SLM);
		fftw_free(gSimFFT_focal);
		fftw_destroy_plan(gSimFFT_plan);
	}
	
	// initialise the Fourier transform structures for the IFTAs
	gSimFFT_SLM  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)  * gSimXsize * gSimYsize);
	gSimFFT_focal = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gSimXsize * gSimYsize);
	gSimFFT_plan =        fftw_plan_dft_2d(gSimYsize, gSimXsize, gSimFFT_SLM, gSimFFT_focal, FFTW_FORWARD,  FFTW_ESTIMATE);
	gSimFFTinverse_plan = fftw_plan_dft_2d(gSimYsize, gSimXsize, gSimFFT_focal, gSimFFT_SLM, FFTW_BACKWARD, FFTW_ESTIMATE);
}


void calcHelmholtzPropagationPhase(double* HelmholtzPhase, double z, double Lx, double Ly, int NxSLM, int NySLM)
{
	// loop over all the sampling points
	for (int k = 0; k < NxSLM; k++)
	for (int l = 0; l < NySLM; l++)
	{
		// compute the corresponding frequency
		double fx;
		double fy;
		
		/*if (k < NxSLM / 2)
			fx = ((double) k) / Lx;
		else
			fx = -((double) (k - NxSLM)) / Lx;
		
		if (l < NySLM / 2)
			fy = ((double) l) / Ly;
		else
			fy = -((double) (l - NySLM)) / Ly;//*/
		
		fx = ((double) k - NxSLM / 2) / Lx;
		fy = ((double) l - NySLM / 2) / Ly;
		
		// compute the phase factor due to the Helmholtz propagation over a distance z
		double lambda = SLM_getWavelength();
		double lfxsq = lambda * lambda * fx * fx;
		double lfysq = lambda * lambda * fy * fy;
		HelmholtzPhase[k + l * NxSLM] = (2 * PI * z / lambda) * sqrt(1.0 - lfxsq - lfysq);
	}
}


void calcHelmholtzLensPhase(double* HelmholtzLensPhase, double f, double Lx, double Ly, int Nx, int Ny)
{
	// calculate the quadratic lens phase
	double qamplx = -PI / (gWavelength * f);
	double qamply = -PI / (gWavelength * f);
	
	// x and y coordinates
	double x, y;
	
	// loop over all pixels of the SLM, and calculate the lens phase there
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)          
	{
		x = k * Lx / ((double) (Nx - 1)) - Lx / 2.0;
		y = l * Ly / ((double) (Ny - 1)) - Ly / 2.0;
			
		HelmholtzLensPhase[k + l * Nx] = qamplx * x * x + qamply * y * y;
	}
}


/// HIFN multiplies a given complex valued array by exp(i * phase)
void multiplyByPhase(fftw_complex* amp, double* phase, int Nx, int Ny)
{
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		double cosphi = cos(phase[k + l * Nx]);
		double sinphi = sin(phase[k + l * Nx]);
		double a = amp[k + l * Nx][0];
		double b = amp[k + l * Nx][1];
		amp[k + l * Nx][0] = a * cosphi - b * sinphi;
		amp[k + l * Nx][1] = b * cosphi + a * sinphi;			
	}	
	
}


// swaps the quadrants of a 2D array:
//
//  I   II	     III  IV
//		    -->
// IV  III   	  II   I
//
void fftshift(fftw_complex* f, int Nx, int Ny)
{
	//fftw_complex* tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny / 4);
	
	int MidX = Nx / 2;
	int MidY = Ny / 2;
	for (int k = 0; k < MidX; k++)
	for (int l = 0; l < MidY; l++)
	{
		double retmp;
		double imtmp;
		
		// swap quadrant I with III
		retmp = f[k + l * Nx][0];
		imtmp = f[k + l * Nx][1];
		f[k + l * Nx][0] = f[(k + MidX) + (l + MidY) * Nx][0];
		f[k + l * Nx][1] = f[(k + MidX) + (l + MidY) * Nx][1];
		f[(k + MidX) + (l + MidY) * Nx][0] = retmp;
		f[(k + MidX) + (l + MidY) * Nx][1] = imtmp;
		
		// swap quadrant II with IV
		retmp = f[k + (l + MidY) * Nx][0];
		imtmp = f[k + (l + MidY) * Nx][1];
		f[k + (l + MidY) * Nx][0] = f[(k + MidX) + l * Nx][0];
		f[k + (l + MidY) * Nx][1] = f[(k + MidX) + l * Nx][1];
		f[(k + MidX) + l * Nx][0] = retmp;
		f[(k + MidX) + l * Nx][1] = imtmp;
	}
}


/// HIFN calculates the intensity pattern in a section of the focal plane of a lens, for the current SLM pattern,
///      and wavelength and focal length settings
void SLM_calculateIntensity(fftw_complex* FFT_focal, fftw_complex* FFT_SLM, fftw_plan* FFT_plan, fftw_plan* FFT_plan_inverse,
	double x0, double y0, double Sx, double Sy, int Nx, int Ny)
{
	// conversion factor from 2 * pi to 256
	double convf = 256.0 / (2 * PI);
		
	// calculate the physical size of the area, in the SLM plane, corresponding to the specified
	// area in the focal plane (by Sx and Sy)
	double Lx = gWavelength * gFocalLength * Nx / Sx;
	double Ly = gWavelength * gFocalLength * Ny / Sy;
	
	// calculate the relative size of the SLM within this area
	double Rx = LxSLM / Lx;
	double Ry = LySLM / Ly;
	
	// calculate the number of pixels it would constitute if 
	// we divide up the area in the SLM plane into Nx-by-Ny pixels
	int NxSLM = Rx * Nx;
	int NySLM = Ry * Ny;
	
	// did the number of sampling points on the SLM change?
	if ((NxSLM != NxSLMprev) || (NySLM != NySLMprev))
	{
		// yes, reallocate the reduced intensity fields
		ReducedSLM       = (double*) realloc(ReducedSLM,       NxSLM * NySLM * sizeof(double));
		ReducedIntensity = (double*) realloc(ReducedIntensity, NxSLM * NySLM * sizeof(double));
		
		// update the number of sampling points
		NxSLMprev = NxSLM;
		NySLMprev = NySLM;
	}
	
	// convert the SLM pixel values to double
	for (int k = 0; k < SLM_getCanvasXres() * SLM_getCanvasYres(); k++)
		gSimSLMsampling[k] = ((double) gSLMPixels[k]) / convf;
	
	// resample the SLM phase pattern and input intensity to this size
	SLM_resamplePhaseInPlace(gSimSLMsampling, 	     ReducedSLM, SLM_getCanvasXres(), SLM_getCanvasYres(), NxSLM, NySLM);
	SLM_resampleBitmapInPlace(gInputAmplitude, ReducedIntensity, gXsize,              gYsize,              NxSLM, NySLM);
	
	// set the sampling points in the SLM plane to zero
	memset(FFT_SLM, 0, Nx * Ny * sizeof(fftw_complex));
	
	// embed the resampled intensity and phase within a rectangle of zero intensity
	int OffsetX = (Nx - NxSLM) / 2;
	int OffsetY = (Ny - NySLM) / 2;
	for (int k = 0; k < NxSLM; k++)
	for (int l = 0; l < NySLM; l++)
	{
		// set the complex amplitude at the SLM surface
		unsigned char phase255 = (unsigned char) (convf * ReducedSLM[l * NxSLM + k]);
		double dphase255 = (double) phase255;
		double amplmod = (1.0 - gAmplitudeModulation * (dphase255 - 128.0) * (dphase255 - 128.0) / (128.0 * 128.0));
		FFT_SLM[(OffsetY + l) * Nx + OffsetX + k][0] = ReducedIntensity[l * NxSLM + k] * cos(ReducedSLM[l * NxSLM + k]) * SLM_AmplitudeModulation(ReducedSLM[l * NxSLM + k]);
		FFT_SLM[(OffsetY + l) * Nx + OffsetX + k][1] = ReducedIntensity[l * NxSLM + k] * sin(ReducedSLM[l * NxSLM + k]) * SLM_AmplitudeModulation(ReducedSLM[l * NxSLM + k]);
	}
	
	// calculate the complex amplitude at the focal plane
	if (gHelmholtzPropagation == 1)
	{
		// should we update the Helmholtz propagation phase factors?
		if ((Nx != Nxprev) || (Ny != Nyprev) || (gSLMz != gSLMzPrev) || (SLM_getLensFocalLength() != gLensFocalLengthPrev))
		{
			gHelmholtzPhase1 = (double*) realloc(gHelmholtzPhase1, Nx * Ny * sizeof(double));
			gHelmholtzPhase2 = (double*) realloc(gHelmholtzPhase2, Nx * Ny * sizeof(double));
			gHelmholtzLens   = (double*) realloc(gHelmholtzLens,   Nx * Ny * sizeof(double));
		
			// update the phase factors
			calcHelmholtzPropagationPhase(gHelmholtzPhase1, gSLMz,                    Lx, Ly, Nx, Ny);
			calcHelmholtzPropagationPhase(gHelmholtzPhase2, SLM_getLensFocalLength(), Lx, Ly, Nx, Ny);
			calcHelmholtzLensPhase(gHelmholtzLens, SLM_getLensFocalLength(), Lx, Ly, Nx, Ny);
			plotField(gHelmholtzLens, Nx, Ny, gDebugPanel, gDebugCanvas1);
		
			// store the current values 
			gSLMzPrev = gSLMz;
			gLensFocalLengthPrev = SLM_getLensFocalLength();
			Nxprev = Nx;
			Nyprev = Ny;
		}			
		
		// compute the fourier transform of the SLM plane, the data ends up in the focal plane array
		
		fftshift(FFT_SLM, Nx, Ny);
		fftw_execute(*FFT_plan);
		fftshift(FFT_focal, Nx, Ny);
		
		// multiply with the Helmholtz phase associated with propagation along a distance z
		multiplyByPhase(FFT_focal, gHelmholtzPhase1, Nx, Ny);
		
		plotFFTField(FFT_focal, Nx, Ny, gDebugPanel, gDebugCanvas2);
		plotFFTFieldPhase(FFT_focal, Nx, Ny, gDebugPanel, gDebugCanvas3);
		
		// compute the inverse fourier transform, to obtain the light field at the lens plane
		fftshift(FFT_focal, Nx, Ny);
		fftw_execute(*FFT_plan_inverse);//*/
		fftshift(FFT_SLM, Nx, Ny);
		plotFFTFieldPhase(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas4);
		plotFFTField(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas5);		
		
		fftshift(FFT_focal, Nx, Ny);
		
		/*// multiply with the lens phase
		multiplyByPhase(FFT_SLM, gHelmholtzLens, Nx, Ny);
		
		// aperture
		for (int k = 0; k < Nx; k++)
		for (int l = 0; l < Ny; l++)          
		{
			double x = k * Lx / ((double) (Nx - 1)) - Lx / 2.0;
			double y = l * Ly / ((double) (Ny - 1)) - Ly / 2.0;
			if (x * x + y * y > 0.005)
			{
				FFT_SLM[k + Nx * l][0] = 0;
				FFT_SLM[k + Nx * l][1] = 0;
			}
		}
		
		plotField(gHelmholtzLens, Nx, Ny, gDebugPanel, gDebugCanvas6);
		plotFFTFieldPhase(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas4);
		plotFFTField(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas5);		
		
		
		// compute fourier transform
		fftshift(FFT_SLM, Nx, Ny);
		//plotFFTField(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas4);		
		fftw_execute(*FFT_plan);
		//plotFFTField(FFT_focal, Nx, Ny, gDebugPanel, gDebugCanvas5);
		fftshift(FFT_focal, Nx, Ny);
		//plotFFTField(FFT_focal, Nx, Ny, gDebugPanel, gDebugCanvas6);
		
		
		// multiply with the Helmholtz phase associated with propagation along a distance f
		multiplyByPhase(FFT_focal, gHelmholtzPhase2, Nx, Ny);
		
		// compute the inverse fourier transform to obtain the light field in the focal plane
		fftshift(FFT_focal, Nx, Ny);
		fftw_execute(*FFT_plan_inverse);
		fftshift(FFT_SLM, Nx, Ny);	   //*/
		
		
		/*
		// x and y coordinates
		double x, y;
	
		// loop over all pixels of the SLM, and calculate the lens phase there
		for (int k = 0; k < Nx; k++)
		for (int l = 0; l < Ny; l++)          
		{
			x = k * Lx / ((double) (Nx - 1)) - Lx / 2.0;
			y = l * Ly / ((double) (Ny - 1)) - Ly / 2.0;
			if (x * x + y * y < (Ly * Ly) / 16)
			{
				FFT_SLM[k + Nx * l][0] = cos(gHelmholtzLens[k + l * Nx]);
				FFT_SLM[k + Nx * l][1] = sin(gHelmholtzLens[k + l * Nx]);
			}
			else
			{
				FFT_SLM[k + Nx * l][0] = 0;
				FFT_SLM[k + Nx * l][1] = 0;
			}
		}
		
		plotFFTField(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas3);
		plotFFTFieldPhase(FFT_SLM, Nx, Ny, gDebugPanel, gDebugCanvas4);
		fftshift(FFT_SLM, Nx, Ny);
		fftw_execute(*FFT_plan);
		fftshift(FFT_focal, Nx, Ny);
		plotFFTField(FFT_focal, Nx, Ny, gDebugPanel, gDebugCanvas1);
		plotFFTFieldPhase(FFT_focal, Nx, Ny, gDebugPanel, gDebugCanvas2);
		*/
		
		
	}
	else
	{
		// use Fresnel approximation, the propagation to the focal plane is simply a Fourier transform
		fftw_execute(*FFT_plan);
	}
}


/// HIFN simulates the intensity pattern in the focal plane of a lens, for the current SLM pattern,
///      and renders it onto a canvas
void SLM_simulateIntensity(int Panel, int Canvas)
{
	// calculate the size of the unzoomed sampling area in the focal plane
	double Sx0 = gFocalLength * gWavelength * gXsize / LxSLM;
	double Sy0 = gFocalLength * gWavelength * gYsize / LySLM;
	
	// compute the simulated intensity
	// (note that we don't zoom if we don't supersample, we'll zoom the bitmap instead)
	if (gSuperSampling)
		SLM_calculateIntensity(gSimFFT_focal, gSimFFT_SLM, &gSimFFT_plan, &gSimFFTinverse_plan, 0, 0, 
			Sx0 / gSimZoom, Sy0 / gSimZoom, gSimXsize, gSimYsize);
	else
		SLM_calculateIntensity(gSimFFT_focal, gSimFFT_SLM, &gSimFFT_plan, &gSimFFTinverse_plan, 0, 0, 
			Sx0, Sy0, gSimXsize, gSimYsize);		
	
	
	// get the canvas size and store it (the easiest way of keeping the 
	// simcanvassize variables up to date, not the cleanest)
	GetCtrlAttribute(Panel, Canvas, ATTR_WIDTH,  &gSimCanvasXsize);
	GetCtrlAttribute(Panel, Canvas, ATTR_HEIGHT, &gSimCanvasYsize);
	
	// now we can start drawing the simulation result
	
	// conversion factor from 2 * pi to 256
	double convf = 256.0 / (2 * PI);
	
	// clear the pixels
	//memset(gSimPixels, 0, gSimXsize * gSimYsize * sizeof(unsigned char));
	
	// point to the correct simulation result (different for Helmholtz propagation)
	fftw_complex* simresult;
	/*if (gHelmholtzPropagation == 0)
		simresult = gSimFFT_focal;
	else
		simresult = gSimFFT_focal;*/
	// OVERRIDE
	simresult = gSimFFT_focal;
	
	// check if we want to plot the intensity or the phase
	if (gSimPhase == 0)
	{
		// we want to show the intensity
		
		// find the maximum value in the output intensity, to fix the color scale
		double fmax = 0.0;
		for (int kx = 10; kx < gSimXsize - 10; kx++)
		for (int ky = 10; ky < gSimYsize - 10; ky++)
		{
			double ftmp = (simresult[ky * gSimXsize + kx][0] * simresult[ky * gSimXsize + kx][0] + 
				           simresult[ky * gSimXsize + kx][1] * simresult[ky * gSimXsize + kx][1]);
			if (ftmp > fmax)
				fmax = ftmp;
		}
	
		// set the simulation pixels
		for (int k = 0; k < gSimXsize * gSimYsize; k++)
		{
			// calculate the color of the pixel
			double pixelval = (256.0 * (simresult[k][0] * simresult[k][0] + 
				                        simresult[k][1] * simresult[k][1]) / fmax) / gSimSaturation;
		
			// convert to unsigned char, with saturation
			gSimPixels[k] = (pixelval > 255.0 ? 255 : (unsigned char) pixelval);
		}
	}
	else
	{
		// we want to show the phase
		double xmul = 1.0;
		for (int k = 0; k < gSimXsize; k++)
		{
			// we have to correct for the mirroring of quadrants, 
			// by multiplying each pixel with alternating signs
			xmul *= -1.0;
			double ymul = 1.0;
			
			// set each pixel
			for (int l = 0; l < gSimYsize; l++)
			{
				// again, keep track of alternating sign, now in the y direction
				ymul *= -1.0;
				
				// set the simulation pixel to the phase
				gSimPixels[gSimXsize * l + k] = (atan2(xmul * ymul * simresult[gSimXsize * l + k][1], xmul * ymul * simresult[gSimXsize * l + k][0]) + PI) * convf; 	
			}
		}
	}
		
	// don't immediately update the canvas while we're drawing on it
	CanvasStartBatchDraw (Panel, Canvas);   
				
	// update the bitmap 
	SetBitmapData(gSimBitmap, -1, 8, gSimColorMap, gSimPixels, NULL);
	
	// .. and copy the bitmap to the simulation canvas
	if (gSuperSampling)
		CanvasDrawBitmap(Panel, Canvas, gSimBitmap, VAL_ENTIRE_OBJECT, VAL_ENTIRE_OBJECT);
	else
		CanvasDrawBitmap(Panel, Canvas, gSimBitmap, MakeRect (0, 0, gSimYsize / gSimZoom, gSimXsize / gSimZoom), VAL_ENTIRE_OBJECT);
	
	// check if we have to draw grid lines
	if (gSimShowGrid)
	{
		// compute how many grid lines fit on the canvas
		int nxgrid = (int) ((Sx0 / gSimZoom) / gSimGridSpacing);
		int nygrid = (int) ((Sy0 / gSimZoom) / gSimGridSpacing);
		
		// calculate the spacing between gridlines in pixels
		double gxspacing = ((double) gSimXsize * gSimGridSpacing) / (Sx0 / gSimZoom);
		double gyspacing = ((double) gSimYsize * gSimGridSpacing) / (Sy0 / gSimZoom);

		// draw the X gridlines
		for (int k = 1; k < nxgrid; k++)
			CanvasDrawLine(Panel, Canvas, MakePoint(k * gxspacing, 0), MakePoint(k * gxspacing, gSimYsize));
		
		// draw the Y gridlines
		for (int l = 1; l < nygrid; l++)
			CanvasDrawLine(Panel, Canvas, MakePoint(0, l * gyspacing), MakePoint(gSimXsize, l * gyspacing));
	}
	
	// end the batch drawing and allow the canvas to be updated
	CanvasEndBatchDraw (Panel, Canvas);
}


/// HIFN set the distance SLM - Lens, for use in the Helmholtz propagator
void SLM_setSimZPos(double zpos)
{
	gSLMz = zpos;	
}


/// HIFN toggle the simulation
void SLM_toggleSim(unsigned int sim, unsigned int simphase, unsigned int simpixel)
{
	// update the indicator variables
	gSimulate = sim;	
	gSimPhase = simphase;
	gSimulateSLMpixels = simpixel;
}


/// HIFN sets the relative saturation level of the simulation (relative to the maximum)
void SLM_setSimSaturation(double simsat)
{
	gSimSaturation = simsat;	
}


/// HIFN sets the simulation zoom factor
void SLM_setSimZoom(double simzoom)
{
	// set the zoom, making sure that we do not exceed the minimum zoom
	gSimZoom = (simzoom > gMinSimZoom ? simzoom : gMinSimZoom);
}


/// HIFN sets the simulation's grid spacing
void SLM_setSimGridSpacing(double simgridspacing)
{
	gSimGridSpacing = simgridspacing;
}


/// HIFN toggles the grid lines
void SLM_setSimShowGrid(int simshowgrid)
{
	gSimShowGrid = simshowgrid;
}


/// HIFN toggles whether the simulation supersamples or not
void SLM_setSimSuperSampling(int ssample)
{
	gSuperSampling = ssample;
	
	// reintialise the simulation
	if (gSuperSampling)
		SLM_initialiseSimulation(gSimCanvasXsize, gSimCanvasYsize);
	else
		SLM_initialiseSimulation(gSimXsize, gSimYsize);		
}
	

/// HIFN toggles whether the simulation should use Helmholtz propagation instead of the Fresnel approximation
void SLM_setSimHelmholtz(int helmholtz)
{
	gHelmholtzPropagation = helmholtz;
}

