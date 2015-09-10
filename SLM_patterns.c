//==============================================================================
//
// Title:       SLM_patterns.c
// Purpose:     Library of functions that generate patterns to be displayed on 
//              a Spatial Light Modulator (SLM).
//
// Created on:  1-9-2011 at 19:46:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#include "SLM_internal.h"

// the mode of the phase grating
int gPhaseGratingGaussianPeaks = 0;

extern int gBitmap;

double* gIFTASimulatedSpotIntensities;
double gIFTASimulatedDiffractionEfficiency;

static double* gUnresampledSignal;
int gUnresampledSignalXsize;
int gUnresampledSignalYsize;


//==============================================================================
// SLM spot array generation functions


/// HIFN gets the current values from the SLM spot parameter input controls, 
///      and updates the SLM pattern accordingly
void SLM_setSpotPattern(unsigned int Nx, unsigned int Ny, double xspacing, double yspacing, 
	double randomampl, double spotlens, unsigned int linspotlens)
{
	// Calculate the number of pixels per spot
	int xpixels = (int) ( ((double) gXsize) / ((double) Nx));
	int ypixels = (int) ( ((double) gYsize) / ((double) Ny));
	
	// loop over all spots, and set the SLM spot pattern
	for (int m = 0; m < Nx; m++)
	for (int n = 0; n < Ny; n++)
	{
		// random fluctiations in position
		double xrand = 0.0, yrand = 0.0;
		
		// add random fluctiation to spot positions
		if (randomampl > 0.0)
		{
			xrand = (double) Random(-0.5, 0.5) * randomampl;
			yrand = (double) Random(-0.5, 0.5) * randomampl;
			
			// create hexagonal structure by displacing odd rows
			//xrand = ((n & 1) ? 0 : 0.5);
		}
		
		// loop over all the pixels in each spot
		for (int mx = 0; mx < xpixels; mx++)
		for (int ny = 0; ny < ypixels; ny++)
		{
			// the scaled coordinates relative to the entire SLM
			double x = ((double) mx) / ((double) gXsize - 1.0);
			double y = ((double) ny) / ((double) gYsize - 1.0);
			
			// the scaled coordinates relative to the current spot
			double sx = ((double) mx) / ((double) xpixels - 1.0) - 0.5;
			double sy = ((double) ny) / ((double) ypixels - 1.0) - 0.5;
			
			// calculate the spotlens factor
			double spotlensfactor;
			if (linspotlens == 1)
				spotlensfactor = spotlens * sqrt(sx * sx + sy * sy); //spotlens * (20.0 * sx * sx + 50.0 * sy * sin(6.3 * sx))
			else								   
				spotlensfactor = spotlens * (sx * sx + sy * sy);
			
			// calculate the value of this pixel	
			gSLMphase[(n * ypixels * gXsize) + ny * gXsize + m * xpixels + mx] = 
				((xspacing / gFocalUnitX) * (m + xrand - ((double) Nx - 1.0) / 2.0) * x + 
				 (yspacing / gFocalUnitY) * (n + yrand - ((double) Ny - 1.0) / 2.0) * y +
				 spotlensfactor) * 2 * PI;
		}
	}
	
	// set the pattern indicator
	gCurrentPattern = SLM_BEAMSPLIT_PROJECT;
	
	// finally, we set the signal values
	
	// clear the old signal
	memset(gSignal, 0, gXsize * gYsize * sizeof(double));
	
	// get the number of focal units separating the spots
	double NFx = (xspacing / SLM_getFocalUnitX());
	double NFy = (yspacing / SLM_getFocalUnitY());
	
	// compute the window area in the focal plane
	gWindowX = (int) -(((double) Nx) * NFx / 2.0);
	gWindowY = (int) -(((double) Ny) * NFy / 2.0);
	gWindowXsize = (int) (((double) Nx) * NFx);
	gWindowYsize = (int) (((double) Ny) * NFy);
	
	// the width of each spot in focal units is proportional to the number of spots 
	// (inversely proportional to the area of the SLM that contributes light)
	double sx = Nx;
	double sy = Ny;
	
	// loop over all the spots, and add a gaussian to the signal
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)   
	{
		// the x and y position of the current spot in the focal plane, in focal units
		int xo = gWindowX + (((double) k) + 0.5) * NFx;
		int yo = gWindowY + (((double) l) + 0.5) * NFy;
		
		// loop over the pixels constituting this spot
		for (int px = -((int) (NFx / 2.0)); px < ((int) (NFx / 2.0)); px++)
		for (int py = -((int) (NFy / 2.0)); py < ((int) (NFy / 2.0)); py++)
		{
			// the x and y coordinate of the current pixel, in focal units, and wrapped around
			int xc =(xo + px + 5 * gXsize) % gXsize;
			int yc =(yo + py + 5 * gYsize) % gYsize;
			gSignal[xc + yc * gXsize] = exp(-1.0 * pow((double) px, 2.0) / pow(sx, 2.0) - 1.0 * pow((double) py, 2.0) / pow(sy, 2.0));
		}
	}
	
	// make the signal window a little larger
	// compute the window area in the focal plane
	gWindowX = (int) -(((double) Nx + 1.0) * NFx / 2.0);
	gWindowY = (int) -(((double) Ny + 1.0) * NFy / 2.0);
	gWindowXsize = (int) (((double) Nx + 2.0) * NFx);
	gWindowYsize = (int) (((double) Ny + 2.0) * NFy);
	
	
}


/// HIFN sets the SLM pattern to a phase grating
void SLM_setPhaseGrating(double a1, double b1, double a2, double b2, double angle)
{
	for (int k = 0; k < gXsize; k++)
	for (int l = 0; l < gYsize; l++)
	{
		// calculate the relative position on the screen
		double x = ((double) k) / ((double) gXsize);		
		double y = ((double) l) / ((double) gYsize);
		
		// calculate the sine and cosine of the phase grating angle
		double ca = cos(angle);
		double sa = sin(angle);
		
		if (gPhaseGratingGaussianPeaks == 0)
		{
			// calculate the first sinusoidal wave
			double s1 = a1 * sin(2 * PI * 4 * x * b1);
		
			// calculate the second sinusoidal wave
			double s2 = a2 * sin(2 * PI * 4 * (ca * x + sa * y) * b2);
		
			// set the phase grating pattern
			gSLMphase[l * gXsize + k] = s1 + s2;
		}
		else
		{
			// set the gaussian peak grating pattern
			gSLMphase[l * gXsize + k] = a1 * exp(-b1 * ((x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5)));
		}
	}
	
	// set the pattern indicator
	gCurrentPattern = SLM_BEAMSPLIT_PHASEGRATING;
}


/// HIFN toggles the phase grating from sinusoids to gaussian peaks
void SLM_togglePhaseGratingGaussianPeaks(int gptoggle)
{
	gPhaseGratingGaussianPeaks = gptoggle;
}


/// HIFN generates an array of spots
void SLM_generateSpotArray(int Panel, int Canvas, int SimPanel, int SimCanvas, int OutputPanel, int OutputTextBox,
	                       int numxspots, int numyspots, int spotxspacing, int spotyspacing, int sigxoffset, 
						   int sigyoffset, int spotxoffset, int spotyoffset, int phaseconstraint, double phasestep)
{
	// intialise parameters
	int Nx = gXsize; 
	int Ny = gYsize;
	
	// indicator variable for distinguishing between a 2D and a 1D calculation
	unsigned char calc2D = 0;
	
	// create arrays for the fourier transforms g and G used in the IFTA
	fftw_complex* g1Dx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
	fftw_complex* g1Dy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ny);
	fftw_complex* G1Dx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
	fftw_complex* G1Dy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ny);
	fftw_complex* g2D  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
	fftw_complex* G2D  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
	
	// create arrays for the signal amplitude (f) and phase
	double* f1Dx = (double*) calloc(Nx, sizeof(double));
	double* f1Dy = (double*) calloc(Ny, sizeof(double));
	double* f2D  = (double*) calloc(Nx * Ny, sizeof(double));
	double* phase1Dx = (double*) calloc(Nx, sizeof(double));
	double* phase1Dy = (double*) calloc(Ny, sizeof(double));
	double* phase2D  = (double*) calloc(Nx * Ny, sizeof(double));
	
	// create arrays for the signal and phase masks
	unsigned char* f1Dxmask = (unsigned char*) calloc(Nx, sizeof(unsigned char));
	unsigned char* f1Dymask = (unsigned char*) calloc(Ny, sizeof(unsigned char));
	unsigned char* f2Dmask  = (unsigned char*) calloc(Nx * Ny, sizeof(unsigned char));
	unsigned char* phase1Dxmask = (unsigned char*) calloc(Nx, sizeof(unsigned char));
	unsigned char* phase1Dymask = (unsigned char*) calloc(Ny, sizeof(unsigned char));
	unsigned char* phase2Dmask  = (unsigned char*) calloc(Nx * Ny, sizeof(unsigned char));
	
	// create FFT plans for the IFTA
	fftw_plan FFTplan_gG_1Dx = fftw_plan_dft_1d(Nx, g1Dx, G1Dx,   FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan FFTplan_gG_1Dy = fftw_plan_dft_1d(Ny, g1Dy, G1Dy,   FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan FFTplan_Gg_1Dx = fftw_plan_dft_1d(Nx, G1Dx, g1Dx,   FFTW_FORWARD,  FFTW_ESTIMATE);
	fftw_plan FFTplan_Gg_1Dy = fftw_plan_dft_1d(Ny, G1Dy, g1Dy,   FFTW_FORWARD,  FFTW_ESTIMATE);
	fftw_plan FFTplan_gG_2D  = fftw_plan_dft_2d(Ny, Nx, g2D, G2D, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan FFTplan_Gg_2D  = fftw_plan_dft_2d(Ny, Nx, G2D, g2D, FFTW_FORWARD,  FFTW_ESTIMATE);
	
	// calculate the size of the spot area in pixels
	int spotareax = numxspots + (numxspots - 1) * spotxspacing;
	int spotareay = numyspots + (numyspots - 1) * spotyspacing;
	
	// compute the x and y indices of the start of the actual spot area
	int xspot1 = sigxoffset + spotxoffset;
	int yspot1 = sigyoffset + spotyoffset;
	
	// now we can also calculate the indices of the end of the signal area
	int sigxend = sigxoffset + 2 * spotxoffset + spotareax;
	int sigyend = sigyoffset + 2 * spotyoffset + spotareay;
	
	
	// fill the signal amplitude and phase arrays, and the phase mask
	for (int k = 0; k < numxspots; k++)
	{
		        f1Dx[xspot1 + k * (spotxspacing + 1)] = 1.0;
		    phase1Dx[xspot1 + k * (spotxspacing + 1)] = k * phasestep;
		phase1Dxmask[xspot1 + k * (spotxspacing + 1)] = 1;
	}
	for (int l = 0; l < numyspots; l++)
	{
		        f1Dy[yspot1 + l * (spotyspacing + 1)] = 1.0;
		    phase1Dy[yspot1 + l * (spotyspacing + 1)] = l * phasestep;
		phase1Dymask[yspot1 + l * (spotyspacing + 1)] = 1;
	}
	
	// create the 2D signal array and phase mask, from the 1D arrays
	memset(gSignal, 0, Nx * Ny * sizeof(double));
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		        f2D[l * Nx + k] = f1Dx[k] * f1Dy[l];
			gSignal[l * Nx + k] = f1Dx[k] * f1Dy[l];	
		    phase2D[l * Nx + k] = phase1Dx[k] + phase1Dy[l];
		phase2Dmask[l * Nx + k] = phase1Dxmask[k] * phase1Dymask[l];
	}
	
	// store the signal window coordinates
	gWindowX = sigxoffset;
	gWindowY = sigyoffset;
	gWindowXsize = sigxend - sigxoffset;
	gWindowYsize = sigyend - sigyoffset;
	
	// fill the amplitude mask arrays
	for (int k = sigxoffset; k < sigxend; k++)
		f1Dxmask[k] = 1;
	for (int l = sigyoffset; l < sigyend; l++)
		f1Dymask[l] = 1;
	
	// create the 2D mask from the 1D masks, assuming rectangular masking area(s)
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
		f2Dmask[l * Nx + k] = (f1Dxmask[k] && f1Dymask[l]);
	
	// conversion factor from 2 * pi to 256
	double convf = 256.0 / (2 * PI);
	
	// arrays to store the optimal phase profiles after phase freedom (PF) IFTA
	fftw_complex* G1Dy_opt = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Ny);
	fftw_complex* G1Dx_opt = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
	memcpy(G1Dx_opt, G1Dx, Nx * sizeof(fftw_complex));
	memcpy(G1Dy_opt, G1Dy, Ny * sizeof(fftw_complex));
	
	// check if the phase is to be contrained
	if (phaseconstraint)
	{
		// the phase and amplitude in the signal window are both specified
		// we need to propagate this signal back to the SLM plane to start the
		// AF iterations
		
		// first fill the 'g' array, this is the signal in the focal plane
		for (int k = 0; k < Nx; k++)
			g1Dx[k][0] = f1Dx[k] * cos(phase1Dx[k]) + Random(0, 0.1);
		for (int l = 0; l < Ny; l++)
			g1Dy[l][0] = f1Dy[l] * cos(phase1Dy[l]) + Random(0, 0.1);
		
		// now propagate the signal back to the SLM plane
		fftw_execute(FFTplan_gG_1Dx);
		fftw_execute(FFTplan_gG_1Dy);
			
	}
	else
	{
		// no, the phase is not constrained
		
		// perform IFTA with only phase freedom (PF), we try it a few times with random initial guesses
		// and pick the phases which yield the highest diffraction efficiencies
		int NumTrials = 100;
		double etamax_x = 0.0, etamax_y = 0.0;
		for (int a = 0; a < NumTrials; a++)
		{
			// initialise G1Dx with a random phase
			for (int k = 0; k < Nx; k++)
			{
				double phase = (double) Random(0.0, 2 * PI);
				G1Dx[k][0] = gInputAmplitudeX[k] * cos(((double) phase) / convf);
				G1Dx[k][1] = gInputAmplitudeX[k] * sin(((double) phase) / convf);
			}
		
			// initialise G1Dy with a random phase
			for (int l = 0; l < Ny; l++)
			{
				double phase = (double) Random(0.0, 2 * PI);
				G1Dy[l][0] = gInputAmplitudeY[l] * cos(((double) phase) / convf);
				G1Dy[l][1] = gInputAmplitudeY[l] * sin(((double) phase) / convf);
			}
			
			// perform the first few iterations of the IFTA
			int NumPF = 30;
			IFTA_optimize(f1Dx, NULL, g1Dx, NULL, G1Dx, &FFTplan_Gg_1Dx, &FFTplan_gG_1Dx, Nx, NULL, NULL, NumPF, 0);
			IFTA_optimize(f1Dy, NULL, g1Dy, NULL, G1Dy, &FFTplan_Gg_1Dy, &FFTplan_gG_1Dy, Ny, NULL, NULL, NumPF, 0);
		
			// propagate g1Dx and g1Dy back to the input plane
			fftw_execute(FFTplan_gG_1Dx);
			fftw_execute(FFTplan_gG_1Dy);
		
			// compute the diffraction efficiency
			double etax = calcEta(NULL, G1Dx, Nx);
			double etay = calcEta(NULL, G1Dy, Ny);
		
			// check if the previous maximum is exceeded
			if (etax > etamax_x)
			{
				memcpy(G1Dx_opt, G1Dx, Nx * sizeof(fftw_complex));
				etamax_x = etax;
			}
			if (etay > etamax_y)
			{
				memcpy(G1Dy_opt, G1Dy, Ny * sizeof(fftw_complex));
				etamax_y = etay;
			}
		}
	
		// apply constraints to the optimal G's and copy them back to G again 
		// (so we can use the old arrays and don't need an extra FFT plan)
		SLM_applyConstraints(G1Dx_opt, Nx, NULL);
		SLM_applyConstraints(G1Dy_opt, Ny, NULL);
		memcpy(G1Dx, G1Dx_opt, Nx * sizeof(fftw_complex));
		memcpy(G1Dy, G1Dy_opt, Ny * sizeof(fftw_complex));
	}
	
	// now allow amplitude freedom, and perform a further set of IFTA iterations
	// we have to distinguish whether we do the calculation in 2D or 1D:
	if (calc2D)
	{
		// generate the 2D problem from the product of the 1D arrays
		for (int k = 0; k < Nx; k++)
		for (int l = 0; l < Ny; l++)
		{
			// compute the real part
			G2D[l * Nx + k][0] = G1Dx[k][0] * G1Dy[l][0] - G1Dx[k][1] * G1Dy[l][1];
		
			// compute the imaginary part
			G2D[l * Nx + k][1] = G1Dx[k][0] * G1Dy[l][1] + G1Dx[k][1] * G1Dy[l][0];
		}
		
		// perform IFTA on the 2D array
		int NumAF = 20;
		IFTA_optimize(f2D, NULL, g2D, NULL, G2D, &FFTplan_Gg_2D, &FFTplan_gG_2D, Nx * Ny, f2Dmask, NULL, NumAF, 0);
	} 
	else
	{
		// perform IFTA for the AF stage, for each of the 1D arrays
		int NumAF = 100;
		IFTA_optimize(f1Dx, NULL, g1Dx, NULL, G1Dx, &FFTplan_Gg_1Dx, &FFTplan_gG_1Dx, Nx, f1Dxmask, NULL, NumAF, 0);	
		IFTA_optimize(f1Dy, NULL, g1Dy, NULL, G1Dy, &FFTplan_Gg_1Dy, &FFTplan_gG_1Dy, Ny, f1Dymask, NULL, NumAF, 0);
		
		// if we have a phase constraint, we follow up with a few 2d iterations to see if we can improve the pattern
		if (phaseconstraint)
		{
			// generate the 2D problem from the product of the 1D arrays
			for (int k = 0; k < Nx; k++)
			for (int l = 0; l < Ny; l++)
			{
				// compute the real part
				G2D[l * Nx + k][0] = G1Dx[k][0] * G1Dy[l][0] - G1Dx[k][1] * G1Dy[l][1];
		
				// compute the imaginary part
				G2D[l * Nx + k][1] = G1Dx[k][0] * G1Dy[l][1] + G1Dx[k][1] * G1Dy[l][0];
			}
		
			// perform IFTA on the 2D array
			int NumAFphase = 30;
			IFTA_optimize(f2D, phase2D, g2D, gInputAmplitude, G2D, &FFTplan_Gg_2D, &FFTplan_gG_2D, Nx * Ny, f2Dmask, phase2Dmask, NumAF, 30);
		}
	}
	
	// the array G now holds the complex valued field 
	// just after the SLM that produces a good approximation to the signal 
	// amplitude f in the focal plane, so, we need to set the SLM phase 
	// pattern to the phase of G
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		// calculate the real and imaginary part of G
		double ImG, ReG;
		if (calc2D)
		{
			// 2D, we can simply copy them from the 2D array
			ImG = G2D[l * Nx + k][1];
			ReG = G2D[l * Nx + k][0];
		} 
		else
		{
			// 1D, we have to compute G as the product of the two 1D arrays
			ImG = G1Dx[k][0] * G1Dy[l][1] + G1Dx[k][1] * G1Dy[l][0];
			ReG = G1Dx[k][0] * G1Dy[l][0] - G1Dx[k][1] * G1Dy[l][1];
			
			// while we're at it, we might as well fill the 2D array too
			G2D[l * Nx + k][0] = ReG;
			G2D[l * Nx + k][1] = ImG;
		}
		
		// the phase is the arctangent of Im(G) / Re(G)
		//gSLMpixels[l * Nx + k] = convf * atan2(ImG, ReG);	
		gSLMphase[ l * Nx + k] = atan2(ImG, ReG);
	}
	
	// compute the diffraction efficiency, by propagating G/g? back to the input plane ...
	// (note that we have filled G2D, even if the calculation was done in 1D)
	//fftw_execute(FFTplan_gG_2D);
		
	// ... and compute the efficiency eta from this back-propagated field
	double eta = calcEta(NULL, G2D, Nx * Ny);
	
	// TODO: write output to textbox!
	
	
	// TODO: destroy FFT plans?
	
	// propagate the light field to the output plane
	fftw_execute(FFTplan_Gg_2D);
	
	// compute the total energy in the output plane, and the energy in the window area
	double Etotal = 0.0;
	double Ewindow = 0.0;
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		Etotal += g2D[k + l * Nx][0] * g2D[k + l * Nx][0] + g2D[k + l * Nx][1] * g2D[k + l * Nx][1];
		
		if (gSignal[k + l * Nx] > 0.0)
			Ewindow += g2D[k + l * Nx][0] * g2D[k + l * Nx][0] + g2D[k + l * Nx][1] * g2D[k + l * Nx][1];
	}
	
	// output diffraction efficiencies
	char msg[1024];
	sprintf(msg, "%i x %i spot array, eta = %.1f [%.1f]", numxspots, numyspots, 100.0 * Ewindow / Etotal, Etotal);
	InsertTextBoxLine (OutputPanel, OutputTextBox, -1, msg);
	
	plotFFTField(g2D, Nx, Ny, gDebugPanel, gDebugCanvas1);
	plotFFTFieldPhase(g2D, Nx, Ny, gDebugPanel, gDebugCanvas2);
	
	// set the pattern indicator
	gCurrentPattern = SLM_BEAMSPLIT_IFTA;
	
	// store the theoretical spot intensities
	gIFTASimulatedSpotIntensities = (double*) realloc(gIFTASimulatedSpotIntensities, numxspots * numyspots * sizeof(double));
	for (int k = 0; k < numxspots; k++)
	for (int l = 0; l < numyspots; l++)
	{	
		int xc = xspot1 + k * (spotxspacing + 1);
		int yc = yspot1 + l * (spotyspacing + 1);
		gIFTASimulatedSpotIntensities[k + l * numxspots] = g2D[xc + yc * Nx][0] * g2D[xc + yc * Nx][0] + g2D[xc + yc * Nx][1] * g2D[xc + yc * Nx][1];
	}
	
	// store the diffraction efficiency
	gIFTASimulatedDiffractionEfficiency = Ewindow / Etotal;
	
	//Use f2D as the unresampled desired signal
	gUnresampledSignal = (double*) realloc(gUnresampledSignal, Nx * Ny * sizeof(double));
	memcpy(gUnresampledSignal, f2D, Nx * Ny * sizeof(double));
	gUnresampledSignalXsize = Nx;
	gUnresampledSignalYsize = Ny;
	
	// free all the intermediate arrays
	fftw_free(g1Dx);
	fftw_free(g1Dy);
	fftw_free(G1Dx);
	fftw_free(G1Dy);
	fftw_free(g2D);
	fftw_free(G2D);
	fftw_free(G1Dy_opt);
	fftw_free(G1Dx_opt);
	free(f1Dx);
	free(f1Dy);
	free(f2D);
	free(f1Dxmask);
	free(f1Dymask);
	free(f2Dmask);
}



//==============================================================================
// SLM beam shaping functions


/// HIFN calculates the phase required for geometric beam shaping of a
/// 	 separabale input (ii = iix * iiy) and output (is = isx * isy) intensity (intensity!)
double* SLM_calcBeamShapingPhase(
	double* iix, double *iiy, int Nxi, int Nyi, double ixmin, double ixmax, double iymin, double iymax, 
	double* isx, double *isy, int Nxs, int Nys, double sxmin, double sxmax, double symin, double symax, int center)
{
	// allocate memory for the output, i.e. the phase to be set at the SLM to shape the beam
	double* iphase = (double*) malloc(Nxi * Nyi * sizeof(double));
	
	// before we start, we first make sure that the intensities integrate to 1
	SLM_normalise(isx, Nxs, sxmin, sxmax);
	SLM_normalise(isy, Nys, symin, symax);
	SLM_normalise(iix, Nxi, ixmin, ixmax);
	SLM_normalise(iiy, Nyi, iymin, iymax);
	
	// compute the coordinate mapping functions gx and gy
	double* gx = SLM_calcCoordinateMapping(iix, isx, ixmin, ixmax, Nxi, sxmin, sxmax, Nxs);
	double* gy = SLM_calcCoordinateMapping(iiy, isy, iymin, iymax, Nyi, symin, symax, Nys);

	// calculate partial integrals of gx and gy
	double* gxI = SLM_cumulativeTrapz(gx, ixmin, ixmax, Nxi);
	double* gyI = SLM_cumulativeTrapz(gy, iymin, iymax, Nyi);
	
	// for each of the grid points on the SLM, calculate the phase as the
	// (interpolated) integrals gxI, and gyI
	for (int k = 0; k < Nxi; k++)
	for (int l = 0; l < Nyi; l++)
	{
		// the phase at SLM gridpoint (k, l) is the sum of the
		// interpolated values of gxI and gyI
		iphase[l * Nxi + k] = 2 * PI * (gxI[k] + gyI[l]) / (gFocalLength * gWavelength);
		
		// do we want center the pattern?
		if (center)
		{
			// add a linear phase to center the signal area at the origin
			double htrans = -(sxmax - sxmin) / 2.0;
			double vtrans = -(symax - symin) / 2.0;
			iphase[l * Nxi + k] += 2 * PI * ( ((double) k) / ((double) Nxi - 1) ) * (htrans / gFocalUnitX)
								+  2 * PI * ( ((double) l) / ((double) Nyi - 1)) * (vtrans / gFocalUnitY);
		}
	}
	
	// clean up intermediate arrays
	free(gx);
	free(gy);
	free(gxI);
	free(gyI);
	
	// return the calculated phase
	return iphase;
}


/// HIFN calculates the phase required for geometric beam shaping of the intensity 
///      incident onto the SLM (gInputAmplitude) to output intensity (is = isx * isy)
void SLM_generateBeamShape(double* isx, double* isy, int Nxs, int Nys, 
	double sxmin, double sxmax, double symin, double symax, int center)
{
	// the input area is given by that of the SLM
	double ixmin = 0.0;
	double ixmax = LxSLM;
	double iymin = 0.0;
	double iymax = LySLM;
	
	// create arrays of the input intensity, from the square root which is stored in the 
	// global gInputAmplitude arrays
	double* iix = (double*) malloc(gXsize * sizeof(double));
	double* iiy = (double*) malloc(gYsize * sizeof(double));
	for(int k = 0; k < gXsize; k++)
		iix[k] = gInputAmplitudeX[k] * gInputAmplitudeX[k];
	for(int l = 0; l < gYsize; l++)
		iiy[l] = gInputAmplitudeY[l] * gInputAmplitudeY[l];
	
	// sxmin, sxmax, symin and symax hold the corner coordinates of the signal window
	// in units of meters, we convert this to pixels in the focal plane and store them
	gWindowX = sxmin / gFocalUnitX;
	gWindowY = symin / gFocalUnitY;
	gWindowXsize = (sxmax - sxmin) / gFocalUnitX;
	gWindowYsize = (symax - symin) / gFocalUnitY;
	
	// next we first resample the signal to the window dimensions
	double* tmpsignal = (double*) malloc(Nxs * Nys * sizeof(double));
	for (int k = 0; k < Nxs; k++)
	for (int l = 0; l < Nys; l++)
		tmpsignal[k + l * Nxs] = isx[k] * isy[l];
	double* tmpsignal2 = SLM_resampleBilinear(tmpsignal, Nxs, Nys, 
		gWindowXsize, gWindowYsize);
	free(tmpsignal);
	
	// now copy the window to the gSignal array (set the remainder to 0)
	memset(gSignal, 0, gXsize * gYsize * sizeof(double));
	for (int k = 0; k < gWindowXsize; k++)
	for (int l = 0; l < gWindowYsize; l++)
	{
		gSignal[k + gWindowX + (l + gWindowY) * gXsize] = 
			tmpsignal2[k + l * gWindowXsize];
	}
	free(tmpsignal2);
	tmpsignal2 = NULL;
	
	// calculate the phase of the input field
	double* iphase = SLM_calcBeamShapingPhase(iix, iiy, gXsize, gYsize, ixmin, ixmax, iymin, iymax, 
											  isx, isy, Nxs,    Nys,    sxmin, sxmax, symin, symax, center);
	
	// copy the phase to the SLM
	memcpy(gSLMphase, iphase, gXsize * gYsize * sizeof(double));
	
	// free arrays we will no longer use
	free(iphase);
	free(iix);
	free(iiy);
	
	// set the pattern indicator
	gCurrentPattern = SLM_BEAMSHAPE_STD;
}


/// HIFN Generates the phase pattern for an arbitrary signal intensity ioriginal
void SLM_generatePhase(double* ioriginal, int Nxo, int Nyo, int WindowX, int WindowY, 
	int WindowXsize, int WindowYsize, int numPF, int numAF, double SigAmp, int keepcurrentphase, int softop, int MRAF, int storeUnresampledSignal)
{
	
	//CVIProfSetCurrentThreadProfiling (1);
	
	// store the original, unresampled signal
	if (storeUnresampledSignal)
	{
		gUnresampledSignal = (double*) realloc(gUnresampledSignal, Nxo * Nyo * sizeof(double));
		memcpy(gUnresampledSignal, ioriginal, Nxo * Nyo * sizeof(double));
		gUnresampledSignalXsize = Nxo;
		gUnresampledSignalYsize = Nyo;
	}
	
	// store the window boundaries of the signal
	gWindowX = WindowX;
	gWindowY = WindowY;
	gWindowXsize = WindowXsize;
	gWindowYsize = WindowYsize;
	
	// before we can start the IFTA, we have to resample the signal such that it fits in the
	// signal window in the focal plane
	
	// resample the signal to fit the signal window
	double* isignalw = SLM_resampleBitmap(ioriginal, Nxo, Nyo, WindowXsize, WindowYsize);		
	
	// create a mask for indicating the signal area
	unsigned char* signalmask = (unsigned char*) calloc(gXsize * gYsize, sizeof(unsigned char));
	
	// now we create the total signal as an array of zeros with the windowed signal embedded in it.
	// simultaneously, we set the contents of the gSignal array
	double* isignal = (double*) calloc(gXsize * gYsize, sizeof(double));
	memset(gSignal, 0, gXsize * gYsize * sizeof(double));
	for (int k = 0; k < WindowXsize; k++)
	for (int l = 0; l < WindowYsize; l++)
	{
		// write the windowed signal, and take the square root (!)
		isignal[   (l + WindowY) * gXsize + (k + WindowX)] = sqrt(isignalw[l * WindowXsize + k]);
		signalmask[(l + WindowY) * gXsize + (k + WindowX)] = 1;
		
		// copy the intensity to the gSignal array (so this is the square of isignal!)
		gSignal[   (l + WindowY) * gXsize + (k + WindowX)] = (isignalw[l * WindowXsize + k]);
	}
	
	// normalise the windowed signal intensity (!) to 1
	SLM_normalise_sq_2D(isignal, gXsize, gYsize, 0.0, gXsize * SLM_getFocalUnitX(), 0.0, gYsize * SLM_getFocalUnitY());
	
	// multiply the signal by the SigAmp factor, such that the total intensity in isignal is now 'SigAmp' (as opposed to '1')
	for (int k = 0; k < gXsize * gYsize; k++)
		isignal[k] = isignal[k] * SigAmp;
	
	// check if we want to use the current SLM phase pattern as an initial guess,
	// or whether we have to assemble a new initial guess
	if (keepcurrentphase)
	{
		// fill the G-array	from the gSLMphase variable
		SLM_setG(gInputAmplitude, gSLMphase, gXsize * gYsize);
	}
	else
	{
		// we want to assemble a new initial guess, which is done using an inverse beam shaping problem
		
		// create a separable intensity pattern ix * iy that approximates the signal
		// note that the beam shaping works with intensities, so we apply SepOp to the ioriginal array, which is an intensity
		double *isx = (double*) calloc(gXsize, sizeof(double));
		double *isy = (double*) calloc(gYsize, sizeof(double));
		double *tmpisignal = SLM_resampleBitmap(ioriginal, Nxo, Nyo, gXsize, gYsize);
		SLM_SepOp(tmpisignal, gXsize, gYsize, &isx, &isy);
		free(tmpisignal);
		tmpisignal = NULL;
	
		// create arrays for the input intensity
		double *iix = (double*) calloc(gXsize, sizeof(double));
		double *iiy = (double*) calloc(gYsize, sizeof(double));
	
		// fill the input intensity arrays with the SQUARE of gInputAmplitude 
		for (int k = 0; k < gXsize; k++)
			iix[k] = gInputAmplitudeX[k] * gInputAmplitudeX[k];
		for (int l = 0; l < gYsize; l++)
			iiy[l] = gInputAmplitudeY[l] * gInputAmplitudeY[l];
		
		// compute the physical positions of the corner coordinates of the window
		double sxmin = gWindowX * SLM_getFocalUnitX();
		double sxmax = sxmin + gWindowXsize * SLM_getFocalUnitX();
		double symin = gWindowY * SLM_getFocalUnitY();
		double symax = symin + gWindowYsize * SLM_getFocalUnitY();
	
		// apply the geometric beam shaping now INVERSELY, that is, we find the phase pattern for the
		// (separable approximation of) the signal intensity, that matches the input intensity as closely as possible
		// when propagated backwards
		double* iphase = SLM_calcBeamShapingPhase(isx, isy, gXsize, gYsize, sxmin, sxmax, symin, symax, 
								 iix, iiy, gXsize, gYsize, 0.0, LxSLM, 0.0, LySLM, 0);
		
		//plotField(iphase, gXsize, gYsize, gDebugPanel, gDebugCanvas2);
	
		// we have to correct the sign of iphase, as it is an inverse beam shaping problem that was solved
		for (int k = 0; k < gXsize * gYsize; k++)
			iphase[k] *= -1.0;
	
		// iphase now holds a suitable initial phase for the signal (!) field
		
		// resample the initial signal phase to fit the signal window
		double* restrict iphasew  = SLM_reduceBitmap(iphase, gXsize, gYsize, WindowXsize, WindowYsize);
	
		// now we write the windowed signal phase over the old signal phase
		memset(iphase,  0, gXsize * gYsize * sizeof(double));
		for (int k = 0; k < WindowXsize; k++)
		for (int l = 0; l < WindowYsize; l++)
			iphase[(l + WindowY) * gXsize + (k + WindowX)] = iphasew[l * WindowXsize + k];
		
		// the phase is ready, we can now continue with setting up the IFTA
	
		// compose the signal (g) from its amplitude and phase
		for (int k = 0; k < gXsize; k++)
		for (int l = 0; l < gYsize; l++)
		{
			gFFTout[l * gXsize + k][0] = isignal[l * gXsize + k] * cos(iphase[l * gXsize + k]);
			gFFTout[l * gXsize + k][1] = isignal[l * gXsize + k] * sin(iphase[l * gXsize + k]);
		}
	
		// propagate the signal (g) backwards to obtain G
		fftw_execute(gFFTplan_gG);
		
		// clean up arrays we no longer need
		free(iix);
		free(iiy);
		free(isx);
		free(isy);
		free(iphase);
		free(iphasew);
	}
	
	// the G-array now holds our initial guess, all set, we can start with the IFTA!
	
	// DEBUG: rms values during the iterations
	double* rms = (double*) malloc((numPF + numAF) * sizeof(double));
	
	// first, perform the Phase Freedom (PF) iterations
	IFTA_optimize_rms(isignal, NULL, gFFTout, gInputAmplitude, gFFTin, 
		&gFFTplan_Gg, &gFFTplan_gG, gXsize * gYsize, NULL, NULL, numPF, softop, MRAF, rms, gXsize, gYsize, WindowX, WindowX + WindowXsize, WindowY, WindowY + WindowYsize);
	
	// next, perform the Amplitude Freedom (AF) iterations (using the signal mask)
	IFTA_optimize_rms(isignal, NULL, gFFTout, gInputAmplitude, gFFTin, 
		&gFFTplan_Gg, &gFFTplan_gG, gXsize * gYsize, signalmask, NULL, numAF, softop, MRAF, rms + numPF, gXsize, gYsize, WindowX, WindowX + WindowXsize, WindowY, WindowY + WindowYsize);	 
	
	// DEBUG plot the rms
	//DeleteGraphPlot (gDebugPanel, gDebugGraph1, -1, VAL_DELAYED_DRAW);
	//PlotY(gDebugPanel, gDebugGraph1, rms, numPF + numAF, VAL_DOUBLE,
	//				VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, 0x00FF0000);
	//RefreshGraph (gDebugPanel, gDebugGraph1);
			
	//GeneticOptimize(numPF, DebugPanel, dbc1, dbc2, dbc3, dbc4, dbc5, dbc6, debuggraph1);
	
	// conversion factor from 2 * pi to 256
	double convf = 256.0 / (2 * PI); 
	
	// the array G now holds the complex valued field
	// just after the SLM that produces a good approximation to the signal 
	// amplitude f in the focal plane, so, we need to set the SLM phase 
	// pattern to the phase of G
	for (int k = 0; k < gXsize; k++)
	for (int l = 0; l < gYsize; l++)
	{
		double ImG, ReG;
		
		// calculate the real and imaginary part of G
		ImG = gFFTin[l * gXsize + k][1];
		ReG = gFFTin[l * gXsize + k][0];
		
		// the phase is the arctangent of Im(G) / Re(G)
		gSLMphase[ l * gXsize + k] = atan2(ImG, ReG);
		//gSLMpixels[l * gXsize + k] = convf * gSLMphase[l * gXsize + k];	
	}
	
	// clean up
	free(isignal);
	free(isignalw);
	free(signalmask);
	free(rms);
	
	// set the pattern indicator
	gCurrentPattern = SLM_BEAMSHAPE_ARB;
}



/// HIFN get the original, unresampled signal (in the case of beam shaping from file)
double* SLM_getUnresampledSignal(void)
{
	return gUnresampledSignal;
}

/*
double SLM_UpdateSpotarray(gSLMSignal, pnlSLMpixels, pnlSimPanel, SimPanel_CANVAS, panel, OutputTextBox)
{

*/

	
