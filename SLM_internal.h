//==============================================================================
//
// Title:       SLM_internal.h
// Purpose:     Internal header file for all SLM source code files, contains
//              variables and functions that should not be visible outside
//              the SLM code. Externally visible variables and functions are
//              included in SLM.h
//
// Created on:  29-10-2011 at 19:15:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#include <ansi_c.h>
#include <cvirte.h>     
#include <userint.h>
#include <math.h>
#include "toolbox.h"
#include "FFTW/fftw3.h"
#include "SLM.h"


#ifndef PI
#define PI 3.14159265358979323846 
#endif

//==============================================================================
// Global variables that are shared by multiple library source files

// Number of pixels in each dimension on the SLM pixel canvas
extern int gXsize;
extern int gYsize;

// the current pattern mode
extern int gCurrentPattern;

// the SLM base phase pattern
extern double* gSLMphase;

// the total SLM phase pattern, without having taken the modulus
extern double* gSLMtotal;

// the lens focal length (in m)
extern double gFocalLength;

// the wavelength (in m)
extern double gWavelength;

// the bias phase offset
extern unsigned char gBias;

// the current settings for translation, aberration and lens phase (focal length in m) and lens position
extern double gLensXphase, gLensYphase, gHorizTrans, gVertTrans, gAberration, gAberrationAnisotropy, gLensX, gLensY;

// SLM pixel data in which we can directly write
extern unsigned char* gSLMPixels;

// the SLM aberration correction phase pattern, that will be SUBTRACTED from the global SLM phase
extern double* gSLMaberration;

// the (Gaussian) input amplitude (NOTE: this is the square root of the intensity)
extern double* gInputAmplitude;

// the separable version of the input amplitude pattern 
extern double* gInputAmplitudeX;
extern double* gInputAmplitudeY;

// parameters of the Gaussian input intensity
extern double gIxcenter, gIycenter, gIxsigma, gIysigma;

// the desired signal intensity (set by most pattern generators, and used for feedback purposes)
// each value corresponds to a focal unit pixel in the SLM plane
extern double* gSignal;

// in case of beam shaping from file, we also store the original, unresampled signal (used in the feedback procedure)
//extern double* gUnresampledSignal;
extern int gUnresampledSignalXsize;
extern int gUnresampledSignalYsize;

// variables indicating the window position and size within the gSignal array
extern int gWindowX, gWindowY, gWindowXsize, gWindowYsize;

// the focal unit, i.e. the elementary pixel size in the focal plane
extern double gFocalUnitX, gFocalUnitY;

// the zoom factor of the pattern 
extern double gZoom;

// subsampling factor
extern int gSubSampleFactor;

// variables for the internal FFTs used for IFTA
extern fftw_complex *gFFTin, *gFFTout;
extern fftw_plan gFFTplan_Gg, gFFTplan_gG;

// the amplitude modulation of the SLM
extern double gAmplitudeModulation;

// the simulated spot intensities for a spot array
extern double* gIFTASimulatedSpotIntensities;

// the simulated spot array diffraction efficiency
extern double gIFTASimulatedDiffractionEfficiency;

// debug plotting controls
extern int gDebugPanel, gDebugCanvas1, gDebugCanvas2, gDebugCanvas3, gDebugCanvas4,  gDebugCanvas5, gDebugCanvas6, gDebugGraph1, gDebugGraph2, gDebugGraph3;  


// =================================================================================
// Support Functions (implemented in SLM_support.c)

/// HIFN initialises to complex valued array G representing the light field
///      in the SLM plane to given values
void SLM_setG(double* ampl, double* phase, int NG);

/// HIFN applies the constraints of unit amplitude in the SLM plane
void SLM_applyConstraints(fftw_complex* G, int NG, double* constraint);

/// HIFN performs iterative fourier transform algorithm (IFTA) to find the optimal phase pattern
void IFTA_optimize(double* f, double* phasef, fftw_complex* g, double* constraintG, fftw_complex* G, fftw_plan* FFTplan_Gg, 
						  fftw_plan* FFTplan_gG, int Nf, unsigned char* signalmask, unsigned char* phasemask,
						  int Nit, int UseSoftOp);

/// HIFN performs iterative fourier transform algorithm (IFTA) to find the optimal phase pattern in the SLM plane, and stores
///      the rms deviation of the signal with the desired signal at each iteration
///      THIS FUNCTION IS ONLY USED TO GENERATE SOME FIGURES FOR THE PAPER (ugly and useless)
void IFTA_optimize_rms(double* f, double* phasef, fftw_complex* g, double* constraintG, fftw_complex* G, fftw_plan* FFTplan_Gg, 
						  fftw_plan* FFTplan_gG, int Nf, unsigned char* signalmask, unsigned char* phasemask,
						  int Nit, int UseSoftOp, int MRAF, double* rms, int Nx, int Ny, int sigxmin, int sigxmax, int sigymin, int sigymax);

/// HIFN computes the diffraction efficiency eta, as defined in the book by Turunen and Wyrowski
double calcEta(double* Ui, fftw_complex* Gh, int Nf);

/// HIFN calculates the integral over the intensity up to each point in the grid specified by xmin, xmax
double* SLM_cumulativeTrapz(double* intensity, double xmin, double xmax, int N);

/// HIFN inverts the intensity integral I(x), i.e. finds the corresponding x at which I(x) = I
double SLM_invertIntensityIntegral(double* Ix, double I, double xmin, double xmax, int N);

/// HIFN calculates the coordinate mapping function g
double* SLM_calcCoordinateMapping(double* ii, double* is, 
	double ixmin, double ixmax, int Ni,
	double sxmin, double sxmax, int Ns);

// normalises the integral of i (over the range xmin - xmax) to 1
void SLM_normalise(double* i, int N, double xmin, double xmax);

/// HIFN Approximates an arbitrary 2D intensity pattern by a separable one, by applying the SepOp operator (see Wyrowski book)
void SLM_SepOp(double* intensity, int Nx, int Ny, double** ix, double** iy);

/// HIFN resamples a 2D input intensity given on one grid, onto another grid, by means of bilinear interpolation
double* SLM_resampleBilinear(double* input, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN resamples a 2D input intensity given on one grid, onto another grid, by means of bilinear interpolation
void SLM_resampleBilinearInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo);             

/// HIFN Resample a phase, which takes some extra effort to avoid trouble at the branch cut at +/- pi.
double* SLM_resamplePhase(double* input, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN Resample a phase, which takes some extra effort to avoid trouble at the branch cut at +/- pi.
void SLM_resamplePhaseInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN function for resampling a 2D input on one grid, onto another, assuming the output grid is SMALLER than the input
double* SLM_reduceBitmap(double* input, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN Shrinks a 2D pixel array to a smaller number of pixels
void SLM_reduceBitmapInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN General resampling method, uses bilinear or reduce when appropriate
double* SLM_resampleBitmap(double* input, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN General resampling method, uses bilinear or reduce when appropriate, output array is already allocated
void SLM_resampleBitmapInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo);

/// HIFN the amplitude modulation of the SLM as a function of the phase
double SLM_AmplitudeModulation(double phase);

/// HIFN debug function to plot an array of doubles
void plotField(double* f, int Nx, int Ny, int Panel, int Canvas);

/// HIFN Debug function to plot the absolute value of a complex-valued array
void plotFFTField(fftw_complex* fftfield, int Nx, int Ny, int Panel, int Canvas);

/// HIFN Debug function to plot the phase of a complex-valued array
void plotFFTFieldPhase(fftw_complex* fftfield, int Nx, int Ny, int Panel, int Canvas);

/// HIFN generates a nice RGB colormap
int* SLM_CreateColorMap(int satmax);

/// HIFN matches two images with a cross correlation
void SLM_crossCorrelate(double* im1, double* im2, int Nx, int Ny, int *xpos, int *ypos);

/// HIFN compute the dot product of two arrays of doubles
//double DotProduct(double* A, double* B, int N);

/// HIFN calculate the total squared difference between two arrays of doubles
double calcSquaredDifference(double* f1, double *f2, int N);

/// HIFN calculate the quality of the intensity resulting from a given phase pattern
double PhaseQuality(unsigned char * restrict phase);

/// HIFN Calculates the RMS difference between the intensities represented by two arrays that represent amplitudes (so they will be squared),
///      one of which is complex valued, and the difference is only taken within a certain window
double calcRMSdiff(fftw_complex* f, double* v, int Nx, int Ny, int SigXmin, int SigXmax, int SigYmin, int SigYmax);


int GeneticOptimize(int NumEval, int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6, int debuggraph1);

// normalises the integral of i^2 (over the range xmin - xmax) to 1
void SLM_normalise_sq(double* i, int N, double xmin, double xmax);

// calculates the 2D integral of i (over the range xmin - xmax, ymin - ymax)
double SLM_integrate_2D(double* i, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax);

// normalises the 2D integral of i (over the range xmin - xmax, ymin - ymax) to 1
void SLM_normalise_2D(double* i, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax);

// normalises the 2D integral of i^2 (over the range xmin - xmax, ymin - ymax) to 1
void SLM_normalise_sq_2D(double* i, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax);

// normalises the 2D integral of |i|^2 for complex valued arrays (over the range xmin - xmax, ymin - ymax) to 1
void SLM_normalise_abssq_2D(fftw_complex* f, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax);


