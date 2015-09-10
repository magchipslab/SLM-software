//==============================================================================
//
// Title:       SLM.h
// Purpose:     Library of functions that generate patterns to be displayed on 
//              a Spatial Light Modulator (SLM), and simulate the resulting 
//              diffraction pattern intensity.
//
// Created on:  9-9-2011 at 19:46:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================


// includes needed for the MATfile type
#include "libmat/mat.h" 
#include "libmat/matrix.h"

// the physical dimensions of the SLM (in m)
extern const double LxSLM;
extern const double LySLM;

// possible pattern modes on the SLM
#define SLM_BEAMSHAPE_STD			0
#define SLM_BEAMSHAPE_ARB			1
#define SLM_BEAMSPLIT_PROJECT		2
#define SLM_BEAMSPLIT_IFTA			3
#define SLM_BEAMSPLIT_PHASEGRATING  4
#define SLM_SHACKHARTMANN			5
#define SLM_BEAMSHAPE_CHAR			6

/// HIFN Get the current pattern on the SLM
int SLM_getCurrentPattern(void);

/// HIFN get the SLM x resolution
int SLM_getXres(void);

/// HIFN get the SLM y resolution
int SLM_getYres(void);

/// HIFN get the SLM canvas x resolution
int SLM_getCanvasXres(void);

/// HIFN get the SLM canvas y resolution
int SLM_getCanvasYres(void);



/// HIFN get the current base phase pattern of the SLM
double* SLM_getPhase(void);

/// HIFN get the current zoom of the SLM pattern
double SLM_getZoom(void);

/// HIFN get the current amplitude modulation
double SLM_getAmplitudeModulation(void);

/// HIFN copies the SLM phase from a specified memory location
void SLM_setPhase(double* src);

/// HIFN clear the SLM base phase pattern (i.e. not the lens, translation, corrections.. etc.)
void SLM_Clear(void);

/// HIFN (Re)initialise the SLM with a certain size; i.e. allocate all memory buffers etc.
void SLM_initialise(int Xsize, int Ysize, int SimXsize, int SimYsize, int SubSampleFactor);

/// HIFN Render the pixels of the SLM onto a canvas
void SLM_update(int PanelID, int CanvasID, int SimPanelID, int SimCanvasID, int pattern_updated);

/// HIFN simulates the intensity pattern in the focal plane of a lens, for the current SLM pattern,
///      and renders it onto a canvas
void SLM_simulateIntensity(int Panel, int Canvas);

/// HIFN adjusts the horizontal translation of the intensity pattern
void SLM_setHorizTrans(double htrans);

/// HIFN adjusts the horizontal translation of the intensity pattern
void SLM_setVertTrans(double vtrans);

/// HIFN sets the SLM lens phase pattern
void SLM_setLensPhase(double xstrength, double ystrength);

// HIFN sets the SLM aberration correction phase pattern
void SLM_setAberrationCorrection(double* gZernikeCoefficients, int Nz, int mx, int my);

/// HIFN toggle the simulation
void SLM_toggleSim(unsigned int sim, unsigned int simphase, unsigned int simpixel);

/// HIFN set the simulation saturation, relative to the maximum intensity
void SLM_setSimSaturation(double simsat);

/// HIFN set the simulation lens focal length (in m)
void SLM_setLensFocalLength(double f);

/// HIFN set the magnification of the focal length
void SLM_setMagnification(double h);

/// HIFN loads a factory correction
void LoadFactoryCorrection(char FactoryCorrectionfilename[], int FactoryCorrectionCheck);

/// HIFN loads a SH surface correction
void LoadSHcorrection(char filename[], int check);

/// HIFN get the lens focal length (in m)
double SLM_getLensFocalLength(void);

/// HIFN set the simulation wavelength (in m)
void SLM_setWavelength(double lambda);
		  
/// HIFN get the wavelength (in m)
double SLM_getWavelength(void);

/// HIFN set the Gaussian input intensity (parameters are all lengths, in m)
void SLM_setInputIntensityGaussian(double x, double y, double sigmax, double sigmay, int subsamplefactor);

/// HIFN set an annular Gaussian input intensity (parameters are all lengths, in m)
void SLM_setInputIntensityAnnular(double x, double y, double sigmax, double sigmay, int subsamplefactor);

/// HIFN set an arbitrary input intensity (parameters are all lengths, in m)
void SLM_setInputIntensityFromCameraData(double x, double y, double sigmax, double sigmay, int subsamplefactor);

/// HIFN sets the SLM pattern to a spot projector
void SLM_setSpotPattern(unsigned int Nx, unsigned int Ny, double xspacing, double yspacing, 
	double randomampl, double spotlens, unsigned int linspotlens);

/// HIFN adjusts the SLM spot mask that controls the number of pixels contributing to a single spot
void SLM_setSpotMask(unsigned int Nx, unsigned int Ny, double size);

/// HIFN sets the SLM pattern to a phase grating
void SLM_setPhaseGrating(double a1, double b1, double a2, double b2, double angle);

/// HIFN changes the mode of the phase grating to a series of Gaussian peaks
void SLM_togglePhaseGratingGaussianPeaks(int gptoggle);

/// HIFN set the SLM zoom factor
void SLM_setZoom(double zoom);

/// HIFN generate a beam splitting phase pattern
void SLM_generateSpotArray(int Panel, int Canvas, int sSimPanel, int sSimCanvas, int OutputPanel, int OutputTextBox,
	                       int numxspots, int numyspots, int spotxspacing, int spotyspacing, int sigxoffset, 
						   int sigyoffset, int spotxoffset, int spotyoffset, int constrainphase, double phasestep);

/// HIFN generate a geometric beam shaping phase pattern
void SLM_generateBeamShape(double* isx, double *isy, int Nxs, int Nys, 
	double sxmin, double sxmax, double symin, double symax, int center);

/// HIFN generate phase pattern for arbitrary 2D input (typically from a picture)
void SLM_generatePhase(double* ioriginal, int Nxo, int Nyo, int WindowX, int WindowY, 
	int WindowXsize, int WindowYsize, int numPF, int numAF, double SigAmp, int keepcurrentphase, int softop, int MRAF, int storeunresampledisgnal);

/// HIFN Finds the pixel coordinates and spacing of a spot pattern in an image
void SLM_findSpotPattern(double* image, int Nxpixels, int Nypixels, int Nxspots, int Nyspots, double* xpos, double* ypos, double* xspacing, double* yspacing);

/// HIFN Initialise the simulation
void SLM_initialiseSimulation(int Xsize, int Ysize);

/// HIFN sets the simulation zoom factor
void SLM_setSimZoom(double simzoom);

/// HIFN sets the simulation's grid spacing
void SLM_setSimGridSpacing(double simgridspacing);

/// HIFN toggles the grid lines
void SLM_setSimShowGrid(int simshowgrid);

/// HIFN sets the amplitude modulation
void SLM_setAmplitudeModulation(double amplmod);

/// HIFN set the distance SLM - Lens, for the simulation with Helmholtz propagation
void SLM_setSimZPos(double zpos);

/// HIFN sets the bias phase (0..255)
void SLM_setBias(int bias);

/// HIFN get the bias phase
int SLM_getBias();

/// HIFN toggles whether the simulation supersamples or not
void SLM_setSimSuperSampling(int ssample);

/// HIFN toggles whether the simulation should use Helmholtz propagation
void SLM_setSimHelmholtz(int helmholtz);

/// HIFN sets the coordinates of the focusing lens
void SLM_setLensCoordinates(double lensX, double lensY);

/// HIFN writes the basic SLM settings to a mat-file
void SLM_WriteSettingsToFile(MATFile *pmat);

/// HIFN Convenience function to write a scalar double value to a matlab .mat file
void writeMatDoubleScalar(MATFile *matfile, char varname[], double varvalue);

/// HIFN Convenience function to write an array of double values to a matlab .mat file
void writeMatDoubleArray(MATFile *matfile, const char varname[], double* data, int M, int N);

// /// HIFN Write an array of integers as an M x N matrix to a matlab .mat file
//void writeMatIntArray(MATFile *matfile, const char varname[], int* data, int M, int N);

/// HIFN Generates an array of logarithmically spaced values
double* logspace(double xmin, double xmax, double Nx);

/// HIFN writes an array of integers to a mat-file
void writeMatIntArray(MATFile *matfile, char varname[], int* varvalues, int Nv);

/// HIFN Write an array of unsigned chars as an M x N matrix to a matlab .mat file
void writeMatUnsignedCharArray(MATFile *matfile, const char varname[], unsigned char* data, int M, int N);

/// HIFN Convenience function to read a scalar double value from a matlab .mat file
double readMatDoubleScalar(MATFile *matfile, char varname[]);

/// HIFN Read an array of doubles as an M x N matrix from a matlab .mat file
double* readMatDoubleArray(MATFile *matfile, const char varname[], int *M, int *N);

/// HIFN computes the focal unit, and updates gFocalUnitX and gFocalUnitY
void SLM_updateFocalUnit(void);

/// HIFN get the focal unit in the X direction (in meters)
double SLM_getFocalUnitX(void);

/// HIFN get the focal unit in the Y direction (in meters)
double SLM_getFocalUnitY(void);

/// HIFN get the signal window area position and dimensions
void SLM_getSignalWindowParameters(int *wx, int *wy, int *wxsize, int* wysize);

/// HIFN get the signal window contents
double* SLM_getSignalWindow(void);

/// HIFN get the original, unresampled signal (in the case of beam shaping from file)
double* SLM_getUnresampledSignal(void);

/// HIFN get the physical size of the SLM in the x-direction (in m)
double SLM_getLxSLM(void);

/// HIFN get the physical size of the SLM in the y-direction (in m)
double SLM_getLySLM(void);

/// HIFN Projects a phase field onto a subspace spanned by the lowest few Zernike polynomials
///      Array 'cj' holds the coefficients of the Zernike polynomials, after the projection
void ProjectPhaseOntoZernikeSpace(double* phase, double* cj, int Nx, int Ny, int NZP);

/// HIFN Unwraps a phase, i.e. removes jumps of 2 pi, of a 2D array
void SLM_unwrapPhase(double *phase, int Nx, int Ny);

/// HIFN Performs a weighted average of two phase angles, taking the branch cut at 2 pi into account
double SLM_averagePhase(double t1, double t2, double a);

/// HIFN Initialise the data structures for the phase error characterisation
void SLM_initPhaseErrorCharacterisation(int NxPE, int NyPE, int NZP, double noisecutoff, int prealloc, 
	int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6);

/// HIFN Performs 'N' IFTA iterations, where the projection step in the SLM plane consists of 
///      projecting onto a subspace of a fixed SLM pattern + a linear combination of the 
///      lowest few Zernike rectangle polynomials.
///      Function returns the coefficients of the Zernike polynomials.
double* SLM_IFTA_PE(double* cj0, int N, int UseSoftOp, int NPmax, int prjinterval, double beta, double* RefAmpMatch, 
	int debugpanel, int dbc1, int dbc2, int dbc3, int dbc4, int dbc5, int dbc6, int debuggraph, int debuggraph2, int debuggraph3);

/// HIFN get the coefficients of rectangular Zernike polynomials describing the aberration correction of the SLM
double* SLM_getZernikeAberrationCorrection();

/// HIFN get the SLM base phase pattern used in the phase error characterisation
double* SLM_getCharacterisationBasePhase();

/// HIFN get the reference amplitude and its parameters used in the phase error characterisation
void SLM_getCharacterisationReferenceAmplitude(double** refampl, int* wx, int* wy, int* wxpos, int* wypos);

/// HIFN sets a pattern used for Shack-Hartmann aberration correction
void SLM_SetShackHartmannPattern(double spotx, double spoty, double spotdiameter, double strayx, double strayy);

/// HFIN store the raw data for the input intensity, when a measurement from a camera is used
void SLM_setInputIntensityCameraData(unsigned char* camdata, int camxres, int camyres, int camxmin, int camymin, int camxsize, int camysize);

/// HIFN assigns the debug plotting controls
void SLM_setDebugPlottingControls(int dbgpanel, int dbgc1, int dbgc2, int dbgc3, int dbgc4, int dbgc5, int dbgc6, int dbgg1, int dbgg2, int dbgg3);



													   
