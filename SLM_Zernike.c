//==============================================================================
//
// Title:       SLM_Zernike.c
// Purpose:     Library of rectangular Zernike polynomials,
//              see: V. N. Mahajan and G.-m. Dai, J. Opt. Soc. Am, 24, p2994 (2007)
//
// Created on:  02-04-2012 at 14:46:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#include <ansi_c.h>
#include <analysis.h>
#include "SLM_Zernike.h"


// array of precalculated Zernike polynomials
static double** gZernikePolynomials;

// dimensions of the precalculated Zernike polynomials
static int gPrecalcX = 0, gPrecalcY = 0;

// number of precalculated Zernike polynomials
static int gPrecalcNZ = 0;

void CalculateZernikeJ(double* ZJ, int j, int Nx, int Ny);

/// HIFN Precalculate arrays for the orthogonal Zernike rectangle polynomials
void PreCalculateZernikePolynomials(int Nx, int Ny, int NZP)
{
	// first check whether we have to free up the memory first
	if (((gPrecalcNZ > 0) && (gPrecalcNZ != NZP)) || (gPrecalcX != Nx) || (gPrecalcY != Ny))
	{
		// yes, first free the memory for the polynomials
		for (int j = 1; j <= gPrecalcNZ; j++)
			free(gZernikePolynomials[j - 1]);
		
		// next, free the main pointer itself
		free(gZernikePolynomials);
	}
	
	// store the current dimensions and nr. of polynomials
	gPrecalcX = Nx;
	gPrecalcY = Ny;
	gPrecalcNZ = NZP;
	
	if (NZP > 0)
	{
		// (re)allocate the array of pointers to Zernike polynomials
		gZernikePolynomials = (double**) malloc(NZP * sizeof(double*));
	
		// allocate and sample the individual Zernike polynomials
		for (int j = 1; j <= NZP; j++)
		{
			// allocate memory
			gZernikePolynomials[j - 1] = (double*) malloc(Nx * Ny * sizeof(double));	
		
			// sample the polynomial
			CalculateZernikeJ(gZernikePolynomials[j - 1], j, Nx, Ny);
		}
	}
}


/// HIFN Sample the j-th Zernike rectangle polynomial on a Nx-by-Ny grid covering the SLM rectangle
double* SampleZernikeJ(int j, int Nx, int Ny)
{
	// allocate memory to hold the array
	double* ZJ = (double*) malloc(Nx * Ny * sizeof(double));
	
	// fill the array with sampling points of the j-th polynomial
	FillArrayZernikeJ(ZJ, j, Nx, Ny);
	
	return ZJ;
}


/// HIFN Sample the j-th Zernike rectangle polynomial on a Nx-by-Ny grid covering the SLM rectangle,
///      filling a pre-allocated array for holding the sampling points, uses precalculated Zernike polynomials if available
void FillArrayZernikeJ(double* ZJ, int j, int Nx, int Ny)
{	
	// is there a precalculated polynomial available?
	if ((gPrecalcX == Nx) && (gPrecalcY == Ny) && (j <= gPrecalcNZ))
	{
		// yes, just copy the memory
		memcpy(ZJ, gZernikePolynomials[j - 1], Nx * Ny * sizeof(double));
	}
	else
	{
		// no, calculate the values on the spot
		CalculateZernikeJ(ZJ, j, Nx, Ny);
	}
}


/// HIFN Calculate the j-th Zernike rectangle polynomial on a Nx-by-Ny grid covering the SLM rectangle,
///      filling a pre-allocated array ZJ for holding the sampling points.
void CalculateZernikeJ(double* ZJ, int j, int Nx, int Ny)
{
	// get the rectangle sides
	double LZx = 2.0 * LxSLM / sqrt(LxSLM * LxSLM + LySLM * LySLM);
	double LZy = 2.0 * LySLM / sqrt(LxSLM * LxSLM + LySLM * LySLM);

	// get a function pointer to the j-th Zernike rectangle polynomial
	dFdd rzj = RectZernikePolynomial(j);

	// loop over all the grid points
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		// determine the x and y coordinates
		double x = k * LZx / ((double) (Nx - 1)) - 0.5 * LZx;
		double y = l * LZy / ((double) (Ny - 1)) - 0.5 * LZy; 
	
		// sample the J-th Zernike polynomial
		ZJ[k + l * Nx] = rzj(x, y);
	}
}
	

/// HIFN Projects a phase field onto a subspace spanned by the lowest few Zernike polynomials
///      Array 'cj' holds the coefficients of the Zernike polynomials, after the projection
void ProjectPhaseOntoZernikeSpace(double* phase, double* cj, int Nx, int Ny, int NZP)
{
	// make sure the cj-array can hold all the coefficients
	//cj = (double*) realloc(cj, NZP * sizeof(double));
	
	// allocate memory to temporarily hold a Zernike polynomial and the product of it with the phase
	double* ZJ = (double*) malloc(Nx * Ny * sizeof(double));
	
	// allocate temporary array holding the projected phase
	double* prjphase = (double*) calloc(Nx * Ny, sizeof(double));
	
	// check if there are precalculated Zernike polynomials available
	int precalc;
	if ((gPrecalcX == Nx) && (gPrecalcY == Ny) && (gPrecalcNZ == NZP))
		precalc = 1;
	else
		precalc = 0;

	// loop over all the Zernike polynomials
	for (int j = 1; j <= NZP; j++)
	{
		// calculate the pointwise multiplication of the phase with the j-th Zernike polynomial
		if (precalc)
		{
			for (int k = 0; k < Nx * Ny; k++)
				ZJ[k] = gZernikePolynomials[j - 1][k] * phase[k];	
		}
		else
		{
			// create a sampling of the j-th polynomial on the SLM rectangle
			FillArrayZernikeJ(ZJ, j, Nx, Ny);
		
			// perform pointwise multiplication with the phase array
			for (int k = 0; k < Nx * Ny; k++)
				ZJ[k] *= phase[k];
		}
		
		// complete the dot-product by integrating over the SLM rectangle
		cj[j - 1] = ZernikeRectangleIntegral(ZJ, Nx, Ny);
		
		// add the contribution of this Zernike polynomial to the projected phase
		if (precalc)
			for (int k = 0; k < Nx * Ny; k++)
				prjphase[k] += cj[j - 1] * gZernikePolynomials[j - 1][k];			
		else
			for (int k = 0; k < Nx * Ny; k++)
				prjphase[k] += cj[j - 1] * ZJ[k] / phase[k];
	}
	
	// copy the projected phase over the original phase
	memcpy(phase, prjphase, Nx * Ny * sizeof(double));
	
	// clean up the temporary arrays
	free(prjphase);
	free(ZJ);
}


/// HIFN Numerically integrates over the specific rectangle on which the Zernike polynomials are defined
double ZernikeRectangleIntegral(double* f, int Nx, int Ny)
{
	// allocate array for holding the x-integrals, for each value of y
	double* xintegrals = (double*) malloc(Ny * sizeof(double));
	
	// get the rectangle sides
	double LZx = 2.0 * LxSLM / sqrt(LxSLM * LxSLM + LySLM * LySLM);
	double LZy = 2.0 * LySLM / sqrt(LxSLM * LxSLM + LySLM * LySLM);

	// compute the x and y spacings
	double dx = (LZx) / ((double) (Nx - 1));
	double dy = (LZy) / ((double) (Ny - 1));
	
	// loop over all the y points
	for (int l = 0; l < Ny; l++)
	{
		// integrate over x
		NumericIntegration(&(f[l * Nx]), Nx, dx, 2, &(xintegrals[l]));
	}
	
	// now we integrate over y
	double result;
	NumericIntegration(xintegrals, Ny, dy, 2, &result);
	
	// clean up
	free(xintegrals);
	
	// return the result
	return result;
}
		
		

/// HIFN Evaluate the j-th Zernike Circle Polynomial (for j up to, and including, 91)
double ZernikeJ(int j, double x, double y)
{
	double Z;
	
	// evaluate the j-th Zernike circle polynomial, code generated by Mathematica (ZernikeCirclePolynomials.nb)
	switch (j)
	{
		case 1:
			 Z = 1;
			 break;
		case 2:
			 Z = 2 * x;
			 break;
		case 3:
			 Z = 2 * y;
			 break;
		case 4:
			 Z = -sqrt(3) + 2 * sqrt(3) * pow(x, 2) + 2 * sqrt(3) * pow(y, 2);
			 break;
		case 5:
			 Z = 2 * sqrt(6) * x * y;
			 break;
		case 6:
			 Z = sqrt(6) * pow(x, 2) - sqrt(6) * pow(y, 2);
			 break;
		case 7:
			 Z = 6 * sqrt(2) * pow(x, 2) * y + y * (-4 * sqrt(2) + 6 * sqrt(2) * pow(y, 2));
			 break;
		case 8:
			 Z = x * (-4 * sqrt(2) + 6 * sqrt(2) * pow(x, 2) + 6 * sqrt(2) * pow(y, 2));
			 break;
		case 9:
			 Z = 6 * sqrt(2) * pow(x, 2) * y - 2 * sqrt(2) * pow(y, 3);
			 break;
		case 10:
			 Z = x * (2 * sqrt(2) * pow(x, 2) - 6 * sqrt(2) * pow(y, 2));
			 break;
		case 11:
			 Z = sqrt(5) + pow(y, 2) * (-6 * sqrt(5) + 6 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (-6 * sqrt(5) + 6 * sqrt(5) * pow(x, 2) + 12 * sqrt(5) * pow(y, 2));
			 break;
		case 12:
			 Z = pow(x, 2) * (-3 * sqrt(10) + 4 * sqrt(10) * pow(x, 2)) + pow(y, 2) * (3 * sqrt(10) - 4 * sqrt(10) * pow(y, 2));
			 break;
		case 13:
			 Z = x * (8 * sqrt(10) * pow(x, 2) * y + y * (-6 * sqrt(10) + 8 * sqrt(10) * pow(y, 2)));
			 break;
		case 14:
			 Z = sqrt(10) * pow(y, 4) + pow(x, 2) * (sqrt(10) * pow(x, 2) - 6 * sqrt(10) * pow(y, 2));
			 break;
		case 15:
			 Z = x * (4 * sqrt(10) * pow(x, 2) * y - 4 * sqrt(10) * pow(y, 3));
			 break;
		case 16:
			 Z = x * (6 * sqrt(3) + pow(y, 2) * (-24 * sqrt(3) + 20 * sqrt(3) * pow(y, 2)) + pow(x, 2) * (-24 * sqrt(3) + 20 * sqrt(3) * pow(x, 2) + 40 * sqrt(3) * pow(y, 2)));
			 break;
		case 17:
			 Z = y * (6 * sqrt(3) + pow(y, 2) * (-24 * sqrt(3) + 20 * sqrt(3) * pow(y, 2))) + pow(x, 2) * (20 * sqrt(3) * pow(x, 2) * y + y * (-24 * sqrt(3) + 40 * sqrt(3) * pow(y, 2)));
			 break;
		case 18:
			 Z = x * (pow(y, 2) * (24 * sqrt(3) - 30 * sqrt(3) * pow(y, 2)) + pow(x, 2) * (-8 * sqrt(3) + 10 * sqrt(3) * pow(x, 2) - 20 * sqrt(3) * pow(y, 2)));
			 break;
		case 19:
			 Z = pow(y, 3) * (8 * sqrt(3) - 10 * sqrt(3) * pow(y, 2)) + pow(x, 2) * (30 * sqrt(3) * pow(x, 2) * y + y * (-24 * sqrt(3) + 20 * sqrt(3) * pow(y, 2)));
			 break;
		case 20:
			 Z = x * (10 * sqrt(3) * pow(y, 4) + pow(x, 2) * (2 * sqrt(3) * pow(x, 2) - 20 * sqrt(3) * pow(y, 2)));
			 break;
		case 21:
			 Z = 2 * sqrt(3) * pow(y, 5) + pow(x, 2) * (10 * sqrt(3) * pow(x, 2) * y - 20 * sqrt(3) * pow(y, 3));
			 break;
		case 22:
			 Z = -sqrt(7) + pow(y, 2) * (12 * sqrt(7) + pow(y, 2) * (-30 * sqrt(7) + 20 * sqrt(7) * pow(y, 2))) + pow(x, 2) * (12 * sqrt(7) + pow(y, 2) * (-60 * sqrt(7) + 60 * sqrt(7) * pow(y, 2)) + pow(x, 2) * (-30 * sqrt(7) + 20 * sqrt(7) * pow(x, 2) + 60 * sqrt(7) * pow(y, 2)));
			 break;
		case 23:
			 Z = x * (y * (12 * sqrt(14) + pow(y, 2) * (-40 * sqrt(14) + 30 * sqrt(14) * pow(y, 2))) + pow(x, 2) * (30 * sqrt(14) * pow(x, 2) * y + y * (-40 * sqrt(14) + 60 * sqrt(14) * pow(y, 2))));
			 break;
		case 24:
			 Z = pow(y, 2) * (-6 * sqrt(14) + pow(y, 2) * (20 * sqrt(14) - 15 * sqrt(14) * pow(y, 2))) + pow(x, 2) * (6 * sqrt(14) - 15 * sqrt(14) * pow(y, 4) + pow(x, 2) * (-20 * sqrt(14) + 15 * sqrt(14) * pow(x, 2) + 15 * sqrt(14) * pow(y, 2)));
			 break;
		case 25:
			 Z = x * (pow(x, 2) * (-20 * sqrt(14) * y + 24 * sqrt(14) * pow(x, 2) * y) + pow(y, 3) * (20 * sqrt(14) - 24 * sqrt(14) * pow(y, 2)));
			 break;
		case 26:
			 Z = pow(y, 4) * (-5 * sqrt(14) + 6 * sqrt(14) * pow(y, 2)) + pow(x, 2) * (pow(y, 2) * (30 * sqrt(14) - 30 * sqrt(14) * pow(y, 2)) + pow(x, 2) * (-5 * sqrt(14) + 6 * sqrt(14) * pow(x, 2) - 30 * sqrt(14) * pow(y, 2)));
			 break;
		case 27:
			 Z = x * (6 * sqrt(14) * pow(y, 5) + pow(x, 2) * (6 * sqrt(14) * pow(x, 2) * y - 20 * sqrt(14) * pow(y, 3)));
			 break;
		case 28:
			 Z = -(sqrt(14) * pow(y, 6)) + pow(x, 2) * (15 * sqrt(14) * pow(y, 4) + pow(x, 2) * (sqrt(14) * pow(x, 2) - 15 * sqrt(14) * pow(y, 2)));
			 break;
		case 29:
			 Z = y * (-16 + pow(x, 2) * (120 + pow(x, 2) * (-240 + 140 * pow(x, 2))) + pow(y, 2) * (120 + pow(x, 2) * (-480 + 420 * pow(x, 2)) + pow(y, 2) * (-240 + 420 * pow(x, 2) + 140 * pow(y, 2))));
			 break;
		case 30:
			 Z = x * (-16 + pow(y, 2) * (120 + pow(y, 2) * (-240 + 140 * pow(y, 2))) + pow(x, 2) * (120 + pow(y, 2) * (-480 + 420 * pow(y, 2)) + pow(x, 2) * (-240 + 140 * pow(x, 2) + 420 * pow(y, 2))));
			 break;
		case 31:
			 Z = pow(y, 3) * (-40 + pow(y, 2) * (120 - 84 * pow(y, 2))) + pow(x, 2) * (y * (120 + pow(y, 2) * (-240 + 84 * pow(y, 2))) + pow(x, 2) * (252 * pow(x, 2) * y + y * (-360 + 420 * pow(y, 2))));
			 break;
		case 32:
			 Z = x * (pow(y, 2) * (-120 + pow(y, 2) * (360 - 252 * pow(y, 2))) + pow(x, 2) * (40 + pow(y, 2) * (240 - 420 * pow(y, 2)) + pow(x, 2) * (-120 + 84 * pow(x, 2) - 84 * pow(y, 2))));
			 break;
		case 33:
			 Z = pow(y, 5) * (-24 + 28 * pow(y, 2)) + pow(x, 2) * (pow(y, 3) * (240 - 252 * pow(y, 2)) + pow(x, 2) * (140 * pow(x, 2) * y + y * (-120 - 140 * pow(y, 2))));
			 break;
		case 34:
			 Z = x * (pow(y, 4) * (-120 + 140 * pow(y, 2)) + pow(x, 2) * (pow(x, 2) * (-24 + 28 * pow(x, 2) - 252 * pow(y, 2)) + pow(y, 2) * (240 - 140 * pow(y, 2))));
			 break;
		case 35:
			 Z = -4 * pow(y, 7) + pow(x, 2) * (84 * pow(y, 5) + pow(x, 2) * (28 * pow(x, 2) * y - 140 * pow(y, 3)));
			 break;
		case 36:
			 Z = x * (-28 * pow(y, 6) + pow(x, 2) * (140 * pow(y, 4) + pow(x, 2) * (4 * pow(x, 2) - 84 * pow(y, 2))));
			 break;
		case 37:
			 Z = 3 + pow(y, 2) * (-60 + pow(y, 2) * (270 + pow(y, 2) * (-420 + 210 * pow(y, 2)))) + pow(x, 2) * (-60 + pow(y, 2) * (540 + pow(y, 2) * (-1260 + 840 * pow(y, 2))) + pow(x, 2) * (270 + pow(x, 2) * (-420 + 210 * pow(x, 2) + 840 * pow(y, 2)) + pow(y, 2) * (-1260 + 1260 * pow(y, 2))));
			 break;
		case 38:
			 Z = pow(y, 2) * (30 * sqrt(2) + pow(y, 2) * (-180 * sqrt(2) + pow(y, 2) * (315 * sqrt(2) - 168 * sqrt(2) * pow(y, 2)))) + pow(x, 2) * (-30 * sqrt(2) + pow(y, 4) * (315 * sqrt(2) - 336 * sqrt(2) * pow(y, 2)) + pow(x, 2) * (180 * sqrt(2) - 315 * sqrt(2) * pow(y, 2) + pow(x, 2) * (-315 * sqrt(2) + 168 * sqrt(2) * pow(x, 2) + 336 * sqrt(2) * pow(y, 2))));
			 break;
		case 39:
			 Z = x * (y * (-60 * sqrt(2) + pow(y, 2) * (360 * sqrt(2) + pow(y, 2) * (-630 * sqrt(2) + 336 * sqrt(2) * pow(y, 2)))) + pow(x, 2) * (y * (360 * sqrt(2) + pow(y, 2) * (-1260 * sqrt(2) + 1008 * sqrt(2) * pow(y, 2))) + pow(x, 2) * (336 * sqrt(2) * pow(x, 2) * y + y * (-630 * sqrt(2) + 1008 * sqrt(2) * pow(y, 2)))));
			 break;
		case 40:
			 Z = pow(y, 4) * (45 * sqrt(2) + pow(y, 2) * (-126 * sqrt(2) + 84 * sqrt(2) * pow(y, 2))) + pow(x, 2) * (pow(y, 2) * (-270 * sqrt(2) + pow(y, 2) * (630 * sqrt(2) - 336 * sqrt(2) * pow(y, 2))) + pow(x, 2) * (45 * sqrt(2) + pow(y, 2) * (630 * sqrt(2) - 840 * sqrt(2) * pow(y, 2)) + pow(x, 2) * (-126 * sqrt(2) + 84 * sqrt(2) * pow(x, 2) - 336 * sqrt(2) * pow(y, 2))));
			 break;
		case 41:
			 Z = x * (pow(y, 3) * (-180 * sqrt(2) + pow(y, 2) * (504 * sqrt(2) - 336 * sqrt(2) * pow(y, 2))) + pow(x, 2) * (y * (180 * sqrt(2) - 336 * sqrt(2) * pow(y, 4)) + pow(x, 2) * (336 * sqrt(2) * pow(x, 2) * y + y * (-504 * sqrt(2) + 336 * sqrt(2) * pow(y, 2)))));
			 break;
		case 42:
			 Z = pow(y, 6) * (21 * sqrt(2) - 24 * sqrt(2) * pow(y, 2)) + pow(x, 2) * (pow(y, 4) * (-315 * sqrt(2) + 336 * sqrt(2) * pow(y, 2)) + pow(x, 2) * (315 * sqrt(2) * pow(y, 2) + pow(x, 2) * (-21 * sqrt(2) + 24 * sqrt(2) * pow(x, 2) - 336 * sqrt(2) * pow(y, 2))));
			 break;
		case 43:
			 Z = x * (pow(y, 5) * (-126 * sqrt(2) + 144 * sqrt(2) * pow(y, 2)) + pow(x, 2) * (pow(y, 3) * (420 * sqrt(2) - 336 * sqrt(2) * pow(y, 2)) + pow(x, 2) * (144 * sqrt(2) * pow(x, 2) * y + y * (-126 * sqrt(2) - 336 * sqrt(2) * pow(y, 2)))));
			 break;
		case 44:
			 Z = 3 * sqrt(2) * pow(y, 8) + pow(x, 2) * (-84 * sqrt(2) * pow(y, 6) + pow(x, 2) * (210 * sqrt(2) * pow(y, 4) + pow(x, 2) * (3 * sqrt(2) * pow(x, 2) - 84 * sqrt(2) * pow(y, 2))));
			 break;
		case 45:
			 Z = x * (-24 * sqrt(2) * pow(y, 7) + pow(x, 2) * (168 * sqrt(2) * pow(y, 5) + pow(x, 2) * (24 * sqrt(2) * pow(x, 2) * y - 168 * sqrt(2) * pow(y, 3))));
			 break;
		case 46:
			 Z = x * (10 * sqrt(5) + pow(y, 2) * (-120 * sqrt(5) + pow(y, 2) * (420 * sqrt(5) + pow(y, 2) * (-560 * sqrt(5) + 252 * sqrt(5) * pow(y, 2)))) + pow(x, 2) * (-120 * sqrt(5) + pow(y, 2) * (840 * sqrt(5) + pow(y, 2) * (-1680 * sqrt(5) + 1008 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (420 * sqrt(5) + pow(x, 2) * (-560 * sqrt(5) + 252 * sqrt(5) * pow(x, 2) + 1008 * sqrt(5) * pow(y, 2)) + pow(y, 2) * (-1680 * sqrt(5) + 1512 * sqrt(5) * pow(y, 2)))));
			 break;
		case 47:
			 Z = y * (10 * sqrt(5) + pow(y, 2) * (-120 * sqrt(5) + pow(y, 2) * (420 * sqrt(5) + pow(y, 2) * (-560 * sqrt(5) + 252 * sqrt(5) * pow(y, 2))))) + pow(x, 2) * (y * (-120 * sqrt(5) + pow(y, 2) * (840 * sqrt(5) + pow(y, 2) * (-1680 * sqrt(5) + 1008 * sqrt(5) * pow(y, 2)))) + pow(x, 2) * (pow(x, 2) * (252 * sqrt(5) * pow(x, 2) * y + y * (-560 * sqrt(5) + 1008 * sqrt(5) * pow(y, 2))) + y * (420 * sqrt(5) + pow(y, 2) * (-1680 * sqrt(5) + 1512 * sqrt(5) * pow(y, 2)))));
			 break;
		case 48:
			 Z = x * (pow(x, 2) * (-40 * sqrt(5) + pow(y, 2) * (-420 * sqrt(5) + pow(y, 2) * (1680 * sqrt(5) - 1344 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (210 * sqrt(5) + pow(x, 2) * (-336 * sqrt(5) + 168 * sqrt(5) * pow(x, 2)) + pow(y, 2) * (336 * sqrt(5) - 1008 * sqrt(5) * pow(y, 2)))) + pow(y, 2) * (120 * sqrt(5) + pow(y, 2) * (-630 * sqrt(5) + pow(y, 2) * (1008 * sqrt(5) - 504 * sqrt(5) * pow(y, 2)))));
			 break;
		case 49:
			 Z = pow(y, 3) * (40 * sqrt(5) + pow(y, 2) * (-210 * sqrt(5) + pow(y, 2) * (336 * sqrt(5) - 168 * sqrt(5) * pow(y, 2)))) + pow(x, 2) * (y * (-120 * sqrt(5) + pow(y, 2) * (420 * sqrt(5) - 336 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (y * (630 * sqrt(5) + pow(y, 2) * (-1680 * sqrt(5) + 1008 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (504 * sqrt(5) * pow(x, 2) * y + y * (-1008 * sqrt(5) + 1344 * sqrt(5) * pow(y, 2)))));
			 break;
		case 50:
			 Z = x * (pow(y, 4) * (210 * sqrt(5) + pow(y, 2) * (-560 * sqrt(5) + 360 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (pow(y, 2) * (-420 * sqrt(5) + 560 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (42 * sqrt(5) + pow(y, 2) * (1008 * sqrt(5) - 1008 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (-112 * sqrt(5) + 72 * sqrt(5) * pow(x, 2) - 576 * sqrt(5) * pow(y, 2)))));
			 break;
		case 51:
			 Z = pow(y, 5) * (42 * sqrt(5) + pow(y, 2) * (-112 * sqrt(5) + 72 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (pow(y, 3) * (-420 * sqrt(5) + pow(y, 2) * (1008 * sqrt(5) - 576 * sqrt(5) * pow(y, 2))) + pow(x, 2) * (pow(x, 2) * (-560 * sqrt(5) * y + 360 * sqrt(5) * pow(x, 2) * y) + y * (210 * sqrt(5) + pow(y, 2) * (560 * sqrt(5) - 1008 * sqrt(5) * pow(y, 2)))));
			 break;
		case 52:
			 Z = x * (pow(y, 6) * (112 * sqrt(5) - 126 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (pow(y, 4) * (-560 * sqrt(5) + 504 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (pow(x, 2) * (-16 * sqrt(5) + 18 * sqrt(5) * pow(x, 2) - 360 * sqrt(5) * pow(y, 2)) + pow(y, 2) * (336 * sqrt(5) + 252 * sqrt(5) * pow(y, 2)))));
			 break;
		case 53:
			 Z = pow(y, 7) * (16 * sqrt(5) - 18 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (pow(y, 5) * (-336 * sqrt(5) + 360 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (pow(y, 3) * (560 * sqrt(5) - 252 * sqrt(5) * pow(y, 2)) + pow(x, 2) * (126 * sqrt(5) * pow(x, 2) * y + y * (-112 * sqrt(5) - 504 * sqrt(5) * pow(y, 2)))));
			 break;
		case 54:
			 Z = x * (18 * sqrt(5) * pow(y, 8) + pow(x, 2) * (-168 * sqrt(5) * pow(y, 6) + pow(x, 2) * (252 * sqrt(5) * pow(y, 4) + pow(x, 2) * (2 * sqrt(5) * pow(x, 2) - 72 * sqrt(5) * pow(y, 2)))));
			 break;
		case 55:
			 Z = 2 * sqrt(5) * pow(y, 9) + pow(x, 2) * (-72 * sqrt(5) * pow(y, 7) + pow(x, 2) * (252 * sqrt(5) * pow(y, 5) + pow(x, 2) * (18 * sqrt(5) * pow(x, 2) * y - 168 * sqrt(5) * pow(y, 3))));
			 break;
		case 56:
			 Z = -sqrt(11) + pow(y, 2) * (30 * sqrt(11) + pow(y, 2) * (-210 * sqrt(11) + pow(y, 2) * (560 * sqrt(11) + pow(y, 2) * (-630 * sqrt(11) + 252 * sqrt(11) * pow(y, 2))))) + pow(x, 2) * (30 * sqrt(11) + pow(y, 2) * (-420 * sqrt(11) + pow(y, 2) * (1680 * sqrt(11) + pow(y, 2) * (-2520 * sqrt(11) + 1260 * sqrt(11) * pow(y, 2)))) + pow(x, 2) * (-210 * sqrt(11) + pow(y, 2) * (1680 * sqrt(11) + pow(y, 2) * (-3780 * sqrt(11) + 2520 * sqrt(11) * pow(y, 2))) + pow(x, 2) * (560 * sqrt(11) + pow(x, 2) * (-630 * sqrt(11) + 252 * sqrt(11) * pow(x, 2) + 1260 * sqrt(11) * pow(y, 2)) + pow(y, 2) * (-2520 * sqrt(11) + 2520 * sqrt(11) * pow(y, 2)))));
			 break;
		case 57:
			 Z = x * (y * (30 * sqrt(22) + pow(y, 2) * (-280 * sqrt(22) + pow(y, 2) * (840 * sqrt(22) + pow(y, 2) * (-1008 * sqrt(22) + 420 * sqrt(22) * pow(y, 2))))) + pow(x, 2) * (y * (-280 * sqrt(22) + pow(y, 2) * (1680 * sqrt(22) + pow(y, 2) * (-3024 * sqrt(22) + 1680 * sqrt(22) * pow(y, 2)))) + pow(x, 2) * (pow(x, 2) * (420 * sqrt(22) * pow(x, 2) * y + y * (-1008 * sqrt(22) + 1680 * sqrt(22) * pow(y, 2))) + y * (840 * sqrt(22) + pow(y, 2) * (-3024 * sqrt(22) + 2520 * sqrt(22) * pow(y, 2))))));
			 break;
		case 58:
			 Z = pow(y, 2) * (-15 * sqrt(22) + pow(y, 2) * (140 * sqrt(22) + pow(y, 2) * (-420 * sqrt(22) + pow(y, 2) * (504 * sqrt(22) - 210 * sqrt(22) * pow(y, 2))))) + pow(x, 2) * (15 * sqrt(22) + pow(y, 4) * (-420 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 630 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (-140 * sqrt(22) + pow(y, 2) * (420 * sqrt(22) - 420 * sqrt(22) * pow(y, 4)) + pow(x, 2) * (420 * sqrt(22) + pow(y, 2) * (-1008 * sqrt(22) + 420 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (-504 * sqrt(22) + 210 * sqrt(22) * pow(x, 2) + 630 * sqrt(22) * pow(y, 2)))));
			 break;
		case 59:
			 Z = x * (pow(y, 3) * (140 * sqrt(22) + pow(y, 2) * (-672 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 480 * sqrt(22) * pow(y, 2)))) + pow(x, 2) * (y * (-140 * sqrt(22) + pow(y, 4) * (1008 * sqrt(22) - 960 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (y * (672 * sqrt(22) - 1008 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (480 * sqrt(22) * pow(x, 2) * y + y * (-1008 * sqrt(22) + 960 * sqrt(22) * pow(y, 2))))));
			 break;
		case 60:
			 Z = pow(y, 4) * (-35 * sqrt(22) + pow(y, 2) * (168 * sqrt(22) + pow(y, 2) * (-252 * sqrt(22) + 120 * sqrt(22) * pow(y, 2)))) + pow(x, 2) * (pow(y, 2) * (210 * sqrt(22) + pow(y, 2) * (-840 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 360 * sqrt(22) * pow(y, 2)))) + pow(x, 2) * (-35 * sqrt(22) + pow(y, 2) * (-840 * sqrt(22) + pow(y, 2) * (2520 * sqrt(22) - 1680 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (168 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 1680 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (-252 * sqrt(22) + 120 * sqrt(22) * pow(x, 2) - 360 * sqrt(22) * pow(y, 2)))));
			 break;
		case 61:
			 Z = x * (pow(y, 5) * (168 * sqrt(22) + pow(y, 2) * (-432 * sqrt(22) + 270 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (pow(y, 3) * (-560 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 360 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (y * (168 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 1260 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (270 * sqrt(22) * pow(x, 2) * y + y * (-432 * sqrt(22) - 360 * sqrt(22) * pow(y, 2))))));
			 break;
		case 62:
			 Z = pow(y, 6) * (-28 * sqrt(22) + pow(y, 2) * (72 * sqrt(22) - 45 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (pow(y, 4) * (420 * sqrt(22) + pow(y, 2) * (-1008 * sqrt(22) + 585 * sqrt(22) * pow(y, 2))) + pow(x, 2) * (pow(y, 2) * (-420 * sqrt(22) + 630 * sqrt(22) * pow(y, 4)) + pow(x, 2) * (28 * sqrt(22) + pow(y, 2) * (1008 * sqrt(22) - 630 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (-72 * sqrt(22) + 45 * sqrt(22) * pow(x, 2) - 585 * sqrt(22) * pow(y, 2)))));
			 break;
		case 63:
			 Z = x * (pow(y, 7) * (72 * sqrt(22) - 80 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (pow(y, 5) * (-504 * sqrt(22) + 480 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (504 * sqrt(22) * pow(y, 3) + pow(x, 2) * (80 * sqrt(22) * pow(x, 2) * y + y * (-72 * sqrt(22) - 480 * sqrt(22) * pow(y, 2))))));
			 break;
		case 64:
			 Z = pow(y, 8) * (-9 * sqrt(22) + 10 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (pow(y, 6) * (252 * sqrt(22) - 270 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (pow(y, 4) * (-630 * sqrt(22) + 420 * sqrt(22) * pow(y, 2)) + pow(x, 2) * (pow(x, 2) * (-9 * sqrt(22) + 10 * sqrt(22) * pow(x, 2) - 270 * sqrt(22) * pow(y, 2)) + pow(y, 2) * (252 * sqrt(22) + 420 * sqrt(22) * pow(y, 2)))));
			 break;
		case 65:
			 Z = x * (10 * sqrt(22) * pow(y, 9) + pow(x, 2) * (-120 * sqrt(22) * pow(y, 7) + pow(x, 2) * (252 * sqrt(22) * pow(y, 5) + pow(x, 2) * (10 * sqrt(22) * pow(x, 2) * y - 120 * sqrt(22) * pow(y, 3)))));
			 break;
		case 66:
			 Z = -(sqrt(22) * pow(y, 10)) + pow(x, 2) * (45 * sqrt(22) * pow(y, 8) + pow(x, 2) * (-210 * sqrt(22) * pow(y, 6) + pow(x, 2) * (210 * sqrt(22) * pow(y, 4) + pow(x, 2) * (sqrt(22) * pow(x, 2) - 45 * sqrt(22) * pow(y, 2)))));
			 break;
		case 67:
			 Z = y * (-12 * sqrt(6) + pow(y, 2) * (210 * sqrt(6) + pow(y, 2) * (-1120 * sqrt(6) + pow(y, 2) * (2520 * sqrt(6) + pow(y, 2) * (-2520 * sqrt(6) + 924 * sqrt(6) * pow(y, 2)))))) + pow(x, 2) * (y * (210 * sqrt(6) + pow(y, 2) * (-2240 * sqrt(6) + pow(y, 2) * (7560 * sqrt(6) + pow(y, 2) * (-10080 * sqrt(6) + 4620 * sqrt(6) * pow(y, 2))))) + pow(x, 2) * (y * (-1120 * sqrt(6) + pow(y, 2) * (7560 * sqrt(6) + pow(y, 2) * (-15120 * sqrt(6) + 9240 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (pow(x, 2) * (924 * sqrt(6) * pow(x, 2) * y + y * (-2520 * sqrt(6) + 4620 * sqrt(6) * pow(y, 2))) + y * (2520 * sqrt(6) + pow(y, 2) * (-10080 * sqrt(6) + 9240 * sqrt(6) * pow(y, 2))))));
			 break;
		case 68:
			 Z = x * (-12 * sqrt(6) + pow(y, 2) * (210 * sqrt(6) + pow(y, 2) * (-1120 * sqrt(6) + pow(y, 2) * (2520 * sqrt(6) + pow(y, 2) * (-2520 * sqrt(6) + 924 * sqrt(6) * pow(y, 2))))) + pow(x, 2) * (210 * sqrt(6) + pow(y, 2) * (-2240 * sqrt(6) + pow(y, 2) * (7560 * sqrt(6) + pow(y, 2) * (-10080 * sqrt(6) + 4620 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (-1120 * sqrt(6) + pow(y, 2) * (7560 * sqrt(6) + pow(y, 2) * (-15120 * sqrt(6) + 9240 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (2520 * sqrt(6) + pow(x, 2) * (-2520 * sqrt(6) + 924 * sqrt(6) * pow(x, 2) + 4620 * sqrt(6) * pow(y, 2)) + pow(y, 2) * (-10080 * sqrt(6) + 9240 * sqrt(6) * pow(y, 2))))));
			 break;
		case 69:
			 Z = pow(y, 3) * (-70 * sqrt(6) + pow(y, 2) * (560 * sqrt(6) + pow(y, 2) * (-1512 * sqrt(6) + pow(y, 2) * (1680 * sqrt(6) - 660 * sqrt(6) * pow(y, 2))))) + pow(x, 2) * (y * (210 * sqrt(6) + pow(y, 2) * (-1120 * sqrt(6) + pow(y, 2) * (1512 * sqrt(6) - 660 * sqrt(6) * pow(y, 4)))) + pow(x, 2) * (y * (-1680 * sqrt(6) + pow(y, 2) * (7560 * sqrt(6) + pow(y, 2) * (-10080 * sqrt(6) + 3960 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (pow(x, 2) * (1980 * sqrt(6) * pow(x, 2) * y + y * (-5040 * sqrt(6) + 7260 * sqrt(6) * pow(y, 2))) + y * (4536 * sqrt(6) + pow(y, 2) * (-13440 * sqrt(6) + 9240 * sqrt(6) * pow(y, 2))))));
			 break;
		case 70:
			 Z = x * (pow(y, 2) * (-210 * sqrt(6) + pow(y, 2) * (1680 * sqrt(6) + pow(y, 2) * (-4536 * sqrt(6) + pow(y, 2) * (5040 * sqrt(6) - 1980 * sqrt(6) * pow(y, 2))))) + pow(x, 2) * (70 * sqrt(6) + pow(y, 2) * (1120 * sqrt(6) + pow(y, 2) * (-7560 * sqrt(6) + pow(y, 2) * (13440 * sqrt(6) - 7260 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (-560 * sqrt(6) + pow(y, 2) * (-1512 * sqrt(6) + pow(y, 2) * (10080 * sqrt(6) - 9240 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (1512 * sqrt(6) - 3960 * sqrt(6) * pow(y, 4) + pow(x, 2) * (-1680 * sqrt(6) + 660 * sqrt(6) * pow(x, 2) + 660 * sqrt(6) * pow(y, 2))))));
			 break;
		case 71:
			 Z = pow(y, 5) * (-112 * sqrt(6) + pow(y, 2) * (504 * sqrt(6) + pow(y, 2) * (-720 * sqrt(6) + 330 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (pow(y, 3) * (1120 * sqrt(6) + pow(y, 2) * (-4536 * sqrt(6) + pow(y, 2) * (5760 * sqrt(6) - 2310 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (y * (-560 * sqrt(6) + pow(y, 2) * (-2520 * sqrt(6) + pow(y, 2) * (10080 * sqrt(6) - 7260 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (y * (2520 * sqrt(6) - 4620 * sqrt(6) * pow(y, 4)) + pow(x, 2) * (1650 * sqrt(6) * pow(x, 2) * y + y * (-3600 * sqrt(6) + 1650 * sqrt(6) * pow(y, 2))))));
			 break;
		case 72:
			 Z = x * (pow(y, 4) * (-560 * sqrt(6) + pow(y, 2) * (2520 * sqrt(6) + pow(y, 2) * (-3600 * sqrt(6) + 1650 * sqrt(6) * pow(y, 2)))) + pow(x, 2) * (pow(y, 2) * (1120 * sqrt(6) + pow(y, 2) * (-2520 * sqrt(6) + 1650 * sqrt(6) * pow(y, 4))) + pow(x, 2) * (-112 * sqrt(6) + pow(y, 2) * (-4536 * sqrt(6) + pow(y, 2) * (10080 * sqrt(6) - 4620 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (504 * sqrt(6) + pow(y, 2) * (5760 * sqrt(6) - 7260 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (-720 * sqrt(6) + 330 * sqrt(6) * pow(x, 2) - 2310 * sqrt(6) * pow(y, 2))))));
			 break;
		case 73:
			 Z = pow(y, 7) * (-72 * sqrt(6) + pow(y, 2) * (180 * sqrt(6) - 110 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (pow(y, 5) * (1512 * sqrt(6) + pow(y, 2) * (-3600 * sqrt(6) + 2090 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (pow(y, 3) * (-2520 * sqrt(6) + pow(y, 2) * (2520 * sqrt(6) + 660 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (y * (504 * sqrt(6) + pow(y, 2) * (5040 * sqrt(6) - 4620 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (770 * sqrt(6) * pow(x, 2) * y + y * (-1260 * sqrt(6) - 2310 * sqrt(6) * pow(y, 2))))));
			 break;
		case 74:
			 Z = x * (pow(y, 6) * (-504 * sqrt(6) + pow(y, 2) * (1260 * sqrt(6) - 770 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (pow(y, 4) * (2520 * sqrt(6) + pow(y, 2) * (-5040 * sqrt(6) + 2310 * sqrt(6) * pow(y, 2))) + pow(x, 2) * (pow(x, 2) * (72 * sqrt(6) + pow(x, 2) * (-180 * sqrt(6) + 110 * sqrt(6) * pow(x, 2) - 2090 * sqrt(6) * pow(y, 2)) + pow(y, 2) * (3600 * sqrt(6) - 660 * sqrt(6) * pow(y, 2))) + pow(y, 2) * (-1512 * sqrt(6) + pow(y, 2) * (-2520 * sqrt(6) + 4620 * sqrt(6) * pow(y, 2))))));
			 break;
		case 75:
			 Z = pow(y, 9) * (-20 * sqrt(6) + 22 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (pow(y, 7) * (720 * sqrt(6) - 770 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (pow(y, 5) * (-2520 * sqrt(6) + 1980 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (pow(y, 3) * (1680 * sqrt(6) + 924 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (198 * sqrt(6) * pow(x, 2) * y + y * (-180 * sqrt(6) - 1650 * sqrt(6) * pow(y, 2))))));
			 break;
		case 76:
			 Z = x * (pow(y, 8) * (-180 * sqrt(6) + 198 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (pow(y, 6) * (1680 * sqrt(6) - 1650 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (pow(y, 4) * (-2520 * sqrt(6) + 924 * sqrt(6) * pow(y, 2)) + pow(x, 2) * (pow(x, 2) * (-20 * sqrt(6) + 22 * sqrt(6) * pow(x, 2) - 770 * sqrt(6) * pow(y, 2)) + pow(y, 2) * (720 * sqrt(6) + 1980 * sqrt(6) * pow(y, 2))))));
			 break;
		case 77:
			 Z = -2 * sqrt(6) * pow(y, 11) + pow(x, 2) * (110 * sqrt(6) * pow(y, 9) + pow(x, 2) * (-660 * sqrt(6) * pow(y, 7) + pow(x, 2) * (924 * sqrt(6) * pow(y, 5) + pow(x, 2) * (22 * sqrt(6) * pow(x, 2) * y - 330 * sqrt(6) * pow(y, 3)))));
			 break;
		case 78:
			 Z = x * (-22 * sqrt(6) * pow(y, 10) + pow(x, 2) * (330 * sqrt(6) * pow(y, 8) + pow(x, 2) * (-924 * sqrt(6) * pow(y, 6) + pow(x, 2) * (660 * sqrt(6) * pow(y, 4) + pow(x, 2) * (2 * sqrt(6) * pow(x, 2) - 110 * sqrt(6) * pow(y, 2))))));
			 break;
		case 79:
			 Z = sqrt(13) + pow(y, 2) * (-42 * sqrt(13) + pow(y, 2) * (420 * sqrt(13) + pow(y, 2) * (-1680 * sqrt(13) + pow(y, 2) * (3150 * sqrt(13) + pow(y, 2) * (-2772 * sqrt(13) + 924 * sqrt(13) * pow(y, 2)))))) + pow(x, 2) * (-42 * sqrt(13) + pow(y, 2) * (840 * sqrt(13) + pow(y, 2) * (-5040 * sqrt(13) + pow(y, 2) * (12600 * sqrt(13) + pow(y, 2) * (-13860 * sqrt(13) + 5544 * sqrt(13) * pow(y, 2))))) + pow(x, 2) * (420 * sqrt(13) + pow(y, 2) * (-5040 * sqrt(13) + pow(y, 2) * (18900 * sqrt(13) + pow(y, 2) * (-27720 * sqrt(13) + 13860 * sqrt(13) * pow(y, 2)))) + pow(x, 2) * (-1680 * sqrt(13) + pow(x, 2) * (3150 * sqrt(13) + pow(x, 2) * (-2772 * sqrt(13) + 924 * sqrt(13) * pow(x, 2) + 5544 * sqrt(13) * pow(y, 2)) + pow(y, 2) * (-13860 * sqrt(13) + 13860 * sqrt(13) * pow(y, 2))) + pow(y, 2) * (12600 * sqrt(13) + pow(y, 2) * (-27720 * sqrt(13) + 18480 * sqrt(13) * pow(y, 2))))));
			 break;
		case 80:
			 Z = pow(y, 2) * (21 * sqrt(26) + pow(y, 2) * (-280 * sqrt(26) + pow(y, 2) * (1260 * sqrt(26) + pow(y, 2) * (-2520 * sqrt(26) + pow(y, 2) * (2310 * sqrt(26) - 792 * sqrt(26) * pow(y, 2)))))) + pow(x, 2) * (-21 * sqrt(26) + pow(y, 4) * (1260 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (6930 * sqrt(26) - 3168 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (280 * sqrt(26) + pow(y, 2) * (-1260 * sqrt(26) + pow(y, 4) * (4620 * sqrt(26) - 3960 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (-1260 * sqrt(26) + pow(y, 2) * (5040 * sqrt(26) - 4620 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (2520 * sqrt(26) + pow(x, 2) * (-2310 * sqrt(26) + 792 * sqrt(26) * pow(x, 2) + 3168 * sqrt(26) * pow(y, 2)) + pow(y, 2) * (-6930 * sqrt(26) + 3960 * sqrt(26) * pow(y, 2))))));
			 break;
		case 81:
			 Z = x * (y * (-42 * sqrt(26) + pow(y, 2) * (560 * sqrt(26) + pow(y, 2) * (-2520 * sqrt(26) + pow(y, 2) * (5040 * sqrt(26) + pow(y, 2) * (-4620 * sqrt(26) + 1584 * sqrt(26) * pow(y, 2)))))) + pow(x, 2) * (y * (560 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (15120 * sqrt(26) + pow(y, 2) * (-18480 * sqrt(26) + 7920 * sqrt(26) * pow(y, 2))))) + pow(x, 2) * (y * (-2520 * sqrt(26) + pow(y, 2) * (15120 * sqrt(26) + pow(y, 2) * (-27720 * sqrt(26) + 15840 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (pow(x, 2) * (1584 * sqrt(26) * pow(x, 2) * y + y * (-4620 * sqrt(26) + 7920 * sqrt(26) * pow(y, 2))) + y * (5040 * sqrt(26) + pow(y, 2) * (-18480 * sqrt(26) + 15840 * sqrt(26) * pow(y, 2)))))));
			 break;
		case 82:
			 Z = pow(y, 4) * (70 * sqrt(26) + pow(y, 2) * (-504 * sqrt(26) + pow(y, 2) * (1260 * sqrt(26) + pow(y, 2) * (-1320 * sqrt(26) + 495 * sqrt(26) * pow(y, 2))))) + pow(x, 2) * (pow(y, 2) * (-420 * sqrt(26) + pow(y, 2) * (2520 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (3960 * sqrt(26) - 990 * sqrt(26) * pow(y, 2))))) + pow(x, 2) * (70 * sqrt(26) + pow(y, 2) * (2520 * sqrt(26) + pow(y, 2) * (-12600 * sqrt(26) + pow(y, 2) * (18480 * sqrt(26) - 8415 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (-504 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (18480 * sqrt(26) - 13860 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (1260 * sqrt(26) + pow(y, 2) * (3960 * sqrt(26) - 8415 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (-1320 * sqrt(26) + 495 * sqrt(26) * pow(x, 2) - 990 * sqrt(26) * pow(y, 2))))));
			 break;
		case 83:
			 Z = x * (pow(y, 3) * (-280 * sqrt(26) + pow(y, 2) * (2016 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (5280 * sqrt(26) - 1980 * sqrt(26) * pow(y, 2))))) + pow(x, 2) * (y * (280 * sqrt(26) + pow(y, 4) * (-5040 * sqrt(26) + pow(y, 2) * (10560 * sqrt(26) - 5940 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (y * (-2016 * sqrt(26) + pow(y, 2) * (5040 * sqrt(26) - 3960 * sqrt(26) * pow(y, 4))) + pow(x, 2) * (y * (5040 * sqrt(26) + pow(y, 2) * (-10560 * sqrt(26) + 3960 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (1980 * sqrt(26) * pow(x, 2) * y + y * (-5280 * sqrt(26) + 5940 * sqrt(26) * pow(y, 2)))))));
			 break;
		case 84:
			 Z = pow(y, 6) * (84 * sqrt(26) + pow(y, 2) * (-360 * sqrt(26) + pow(y, 2) * (495 * sqrt(26) - 220 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (pow(y, 4) * (-1260 * sqrt(26) + pow(y, 2) * (5040 * sqrt(26) + pow(y, 2) * (-6435 * sqrt(26) + 2640 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (pow(y, 2) * (1260 * sqrt(26) + pow(y, 4) * (-6930 * sqrt(26) + 5940 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (-84 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + 6930 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (360 * sqrt(26) + pow(y, 2) * (6435 * sqrt(26) - 5940 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (-495 * sqrt(26) + 220 * sqrt(26) * pow(x, 2) - 2640 * sqrt(26) * pow(y, 2))))));
			 break;
		case 85:
			 Z = x * (pow(y, 5) * (-504 * sqrt(26) + pow(y, 2) * (2160 * sqrt(26) + pow(y, 2) * (-2970 * sqrt(26) + 1320 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (pow(y, 3) * (1680 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (3960 * sqrt(26) - 440 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (y * (-504 * sqrt(26) + pow(y, 2) * (-5040 * sqrt(26) + pow(y, 2) * (13860 * sqrt(26) - 7920 * sqrt(26) * pow(y, 2)))) + pow(x, 2) * (y * (2160 * sqrt(26) + pow(y, 2) * (3960 * sqrt(26) - 7920 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (1320 * sqrt(26) * pow(x, 2) * y + y * (-2970 * sqrt(26) - 440 * sqrt(26) * pow(y, 2)))))));
			 break;
		case 86:
			 Z = pow(y, 8) * (45 * sqrt(26) + pow(y, 2) * (-110 * sqrt(26) + 66 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (pow(y, 6) * (-1260 * sqrt(26) + pow(y, 2) * (2970 * sqrt(26) - 1716 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (pow(y, 4) * (3150 * sqrt(26) + pow(y, 2) * (-4620 * sqrt(26) + 990 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (pow(x, 2) * (45 * sqrt(26) + pow(x, 2) * (-110 * sqrt(26) + 66 * sqrt(26) * pow(x, 2) - 1716 * sqrt(26) * pow(y, 2)) + pow(y, 2) * (2970 * sqrt(26) + 990 * sqrt(26) * pow(y, 2))) + pow(y, 2) * (-1260 * sqrt(26) + pow(y, 2) * (-4620 * sqrt(26) + 5544 * sqrt(26) * pow(y, 2))))));
			 break;
		case 87:
			 Z = x * (pow(y, 7) * (-360 * sqrt(26) + pow(y, 2) * (880 * sqrt(26) - 528 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (pow(y, 5) * (2520 * sqrt(26) + pow(y, 2) * (-5280 * sqrt(26) + 2640 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (pow(y, 3) * (-2520 * sqrt(26) + 3168 * sqrt(26) * pow(y, 4)) + pow(x, 2) * (y * (360 * sqrt(26) + pow(y, 2) * (5280 * sqrt(26) - 3168 * sqrt(26) * pow(y, 2))) + pow(x, 2) * (528 * sqrt(26) * pow(x, 2) * y + y * (-880 * sqrt(26) - 2640 * sqrt(26) * pow(y, 2)))))));
			 break;
		case 88:
			 Z = pow(y, 10) * (11 * sqrt(26) - 12 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (pow(y, 8) * (-495 * sqrt(26) + 528 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (pow(y, 6) * (2310 * sqrt(26) - 1980 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (-2310 * sqrt(26) * pow(y, 4) + pow(x, 2) * (pow(x, 2) * (-11 * sqrt(26) + 12 * sqrt(26) * pow(x, 2) - 528 * sqrt(26) * pow(y, 2)) + pow(y, 2) * (495 * sqrt(26) + 1980 * sqrt(26) * pow(y, 2))))));
			 break;
		case 89:
			 Z = x * (pow(y, 9) * (-110 * sqrt(26) + 120 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (pow(y, 7) * (1320 * sqrt(26) - 1320 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (pow(y, 5) * (-2772 * sqrt(26) + 1584 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (pow(y, 3) * (1320 * sqrt(26) + 1584 * sqrt(26) * pow(y, 2)) + pow(x, 2) * (120 * sqrt(26) * pow(x, 2) * y + y * (-110 * sqrt(26) - 1320 * sqrt(26) * pow(y, 2)))))));
			 break;
		case 90:
			 Z = sqrt(26) * pow(y, 12) + pow(x, 2) * (-66 * sqrt(26) * pow(y, 10) + pow(x, 2) * (495 * sqrt(26) * pow(y, 8) + pow(x, 2) * (-924 * sqrt(26) * pow(y, 6) + pow(x, 2) * (495 * sqrt(26) * pow(y, 4) + pow(x, 2) * (sqrt(26) * pow(x, 2) - 66 * sqrt(26) * pow(y, 2))))));
			 break;
		case 91:
			 Z = x * (-12 * sqrt(26) * pow(y, 11) + pow(x, 2) * (220 * sqrt(26) * pow(y, 9) + pow(x, 2) * (-792 * sqrt(26) * pow(y, 7) + pow(x, 2) * (792 * sqrt(26) * pow(y, 5) + pow(x, 2) * (12 * sqrt(26) * pow(x, 2) * y - 220 * sqrt(26) * pow(y, 3))))));
			 break;
		default:
			 Z = 0.0;
			 break;
	}
	
	return Z;
}



double RectZernike1(double x, double y)
{
	return 0.7759402897989853;
}

double RectZernike2(double x, double y)
{
	return 0.0 + 1.5231637396477693 * x;
}

double RectZernike3(double x, double y)
{
	return 0.0 + 2.855932011839568 * y;
}

double RectZernike4(double x, double y)
{
	return -1.0717762752241962 + 3.215328825672589 * pow(x, 2) + 3.215328825672589 * pow(y, 2);
}

double RectZernike5(double x, double y)
{
	return 0.0 + 5.606168593797669 * x * y;
}

double RectZernike6(double x, double y)
{
	return 0.5970795166473906 + 0.9145824215246474 * pow(x, 2) - 11.303890402755199 * pow(y, 2);
}

double RectZernike7(double x, double y)
{
	return 0.0 + y * (-4.683898944411964 + 11.936920590256237 * pow(x, 2) + 11.936920590256237 * pow(y, 2));
}

double RectZernike8(double x, double y)
{
	return 0.0 + x * (-3.8449101755941104 + 7.107755058081223 * pow(x, 2) + 7.107755058081223 * pow(y, 2));
}

double RectZernike9(double x, double y)
{
	return 0.0 + y * (5.574706310126192 + 2.9822058687518833 * pow(x, 2) - 47.7800928068453 * pow(y, 2));
}

double RectZernike10(double x, double y)
{
	return 0.0 + x * (0.5448360277642879 + 2.301876658170394 * pow(x, 2) - 21.94738879110386 * pow(y, 2));
}

double RectZernike11(double x, double y)
{
	return 1.3195274681326268 + pow(y, 2) * (-9.900792113866308 + 13.967435457706781 * pow(y, 2)) + pow(x, 2) * (-11.382919059782129 + 13.967435457706781 * pow(x, 2) + 27.934870915413562 * pow(y, 2));
}

double RectZernike12(double x, double y)
{
	return -0.5045877549716362 + pow(y, 2) * (19.952467628531956 - 49.74452471448675 * pow(y, 2)) + pow(x, 2) * (-3.18027421588698 + 9.245582767520709 * pow(x, 2) - 40.49894194696604 * pow(y, 2));
}

double RectZernike13(double x, double y)
{
	return 0.0 + x * (26.44937324586838 * pow(x, 2) * y + y * (-15.86962394752103 + 26.44937324586838 * pow(y, 2)));
}

double RectZernike14(double x, double y)
{
	return 0.6848481745752704 + pow(x, 2) * (0.003922502485664836 + 1.316669436900682 * pow(x, 2) - 11.95606447998486 * pow(y, 2)) + pow(y, 2) * (-35.07601196383209 + 201.13465118443742 * pow(y, 2));
}

double RectZernike15(double x, double y)
{
	return 0.0 + x * (7.5233776126249765 * pow(x, 2) * y + y * (8.840863475365907 - 92.98608130429997 * pow(y, 2)));
}

double RectZernike16(double x, double y)
{
	return 0.0 + x * (6.827090793163167 + pow(y, 2) * (-35.28569993872314 + 31.390928361378254 * pow(y, 2)) + pow(x, 2) * (-31.78919849235442 + 31.390928361378254 * pow(x, 2) + 62.78185672275651 * pow(y, 2)));
}

double RectZernike17(double x, double y)
{
	return 0.0 + y * (7.071023809633416 + pow(x, 2) * (-49.75846905161568 + 53.32770868157581 * pow(x, 2)) + pow(y, 2) * (-40.800515317357984 + 106.65541736315161 * pow(x, 2) + 53.32770868157581 * pow(y, 2)));
}

double RectZernike18(double x, double y)
{
	return 0.0 + x * (-0.5142579769575932 + pow(y, 2) * (62.568513388372054 - 110.14153261350096 * pow(y, 2)) + pow(x, 2) * (-11.543134719694734 + 20.95449126780032 * pow(x, 2) - 89.18704134570064 * pow(y, 2)));
}

double RectZernike19(double x, double y)
{
	return 0.0 + y * (-6.8996500407134995 + pow(x, 2) * (2.7536204169073955 + 31.02663667980761 * pow(x, 2)) + pow(y, 2) * (96.89341622570836 - 176.54960413298448 * pow(x, 2) - 207.57624081279207 * pow(y, 2)));
}

double RectZernike20(double x, double y)
{
	return 0.0 + x * (1.0929087596150495 + pow(x, 2) * (-0.6990451809320746 + 3.3805756580762605 * pow(x, 2) - 30.146251031213524 * pow(y, 2)) + pow(y, 2) * (-60.17931780586514 + 391.22650589674697 * pow(y, 2)));
}

double RectZernike21(double x, double y)
{
	return 0.0 + y * (8.44852453647809 + pow(x, 2) * (3.7869794501906426 + 4.223230971634166 * pow(x, 2)) + pow(y, 2) * (-196.63258190927957 - 49.7112959833446 * pow(x, 2) + 851.5561879171631 * pow(y, 2)));
}

double RectZernike22(double x, double y)
{
	return -1.489685334333457 + pow(y, 2) * (19.246482248063117 + pow(y, 2) * (-63.80499990379534 + 59.049818050083445 * pow(y, 2))) + pow(x, 2) * (25.740292225289032 + pow(y, 2) * (-151.84238928818252 + 177.14945415025034 * pow(y, 2)) + pow(x, 2) * (-75.7672877616197 + 59.049818050083445 * pow(x, 2) + 177.14945415025034 * pow(y, 2)));
}

double RectZernike23(double x, double y)
{
	return 0.0 + x * (y * (33.59411078744737 + pow(y, 2) * (-141.54685357801458 + 119.92292045542553 * pow(y, 2))) + pow(x, 2) * (119.92292045542553 * pow(x, 2) * y + y * (-135.60834045356893 + 239.84584091085105 * pow(y, 2))));
}

double RectZernike24(double x, double y)
{
	return 0.46897803744704003 + pow(y, 2) * (-34.60696389287215 + pow(y, 2) * (191.77095292606123 - 252.01210542891937 * pow(y, 2))) + pow(x, 2) * (7.525020472546814 + pow(y, 2) * (174.94370476458073 - 445.7069203634947 * pow(y, 2)) + pow(x, 2) * (-51.91955433314999 + 58.317290494344036 * pow(x, 2) - 135.3775244402313 * pow(y, 2)));
}

double RectZernike25(double x, double y)
{
	return 0.0 + x * (y * (-19.37931568023186 + pow(y, 2) * (295.0125337093933 - 459.7496579098531 * pow(y, 2))) + pow(x, 2) * (70.37695135234335 * pow(x, 2) * y + y * (-9.142995546236307 - 389.37270655750973 * pow(y, 2))));
}

double RectZernike26(double x, double y)
{
	return -0.742808744181275 + pow(x, 2) * (2.8096334988144918 + pow(x, 2) * (-13.165444753862289 + 20.84979906022079 * pow(x, 2) - 121.5133071297899 * pow(y, 2)) + pow(y, 2) * (-56.682214344189674 + 725.8101243118249 * pow(y, 2))) + pow(y, 2) * (48.79360237686487 + pow(y, 2) * (-450.5315373722875 + 868.1732305018356 * pow(y, 2)));
}

double RectZernike27(double x, double y)
{
	return 0.0 + x * (pow(x, 2) * (10.847251470433093 * pow(x, 2) * y + y * (7.277528736077613 - 125.39114893203441 * pow(y, 2))) + y * (13.131119469256987 + pow(y, 2) * (-349.15510244219985 + 1657.0356717561892 * pow(y, 2))));
}

double RectZernike28(double x, double y)
{
	return 0.6654795024917632 + pow(y, 2) * (-71.00244966018636 + pow(y, 2) * (1035.096014442176 - 3607.0515415496975 * pow(y, 2))) + pow(x, 2) * (0.5718641309392751 + pow(x, 2) * (-0.7816539515729346 + 1.9104977960649876 * pow(x, 2) - 16.8879891620154 * pow(y, 2)) + pow(y, 2) * (-28.350401958265138 + 208.7281355712703 * pow(y, 2)));
}



/// HIFN Returns a pointer to a Rectangular Zernike Polynomial, orthogonalised on a rectangle with the
///      same aspect ratio as the SLM, and fitted precisely inside a unit circle
dFdd RectZernikePolynomial(int i)
{
	switch(i)
	{
		case 1:
			return &RectZernike1;
			break;
		case 2:
			return &RectZernike2;
			break;
		case 3:
			return &RectZernike3;
			break;
		case 4:
			return &RectZernike4;
			break;
		case 5:
			return &RectZernike5;
			break;
		case 6:
			return &RectZernike6;
			break;
		case 7:
			return &RectZernike7;
			break;
		case 8:
			return &RectZernike8;
			break;
		case 9:
			return &RectZernike9;
			break;
		case 10:
			return &RectZernike10;
			break;
		case 11:
			return &RectZernike11;
			break;
		case 12:
			return &RectZernike12;
			break;
		case 13:
			return &RectZernike13;
			break;
		case 14:
			return &RectZernike14;
			break;
		case 15:
			return &RectZernike15;
			break;
		case 16:
			return &RectZernike16;
			break;
		case 17:
			return &RectZernike17;
			break;
		case 18:
			return &RectZernike18;
			break;
		case 19:
			return &RectZernike19;
			break;
		case 20:
			return &RectZernike20;
			break;
		case 21:
			return &RectZernike21;
			break;
		case 22:
			return &RectZernike22;
			break;
		case 23:
			return &RectZernike23;
			break;
		case 24:
			return &RectZernike24;
			break;
		case 25:
			return &RectZernike25;
			break;
		case 26:
			return &RectZernike26;
			break;
		case 27:
			return &RectZernike27;
			break;
		case 28:
			return &RectZernike28;
			break;
		default:
			return &RectZernike1;
			break;
	}
	
	return NULL;
}



/// HIFN Generate the matrix of dot products between the first N (max. 28) Zernike Circle polynomials, evaluated on a rectangle 
double* generateCZZ(double a, double b, int N)
{
	// round off N to be at most 28
	N = (N < 28 ? N : 28);
	
	// allocate the array
	double* CZZ = (double*) calloc(N * N, sizeof(double));
	
	// fill the array with the entries calculated by Mathematica (ZernikeCirclePolynomials.nb)
	switch (N)
	{
		case 28:
			CZZ[27 * N + 0] = 4 * sqrt(0.2857142857142857) * a * b * (pow(a, 6) - 7 * pow(a, 4) * pow(b, 2) + 7 * pow(a, 2) * pow(b, 4) - pow(b, 6));
			CZZ[0 * N + 27] = CZZ[27 * N + 0];
			CZZ[27 * N + 3] = (4 * sqrt(0.09523809523809523) * a * b * (14 * pow(a, 8) + 63 * pow(a, 4) * pow(b, 2) + 9 * pow(b, 6) - 14 * pow(b, 8) + 21 * pow(a, 2) * pow(b, 4) * (-3 + 4 * pow(b, 2)) - 3 * pow(a, 6) * (3 + 28 * pow(b, 2))))/3.;
			CZZ[3 * N + 27] = CZZ[27 * N + 3];
			CZZ[27 * N + 5] = (8 * a * b * (35 * pow(a, 8) - 240 * pow(a, 6) * pow(b, 2) + 378 * pow(a, 4) * pow(b, 4) - 240 * pow(a, 2) * pow(b, 6) + 35 * pow(b, 8)))/(15. * sqrt(21));
			CZZ[5 * N + 27] = CZZ[27 * N + 5];
			CZZ[27 * N + 10] = (4 * sqrt(0.05714285714285714) * a * b * (1890 * pow(a, 10) - 770 * pow(a, 8) * (3 + 13 * pow(b, 2)) + 693 * pow(a, 4) * pow(b, 2) * (-5 + 12 * pow(b, 4)) + 385 * pow(a, 2) * pow(b, 4) * (9 - 36 * pow(b, 2) + 26 * pow(b, 4)) - 99 * pow(a, 6) * (-5 - 140 * pow(b, 2) + 84 * pow(b, 4)) - 15 * pow(b, 6) * (33 - 154 * pow(b, 2) + 126 * pow(b, 4))))/99.;
			CZZ[10 * N + 27] = CZZ[27 * N + 10];
			CZZ[27 * N + 11] = (8 * a * b * (1260 * pow(a, 10) + 1386 * pow(a, 4) * pow(b, 4) * (-9 + 4 * pow(b, 2)) + 792 * pow(a, 6) * pow(b, 2) * (10 + 7 * pow(b, 2)) + 105 * pow(b, 8) * (-11 + 12 * pow(b, 2)) - 385 * pow(a, 8) * (3 + 20 * pow(b, 2)) - 220 * pow(a, 2) * pow(b, 6) * (-36 + 35 * pow(b, 2))))/(99. * sqrt(35));
			CZZ[11 * N + 27] = CZZ[27 * N + 11];
			CZZ[27 * N + 13] = (8 * a * b * (315 * pow(a, 10) - 2695 * pow(a, 8) * pow(b, 2) + 10494 * pow(a, 6) * pow(b, 4) - 10494 * pow(a, 4) * pow(b, 6) + 2695 * pow(a, 2) * pow(b, 8) - 315 * pow(b, 10)))/(99. * sqrt(35));
			CZZ[13 * N + 27] = CZZ[27 * N + 13];
			CZZ[27 * N + 21] = (4 * sqrt(2) * a * b * (13860 * pow(a, 12) - 8190 * pow(a, 10) * (3 + 8 * pow(b, 2)) - 2002 * pow(a, 8) * (-6 - 65 * pow(b, 2) + 54 * pow(b, 4)) + 1287 * pow(a, 6) * (-1 - 56 * pow(b, 2) + 84 * pow(b, 4)) + 91 * pow(a, 2) * pow(b, 4) * (-99 + 792 * pow(b, 2) - 1430 * pow(b, 4) + 720 * pow(b, 6)) - 3 * pow(b, 6) * (-429 + 4004 * pow(b, 2) - 8190 * pow(b, 4) + 4620 * pow(b, 6)) + 9009 * pow(a, 4) * (pow(b, 2) - 12 * pow(b, 6) + 12 * pow(b, 8))))/1287.;
			CZZ[21 * N + 27] = CZZ[27 * N + 21];
			CZZ[27 * N + 23] = (840 * pow(a, 13) * b)/13.0 - (560 * pow(a, 11) * b * (2 + 7 * pow(b, 2)))/11.0 - (56 * pow(a, 9) * b * (-6 - 100 * pow(b, 2) + 3 * pow(b, 4)))/9.0 - (56 * pow(a, 5) * pow(b, 5) * (-108 + 120 * pow(b, 2) + 5 * pow(b, 4)))/15.0 + 32 * pow(a, 7) * pow(b, 3) * (-8 - 14 * pow(b, 2) + 15 * pow(b, 4)) + (56 * a * pow(b, 9) * (286 - 780 * pow(b, 2) + 495 * pow(b, 4)))/429.0 - (16 * pow(a, 3) * pow(b, 7) * (1584 - 3850 * pow(b, 2) + 2205 * pow(b, 4)))/99.;
			CZZ[23 * N + 27] = CZZ[27 * N + 23];
			CZZ[27 * N + 25] = (8 * a * b * (4158 * pow(a, 12) - 136422 * pow(a, 6) * pow(b, 4) + 63 * pow(b, 10) * (65 - 66 * pow(b, 2)) - 4095 * pow(a, 10) * (1 + 8 * pow(b, 2)) + 455 * pow(a, 2) * pow(b, 8) * (-77 + 72 * pow(b, 2)) + 1001 * pow(a, 8) * pow(b, 2) * (35 + 102 * pow(b, 2)) - 858 * pow(a, 4) * pow(b, 6) * (-159 + 119 * pow(b, 2))))/1287.;
			CZZ[25 * N + 27] = CZZ[27 * N + 25];
			CZZ[27 * N + 27] = (56 * pow(a, 13) * b)/13.0 - (560 * pow(a, 11) * pow(b, 3))/11.0 + (952 * pow(a, 9) * pow(b, 5))/3.0 - (3616 * pow(a, 7) * pow(b, 7))/7.0 + (952 * pow(a, 5) * pow(b, 9))/3.0 - (560 * pow(a, 3) * pow(b, 11))/11.0 + (56 * a * pow(b, 13))/13.;
		case 27:
			CZZ[26 * N + 4] = (32 * sqrt(0.42857142857142855) * pow(a, 3) * pow(b, 3) * (5 * pow(a, 4) - 14 * pow(a, 2) * pow(b, 2) + 5 * pow(b, 4)))/5.;
			CZZ[4 * N + 26] = CZZ[26 * N + 4];
			CZZ[26 * N + 12] = (32 * pow(a, 3) * pow(b, 3) * (140 * pow(a, 6) - 126 * pow(a, 2) * pow(b, 2) * (-3 + 2 * pow(b, 2)) + 5 * pow(b, 4) * (-27 + 28 * pow(b, 2)) - 9 * pow(a, 4) * (15 + 28 * pow(b, 2))))/(9. * sqrt(35));
			CZZ[12 * N + 26] = CZZ[26 * N + 12];
			CZZ[26 * N + 14] = (64 * pow(a, 3) * pow(b, 3) * (35 * pow(a, 6) - 117 * pow(a, 4) * pow(b, 2) + 117 * pow(a, 2) * pow(b, 4) - 35 * pow(b, 6)))/(9. * sqrt(35));
			CZZ[14 * N + 26] = CZZ[26 * N + 14];
			CZZ[26 * N + 22] = (32 * pow(a, 3) * pow(b, 3) * (4725 * pow(a, 8) - 1540 * pow(a, 6) * (5 + 3 * pow(b, 2)) - 924 * pow(a, 2) * pow(b, 2) * (9 - 15 * pow(b, 2) + 5 * pow(b, 4)) - 990 * pow(a, 4) * (-3 - 14 * pow(b, 2) + 15 * pow(b, 4)) + 5 * pow(b, 4) * (594 - 1540 * pow(b, 2) + 945 * pow(b, 4))))/495.;
			CZZ[22 * N + 26] = CZZ[26 * N + 22];
			CZZ[26 * N + 24] = (64 * pow(a, 3) * pow(b, 3) * (378 * pow(a, 8) + 1287 * pow(a, 4) * pow(b, 2) + 7 * pow(b, 6) * (55 - 54 * pow(b, 2)) - 77 * pow(a, 6) * (5 + 12 * pow(b, 2)) + 33 * pow(a, 2) * pow(b, 4) * (-39 + 28 * pow(b, 2))))/99.;
			CZZ[24 * N + 26] = CZZ[26 * N + 24];
			CZZ[26 * N + 26] = (32 * pow(a, 3) * pow(b, 3) * (441 * pow(a, 8) - 2156 * pow(a, 6) * pow(b, 2) + 3894 * pow(a, 4) * pow(b, 4) - 2156 * pow(a, 2) * pow(b, 6) + 441 * pow(b, 8)))/231.;
		case 26:
			CZZ[25 * N + 0] = (4 * sqrt(0.2857142857142857) * a * b * (18 * pow(a, 6) - 21 * pow(a, 4) * (1 + 2 * pow(b, 2)) + 3 * pow(b, 4) * (-7 + 6 * pow(b, 2)) + pow(a, 2) * (70 * pow(b, 2) - 42 * pow(b, 4))))/3.;
			CZZ[0 * N + 25] = CZZ[25 * N + 0];
			CZZ[25 * N + 3] = (4 * sqrt(0.09523809523809523) * a * b * (140 * pow(a, 8) - 240 * pow(a, 6) * (1 + pow(b, 2)) - 10 * pow(a, 2) * pow(b, 2) * (35 - 56 * pow(b, 2) + 24 * pow(b, 4)) + 5 * pow(b, 4) * (21 - 48 * pow(b, 2) + 28 * pow(b, 4)) - 7 * pow(a, 4) * (-15 - 80 * pow(b, 2) + 72 * pow(b, 4))))/5.;
			CZZ[3 * N + 25] = CZZ[25 * N + 3];
			CZZ[25 * N + 5] = (8 * a * b * (14 * pow(a, 8) + 49 * pow(a, 4) * pow(b, 2) + 15 * pow(b, 6) - 14 * pow(b, 8) - 3 * pow(a, 6) * (5 + 12 * pow(b, 2)) + pow(a, 2) * pow(b, 4) * (-49 + 36 * pow(b, 2))))/sqrt(21);
			CZZ[5 * N + 25] = CZZ[25 * N + 5];
			CZZ[25 * N + 10] = (4 * sqrt(0.05714285714285714) * a * b * (3780 * pow(a, 10) - 770 * pow(a, 8) * (11 + 6 * pow(b, 2)) - 132 * pow(a, 6) * (-45 - 110 * pow(b, 2) + 126 * pow(b, 4)) - 110 * pow(a, 2) * pow(b, 2) * (-35 + 126 * pow(b, 2) - 132 * pow(b, 4) + 42 * pow(b, 6)) - 231 * pow(a, 4) * (5 + 60 * pow(b, 2) - 132 * pow(b, 4) + 72 * pow(b, 6)) + 5 * pow(b, 4) * (-231 + 1188 * pow(b, 2) - 1694 * pow(b, 4) + 756 * pow(b, 6))))/33.;
			CZZ[10 * N + 25] = CZZ[25 * N + 10];
			CZZ[25 * N + 11] = (8 * a * b * (7560 * pow(a, 10) - 770 * pow(a, 8) * (19 + 20 * pow(b, 2)) + 99 * pow(a, 4) * pow(b, 2) * (-245 + 144 * pow(b, 4)) - 99 * pow(a, 6) * (-75 - 380 * pow(b, 2) + 144 * pow(b, 4)) + 55 * pow(a, 2) * pow(b, 4) * (441 - 684 * pow(b, 2) + 280 * pow(b, 4)) - 5 * pow(b, 6) * (1485 - 2926 * pow(b, 2) + 1512 * pow(b, 4))))/(99. * sqrt(35));
			CZZ[11 * N + 25] = CZZ[25 * N + 11];
			CZZ[25 * N + 13] = (8 * a * b * (1890 * pow(a, 10) - 385 * pow(a, 8) * (5 + 22 * pow(b, 2)) + 396 * pow(a, 6) * pow(b, 2) * (25 + 39 * pow(b, 2)) + 35 * pow(b, 8) * (-55 + 54 * pow(b, 2)) - 110 * pow(a, 2) * pow(b, 6) * (-90 + 77 * pow(b, 2)) + 198 * pow(a, 4) * pow(b, 4) * (-133 + 78 * pow(b, 2))))/(99. * sqrt(35));
			CZZ[13 * N + 25] = CZZ[25 * N + 13];
			CZZ[25 * N + 21] = (4 * sqrt(2) * a * b * (415800 * pow(a, 12) - 163800 * pow(a, 10) * (7 + 2 * pow(b, 2)) - 10010 * pow(a, 8) * (-111 - 140 * pow(b, 2) + 204 * pow(b, 4)) - 12870 * pow(a, 6) * (33 + 148 * pow(b, 2) - 392 * pow(b, 4) + 240 * pow(b, 6)) - 3003 * pow(a, 4) * (-15 - 330 * pow(b, 2) + 1332 * pow(b, 4) - 1680 * pow(b, 6) + 680 * pow(b, 8)) - 130 * pow(a, 2) * pow(b, 2) * (1155 - 7623 * pow(b, 2) + 14652 * pow(b, 4) - 10780 * pow(b, 6) + 2520 * pow(b, 8)) + 15 * pow(b, 4) * (3003 - 28314 * pow(b, 2) + 74074 * pow(b, 4) - 76440 * pow(b, 6) + 27720 * pow(b, 8))))/6435.;
			CZZ[21 * N + 25] = CZZ[25 * N + 21];
			CZZ[25 * N + 23] = (5040 * pow(a, 13) * b)/13.0 - (840 * pow(a, 11) * b * (13 + 8 * pow(b, 2)))/11.0 + 16 * pow(a, 7) * b * (-15 - 136 * pow(b, 2) + 117 * pow(b, 4)) - (56 * pow(a, 9) * b * (-136 - 325 * pow(b, 2) + 198 * pow(b, 4)))/9.0 + 16 * pow(a, 5) * pow(b, 3) * (49 - 117 * pow(b, 4) + 77 * pow(b, 6)) + (8 * pow(a, 3) * pow(b, 5) * (-9702 + 26928 * pow(b, 2) - 25025 * pow(b, 4) + 7560 * pow(b, 6)))/99.0 + a * (240 * pow(b, 7) - (7616 * pow(b, 9))/9.0 + (10920 * pow(b, 11))/11.0 - (5040 * pow(b, 13))/13.);
			CZZ[23 * N + 25] = CZZ[25 * N + 23];
			CZZ[25 * N + 25] = (2016 * pow(a, 13) * b)/13.0 - (3360 * pow(a, 11) * (b + 2 * pow(b, 3)))/11.0 + 16 * pow(a, 5) * pow(b, 5) * (133 - 156 * pow(b, 2) + 42 * pow(b, 4)) + (56 * pow(a, 9) * b * (25 + 220 * pow(b, 2) + 108 * pow(b, 4)))/9.0 - (160 * pow(a, 3) * pow(b, 7) * (495 - 847 * pow(b, 2) + 378 * pow(b, 4)))/99.0 + (32 * pow(a, 7) * pow(b, 3) * (-175 - 546 * pow(b, 2) + 468 * pow(b, 4)))/7.0 + (56 * a * pow(b, 9) * (3575 - 7020 * pow(b, 2) + 3564 * pow(b, 4)))/1287.;
		case 25:
			CZZ[24 * N + 4] = (64 * pow(a, 3) * pow(b, 3) * (-7 * pow(a, 2) + 6 * pow(a, 4) + 7 * pow(b, 2) - 6 * pow(b, 4)))/sqrt(21);
			CZZ[4 * N + 24] = CZZ[24 * N + 4];
			CZZ[24 * N + 12] = (64 * pow(a, 3) * (a - b) * pow(b, 3) * (a + b) * (315 + 280 * pow(a, 4) - 570 * pow(b, 2) + 280 * pow(b, 4) + pow(a, 2) * (-570 + 496 * pow(b, 2))))/(9. * sqrt(35));
			CZZ[12 * N + 24] = CZZ[24 * N + 12];
			CZZ[24 * N + 14] = (128 * pow(a, 3) * pow(b, 3) * (70 * pow(a, 6) - 75 * pow(b, 4) + 70 * pow(b, 6) - 18 * pow(a, 2) * pow(b, 2) * (-7 + 3 * pow(b, 2)) - 3 * pow(a, 4) * (25 + 18 * pow(b, 2))))/(9. * sqrt(35));
			CZZ[14 * N + 24] = CZZ[24 * N + 14];
			CZZ[24 * N + 22] = (64 * pow(a, 3) * pow(b, 3) * (1890 * pow(a, 8) + 1386 * pow(b, 2) - 4488 * pow(b, 4) + 5005 * pow(b, 6) - 1890 * pow(b, 8) + pow(a, 4) * (4488 - 3861 * pow(b, 2)) + 77 * pow(a, 6) * (-65 + 36 * pow(b, 2)) - 99 * pow(a, 2) * (14 - 39 * pow(b, 4) + 28 * pow(b, 6))))/99.;
			CZZ[22 * N + 24] = CZZ[24 * N + 22];
			CZZ[24 * N + 24] = (128 * pow(a, 3) * pow(b, 3) * (-10780 * pow(a, 6) + 5292 * pow(a, 8) + 1386 * pow(a, 2) * pow(b, 2) * (-7 + 6 * pow(b, 2)) + pow(a, 4) * (5775 + 8316 * pow(b, 2) - 7128 * pow(b, 4)) + 7 * pow(b, 4) * (825 - 1540 * pow(b, 2) + 756 * pow(b, 4))))/693.;
		case 24:
			CZZ[23 * N + 0] = 4 * sqrt(0.2857142857142857) * a * b * (15 * pow(a, 6) - 14 * pow(b, 2) + 28 * pow(b, 4) - 15 * pow(b, 6) + 7 * pow(a, 4) * (-4 + pow(b, 2)) - 7 * pow(a, 2) * (-2 + pow(b, 4)));
			CZZ[0 * N + 23] = CZZ[23 * N + 0];
			CZZ[23 * N + 3] = (4 * sqrt(0.09523809523809523) * a * (a - b) * b * (a + b) * (-210 + 350 * pow(a, 6) + 672 * pow(b, 2) - 825 * pow(b, 4) + 350 * pow(b, 6) + 25 * pow(a, 4) * (-33 + 26 * pow(b, 2)) + 2 * pow(a, 2) * (336 - 605 * pow(b, 2) + 325 * pow(b, 4))))/5.;
			CZZ[3 * N + 23] = CZZ[23 * N + 3];
			CZZ[23 * N + 5] = (8 * a * b * (-300 * pow(a, 6) + 175 * pow(a, 8) + 140 * pow(a, 2) * pow(b, 2) * (-1 + pow(b, 2)) - 14 * pow(a, 4) * (-9 - 10 * pow(b, 2) + 9 * pow(b, 4)) + pow(b, 4) * (126 - 300 * pow(b, 2) + 175 * pow(b, 4))))/(5. * sqrt(21));
			CZZ[5 * N + 23] = CZZ[23 * N + 5];
			CZZ[23 * N + 10] = (4 * sqrt(0.05714285714285714) * a * b * (9450 * pow(a, 10) + 3850 * pow(a, 8) * (-7 + 3 * pow(b, 2)) + 165 * pow(a, 6) * (171 - 140 * pow(b, 2) + 36 * pow(b, 4)) - 33 * pow(a, 4) * (392 - 399 * pow(b, 2) + 180 * pow(b, 6)) + pow(b, 2) * (-2310 + 12936 * pow(b, 2) - 28215 * pow(b, 4) + 26950 * pow(b, 6) - 9450 * pow(b, 8)) - 231 * pow(a, 2) * (-10 + 57 * pow(b, 4) - 100 * pow(b, 6) + 50 * pow(b, 8))))/33.;
			CZZ[10 * N + 23] = CZZ[23 * N + 10];
			CZZ[23 * N + 11] = (8 * a * b * (18900 * pow(a, 10) + 1925 * pow(a, 8) * (-25 + 4 * pow(b, 2)) - 5940 * pow(a, 6) * (-7 + 2 * pow(b, 4)) + 308 * pow(a, 2) * pow(b, 2) * (45 - 63 * pow(b, 2) + 25 * pow(b, 6)) - 198 * pow(a, 4) * (63 + 98 * pow(b, 2) - 175 * pow(b, 4) + 60 * pow(b, 6)) + 7 * pow(b, 4) * (-1782 + 5940 * pow(b, 2) - 6875 * pow(b, 4) + 2700 * pow(b, 6))))/(99. * sqrt(35));
			CZZ[11 * N + 23] = CZZ[23 * N + 11];
			CZZ[23 * N + 13] = (8 * a * b * (4725 * pow(a, 10) - 1925 * pow(a, 8) * (4 + 5 * pow(b, 2)) - 990 * pow(a, 6) * (-3 - 20 * pow(b, 2) + 9 * pow(b, 4)) + 198 * pow(a, 4) * pow(b, 2) * (-49 + 45 * pow(b, 4)) + 11 * pow(a, 2) * pow(b, 4) * (882 - 1800 * pow(b, 2) + 875 * pow(b, 4)) - 5 * pow(b, 6) * (594 - 1540 * pow(b, 2) + 945 * pow(b, 4))))/(99. * sqrt(35));
			CZZ[13 * N + 23] = CZZ[23 * N + 13];
			CZZ[23 * N + 21] = (-4 * sqrt(2) * a * b * (-1039500 * pow(a, 12) - 204750 * pow(a, 10) * (-17 + 8 * pow(b, 2)) - 250250 * pow(a, 8) * (18 - 17 * pow(b, 2) + 6 * pow(b, 4)) + 32175 * pow(a, 6) * (87 - 120 * pow(b, 2) + 68 * pow(b, 4)) + 429 * pow(a, 4) * (-1932 + 3045 * pow(b, 2) - 5100 * pow(b, 6) + 3500 * pow(b, 8)) + 65 * pow(a, 2) * (1386 - 20097 * pow(b, 4) + 59400 * pow(b, 6) - 65450 * pow(b, 8) + 25200 * pow(b, 10)) + 9 * pow(b, 2) * (-10010 + 92092 * pow(b, 2) - 311025 * pow(b, 4) + 500500 * pow(b, 6) - 386750 * pow(b, 8) + 115500 * pow(b, 10))))/6435.;
			CZZ[21 * N + 23] = CZZ[23 * N + 21];
			CZZ[23 * N + 23] = (12600 * pow(a, 13) * b)/13.0 + (8400 * pow(a, 11) * b * (-4 + pow(b, 2)))/11.0 - (280 * pow(a, 9) * b * (-116 + 40 * pow(b, 2) + 9 * pow(b, 4)))/9.0 - (480 * pow(a, 7) * b * (28 - 28 * pow(b, 4) + 15 * pow(b, 6)))/7.0 - (8 * pow(a, 5) * b * (-252 - 560 * pow(b, 2) + 1624 * pow(b, 4) - 1200 * pow(b, 6) + 175 * pow(b, 8)))/5.0 + (112 * pow(a, 3) * pow(b, 3) * (-396 + 792 * pow(b, 2) - 1100 * pow(b, 6) + 675 * pow(b, 8)))/99.0 + a * ((2016 * pow(b, 5))/5.0 - 1920 * pow(b, 7) + (32480 * pow(b, 9))/9.0 - (33600 * pow(b, 11))/11.0 + (12600 * pow(b, 13))/13.);
		case 23:
			CZZ[22 * N + 4] = (32 * pow(a, 3) * pow(b, 3) * (70 + 75 * pow(a, 4) - 140 * pow(b, 2) + 75 * pow(b, 4) + 14 * pow(a, 2) * (-10 + 9 * pow(b, 2))))/(5. * sqrt(21));
			CZZ[4 * N + 22] = CZZ[22 * N + 4];
			CZZ[22 * N + 12] = (32 * pow(a, 3) * pow(b, 3) * (-630 + 700 * pow(a, 6) + 1764 * pow(b, 2) - 1875 * pow(b, 4) + 700 * pow(b, 6) + 15 * pow(a, 4) * (-125 + 108 * pow(b, 2)) + 18 * pow(a, 2) * (98 - 175 * pow(b, 2) + 90 * pow(b, 4))))/(9. * sqrt(35));
			CZZ[12 * N + 22] = CZZ[22 * N + 12];
			CZZ[22 * N + 14] = (64 * pow(a, 3) * (a - b) * pow(b, 3) * (a + b) * (126 - 300 * pow(b, 2) + 5 * (35 * pow(a, 4) + 35 * pow(b, 4) + pow(a, 2) * (-60 + 62 * pow(b, 2)))))/(9. * sqrt(35));
			CZZ[14 * N + 22] = CZZ[22 * N + 14];
			CZZ[22 * N + 22] = (32 * pow(a, 3) * pow(b, 3) * (165375 * pow(a, 8) + 53900 * pow(a, 6) * (-10 + 9 * pow(b, 2)) + 1650 * pow(a, 4) * (406 - 756 * pow(b, 2) + 405 * pow(b, 4)) + 2772 * pow(a, 2) * (-140 + 406 * pow(b, 2) - 450 * pow(b, 4) + 175 * pow(b, 6)) + 35 * (2772 - 11088 * pow(b, 2) + 19140 * pow(b, 4) - 15400 * pow(b, 6) + 4725 * pow(b, 8))))/3465.;
		case 22:
			CZZ[21 * N + 0] = (4 * a * b * (60 * pow(a, 6) + 42 * pow(a, 4) * (-3 + 2 * pow(b, 2)) + 28 * pow(a, 2) * (3 - 5 * pow(b, 2) + 3 * pow(b, 4)) + 3 * (-7 + 28 * pow(b, 2) - 42 * pow(b, 4) + 20 * pow(b, 6))))/(3. * sqrt(7));
			CZZ[0 * N + 21] = CZZ[21 * N + 0];
			CZZ[21 * N + 3] = (4 * a * b * (315 + 1400 * pow(a, 8) - 1470 * pow(b, 2) + 3402 * pow(b, 4) - 3600 * pow(b, 6) + 1400 * pow(b, 8) + 1200 * pow(a, 6) * (-3 + 2 * pow(b, 2)) + 126 * pow(a, 4) * (27 - 40 * pow(b, 2) + 24 * pow(b, 4)) + 30 * pow(a, 2) * (-49 + 126 * pow(b, 2) - 168 * pow(b, 4) + 80 * pow(b, 6))))/(15. * sqrt(21));
			CZZ[3 * N + 21] = CZZ[21 * N + 3];
			CZZ[21 * N + 5] = (4 * sqrt(0.09523809523809523) * a * b * (700 * pow(a, 8) + 105 * pow(b, 2) - 756 * pow(b, 4) + 1350 * pow(b, 6) - 700 * pow(b, 8) + 150 * pow(a, 6) * (-9 + 4 * pow(b, 2)) - 126 * pow(a, 4) * (-6 + 5 * pow(b, 2)) - 15 * pow(a, 2) * (7 - 42 * pow(b, 4) + 40 * pow(b, 6))))/15.;
			CZZ[5 * N + 21] = CZZ[21 * N + 5];
			CZZ[21 * N + 10] = (4 * a * b * (37800 * pow(a, 10) + 38500 * pow(a, 8) * (-3 + 2 * pow(b, 2)) + 7920 * pow(a, 6) * (17 - 25 * pow(b, 2) + 15 * pow(b, 4)) + 396 * pow(a, 4) * (-189 + 476 * pow(b, 2) - 630 * pow(b, 4) + 300 * pow(b, 6)) + 22 * pow(a, 2) * (945 - 3780 * pow(b, 2) + 8568 * pow(b, 4) - 9000 * pow(b, 6) + 3500 * pow(b, 8)) + 3 * (-1155 + 6930 * pow(b, 2) - 24948 * pow(b, 4) + 44880 * pow(b, 6) - 38500 * pow(b, 8) + 12600 * pow(b, 10))))/(99. * sqrt(35));
			CZZ[10 * N + 21] = CZZ[21 * N + 10];
			CZZ[21 * N + 11] = (4 * sqrt(0.05714285714285714) * a * b * (25200 * pow(a, 10) + 7700 * pow(a, 8) * (-9 + 4 * pow(b, 2)) + 990 * pow(a, 6) * (69 - 60 * pow(b, 2) + 16 * pow(b, 4)) - 198 * pow(a, 4) * (140 - 161 * pow(b, 2) + 80 * pow(b, 6)) - 45 * pow(b, 2) * (77 - 616 * pow(b, 2) + 1518 * pow(b, 4) - 1540 * pow(b, 6) + 560 * pow(b, 8)) - 11 * pow(a, 2) * (-315 + 2898 * pow(b, 4) - 5400 * pow(b, 6) + 2800 * pow(b, 8))))/99.;
			CZZ[11 * N + 21] = CZZ[21 * N + 11];
			CZZ[21 * N + 13] = (4 * sqrt(0.05714285714285714) * a * b * (6300 * pow(a, 10) - 3850 * pow(a, 8) * (3 + 2 * pow(b, 2)) - 1980 * pow(a, 6) * (-3 - 10 * pow(b, 2) + 14 * pow(b, 4)) - 693 * pow(a, 4) * (1 + 20 * pow(b, 2) - 60 * pow(b, 4) + 40 * pow(b, 6)) - 110 * pow(a, 2) * pow(b, 2) * (-21 + 126 * pow(b, 2) - 180 * pow(b, 4) + 70 * pow(b, 6)) + 3 * pow(b, 4) * (-231 + 1980 * pow(b, 2) - 3850 * pow(b, 4) + 2100 * pow(b, 6))))/99.;
			CZZ[13 * N + 21] = CZZ[21 * N + 13];
			CZZ[21 * N + 21] = (11200 * pow(a, 13) * b)/13.0 + (11200 * pow(a, 11) * b * (-3 + 2 * pow(b, 2)))/11.0 + (560 * pow(a, 9) * b * (69 - 100 * pow(b, 2) + 60 * pow(b, 4)))/9.0 + (160 * pow(a, 7) * b * (-133 + 322 * pow(b, 2) - 420 * pow(b, 4) + 200 * pow(b, 6)))/7.0 + (16 * pow(a, 5) * b * (1071 - 3990 * pow(b, 2) + 8694 * pow(b, 4) - 9000 * pow(b, 6) + 3500 * pow(b, 8)))/15.0 + (32 * pow(a, 3) * b * (-693 + 3927 * pow(b, 2) - 13167 * pow(b, 4) + 22770 * pow(b, 6) - 19250 * pow(b, 8) + 6300 * pow(b, 10)))/99.0 + a * (28 * b - 224 * pow(b, 3) + (5712 * pow(b, 5))/5.0 - 3040 * pow(b, 7) + (12880 * pow(b, 9))/3.0 - (33600 * pow(b, 11))/11.0 + (11200 * pow(b, 13))/13.);
		case 21:
			CZZ[20 * N + 2] = (16 * a * pow(b, 3) * (7 * pow(a, 4) - 14 * pow(a, 2) * pow(b, 2) + 3 * pow(b, 4)))/(7. * sqrt(3));
			CZZ[2 * N + 20] = CZZ[20 * N + 2];
			CZZ[20 * N + 6] = (16 * sqrt(0.6666666666666666) * a * pow(b, 3) * (75 * pow(a, 6) + 5 * pow(b, 4) * (-6 + 7 * pow(b, 2)) - 7 * pow(a, 4) * (10 + 9 * pow(b, 2)) - 5 * pow(a, 2) * pow(b, 2) * (-28 + 27 * pow(b, 2))))/35.;
			CZZ[6 * N + 20] = CZZ[20 * N + 6];
			CZZ[20 * N + 8] = (16 * sqrt(0.6666666666666666) * a * pow(b, 3) * (225 * pow(a, 6) - 441 * pow(a, 4) * pow(b, 2) + 195 * pow(a, 2) * pow(b, 4) - 35 * pow(b, 6)))/105.;
			CZZ[8 * N + 20] = CZZ[20 * N + 8];
			CZZ[20 * N + 16] = (16 * a * pow(b, 3) * (-29700 * pow(a, 6) + 19250 * pow(a, 8) + 2079 * pow(a, 4) * (5 + 12 * pow(b, 2) - 20 * pow(b, 4)) + 45 * pow(b, 4) * (99 - 308 * pow(b, 2) + 210 * pow(b, 4)) - 110 * pow(a, 2) * pow(b, 2) * (189 - 486 * pow(b, 2) + 280 * pow(b, 4))))/3465.;
			CZZ[16 * N + 20] = CZZ[20 * N + 16];
			CZZ[20 * N + 18] = (16 * a * pow(b, 3) * (9625 * pow(a, 8) + 35 * pow(b, 6) * (44 - 45 * pow(b, 2)) - 9900 * pow(a, 6) * (1 + pow(b, 2)) + 220 * pow(a, 2) * pow(b, 4) * (-39 + 35 * pow(b, 2)) - 198 * pow(a, 4) * pow(b, 2) * (-98 + 55 * pow(b, 2))))/1155.;
			CZZ[18 * N + 20] = CZZ[20 * N + 18];
			CZZ[20 * N + 20] = (16 * a * pow(b, 3) * (1925 * pow(a, 8) - 5940 * pow(a, 6) * pow(b, 2) + 6534 * pow(a, 4) * pow(b, 4) - 1540 * pow(a, 2) * pow(b, 6) + 189 * pow(b, 8)))/693.;
		case 20:
			CZZ[19 * N + 1] = (16 * pow(a, 3) * b * (3 * pow(a, 4) - 14 * pow(a, 2) * pow(b, 2) + 7 * pow(b, 4)))/(7. * sqrt(3));
			CZZ[1 * N + 19] = CZZ[19 * N + 1];
			CZZ[19 * N + 7] = (16 * sqrt(0.6666666666666666) * pow(a, 3) * b * (35 * pow(a, 6) - 70 * pow(b, 4) + 75 * pow(b, 6) - 7 * pow(a, 2) * pow(b, 2) * (-20 + 9 * pow(b, 2)) - 15 * pow(a, 4) * (2 + 9 * pow(b, 2))))/35.;
			CZZ[7 * N + 19] = CZZ[19 * N + 7];
			CZZ[19 * N + 9] = (16 * sqrt(0.6666666666666666) * pow(a, 3) * b * (35 * pow(a, 6) - 195 * pow(a, 4) * pow(b, 2) + 441 * pow(a, 2) * pow(b, 4) - 225 * pow(b, 6)))/105.;
			CZZ[9 * N + 19] = CZZ[19 * N + 9];
			CZZ[19 * N + 15] = (480 * pow(a, 11) * b)/11.0 + (96 * pow(a, 5) * pow(b, 3) * (-5 + 6 * pow(b, 2)))/5.0 - (64 * pow(a, 9) * b * (9 + 20 * pow(b, 2)))/9.0 - (48 * pow(a, 7) * b * (-3 - 36 * pow(b, 2) + 28 * pow(b, 4)))/7.0 + pow(a, 3) * (48 * pow(b, 5) - (960 * pow(b, 7))/7.0 + (800 * pow(b, 9))/9.);
			CZZ[15 * N + 19] = CZZ[19 * N + 15];
			CZZ[19 * N + 17] = (16 * pow(a, 3) * b * (1575 * pow(a, 8) + 275 * pow(b, 6) * (36 - 35 * pow(b, 2)) - 1540 * pow(a, 6) * (1 + 5 * pow(b, 2)) + 396 * pow(a, 2) * pow(b, 4) * (-49 + 25 * pow(b, 2)) + 330 * pow(a, 4) * pow(b, 2) * (26 + 33 * pow(b, 2))))/1155.;
			CZZ[17 * N + 19] = CZZ[19 * N + 17];
			CZZ[19 * N + 19] = (16 * pow(a, 3) * b * (189 * pow(a, 8) - 1540 * pow(a, 6) * pow(b, 2) + 6534 * pow(a, 4) * pow(b, 4) - 5940 * pow(a, 2) * pow(b, 6) + 1925 * pow(b, 8)))/693.;
		case 19:
			CZZ[18 * N + 2] = (16 * a * pow(b, 3) * (105 * pow(a, 4) + 84 * pow(b, 2) - 75 * pow(b, 4) + 70 * pow(a, 2) * (-2 + pow(b, 2))))/(35. * sqrt(3));
			CZZ[2 * N + 18] = CZZ[18 * N + 2];
			CZZ[18 * N + 6] = (16 * sqrt(0.6666666666666666) * a * pow(b, 3) * (225 * pow(a, 6) - 168 * pow(b, 2) + 330 * pow(b, 4) - 175 * pow(b, 6) + 21 * pow(a, 4) * (-22 + 15 * pow(b, 2)) + pow(a, 2) * (280 - 308 * pow(b, 2) + 75 * pow(b, 4))))/35.;
			CZZ[6 * N + 18] = CZZ[18 * N + 6];
			CZZ[18 * N + 8] = (16 * sqrt(0.6666666666666666) * a * pow(b, 3) * (675 * pow(a, 6) + 189 * pow(a, 4) * (-4 + pow(b, 2)) + 5 * pow(b, 4) * (-36 + 35 * pow(b, 2)) + pow(a, 2) * (504 * pow(b, 2) - 375 * pow(b, 4))))/105.;
			CZZ[8 * N + 18] = CZZ[18 * N + 8];
			CZZ[18 * N + 16] = (16 * a * pow(b, 3) * (19250 * pow(a, 8) + 9900 * pow(a, 6) * (-5 + 4 * pow(b, 2)) - 66 * pow(a, 2) * (210 - 441 * pow(b, 2) + 250 * pow(b, 4)) + 99 * pow(a, 4) * (441 - 700 * pow(b, 2) + 300 * pow(b, 4)) + 7 * pow(b, 2) * (1188 - 4455 * pow(b, 2) + 5500 * pow(b, 4) - 2250 * pow(b, 6))))/1155.;
			CZZ[16 * N + 18] = CZZ[18 * N + 16];
			CZZ[18 * N + 18] = (16 * a * pow(b, 3) * (86625 * pow(a, 8) + 89100 * pow(a, 6) * (-2 + pow(b, 2)) - 594 * pow(a, 4) * (-168 + 84 * pow(b, 2) + 25 * pow(b, 4)) - 44 * pow(a, 2) * pow(b, 2) * (1512 - 2250 * pow(b, 2) + 875 * pow(b, 4)) + 15 * pow(b, 4) * (1584 - 3080 * pow(b, 2) + 1575 * pow(b, 4))))/3465.;
		case 18:
			CZZ[17 * N + 1] = (16 * pow(a, 3) * b * (75 * pow(a, 4) + 35 * pow(b, 2) * (4 - 3 * pow(b, 2)) - 14 * pow(a, 2) * (6 + 5 * pow(b, 2))))/(35. * sqrt(3));
			CZZ[1 * N + 17] = CZZ[17 * N + 1];
			CZZ[17 * N + 7] = (16 * sqrt(0.6666666666666666) * pow(a, 3) * b * (175 * pow(a, 6) - 280 * pow(b, 2) + 462 * pow(b, 4) - 225 * pow(b, 6) - 15 * pow(a, 4) * (22 + 5 * pow(b, 2)) - 7 * pow(a, 2) * (-24 - 44 * pow(b, 2) + 45 * pow(b, 4))))/35.;
			CZZ[7 * N + 17] = CZZ[17 * N + 7];
			CZZ[17 * N + 9] = (16 * sqrt(0.6666666666666666) * pow(a, 3) * b * (175 * pow(a, 6) + 63 * pow(a, 2) * pow(b, 2) * (8 + 3 * pow(b, 2)) + 27 * pow(b, 4) * (-28 + 25 * pow(b, 2)) - 15 * pow(a, 4) * (12 + 25 * pow(b, 2))))/105.;
			CZZ[9 * N + 17] = CZZ[17 * N + 9];
			CZZ[17 * N + 15] = (16 * pow(a, 3) * b * (-38500 * pow(a, 6) + 15750 * pow(a, 8) - 165 * pow(a, 4) * (-189 - 100 * pow(b, 2) + 180 * pow(b, 4)) + 11 * pow(b, 2) * (1260 - 3969 * pow(b, 2) + 4500 * pow(b, 4) - 1750 * pow(b, 6)) - 198 * pow(a, 2) * (42 + 147 * pow(b, 2) - 350 * pow(b, 4) + 200 * pow(b, 6))))/1155.;
			CZZ[15 * N + 17] = CZZ[17 * N + 15];
			CZZ[17 * N + 17] = (16 * pow(a, 3) * b * (23625 * pow(a, 8) - 7700 * pow(a, 6) * (6 + 5 * pow(b, 2)) - 990 * pow(a, 4) * (-24 - 100 * pow(b, 2) + 15 * pow(b, 4)) + 1188 * pow(a, 2) * pow(b, 2) * (-56 - 42 * pow(b, 2) + 75 * pow(b, 4)) + 99 * pow(b, 4) * (1008 - 1800 * pow(b, 2) + 875 * pow(b, 4))))/3465.;
		case 17:
			CZZ[16 * N + 2] = (16 * a * pow(b, 3) * (70 * pow(a, 4) + 140 * pow(a, 2) * (-1 + pow(b, 2)) + 3 * (35 - 84 * pow(b, 2) + 50 * pow(b, 4))))/(35. * sqrt(3));
			CZZ[2 * N + 16] = CZZ[16 * N + 2];
			CZZ[16 * N + 6] = (16 * sqrt(0.6666666666666666) * a * pow(b, 3) * (150 * pow(a, 6) + 14 * pow(a, 4) * (-28 + 27 * pow(b, 2)) + pow(a, 2) * (385 - 784 * pow(b, 2) + 450 * pow(b, 4)) + 7 * (-30 + 99 * pow(b, 2) - 120 * pow(b, 4) + 50 * pow(b, 6))))/35.;
			CZZ[6 * N + 16] = CZZ[16 * N + 6];
			CZZ[16 * N + 8] = (16 * sqrt(0.6666666666666666) * a * pow(b, 3) * (450 * pow(a, 6) - 189 * pow(b, 2) + 540 * pow(b, 4) - 350 * pow(b, 6) + 126 * pow(a, 4) * (-6 + 5 * pow(b, 2)) + 3 * pow(a, 2) * (105 - 168 * pow(b, 2) + 50 * pow(b, 4))))/105.;
			CZZ[8 * N + 16] = CZZ[16 * N + 8];
			CZZ[16 * N + 16] = (16 * a * pow(b, 3) * (38500 * pow(a, 8) + 118800 * pow(a, 6) * (-1 + pow(b, 2)) + 1188 * pow(a, 4) * (119 - 252 * pow(b, 2) + 150 * pow(b, 4)) + 88 * pow(a, 2) * (-945 + 3213 * pow(b, 2) - 4050 * pow(b, 4) + 1750 * pow(b, 6)) + 9 * (3465 - 16632 * pow(b, 2) + 33660 * pow(b, 4) - 30800 * pow(b, 6) + 10500 * pow(b, 8))))/3465.;
		case 16:
			CZZ[15 * N + 1] = (16 * b * ((30 * pow(a, 7))/7.0 + (4 * pow(a, 5) * (-9 + 5 * pow(b, 2)))/5.0 + pow(a, 3) * (3 - 4 * pow(b, 2) + 2 * pow(b, 4))))/sqrt(3);
			CZZ[1 * N + 15] = CZZ[15 * N + 1];
			CZZ[15 * N + 7] = (16 * sqrt(0.6666666666666666) * pow(a, 3) * b * (-210 + 350 * pow(a, 6) + 385 * pow(b, 2) - 392 * pow(b, 4) + 150 * pow(b, 6) + 30 * pow(a, 4) * (-28 + 15 * pow(b, 2)) + 7 * pow(a, 2) * (99 - 112 * pow(b, 2) + 54 * pow(b, 4))))/35.;
			CZZ[7 * N + 15] = CZZ[15 * N + 7];
			CZZ[15 * N + 9] = (16 * sqrt(0.6666666666666666) * pow(a, 3) * b * (350 * pow(a, 6) - 30 * pow(a, 4) * (18 + 5 * pow(b, 2)) - 63 * pow(a, 2) * (-3 - 8 * pow(b, 2) + 10 * pow(b, 4)) - 9 * pow(b, 2) * (35 - 84 * pow(b, 2) + 50 * pow(b, 4))))/105.;
			CZZ[9 * N + 15] = CZZ[15 * N + 9];
			CZZ[15 * N + 15] = (16 * pow(a, 3) * b * (94500 * pow(a, 8) + 30800 * pow(a, 6) * (-9 + 5 * pow(b, 2)) + 17820 * pow(a, 4) * (17 - 20 * pow(b, 2) + 10 * pow(b, 4)) + 2376 * pow(a, 2) * (-63 + 119 * pow(b, 2) - 126 * pow(b, 4) + 50 * pow(b, 6)) + 11 * (2835 - 7560 * pow(b, 2) + 12852 * pow(b, 4) - 10800 * pow(b, 6) + 3500 * pow(b, 8))))/3465.;
		case 15:
			CZZ[14 * N + 4] = (64 * pow(a, 3) * pow(b, 3) * (pow(a, 2) - pow(b, 2)))/sqrt(15);
			CZZ[4 * N + 14] = CZZ[14 * N + 4];
			CZZ[14 * N + 12] = (64 * pow(a, 3) * pow(b, 3) * (-21 * pow(a, 2) + 20 * pow(a, 4) + 21 * pow(b, 2) - 20 * pow(b, 4)))/21.;
			CZZ[12 * N + 14] = CZZ[14 * N + 12];
			CZZ[14 * N + 14] = (128 * (25 * pow(a, 7) * pow(b, 3) - 42 * pow(a, 5) * pow(b, 5) + 25 * pow(a, 3) * pow(b, 7)))/105.;
		case 14:
			CZZ[13 * N + 0] = (4 * sqrt(0.4) * a * b * (3 * pow(a, 4) - 10 * pow(a, 2) * pow(b, 2) + 3 * pow(b, 4)))/3.;
			CZZ[0 * N + 13] = CZZ[13 * N + 0];
			CZZ[13 * N + 3] = (4 * sqrt(0.13333333333333333) * a * b * (30 * pow(a, 6) - 70 * pow(a, 2) * pow(b, 2) * (-1 + pow(b, 2)) + 3 * pow(b, 4) * (-7 + 10 * pow(b, 2)) - 7 * pow(a, 4) * (3 + 10 * pow(b, 2))))/7.;
			CZZ[3 * N + 13] = CZZ[13 * N + 3];
			CZZ[13 * N + 5] = (8 * a * b * (15 * pow(a, 6) - 49 * pow(a, 4) * pow(b, 2) + 49 * pow(a, 2) * pow(b, 4) - 15 * pow(b, 6)))/(7. * sqrt(15));
			CZZ[5 * N + 13] = CZZ[13 * N + 5];
			CZZ[13 * N + 10] = (4 * sqrt(2) * a * b * (70 * pow(a, 8) - 30 * pow(a, 6) * (3 + 4 * pow(b, 2)) + pow(a, 4) * (21 + 210 * pow(b, 2) - 252 * pow(b, 4)) - 10 * pow(a, 2) * pow(b, 2) * (7 - 21 * pow(b, 2) + 12 * pow(b, 4)) + pow(b, 4) * (21 - 90 * pow(b, 2) + 70 * pow(b, 4))))/21.;
			CZZ[10 * N + 13] = CZZ[13 * N + 10];
			CZZ[13 * N + 11] = (8 * a * b * (140 * pow(a, 8) + 441 * pow(a, 4) * pow(b, 2) + 5 * pow(b, 6) * (27 - 28 * pow(b, 2)) - 45 * pow(a, 6) * (3 + 8 * pow(b, 2)) + 9 * pow(a, 2) * pow(b, 4) * (-49 + 40 * pow(b, 2))))/63.;
			CZZ[11 * N + 13] = CZZ[13 * N + 11];
			CZZ[13 * N + 13] = (8 * (175 * pow(a, 9) * b - 900 * pow(a, 7) * pow(b, 3) + 2394 * pow(a, 5) * pow(b, 5) - 900 * pow(a, 3) * pow(b, 7) + 175 * a * pow(b, 9)))/315.;
		case 13:
			CZZ[12 * N + 4] = (32 * pow(a, 3) * pow(b, 3) * (-5 + 4 * pow(a, 2) + 4 * pow(b, 2)))/sqrt(15);
			CZZ[4 * N + 12] = CZZ[12 * N + 4];
			CZZ[12 * N + 12] = (32 * pow(a, 3) * pow(b, 3) * (525 + 400 * pow(a, 4) - 840 * pow(b, 2) + 400 * pow(b, 4) + 168 * pow(a, 2) * (-5 + 4 * pow(b, 2))))/105.;
		case 12:
			CZZ[11 * N + 0] = 4 * sqrt(0.4) * a * b * (-5 * pow(a, 2) + 4 * pow(a, 4) + 5 * pow(b, 2) - 4 * pow(b, 4));
			CZZ[0 * N + 11] = CZZ[11 * N + 0];
			CZZ[11 * N + 3] = (4 * sqrt(0.13333333333333333) * a * b * (120 * pow(a, 6) + 14 * pow(a, 4) * (-15 + 4 * pow(b, 2)) - 7 * pow(a, 2) * (-15 + 8 * pow(b, 4)) - 15 * pow(b, 2) * (7 - 14 * pow(b, 2) + 8 * pow(b, 4))))/7.;
			CZZ[3 * N + 11] = CZZ[11 * N + 3];
			CZZ[11 * N + 5] = (8 * a * b * (60 * pow(a, 6) - 63 * pow(b, 4) + 60 * pow(b, 6) - 7 * pow(a, 4) * (9 + 4 * pow(b, 2)) + pow(a, 2) * (70 * pow(b, 2) - 28 * pow(b, 4))))/(7. * sqrt(15));
			CZZ[5 * N + 11] = CZZ[11 * N + 5];
			CZZ[11 * N + 10] = (4 * sqrt(2) * a * b * (280 * pow(a, 8) - 42 * pow(a, 4) * (-11 + 7 * pow(b, 2)) + 30 * pow(a, 6) * (-21 + 8 * pow(b, 2)) - 7 * pow(b, 2) * (-15 + 66 * pow(b, 2) - 90 * pow(b, 4) + 40 * pow(b, 6)) - 3 * pow(a, 2) * (35 - 98 * pow(b, 4) + 80 * pow(b, 6))))/21.;
			CZZ[10 * N + 11] = CZZ[11 * N + 10];
			CZZ[11 * N + 11] = (-960 * pow(a, 7) * b)/7.0 + (640 * pow(a, 9) * b)/9.0 + 16 * pow(a, 3) * pow(b, 3) * (-5 + 4 * pow(b, 2)) + pow(a, 5) * (72 * b + 64 * pow(b, 3) - (256 * pow(b, 5))/5.) + a * (72 * pow(b, 5) - (960 * pow(b, 7))/7.0 + (640 * pow(b, 9))/9.);
		case 11:
			CZZ[10 * N + 0] = (4 * a * b * (18 * pow(a, 4) + 10 * pow(a, 2) * (-3 + 2 * pow(b, 2)) + 3 * (5 - 10 * pow(b, 2) + 6 * pow(b, 4))))/(3. * sqrt(5));
			CZZ[0 * N + 10] = CZZ[10 * N + 0];
			CZZ[10 * N + 3] = (4 * a * b * (-105 + 180 * pow(a, 6) + 280 * pow(b, 2) - 378 * pow(b, 4) + 180 * pow(b, 6) + 126 * pow(a, 4) * (-3 + 2 * pow(b, 2)) + 28 * pow(a, 2) * (10 - 15 * pow(b, 2) + 9 * pow(b, 4))))/(7. * sqrt(15));
			CZZ[3 * N + 10] = CZZ[10 * N + 3];
			CZZ[10 * N + 5] = (4 * sqrt(0.13333333333333333) * a * b * (90 * pow(a, 6) - 35 * pow(b, 2) + 126 * pow(b, 4) - 90 * pow(b, 6) + 42 * pow(a, 4) * (-3 + pow(b, 2)) - 7 * pow(a, 2) * (-5 + 6 * pow(b, 4))))/7.;
			CZZ[5 * N + 10] = CZZ[10 * N + 5];
			CZZ[10 * N + 10] = 80 * pow(a, 9) * b + (480 * pow(a, 7) * b * (-3 + 2 * pow(b, 2)))/7.0 + (96 * pow(a, 5) * b * (10 - 15 * pow(b, 2) + 9 * pow(b, 4)))/5.0 + pow(a, 3) * (-80 * b + (640 * pow(b, 3))/3.0 - 288 * pow(b, 5) + (960 * pow(b, 7))/7.) + a * (20 * b - 80 * pow(b, 3) + 192 * pow(b, 5) - (1440 * pow(b, 7))/7.0 + 80 * pow(b, 9));
		case 10:
			CZZ[9 * N + 1] = (16 * sqrt(2) * pow(a, 3) * b * (3 * pow(a, 2) - 5 * pow(b, 2)))/15.;
			CZZ[1 * N + 9] = CZZ[9 * N + 1];
			CZZ[9 * N + 7] = (32 * pow(a, 3) * b * (45 * pow(a, 4) + 70 * pow(b, 2) - 63 * pow(b, 4) - 42 * pow(a, 2) * (1 + pow(b, 2))))/105.;
			CZZ[7 * N + 9] = CZZ[9 * N + 7];
			CZZ[9 * N + 9] = (32 * (5 * pow(a, 7) * b - 14 * pow(a, 5) * pow(b, 3) + 21 * pow(a, 3) * pow(b, 5)))/35.;
		case 9:
			CZZ[8 * N + 2] = (16 * sqrt(2) * a * pow(b, 3) * (5 * pow(a, 2) - 3 * pow(b, 2)))/15.;
			CZZ[2 * N + 8] = CZZ[8 * N + 2];
			CZZ[8 * N + 6] = (32 * a * pow(b, 3) * (63 * pow(a, 4) + 42 * pow(b, 2) - 45 * pow(b, 4) + 14 * pow(a, 2) * (-5 + 3 * pow(b, 2))))/105.;
			CZZ[6 * N + 8] = CZZ[8 * N + 6];
			CZZ[8 * N + 8] = (32 * (21 * pow(a, 5) * pow(b, 3) - 14 * pow(a, 3) * pow(b, 5) + 5 * a * pow(b, 7)))/35.;
		case 8:
			CZZ[7 * N + 1] = (16 * sqrt(2) * pow(a, 3) * b * (9 * pow(a, 2) + 5 * (-2 + pow(b, 2))))/15.;
			CZZ[1 * N + 7] = CZZ[7 * N + 1];
			CZZ[7 * N + 7] = (32 * pow(a, 3) * b * (135 * pow(a, 4) + 126 * pow(a, 2) * (-2 + pow(b, 2)) + 7 * (20 - 20 * pow(b, 2) + 9 * pow(b, 4))))/105.;
		case 7:
			CZZ[6 * N + 2] = (16 * sqrt(2) * a * pow(b, 3) * (-10 + 5 * pow(a, 2) + 9 * pow(b, 2)))/15.;
			CZZ[2 * N + 6] = CZZ[6 * N + 2];
			CZZ[6 * N + 6] = (32 * a * pow(b, 3) * (140 + 63 * pow(a, 4) - 252 * pow(b, 2) + 135 * pow(b, 4) + 14 * pow(a, 2) * (-10 + 9 * pow(b, 2))))/105.;
		case 6:
			CZZ[5 * N + 0] = 4 * sqrt(0.6666666666666666) * a * b * (pow(a, 2) - pow(b, 2));
			CZZ[0 * N + 5] = CZZ[5 * N + 0];
			CZZ[5 * N + 3] = (4 * sqrt(2) * a * b * (-5 * pow(a, 2) + 6 * pow(a, 4) + 5 * pow(b, 2) - 6 * pow(b, 4)))/5.;
			CZZ[3 * N + 5] = CZZ[5 * N + 3];
			CZZ[5 * N + 5] = (8 * (9 * pow(a, 5) * b - 10 * pow(a, 3) * pow(b, 3) + 9 * a * pow(b, 5)))/15.;
		case 5:
			CZZ[4 * N + 4] = (32 * pow(a, 3) * pow(b, 3))/3.;
		case 4:
			CZZ[3 * N + 0] = (4 * a * b * (-3 + 2 * pow(a, 2) + 2 * pow(b, 2)))/sqrt(3);
			CZZ[0 * N + 3] = CZZ[3 * N + 0];
			CZZ[3 * N + 3] = (4 * a * b * (45 + 36 * pow(a, 4) - 60 * pow(b, 2) + 36 * pow(b, 4) + 20 * pow(a, 2) * (-3 + 2 * pow(b, 2))))/15.;
		case 3:
			CZZ[2 * N + 2] = (16 * a * pow(b, 3))/3.;
		case 2:
			CZZ[1 * N + 1] = (16 * pow(a, 3) * b)/3.;
		default:
			CZZ[0 * N + 0] = 4 * a * b;
	}
	
	// return the array
	return CZZ;
}



void GradRectZernike1(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0;
	(*ygrad) = 0;
}

void GradRectZernike2(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 1.5231637396477693;
	(*ygrad) = 0.;
}

void GradRectZernike3(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.;
	(*ygrad) = 2.855932011839568;
}

void GradRectZernike4(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 6.430657651345178 * x;
	(*ygrad) = 6.430657651345178 * y;
}

void GradRectZernike5(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 5.606168593797669 * y;
	(*ygrad) = 5.606168593797669 * x;
}

void GradRectZernike6(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 1.8291648430492948 * x;
	(*ygrad) = -22.607780805510398 * y;
}

void GradRectZernike7(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.0 + 23.873841180512475 * x * y;
	(*ygrad) = -4.683898944411964 + 11.936920590256237 * pow(x, 2) + 35.810761770768714 * pow(y, 2);
}

void GradRectZernike8(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = -3.8449101755941104 + 21.32326517424367 * pow(x, 2) + 7.107755058081223 * pow(y, 2);
	(*ygrad) = 0.0 + 14.215510116162447 * x * y;
}

void GradRectZernike9(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.0 + 5.964411737503767 * x * y;
	(*ygrad) = 5.574706310126192 + 2.9822058687518833 * pow(x, 2) - 143.34027842053592 * pow(y, 2);
}

void GradRectZernike10(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.5448360277642879 + 6.905629974511182 * pow(x, 2) - 21.94738879110386 * pow(y, 2);
	(*ygrad) = 0.0 - 43.89477758220772 * x * y;
}

void GradRectZernike11(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (-22.765838119564258 + 55.869741830827124 * pow(x, 2) + 55.869741830827124 * pow(y, 2));
	(*ygrad) = y * (-19.801584227732615 + 55.869741830827124 * pow(x, 2) + 55.869741830827124 * pow(y, 2));
}

void GradRectZernike12(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (-6.360548431773964 + 36.982331070082836 * pow(x, 2) - 80.99788389393208 * pow(y, 2));
	(*ygrad) = y * (39.904935257063904 - 80.99788389393208 * pow(x, 2) - 198.978098857947 * pow(y, 2));
}

void GradRectZernike13(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = y * (-15.86962394752103 + 79.34811973760515 * pow(x, 2) + 26.44937324586838 * pow(y, 2));
	(*ygrad) = x * (-15.86962394752103 + 26.44937324586838 * pow(x, 2) + 79.34811973760515 * pow(y, 2));
}

void GradRectZernike14(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (0.007845004971329672 + 5.266677747602728 * pow(x, 2) - 23.91212895996972 * pow(y, 2));
	(*ygrad) = y * (-70.15202392766417 - 23.91212895996972 * pow(x, 2) + 804.5386047377497 * pow(y, 2));
}

void GradRectZernike15(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = y * (8.840863475365907 + 22.57013283787493 * pow(x, 2) - 92.98608130429997 * pow(y, 2));
	(*ygrad) = x * (8.840863475365907 + 7.5233776126249765 * pow(x, 2) - 278.9582439128999 * pow(y, 2));
}

void GradRectZernike16(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 6.827090793163168 + pow(y, 2) * (-35.28569993872314 + 31.390928361378254 * pow(y, 2)) + pow(x, 2) * (-95.36759547706328 + 156.95464180689126 * pow(x, 2) + 188.34557016826952 * pow(y, 2));
	(*ygrad) = 0.0 + x * (125.56371344551302 * pow(x, 2) * y + y * (-70.57139987744628 + 125.56371344551302 * pow(y, 2)));
}

void GradRectZernike17(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.0 + x * (213.31083472630323 * pow(x, 2) * y + y * (-99.51693810323135 + 213.31083472630323 * pow(y, 2)));
	(*ygrad) = 7.071023809633417 + pow(y, 2) * (-122.40154595207395 + 266.638543407879 * pow(y, 2)) + pow(x, 2) * (-49.75846905161568 + 53.32770868157581 * pow(x, 2) + 319.9662520894548 * pow(y, 2));
}

void GradRectZernike18(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = -0.5142579769575932 + pow(x, 2) * (-34.6294041590842 + 104.77245633900162 * pow(x, 2) - 267.56112403710193 * pow(y, 2)) + pow(y, 2) * (62.568513388372054 - 110.14153261350096 * pow(y, 2));
	(*ygrad) = 0.0 + x * (-178.37408269140127 * pow(x, 2) * y + y * (125.13702677674411 - 440.5661304540038 * pow(y, 2)));
}

void GradRectZernike19(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.0 + x * (124.10654671923044 * pow(x, 2) * y + y * (5.507240833814819 - 353.09920826596897 * pow(y, 2)));
	(*ygrad) = -6.8996500407134995 + pow(y, 2) * (290.6802486771252 - 1037.8812040639605 * pow(y, 2)) + pow(x, 2) * (2.7536204169073955 + 31.02663667980761 * pow(x, 2) - 529.6488123989534 * pow(y, 2));
}

void GradRectZernike20(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 1.0929087596150495 + pow(x, 2) * (-2.0971355427962095 + 16.90287829038131 * pow(x, 2) - 90.43875309364057 * pow(y, 2)) + pow(y, 2) * (-60.17931780586517 + 391.22650589674697 * pow(y, 2));
	(*ygrad) = 0.0 + x * (-60.29250206242705 * pow(x, 2) * y + y * (-120.35863561173034 + 1564.9060235869879 * pow(y, 2)));
}

void GradRectZernike21(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = 0.0 + x * (16.892923886536664 * pow(x, 2) * y + y * (7.573958900381285 - 99.4225919666892 * pow(y, 2)));
	(*ygrad) = 8.44852453647809 + pow(x, 2) * (3.7869794501906426 + 4.223230971634166 * pow(x, 2) - 149.13388795003357 * pow(y, 2)) + pow(y, 2) * (-589.8977457278387 + 4257.780939585816 * pow(y, 2));
}

void GradRectZernike22(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (51.480584450578064 + pow(y, 2) * (-303.68477857636503 + 354.29890830050067 * pow(y, 2)) + pow(x, 2) * (-303.0691510464788 + 354.29890830050067 * pow(x, 2) + 708.5978166010013 * pow(y, 2)));
	(*ygrad) = y * (38.492964496126234 + pow(x, 2) * (-303.68477857636503 + 354.29890830050067 * pow(x, 2)) + pow(y, 2) * (-255.21999961518137 + 708.5978166010013 * pow(x, 2) + 354.29890830050067 * pow(y, 2)));
}

void GradRectZernike23(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = y * (33.59411078744737 + pow(x, 2) * (-406.8250213607067 + 599.6146022771277 * pow(x, 2)) + pow(y, 2) * (-141.54685357801458 + 719.5375227325532 * pow(x, 2) + 119.92292045542553 * pow(y, 2)));
	(*ygrad) = x * (33.59411078744737 + pow(y, 2) * (-424.64056073404373 + 599.6146022771277 * pow(y, 2)) + pow(x, 2) * (-135.60834045356893 + 119.92292045542553 * pow(x, 2) + 719.5375227325532 * pow(y, 2)));
}

void GradRectZernike24(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (15.050040945093642 + pow(y, 2) * (349.88740952916146 - 891.4138407269894 * pow(y, 2)) + pow(x, 2) * (-207.67821733259996 + 349.9037429660642 * pow(x, 2) - 541.5100977609252 * pow(y, 2)));
	(*ygrad) = y * (-69.2139277857443 + pow(x, 2) * (349.88740952916146 - 270.7550488804626 * pow(x, 2)) + pow(y, 2) * (767.083811704245 - 1782.8276814539788 * pow(x, 2) - 1512.0726325735163 * pow(y, 2)));
}

void GradRectZernike25(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = y * (-19.37931568023186 + pow(x, 2) * (-27.428986638709034 + 351.88475676171674 * pow(x, 2)) + pow(y, 2) * (295.0125337093933 - 1168.1181196725292 * pow(x, 2) - 459.7496579098531 * pow(y, 2)));
	(*ygrad) = x * (-19.37931568023186 + pow(y, 2) * (885.03760112818 - 2298.7482895492653 * pow(y, 2)) + pow(x, 2) * (-9.142995546236307 + 70.37695135234335 * pow(x, 2) - 1168.1181196725292 * pow(y, 2)));
}

void GradRectZernike26(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (5.61926699762904 + pow(x, 2) * (-52.66177901544893 + 125.09879436132496 * pow(x, 2) - 486.0532285191596 * pow(y, 2)) + pow(y, 2) * (-113.36442868837935 + 1451.6202486236498 * pow(y, 2)));
	(*ygrad) = y * (97.58720475372974 + pow(x, 2) * (-113.36442868837935 - 243.0266142595799 * pow(x, 2)) + pow(y, 2) * (-1802.1261494891503 + 2903.2404972473 * pow(x, 2) + 5209.039383011013 * pow(y, 2)));
}

void GradRectZernike27(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = y * (13.131119469256987 + pow(x, 2) * (21.832586208232897 + 54.23625735216547 * pow(x, 2)) + pow(y, 2) * (-349.1551024421999 - 376.17344679610324 * pow(x, 2) + 1657.0356717561895 * pow(y, 2)));
	(*ygrad) = x * (13.131119469256987 + pow(x, 2) * (7.277528736077613 + 10.847251470433093 * pow(x, 2) - 376.17344679610324 * pow(y, 2)) + pow(y, 2) * (-1047.4653073265995 + 8285.178358780946 * pow(y, 2)));
}

void GradRectZernike28(double x, double y, double *xgrad, double *ygrad)
{
	(*xgrad) = x * (1.1437282618785503 + pow(x, 2) * (-3.126615806290829 + 11.462986776389926 * pow(x, 2) - 67.5519566480616 * pow(y, 2)) + pow(y, 2) * (-56.700803916530276 + 417.4562711425406 * pow(y, 2)));
	(*ygrad) = y * (-142.00489932037272 + pow(x, 2) * (-56.700803916530276 - 33.7759783240308 * pow(x, 2)) + pow(y, 2) * (4140.384057768705 + 834.9125422850811 * pow(x, 2) - 21642.30924929818 * pow(y, 2)));
}


/// HIFN Returns a pointer to the gradient of a Rectangular Zernike Polynomial
vFddpp GradRectZernikePolynomial(int i)
{
	switch(i)
	{
		case 1:
			return &GradRectZernike1;
			break;
		case 2:
			return &GradRectZernike2;
			break;
		case 3:
			return &GradRectZernike3;
			break;
		case 4:
			return &GradRectZernike4;
			break;
		case 5:
			return &GradRectZernike5;
			break;
		case 6:
			return &GradRectZernike6;
			break;
		case 7:
			return &GradRectZernike7;
			break;
		case 8:
			return &GradRectZernike8;
			break;
		case 9:
			return &GradRectZernike9;
			break;
		case 10:
			return &GradRectZernike10;
			break;
		case 11:
			return &GradRectZernike11;
			break;
		case 12:
			return &GradRectZernike12;
			break;
		case 13:
			return &GradRectZernike13;
			break;
		case 14:
			return &GradRectZernike14;
			break;
		case 15:
			return &GradRectZernike15;
			break;
		case 16:
			return &GradRectZernike16;
			break;
		case 17:
			return &GradRectZernike17;
			break;
		case 18:
			return &GradRectZernike18;
			break;
		case 19:
			return &GradRectZernike19;
			break;
		case 20:
			return &GradRectZernike20;
			break;
		case 21:
			return &GradRectZernike21;
			break;
		case 22:
			return &GradRectZernike22;
			break;
		case 23:
			return &GradRectZernike23;
			break;
		case 24:
			return &GradRectZernike24;
			break;
		case 25:
			return &GradRectZernike25;
			break;
		case 26:
			return &GradRectZernike26;
			break;
		case 27:
			return &GradRectZernike27;
			break;
		case 28:
			return &GradRectZernike28;
			break;
		default:
			return &GradRectZernike1;
			break;
	}
}
