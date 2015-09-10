//==============================================================================
//
// Title:       SLM_Zernike.h
// Purpose:     Header file with data structures and methods for using the
//              Zernike circle and Zernike rectangle polynomials
//
// Created on:  02-04-2012 at 14:45:02 by Rick van Bijnen
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

// Pointer to function types
typedef double(*dFdd)(double, double);
typedef void(*vFddpp)(double, double, double*, double*);

/// HIFN Returns a pointer to a Rectangular Zernike Polynomial, 
/// defined on a rectangle with the aspect ratio of the SLM, fitted precisely in a unit circle
dFdd RectZernikePolynomial(int i);

/// HIFN Returns a pointer to the gradient of a Rectangular Zernike Polynomial
vFddpp GradRectZernikePolynomial(int i);


double* gCZZ;
double* gMZR;
double** gZernikeRJ;

double ZernikeJ(int j, double x, double y);
double* generateCZZ(double a, double b, int N);


/// HIFN Precalculate arrays for the orthogonal Zernike rectangle polynomials
void PreCalculateZernikePolynomials(int Nx, int Ny, int NZP);

/// HIFN Sample the j-th Zernike rectangle polynomial on a Nx-by-Ny grid covering the SLM rectangle
double* SampleZernikeJ(int j, int Nx, int Ny);

/// HIFN Numerically integrates over the specific rectangle on which the Zernike polynomials are defined
double ZernikeRectangleIntegral(double* f, int Nx, int Ny);

/// HIFN Sample the j-th Zernike rectangle polynomial on a Nx-by-Ny grid covering the SLM rectangle,
///      with pre-allocated array for holding the sampling points
void FillArrayZernikeJ(double* ZJ, int j, int Nx, int Ny);


