//==============================================================================
//
// Title:       SLM_support.c
// Purpose:     Library of support functions that are used by other functions 
//              of the SLM library.
//
// Created on:  1-9-2011 at 19:46:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#include <userint.h>
#include "SLM Control Panel.h"
#include "SLM_internal.h"


/// HIFN initialises to complex valued array G (gFFTin) representing the 
/// light field in the SLM plane to given values
void SLM_setG(double* ampl, double* phase, int NG)
{
	for (int k = 0; k < NG; k++)
	{
		gFFTin[k][0] = ampl[k] * cos(phase[k]);
		gFFTin[k][1] = ampl[k] * sin(phase[k]);
	}
}


/// HIFN applies amplitude constraints to a given complex valued array G
void SLM_applyConstraints(fftw_complex* G, int NG, double* constraint)
{
	// apply constraints for each of the NG entries in G
	for (int k = 0; k < NG; k++)
	{
		// check if the amplitude is zero 
		if ((G[k][0] == 0.0) && (G[k][1] == 0.0))
		{
			// yes, we manually set G
			if (constraint != NULL)
				G[k][0] = constraint[k];
			else
				G[k][0] = 1.0;
		}
		else
		{
			// no, we have to compute the correct amplitude
		
			// the amplitude of G at point k
			double absG = sqrt(G[k][0] * G[k][0] + G[k][1] * G[k][1]);
			
			// check if there is a specified constraint (if not, assume uniform)
			if (constraint != NULL)
			{	
				// rescale the amplitude of G to the specified constraint
				G[k][0] = constraint[k] * G[k][0] / absG;
				G[k][1] = constraint[k] * G[k][1] / absG;
			}
			else
			{
				// rescale the amplitude of G to 1
				G[k][0] = G[k][0] / absG;
				G[k][1] = G[k][1] / absG;
			}
		}
	}	
}


/// HIFN performs iterative fourier transform algorithm (IFTA) to find the optimal phase pattern in the SLM plane
void IFTA_optimize(double* f, double* phasef, fftw_complex* g, double* constraintG, fftw_complex* G, fftw_plan* FFTplan_Gg, 
						  fftw_plan* FFTplan_gG, int Nf, unsigned char* signalmask, unsigned char* phasemask,
						  int Nit, int UseSoftOp)
{
	IFTA_optimize_rms(f, phasef, g, constraintG, G, FFTplan_Gg, FFTplan_gG, Nf, signalmask, phasemask, Nit, UseSoftOp, 0, NULL, 0, 0, 0, 0, 0, 0);
}



/// HIFN performs iterative fourier transform algorithm (IFTA) to find the optimal phase pattern in the SLM plane, and stores
///      the rms deviation of the signal with the desired signal at each iteration
void IFTA_optimize_rms(double* f, double* phasef, fftw_complex* g, double* constraintG, fftw_complex* G, fftw_plan* FFTplan_Gg, 
						  fftw_plan* FFTplan_gG, int Nf, unsigned char* signalmask, unsigned char* phasemask,
						  int Nit, int UseSoftOp, int MRAF, double* rmsvalues, int Nx, int Ny, int sigxmin, int sigxmax, int sigymin, int sigymax)
{
	// array for holding modulated constraints in case of amplitude modulation
	double* modConstraintG = (double*) malloc(Nf * sizeof(double));		
	
	// calculate the total intensity contained in the signal f
	double sumf = 0;
	for (int k = 0; k < Nf; k++)
		sumf += f[k] * f[k];
	
	// allocate memory for a copy of the values of G, for use with the soft operators
	fftw_complex* oldG = (fftw_complex*) fftw_malloc(Nf * sizeof(fftw_complex));
	
	// if we use soft operators, make a copy of G 
	if (UseSoftOp == 1)
	{
		for (int k = 0; k < Nf; k++)
		{
			oldG[k][0] = G[k][0];
			oldG[k][1] = G[k][1];
		}
	}
	
	// keep track of the time elapsed during this calculation
	double t0 = Timer();
	
	// perform the IFTA iterations
	for (int n = 0; n < Nit; n++)
	{
		// propagate (Fourier transform) G to g
		fftw_execute(*FFTplan_Gg);
		
		// check if we want to store the rms values
		if (rmsvalues != NULL)
		{
			rmsvalues[n] = calcRMSdiff(g, f, Nx, Ny, sigxmin, sigxmax, sigymin, sigymax); 
		}
		
		// compute amplitude correction factors
		double m1;
		double m2;
		double A;
		double mf;
		if (MRAF == 1)
		{
			m1 = 0.4;
			m2 = (1.0 - m1);
			double* Amtmp = (double*) malloc(gXsize * gYsize * sizeof(double));
			for (int k = 0; k < Nf; k++)
			{
				// check if we are in the signal window (no signal mask == the entire focal plane is window)
				if ((signalmask == NULL) || (signalmask[k]))
					Amtmp[k] += (g[k][0] * g[k][0] + g[k][1] * g[k][1]);//m1 * (g[k][0] * g[k][0] + g[k][1] * g[k][1]);	
				else
					Amtmp[k] += (g[k][0] * g[k][0] + g[k][1] * g[k][1]);//m2 * (g[k][0] * g[k][0] + g[k][1] * g[k][1]);
			}
			A = SLM_integrate_2D(Amtmp, gXsize, gYsize, 0.0, gXsize * SLM_getFocalUnitX(), 0.0, gYsize * SLM_getFocalUnitY());;
			free(Amtmp);
		}
		else
		{
			m1 = 0.0;
			for (int k = 0; k < Nf; k++)
				m1 += sqrt(g[k][0] * g[k][0] + g[k][1] * g[k][1]) * f[k];
	
			// divide by the total amplitude contained in the signal f
			m1 = m1 / sumf;
			m2 = 1.0;
			A = 1.0;
		}
	
		// apply constraints in the focal plane
		for (int k = 0; k < Nf; k++)
		{
			// check if we are in the signal window (no signal mask == the entire focal plane is window)
			if ((signalmask == NULL) || (signalmask[k]))
			{			
				// check if the signal amplitude is zero 
				if (f[k] == 0)
				{
					// yes, we can quickly set g to zero here and not worry about the phase
					g[k][0] = 0.0;
					g[k][1] = 0.0;
				}
				else
				{
					// no, we have to set the correct amplitude and take the phase into account
			
					// get the phase of g
					double phaseg;
					if ((phasemask == NULL) || (!phasemask[k]))
						phaseg = atan2(g[k][1], g[k][0]);
					else
						phaseg = phasef[k];
			
					// rescale the amplitude of g to that of f, while maintaining the phase
					g[k][0] = m1 * f[k] * cos(phaseg);
					g[k][1] = m1 * f[k] * sin(phaseg);
				}
			}
			else
			{
				// outside the signal window we only rescale the amplitude
				g[k][0] = m2 * g[k][0] / sqrt(A);
				g[k][1] = m2 * g[k][1] / sqrt(A);
			}
		}
	
		// propagate (inverse Fourier transform) g back to G 
		fftw_execute(*FFTplan_gG);
		
		// check if we want to use soft operators (See [Aagedal97])
		if (UseSoftOp == 1)
		{
			// get the amount of softening
			double beta = ((double) n + 1.0) / ((double) Nit);
			
			// apply softening to the phase
			for (int k = 0; k < Nf; k++)
			{
				G[k][0] = beta * G[k][0] + (1.0 - beta) * oldG[k][0];
				G[k][1] = beta * G[k][1] + (1.0 - beta) * oldG[k][1];
				
				// store the new-old G
				oldG[k][0] = G[k][0];
				oldG[k][1] = G[k][1];
			}
		}
		
		// check if we have to correct for amplitude modulation of the SLM
		if ((constraintG != NULL) && (fabs(gAmplitudeModulation) > 1.0e-6))
		{
			// yes, we correct the input intensity to match the 
			// effective intensity after modulation
			for (int k = 0; k < Nf; k++)
			{
				// get the phase
				double phaseG = atan2(G[k][1], G[k][0]);
				
				// get the modulation factor
				double modulation = SLM_AmplitudeModulation(phaseG);
				
				// modulate the amplitude of G
				modConstraintG[k] = modulation * constraintG[k];
			}
			
			// apply modulated constraints to G in the SLM plane
			SLM_applyConstraints(G, Nf, modConstraintG);
		}
		else
		{
			// apply constraints to G in the SLM plane
			SLM_applyConstraints(G, Nf, constraintG);
		}
	}
	
	// clean up
	fftw_free(oldG);
	free(modConstraintG);
}


/// HIFN computes the diffraction efficiency eta, as defined in the book by Turunen and Wyrowski
double calcEta(double* Ui, fftw_complex* Ush, int Nf)
{
	double sumUiUs = 0.0;
	double sumUsUs = 0.0;
	double sumUiUi = 0.0;

	// compute the numerator and denominator as in the book by Turunen and Wyrowski, p. 172, eq. (6.10)
	for (int k = 0; k < Nf; k++)
	{
		double UsUs = Ush[k][0] * Ush[k][0] + Ush[k][1] * Ush[k][1];
		
		double UiUi;
		if (Ui == NULL)
			UiUi = 1.0;
		else
			UiUi = Ui[k] * Ui[k];
			
		sumUiUs += sqrt(UsUs) * sqrt(UiUi);
		sumUsUs += UsUs;
		sumUiUi += UiUi;
	}
	
	// compute the diffraction efficiency eta
	double eta = (sumUiUs * sumUiUs) / (sumUiUi * sumUsUs);
	return eta;
}


/// HIFN calculates the (1D) integral over the intensity up to each point in the grid specified by xmin, xmax
double* SLM_cumulativeTrapz(double* intensity, double xmin, double xmax, int N)
{
	// compute the grid spacing
	double dx = (xmax - xmin) / ((double) N - 1.0);
	
	// allocate array for storing the integral values
	double* Ix = (double*) calloc(N, sizeof(double));
	Ix[0] = 0.0;
	
	// integrate the intensity up to each point on the grid (i.e. we compute N integrals)
	for (int k = 1; k < N; k++)
	{
		// use trapezoidal rule to compute next value
		Ix[k] = Ix[k - 1] + 0.5 * dx * (intensity[k - 1] + intensity[k]);
	}
	
	// return the array of integrals
	return Ix;
}


/// HIFN inverts the intensity integral I(x), i.e. finds the corresponding x at which I(x) = I
double SLM_invertIntensityIntegral(double* Ix, double I, double xmin, double xmax, int N)
{
	// initialise search boundaries (array indices) for a binary search
	int kmin = 0;
	int kmax = N - 1;

	// perform a binary search until the two entries of the array Ix are found, 
	// which are closest to the value I
	while ((kmax - kmin) > 1)
	{
		// find the midpoint (integer value is rounded down),
		// and the corresponding value of Ix
		int kmid = (kmax + kmin) / 2;
		double Imid = Ix[kmid];
		
		// compare the midpoint value against the boundaries
		if (Imid < I)
		{
			// the midpoint value is lower than I, so we have to raise
			// the lower bound (Ix should be monotonically increasing)
			kmin = kmid;
		}
		else
		{
			// the midpoint value is higher than I, lower the upper bound
			kmax = kmid;
		}
	}
	
	// compute the grid spacing
	double dx = (xmax - xmin) / ((double) N - 1.0);
	
	// interpolate between the two entries at kmin and kmax
	if (fabs(Ix[kmax] - Ix[kmin]) < 1e-10)
		return kmin * dx;
	else
		return (dx / (Ix[kmax] - Ix[kmin])) * (I - Ix[kmin]) + kmin * dx;
}


/// HIFN calculates the coordinate mapping function g used in geometric beam shaping
double* SLM_calcCoordinateMapping(double* ii, double* is, 
	double ixmin, double ixmax, int Ni,
	double sxmin, double sxmax, int Ns)
{
	// first we compute the integral over the input intensity ii, and the signal intensity is
	double* Ii = SLM_cumulativeTrapz(ii, ixmin, ixmax, Ni);
	double* Is = SLM_cumulativeTrapz(is, sxmin, sxmax, Ns);	
	
	// initialise array of mapped coordinates g = g(xi)
	double* g = (double*) malloc(Ni * sizeof(double));
	
	// next, for each value of the input coordinate ix, we find the output coordinate sx,
	// by inverting the signal intensity integral (see (6.14 - 6.18) in Wyrowski)
	for (int k = 0; k < Ni; k++)
	{
		// find the sx at which Is(sx) == Ii(ix), with ix the k-th input gridpoint
		g[k] = SLM_invertIntensityIntegral(Is, Ii[k], sxmin, sxmax, Ns);	
	}
	
	// clean up the intermediate integrals
	free(Ii);
	free(Is);
	
	// return the mapped coordinates
	return g;
}


// normalises the integral of i (over the range xmin - xmax) to 1
void SLM_normalise(double* i, int N, double xmin, double xmax)
{
	// get the total integral over x
	double* tmp = SLM_cumulativeTrapz(i, xmin, xmax, N);
	
	// normalise it such that it integrates to 1
	for (int k = 0; k < N; k++)
		i[k] = i[k] / tmp[N - 1];
	
	// clean up the temporary array
	free(tmp);
}


// normalises the integral of i^2 (over the range xmin - xmax) to 1
void SLM_normalise_sq(double* i, int N, double xmin, double xmax)
{
	// compute the squared array i
	double* isq = (double*) malloc(N * sizeof(double));
	for (int k = 0; k < N; k++)
		isq[k] = i[k] * i[k];
	
	// get the total integral over x
	double* tmp = SLM_cumulativeTrapz(isq, xmin, xmax, N);
	
	// normalise i such that its square integrates to 1
	for (int k = 0; k < N; k++)
		i[k] = i[k] / sqrt(tmp[N - 1]);
	
	// clean up the temporary arrays
	free(tmp);
	free(isq);
}


// calculatees the 2D integral of i (over the range xmin - xmax, ymin - ymax)
double SLM_integrate_2D(double* i, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
	double I;
	
	// allocate array for holding results of integrals over x, for each value of y
	double* tmpy = (double*) malloc(Ny * sizeof(double));
	
	// loop over all y coordinates
	for (int l = 0; l < Ny; l++)
	{
		// compute the cumulative integral over x, and store the final value in the tmpy array
		double* tmpx = SLM_cumulativeTrapz(&(i[l * Nx]), xmin, xmax, Nx);
		tmpy[l] = tmpx[Nx - 1];
		free(tmpx);
	}
	
	// the array tmpy now holds all the integrals over x, as a function of y
	// compute the integral over tmpy to obtain the final total integral I
	double* tmpi = SLM_cumulativeTrapz(tmpy, ymin, ymax, Ny);
	I = tmpi[Ny - 1];													  
	
	// clean up
	free(tmpy);
	free(tmpi);
	
	// return the result
	return I;
}


// normalises the 2D integral of i (over the range xmin - xmax, ymin - ymax) to 1
void SLM_normalise_2D(double* i, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
	// first, we compute the integral over i over the 2D domain
	double I = SLM_integrate_2D(i, Nx, Ny, xmin, xmax, ymin, ymax);
	
	// divide all values of i by the integral
	for (int k = 0; k < Nx * Ny; k++)
		i[k] = i[k] / I;
}


// normalises the 2D integral of i^2 (over the range xmin - xmax, ymin - ymax) to 1
void SLM_normalise_sq_2D(double* i, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
	// first, we compute the square of i
	double* isq = (double*) malloc(Nx * Ny * sizeof(double));
	for (int k = 0; k < Nx * Ny; k++)
		isq[k] = i[k] * i[k];
	
	// compute the integral over isq over the 2D domain
	double I = SLM_integrate_2D(isq, Nx, Ny, xmin, xmax, ymin, ymax);
	
	//printf("2D Integral: %.5g\n", I);
	
	// divide all values of i by the square root of the integral
	for (int k = 0; k < Nx * Ny; k++)
		i[k] = i[k] / sqrt(I);
	
	// clean up
	free(isq);
}

// normalises the 2D integral of |i|^2 for complex valued arrays (over the range xmin - xmax, ymin - ymax) to 1
void SLM_normalise_abssq_2D(fftw_complex* f, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax)
{
	// first, we compute the square of i
	double* fsq = (double*) malloc(Nx * Ny * sizeof(double));
	for (int k = 0; k < Nx * Ny; k++)
		fsq[k] = f[k][0] * f[k][0] + f[k][1] * f[k][1];
	
	// compute the integral over isq over the 2D domain
	double I = SLM_integrate_2D(fsq, Nx, Ny, xmin, xmax, ymin, ymax);
	
	// divide all values of i by the square root of the integral
	for (int k = 0; k < Nx * Ny; k++)
	{
		f[k][0] = f[k][0] / sqrt(I);
		f[k][1] = f[k][1] / sqrt(I);
	}
	
	// clean up
	free(fsq);
}


	

/// HIFN Approximates an arbitrary 2D intensity pattern by a separable one, by applying the SepOp operator (see Wyrowski book)
void SLM_SepOp(double* intensity, int Nx, int Ny, double** ix, double** iy)
{
	// get the total intensity in the input signal
	double S = 0.0;
	for (int k = 0; k < Nx * Ny; k++)
		S += intensity[k];
	
	// calculate the elements of ix
	for (int k = 0; k < Nx; k++)
	{
		// sum the intensities of the k-th column 
		for (int l = 0; l < Ny; l++)
			(*ix)[k] += intensity[l * Nx + k] / S;
	}
	
	// calculate the elements of iy
	for (int l = 0; l < Ny; l++)
	{
		// sum the intensities of the l-th row
		for (int k = 0; k < Nx; k++)
			(*iy)[l] += intensity[l * Nx + k] / S;
	}
}


/// HIFN resamples a 2D input intensity given on one grid, onto another grid, by means of bilinear interpolation
double* SLM_resampleBilinear(double* input, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// first allocate memory for the output array
	double* output = (double*) calloc(Nxo * Nyo, sizeof(double));
	
	// use the in-place resampler
	SLM_resampleBilinearInPlace(input, output, Nxi, Nyi, Nxo, Nyo);
	
	// return the output
	return output;
}


/// HIFN resamples a 2D input intensity given on one grid, onto another grid, by means of bilinear interpolation
void SLM_resampleBilinearInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// check if a resize is actually necessary
	if ((Nxi == Nxo) && (Nyi == Nyo))
	{
		// no, we can just copy the memory
		memcpy(output, input, Nxi * Nyi * sizeof(double));
	}
	else
	{
		// yes, there is work to do
		
		// the grid spacings
		double dxi = 1.0 / ((double) Nxi - 1.0);
		double dyi = 1.0 / ((double) Nyi - 1.0);
		double dxo = 1.0 / ((double) Nxo - 1.0);
		double dyo = 1.0 / ((double) Nyo - 1.0);
	
		// loop over all the pixels in the output
		for (int k = 0; k < Nxo; k++)
		for (int l = 0; l < Nyo; l++)
		{
			// get the indices of the four closest input grid points
			int ix1 = (int) (((double) k * dxo) / dxi);
			int ix2 = (ix1 == Nxi - 1 ? ix1 : ix1 + 1);
			int iy1 = (int) (((double) l * dyo) / dyi);
			int iy2 = (iy1 == Nyi - 1 ? iy1 : iy1 + 1);
		
			// get the pixel values of the four closest points
			double p11 = input[iy1 * Nxi + ix1];
			double p12 = input[iy2 * Nxi + ix1];
			double p21 = input[iy1 * Nxi + ix2];
			double p22 = input[iy2 * Nxi + ix2];
		
			// interpolate in the x direction, at y1 and y2
			double p1 = (k * dxo - ix1 * dxi) * (p21 - p11) / dxi + p11;
			double p2 = (k * dxo - ix1 * dxi) * (p22 - p12) / dxi + p12;	
		
			// interpolate in the y direction and store in output array
			output[l * Nxo + k] = (l * dyo - iy1 * dyi) * (p2 - p1) / dyi + p1;
		}
	}
}


/// HIFN Resample a phase, which takes some extra effort to avoid trouble at the branch cut at +/- pi.
double* SLM_resamplePhase(double* input, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// first allocate memory for the output array
	double* output = (double*) calloc(Nxo * Nyo, sizeof(double));
	
	// use the in-place resampler
	SLM_resamplePhaseInPlace(input, output, Nxi, Nyi, Nxo, Nyo);
	
	// return the output
	return output;
}


/// HIFN Resample a phase, which takes some extra effort to avoid trouble at the branch cut at +/- pi.
void SLM_resamplePhaseInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// to correctly interpolate phases, it is better to interpolate e^(i phi),
	// to avoid trouble when interpolating along the branch cut at pi,
	// so we split the phase into its cosine and sine, interpolate those, and then
	// recombine into a new complex number and take the phase of that
	// Note that this only works correctly when the values of the phase to be interpolated
	// are not too different from each other
	double* cosPhase = (double*) malloc(Nxi * Nyi * sizeof(double));
	double* sinPhase = (double*) malloc(Nxi * Nyi * sizeof(double));
	for (int k = 0; k < Nxi * Nyi; k++)
	{
		cosPhase[k] = cos(input[k]);
		sinPhase[k] = sin(input[k]);
	}
	double* cosPhaseResampled = SLM_resampleBitmap(cosPhase, Nxi, Nyi, Nxo, Nyo);
	double* sinPhaseResampled = SLM_resampleBitmap(sinPhase, Nxi, Nyi, Nxo, Nyo);
	
	// now combine those two in a new output phase
	for (int k = 0; k < Nxo * Nyo; k++)
		output[k] = atan2(sinPhaseResampled[k], cosPhaseResampled[k]);
	
	// clean up the temporary structures
	free(cosPhase);
	free(sinPhase);
	free(cosPhaseResampled);
	free(sinPhaseResampled);
}


/// HIFN Shrinks a 2D pixel array to a smaller number of pixels
double* SLM_reduceBitmap(double* input, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// allocate memory for the output array
	double* output = (double*) malloc(Nxo * Nyo * sizeof(double));
	
	// perform the reduction
	SLM_reduceBitmapInPlace(input, output, Nxi, Nyi, Nxo, Nyo);
	
	// return the output array
	return output;
}


/// HIFN Shrinks a 2D pixel array to a smaller number of pixels
void SLM_reduceBitmapInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// check if a resize is actually necessary
	if ((Nxi == Nxo) && (Nyi == Nyo))
	{
		// no, we can just copy the memory
		memcpy(output, input, Nxi * Nyi * sizeof(double));
	}
	else
	{
		// yes, there is work to do
		
		// find the reduction ratio (how many times should the bitmap be shrunk?)
		double Rx = ((double) Nxi) / ((double) Nxo);
		double Ry = ((double) Nyi) / ((double) Nyo);
	
		// find the largest integer less than, or equal to, the reduction ratio
		int iRx = (int) Rx;
		int iRy = (int) Ry;
	
		// find the size that such an integer reduction bitmap would have
		int iNxo = Nxi / iRx;
		int iNyo = Nyi / iRy;
	
		// allocate memory for an intermediate array, for the integer factor reduced bitmap
		double* tmpoutput = (double*) calloc(iNxo * iNyo, sizeof(double));
	
		// perform the integer factor of the reduction as a simple averaging
		for (int k = 0; k < iNxo; k++)
		for (int l = 0; l < iNyo; l++)
		{
			// first, sum all the pixels of the input bitmap falling within the current pixel
			for (int n = 0; n < iRx; n++)
			for (int m = 0; m < iRy; m++)
				tmpoutput[k + l * iNxo] += input[(k * iRx) + n + ((l * iRy) + m) * Nxi];		
		
			// divide by the number of pixels that we've summed over
			tmpoutput[k + l * iNxo] /= iRx * iRy;
		}
	
		// for the final fractional part of the reduction, we use bilinear resampling
		SLM_resampleBilinearInPlace(tmpoutput, output, iNxo, iNyo, Nxo, Nyo);
	
		// free memory for the intermediate array
		free(tmpoutput);
	}
}


/// HIFN General resampling method, uses bilinear or reduce when appropriate
double* SLM_resampleBitmap(double* input, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// allocate array for the output
	double* output = (double*) malloc(Nxo * Nyo * sizeof(double));	
	
	// perform the resampling
	SLM_resampleBitmapInPlace(input, output, Nxi, Nyi, Nxo, Nyo);
	
	// return the output
	return output;
}


/// HIFN General resampling method, uses bilinear or reduce when appropriate
void SLM_resampleBitmapInPlace(double* input, double* output, int Nxi, int Nyi, int Nxo, int Nyo)
{
	// check whether we should use reduce or bilinear
	if ((Nxo < Nxi) && (Nyo < Nyi))
		SLM_reduceBitmapInPlace(input, output, Nxi, Nyi, Nxo, Nyo);
	else
		SLM_resampleBilinearInPlace(input, output, Nxi, Nyi, Nxo, Nyo);
	
}


/// HIFN Finds the pixel coordinates and spacing of a spot pattern in an image
void SLM_findSpotPattern(double* image, int Nxpixels, int Nypixels, int Nxspots, int Nyspots, double* xpos, double* ypos, double* xspacing, double* yspacing)
{
	// create 1D arrays for Fourier transforming rows and columns of the source image
	fftw_complex* ix = (fftw_complex*) fftw_malloc(Nxpixels * sizeof(fftw_complex));
	fftw_complex* iy = (fftw_complex*) fftw_malloc(Nypixels * sizeof(fftw_complex));
	fftw_complex* fx = (fftw_complex*) fftw_malloc(Nxpixels * sizeof(fftw_complex));
	fftw_complex* fy = (fftw_complex*) fftw_malloc(Nypixels * sizeof(fftw_complex));
	
	// create fourier transform plans
	fftw_plan ixfx = fftw_plan_dft_1d(Nxpixels, ix, fx, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan iyfy = fftw_plan_dft_1d(Nypixels, iy, fy, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// set the imaginary component to zero of the input arrays
	for (int k = 0; k < Nxpixels; k++)
		ix[k][1] = 0.0;
	for (int l = 0; l < Nypixels; l++)
		iy[l][1] = 0.0;
	
	// create arrays for holding the summed absolute value squared of the 1d FFTs
	double* sx = (double*) calloc(Nxpixels - 1, sizeof(double));
	double* sy = (double*) calloc(Nypixels - 1, sizeof(double));
	
	// for each row in the input image, calculate the FFT,
	// and the slopes of the abs. squared result, and store in sx
	for (int l = 0; l < Nypixels; l++)
	{
		// copy the data to ix
		for (int k = 0; k < Nxpixels; k++)
			ix[k][0] = image[l * Nxpixels + k];
		
		// calculate the fft
		fftw_execute(ixfx);
		
		// add the difference of the absolute squared values to sx
		for (int k = 0; k < Nxpixels - 1; k++)
		{
			sx[k] += (fx[k + 1][0] * fx[k + 1][0] + fx[k + 1][1] * fx[k + 1][1]);
			sx[k] -= (fx[k + 0][0] * fx[k + 0][0] + fx[k + 0][1] * fx[k + 0][1]);
		}
	}
	
	// for each column in the input image, calculate the FFT,
	// and add the abs. squared result to sy
	for (int k = 0; k < Nxpixels; k++)
	{
		// copy the data to iy
		for (int l = 0; l < Nypixels; l++)
			iy[l][0] = image[l * Nxpixels + k];
		
		// calculate the fft
		fftw_execute(iyfy);
		
		// add the difference of the absolute squared values to sy
		for (int l = 0; l < Nypixels - 1; l++)
		{
			sy[l] += (fy[l + 1][0] * fy[l + 1][0] + fy[l + 1][1] * fy[l + 1][1]);
			sy[l] -= (fy[l + 0][0] * fy[l + 0][0] + fy[l + 0][1] * fy[l + 0][1]);
		}
	}
	
	// plot in the debuggraph
	//sx[0] = 0;
	//sy[0] = 0;
	//PlotY(DebugPanel, DebugGraph, sx, Nxpixels, VAL_DOUBLE, VAL_THIN_LINE, VAL_NO_POINT, VAL_SOLID, 1, 0x00FF0000);
	//PlotY(DebugPanel, DebugGraph, sy, Nypixels, VAL_DOUBLE, VAL_THIN_LINE, VAL_NO_POINT, VAL_SOLID, 1, 0x000000FF);
	
	// find the position of the maxima in the first half of the abs. sq. diff. of the spectrum
	int ixmax = 0;
	for (int k = 1; k < Nxpixels / 2; k++)
		if (sx[k] > sx[ixmax])
			ixmax = k;
	int iymax = 0;
	for (int l = 1; l < Nypixels / 2; l++)
		if (sy[l] > sy[iymax])
			iymax = l;
		
	// next, find the minimum of the first half, beyond the position of the maximum
	// this way, we skip the part containing a (possible) enormous 0th order peak which
	// makes the search for maxima and minima hard
	int ixmin = ixmax + 1;
	for (int k = ixmax + 2; k < Nxpixels / 2; k++)
		if (sx[k] < sx[ixmin])
			ixmin = k;
	int iymin = iymax + 1;
	for (int l = iymax + 2; l < Nypixels / 2; l++)
		if (sy[l] < sy[iymin])
			iymin = l;
		
	// we have now obtained the up and down slopes of the dominant frequency peak in the spectrum
	// (not counting the 0th order), we take the middle of that to be the actual peak
	double xpeak = ((double) ixmax + ixmin + 1) / 2.0;
	double ypeak = ((double) iymax + iymin + 1) / 2.0;
	
	// from the frequency peak we can compute the corresponding spatial period
	double periodx = ((double) Nxpixels) / xpeak;
	double periody = ((double) Nypixels) / ypeak;
	
	
	// we will now proceed with pattern matching using a convolution
	
	// create 2D arrays for Fourier transforming the observed and expected spot pattern,
	// and an array holding the result of the convolution of those two (calculated with FFTs)
	fftw_complex* camspots =     (fftw_complex*) fftw_malloc(Nxpixels * Nypixels * sizeof(fftw_complex));
	fftw_complex* fcamspots =    (fftw_complex*) fftw_malloc(Nxpixels * Nypixels * sizeof(fftw_complex));
	fftw_complex* expspots =     (fftw_complex*) fftw_malloc(Nxpixels * Nypixels * sizeof(fftw_complex));
	fftw_complex* fexpspots =    (fftw_complex*) fftw_malloc(Nxpixels * Nypixels * sizeof(fftw_complex));
	fftw_complex* convolution =  (fftw_complex*) fftw_malloc(Nxpixels * Nypixels * sizeof(fftw_complex));
	fftw_complex* fconvolution = (fftw_complex*) fftw_malloc(Nxpixels * Nypixels * sizeof(fftw_complex));
	
	// create FFT-plans for performing the convolution
	fftw_plan fftcamspots    = fftw_plan_dft_2d(Nypixels, Nxpixels, camspots, fcamspots,       FFTW_FORWARD,  FFTW_ESTIMATE);
	fftw_plan fftexpspots    = fftw_plan_dft_2d(Nypixels, Nxpixels, expspots, fexpspots,       FFTW_FORWARD,  FFTW_ESTIMATE);
	fftw_plan fftconvolution = fftw_plan_dft_2d(Nypixels, Nxpixels, fconvolution, convolution, FFTW_BACKWARD, FFTW_ESTIMATE);
									
	// number of pixels approximately covered by a spot
	int ssx = ((int) periodx + 15.0);
	int ssy = ((int) periody + 15.0);
	
	// set the fft-array of the expected spot pattern to zero
	for (int k = 0; k < Nxpixels * Nypixels; k++)
	{
		expspots[k][0] = 0.0;
		expspots[k][1] = 0.0;
	}
	
	// we estimate the sigma of a single peak to be a fraction of the period
	double sigma   = periodx / 6.0;
	double sigmasq = sigma * sigma;
	
	// fill the fft-array with the expected spot pattern
	// loop over all spots
	for (int nx = 0; nx < Nxspots; nx++)
	for (int ny = 0; ny < Nyspots; ny++)
	{
		// loop over all pixels that we estimate to be covered by this spot
		for (int kx = 0; kx < ssx; kx++)
		for (int ky = 0; ky < ssy; ky++)
		{
			// find the index in the spotpattern array of the current pixel
			int ikx = ((int) (nx * periodx + periodx / 2.0 - ssx / 2.0 + kx));
			int iky = ((int) (ny * periody + periody / 2.0 - ssy / 2.0 + ky));
			
			// check if we have found a valid pair of indices
			if ((ikx >= 0) && (ikx < Nxpixels) && (iky >= 0) && (iky < Nypixels))
			{
				// calculate the x and y coordinates of the center of this pixel
				// relative to the center of the spot
				double xr = ikx - nx * periodx - periodx / 2.0;
				double yr = iky - ny * periody - periody / 2.0;
				
				// add the intensity of the current spot to this pixel
				expspots[ikx + iky * Nxpixels][0] += exp(-(xr * xr + yr * yr) / (2.0 * sigmasq));
			}
		}
		
	}
	
	plotFFTField(expspots, Nxpixels, Nypixels, gDebugPanel, gDebugCanvas2);
	
	// fill the fft-array with the actual camera image
	for (int k = 0; k < Nxpixels * Nypixels; k++)
	{
		camspots[k][0] = image[k];
		camspots[k][1] = 0.0;
	}
	
	// perform the forward FFTs
	fftw_execute(fftcamspots);
	fftw_execute(fftexpspots);
	
	// perform pointwise multiplication of the fourier transformed spot patterns,
	// where we conjugate the fcamspots to obtain a cross-correlation rather than a convolution
	for (int k = 0; k < Nxpixels * Nypixels; k++)
	{
		fconvolution[k][0] = fcamspots[k][0] *   fexpspots[k][0] - (fcamspots[k][1]) * (-fexpspots[k][1]);
		fconvolution[k][1] = fcamspots[k][0] * (-fexpspots[k][1]) + (fcamspots[k][1]) *   fexpspots[k][0];
	}
	
	// perform the backward FFT
	fftw_execute(fftconvolution);
	
	// find the maximum of the convoluted pattern, this is where the best match between the
	// expected and actual spot pattern occurred
	int ixc = 0;
	int iyc = 0;
	for (int k = 0; k < Nxpixels; k++)
	for (int l = 0; l < Nypixels; l++)
	{
		if (convolution[l * Nxpixels + k][0] > convolution[iyc * Nxpixels + ixc][0])
		{
			ixc = k;
			iyc = l;
		}
	}
	
	plotFFTField(convolution, Nxpixels, Nypixels, gDebugPanel, gDebugCanvas3);
	
	// from the indices, calculate the actual pixel position of the upper left spot,
	// and store that in the output variables 
	(*xpos) = ((double) ixc) + periodx / 2.0;
	(*ypos) = ((double) iyc) + periody / 2.0;
	
	// also copy the lattice spacings to the output variables
	(*xspacing) = periodx;
	(*yspacing) = periody;
	
	// clean up
	fftw_free(ix);
	fftw_free(iy);
	fftw_free(fx);
	fftw_free(fy);
	fftw_destroy_plan(ixfx);
	fftw_destroy_plan(iyfy);
	fftw_free( camspots);
	fftw_free(fcamspots);
	fftw_free( expspots);
	fftw_free(fexpspots);
	fftw_free( convolution);
	fftw_free(fconvolution);
	fftw_destroy_plan(fftcamspots);
	fftw_destroy_plan(fftexpspots);
	fftw_destroy_plan(fftconvolution);
	free(sx);
	free(sy);
	
}


/// HIFN Debug function to plot the absolute value of a complex-valued array
void plotFFTField(fftw_complex* fftfield, int Nx, int Ny, int Panel, int Canvas)
{
	// allocate and fill a temporary array for holding the absolute value of the
	// copmplex valued fftfield
	double* tmp = (double*) malloc(Nx * Ny * sizeof(double));
	for (int k = 0; k < Nx * Ny; k++)
		tmp[k] = fftfield[k][0] * fftfield[k][0] + fftfield[k][1] * fftfield[k][1];
	
	// plot the absolute value
	plotField(tmp, Nx, Ny, Panel, Canvas);
	
	// clean up
	free(tmp);
}


/// HIFN Debug function to plot the phase of a complex-valued array
void plotFFTFieldPhase(fftw_complex* fftfield, int Nx, int Ny, int Panel, int Canvas)
{
	// allocate and fill a temporary array for holding the absolute value of the
	// copmplex valued fftfield
	double* tmp = (double*) malloc(Nx * Ny * sizeof(double));
	for (int k = 0; k < Nx * Ny; k++)
		tmp[k] = atan2(fftfield[k][1], fftfield[k][0]);
	
	// plot the absolute value
	plotField(tmp, Nx, Ny, Panel, Canvas);
	
	// clean up
	free(tmp);
}


/// HIFN debug function to plot an array of doubles
void plotField(double* f, int Nx, int Ny, int Panel, int Canvas)
{
	// fill the colormap, following the pattern AARRGGBB
	/*int ColorMap[256];
	for (int k = 0; k < 256; k++)
		ColorMap[k] = k * 0x00010101;
	ColorMap[255] = 0x0000FF00;*/
	int* ColorMap = SLM_CreateColorMap(255);
	
	// find the maximum value of the field
	double fmax = 0.0;
	for (int k = 0; k < Nx * Ny; k++)
	{
		if (fabs(f[k]) > fmax)
			fmax = fabs(f[k]);
	}
	
	// create an unsigned char version of the field
	unsigned char* fb = (unsigned char*) malloc(Nx * Ny * sizeof(unsigned char));
	for (int k = 0; k < Nx * Ny; k++)
		fb[k] = (unsigned char) (255.0 * fabs(f[k]) / fmax);
	
	// create a bitmap
	int bmp = -1;
	NewBitmap(-1, 8, Nx, Ny, ColorMap, fb, NULL, &bmp);

	// and copy the bitmap to the canvas
	CanvasDrawBitmap(Panel, Canvas, bmp, VAL_ENTIRE_OBJECT , VAL_ENTIRE_OBJECT);
	CanvasUpdate (Panel, Canvas, VAL_ENTIRE_OBJECT);

	
	// clean up 
	free(fb);
	DiscardBitmap(bmp);
}


/// HIFN generates a nice RGB colormap
int* SLM_CreateColorMap(int satmax)
{
	// allocate the map
	int* ColorMap = (int*) malloc(256 * sizeof(int));
	
	// fill the gColorMap, following the pattern AARRGGBB
	for (unsigned long k = 0; k < 256; k++)
	{
		// green
		unsigned long val;
		if (k < 32)
			val = 0;
		else if (k < 96)
			val = (k - 32) * 4;
		else if (k < 160)
			val = 255;
		else if (k < 224)
			val = 256 - (k - 160) * 4;
		else 
			val = 0;
		ColorMap[k] = (val < 256 ? val : 255) << 8;
	}
	for (unsigned long k = 0; k < 256; k++)
	{
		// blue
		unsigned long val;
		if (k < 32)
			val = (k + 32) * 4;
		else if (k < 96)
			val = 255;
		else if (k < 160)
			val = 256 - (k - 96) * 4;
		else 
			val = 0;
		ColorMap[k] += (val < 256 ? val : 255) << 0;
	}
	for (unsigned long k = 0; k < 256; k++)
	{
		// red
		unsigned long val;
		if (k < 96)
			val = 0;
		else if (k < 160)
			val = (k - 96) * 4;
		else if (k < 224)
			val = 255;
		else 
			val = 256 - (k - 224) * 4;
		ColorMap[k] += (val < 256 ? val : 255) << 16;
	}
	
	// set the saturated color levels to be white
	for (int k = satmax; k < 256; k++)
		ColorMap[k] = 0x00FFFFFF;
	
	//ColorMap[0] = 0;
	//for (int k = 1; k < 256; k++)
	//	ColorMap[k] = ColorMap[k - 1] + 0x00010101;
	
	return ColorMap;
}


/// HIFN Convenience function to write a scalar double value to a matlab .mat file
void writeMatDoubleScalar(MATFile *matfile, char varname[], double varvalue)
{
	mxArray* mxTemp = mxCreateDoubleScalar(varvalue);
	matPutVariable(matfile, varname, mxTemp);
	mxDestroyArray(mxTemp);
}


/// HIFN Convenience function to read a scalar double value from a matlab .mat file
double readMatDoubleScalar(MATFile *matfile, char varname[])
{
	// create and (try to) read an mxArray from the file	
	mxArray* mxTemp = matGetVariable(matfile, varname);
	
	// extract pointer to the double value
	double* dtmp = mxGetPr(mxTemp);
	
	// copy first element of the pointer
	double result = dtmp[0];
	
	// clean up the mxArray and double pointer
	mxDestroyArray(mxTemp);
	
	// return the result
	return result;
}


// Write an array of doubles as an M x N matrix to a matlab .mat file
void writeMatDoubleArray(MATFile *matfile, const char varname[], double* data, int M, int N)
{
    // create an array for the data
	mxArray *paTemp;
	paTemp = mxCreateNumericMatrix(M, N, mxDOUBLE_CLASS, mxREAL);

	// copy the data to the array
	memcpy((void*) mxGetData(paTemp), (void*) data, M * N * sizeof(double));

	// write the data to the mat file
	matPutVariable(matfile, varname, paTemp);

	// destroy the mxArray
	mxDestroyArray(paTemp);
}


/// HIFN Read an array of doubles as an M x N matrix from a matlab .mat file
double* readMatDoubleArray(MATFile *matfile, const char varname[], int *M, int *N)
{
    // create an array for reading the data from the mat file
    mxArray* mxTemp = matGetVariable(matfile, varname);

    // get the dimensions of the array
    int numrows = (int) mxGetM(mxTemp);
    int numcols = (int) mxGetN(mxTemp);

    // allocate memory for returning the array data
    double* array = (double*) malloc(numrows * numcols * sizeof(double));

    // copy the data
    memcpy(array, mxGetData(mxTemp), numrows * numcols * sizeof(double));

    // return also the dimensions
    *M = numrows;
    *N = numcols;

    // return the array
    return array;
}


/*// Write an array of integers as an M x N matrix to a matlab .mat file
void writeMatIntArray(MATFile *matfile, const char varname[], int* data, int M, int N)
{
    // convert the data to double
    double* tmp = (double*) malloc(M * N * sizeof(double));
    for (uint k = 0; k < M * N; k++)
        tmp[k] = ((double) data[k]);

    // write the double array to file
    writeMatDoubleArray(matfile, varname, tmp, M, N);

    free(tmp);
}					 */


// Write an array of unsigned chars as an M x N matrix to a matlab .mat file
void writeMatUnsignedCharArray(MATFile *matfile, const char varname[], unsigned char* data, int M, int N)
{
    // convert the data to double
    double* tmp = (double*) malloc(M * N * sizeof(double));
    for (int k = 0; k < M * N; k++)
        tmp[k] = ((double) data[k]);

    // write the double array to file
    writeMatDoubleArray(matfile, varname, tmp, M, N);

    free(tmp);
}


/// HIFN Generates an array of logarithmically spaced values
double* logspace(double xmin, double xmax, double Nx)
{
	// allocate memory for the output result
	double* result = (double*) malloc(Nx * sizeof(double));

	// get the logarithms of the minimum and maximum values
	double lxmin = log(xmin);
	double lxmax = log(xmax);
	
	// get the linearl spacing of the exponents
	double dlx = (lxmax - lxmin) / (Nx - 1);
	
	// fill the array
	for (int k = 0; k < Nx; k++)
		 result[k] = exp(lxmin + k * dlx);

	// return the array
	return result;
}	


/// HIFN writes an array of integers to a mat-file
void writeMatIntArray(MATFile *matfile, char varname[], int* varvalues, int Nv)
{
	// create an mxArray to hold the values
	mxArray *paData;
	paData = mxCreateNumericMatrix(1, Nv, mxINT32_CLASS, mxREAL);

	// copy the data to the mxArray
	memcpy((void*) mxGetData(paData), (void*) varvalues, Nv * sizeof(int));			

	// write the data to the mat file
	matPutVariable(matfile, varname, paData);

	// destroy the mxArray
	mxDestroyArray(paData);	 
}



/// HIFN writes the SLM settings to a .mat file
void SLM_WriteSettingsToFile(MATFile *pmat)
{
	// create an array for the SLM data
	mxArray *paSLMData;
	paSLMData = mxCreateNumericMatrix(SLM_getXres(), SLM_getYres(), mxDOUBLE_CLASS, mxREAL);

	// copy the SLM data to the array
	memcpy((void*) mxGetData(paSLMData), (void*) SLM_getPhase(), SLM_getXres() * SLM_getYres() * sizeof(double));			

	// write the data to the mat file
	matPutVariable(pmat, "SLM_phase", paSLMData);

	// destroy the mxArray
	mxDestroyArray(paSLMData);
	
	// create an array for the signal data
	mxArray *paSignal;
	paSignal = mxCreateNumericMatrix(gXsize, gYsize, mxDOUBLE_CLASS, mxREAL);
	
	// copy the signal data to the array
	memcpy((void*) mxGetData(paSignal), (void*) gSignal, gXsize * gYsize * sizeof(double));
	
	// write the data to the mat file
	matPutVariable(pmat, "SLM_signal", paSignal);
	
	// clean up the mxArray
	mxDestroyArray(paSignal);
	
	// write the aberration correction pattern
	writeMatDoubleArray(pmat, "SLM_aberrationcorrection", gSLMaberration, gXsize, gYsize);
	
	// write SLM settings
	writeMatDoubleScalar(pmat, "SLM_x_size", (double) gXsize);
	writeMatDoubleScalar(pmat, "SLM_y_size", (double) gYsize);
	writeMatDoubleScalar(pmat, "SLM_focallength", gFocalLength);			
	writeMatDoubleScalar(pmat, "SLM_wavelength", gWavelength);
	writeMatDoubleScalar(pmat, "SLM_Lx", LxSLM);			
	writeMatDoubleScalar(pmat, "SLM_Ly", LySLM);
	writeMatDoubleScalar(pmat, "SLM_Bias", (double) gBias);	
	writeMatDoubleScalar(pmat, "SLM_LensXphase", (double) gLensXphase);
	writeMatDoubleScalar(pmat, "SLM_LensYphase", (double) gLensYphase);	         
	writeMatDoubleScalar(pmat, "SLM_LensX", (double) gLensX);	         
	writeMatDoubleScalar(pmat, "SLM_LensY", (double) gLensY);	         
	writeMatDoubleScalar(pmat, "SLM_VertTrans", (double) gVertTrans);	         
	writeMatDoubleScalar(pmat, "SLM_HorizTrans", (double) gHorizTrans);
	writeMatDoubleScalar(pmat, "SLM_AmplitudeModulation", (double) gAmplitudeModulation);
	writeMatDoubleScalar(pmat, "SLM_SubSampleFactor", (double) gSubSampleFactor);
	writeMatDoubleScalar(pmat, "SLM_window_x",      (double) gWindowX);
	writeMatDoubleScalar(pmat, "SLM_window_y",      (double) gWindowY);
	writeMatDoubleScalar(pmat, "SLM_window_x_size", (double) gWindowXsize);
	writeMatDoubleScalar(pmat, "SLM_window_y_size", (double) gWindowYsize);
	writeMatDoubleScalar(pmat, "SLM_focal_x", gFocalUnitX);
	writeMatDoubleScalar(pmat, "SLM_focal_y", gFocalUnitY);
	
	// write input intensity settings
	writeMatDoubleScalar(pmat, "input_x_center", gIxcenter);
	writeMatDoubleScalar(pmat, "input_y_center", gIycenter);
	writeMatDoubleScalar(pmat, "input_x_sigma", gIxsigma);
	writeMatDoubleScalar(pmat, "input_y_sigma", gIysigma);
}


/// HIFN computes the cross correlation of two images, returns the offset (xpos, ypos)
///      with which the second image (im2) needs to be shifted to obtain maximum overlap
void SLM_crossCorrelate(double* im1, double* im2, int Nx, int Ny, int *xpos, int *ypos)
{
	// allocate FFT arrays for holding the image data, and the convolution
	fftw_complex* im1f = (fftw_complex*) fftw_malloc(Nx * Ny * sizeof(fftw_complex));
	fftw_complex* im2f = (fftw_complex*) fftw_malloc(Nx * Ny * sizeof(fftw_complex));
	fftw_complex* conv = (fftw_complex*) fftw_malloc(Nx * Ny * sizeof(fftw_complex));
	
	// create FFT-plans for performing the convolution
	fftw_plan fftim1  = fftw_plan_dft_2d(Ny, Nx, im1f, im1f, FFTW_FORWARD,  FFTW_ESTIMATE);
	fftw_plan fftim2  = fftw_plan_dft_2d(Ny, Nx, im2f, im2f, FFTW_FORWARD,  FFTW_ESTIMATE);
	fftw_plan fftconv = fftw_plan_dft_2d(Ny, Nx, conv, conv, FFTW_BACKWARD, FFTW_ESTIMATE);

	// fill the FFT arrays with data
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		im1f[k + l * Nx][0] = im1[k + l * Nx];
		im1f[k + l * Nx][1] = 0.0;
		im2f[k + l * Nx][0] = im2[k + l * Nx];
		im2f[k + l * Nx][1] = 0.0;
	}
	
	// fourier transform both images
	fftw_execute(fftim1);
	fftw_execute(fftim2);
	
	// perform pointwise multiplication of the fourier transformed images,
	// where we conjugate the first (im1f)
	for (int k = 0; k < Nx * Ny; k++)
	{
		conv[k][0] = im1f[k][0] * im2f[k][0] - (-im1f[k][1]) * im2f[k][1];
		conv[k][1] = im1f[k][0] * im2f[k][1] + (-im1f[k][1]) * im2f[k][0];
	}
	
	// perform the backward fourier transform and thus the cross correlation
	fftw_execute(fftconv);
	
	// find the maximum of the convolution
	*xpos = 0;
	*ypos = 0;
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		if (conv[k + l * Nx][0] > conv[*xpos + *ypos * Nx][0])
		{
			*xpos = k;
			*ypos = l;
		}
	}
	
	// we now have the coordinates to which im2 needs to be shifted for maximum overlap
	// and have stored them already in the output variables	
	
	// all that's left is to clean up the fft arrays and plans
	fftw_free(im1f);
	fftw_free(im2f);
	fftw_free(conv);
	fftw_destroy_plan(fftim1);
	fftw_destroy_plan(fftim2);		
	fftw_destroy_plan(fftconv);
}


/// HIFN Unwraps a phase, i.e. removes jumps of 2 pi, of a 2D array
void SLM_unwrapPhase(double *phase, int Nx, int Ny)
{
	// loop over all rows
	for (int l = 0; l < Ny; l++)
	{
		// the offset to add
		double offset = 0.0;
		
		// loop over all entries in this row
		for (int k = 1; k < Nx; k++)
		{
			// the current index into the phase array
			int index = k + l * Nx;
			
			// the difference between the current entry and the previous
			double diff = phase[index] + offset - phase[index - 1];
			
			// check if there is a jump larger than pi between this entry and the previous
			if (abs(diff) > PI)
			{
				// yes, determine the sign of the jump and adjust the offset
				if (diff < 0.0)
				{
					// we jumped DOWN, so we need to ADD 2 pi
					offset += 2 * PI;
				}
				else
				{
					// we jumped UP, so we need to SUBTRACT 2 pi
					offset -= 2 * PI;
				}
			}
			
			// add the current offset
			phase[index] += offset;
		}
	}
	
	// now we repeat the procedure, but traverse the array along the other dimension
	
	// loop over all columns
	for (int k = 0; k < Nx; k++)
	{
		// the offset to add
		double offset = 0.0;
		
		// loop over all entries in this row
		for (int l = 1; l < Ny; l++)
		{
			// the current index into the phase array
			int index = k + l * Nx;
			
			// the difference between the current entry and the previous
			double diff = phase[index] + offset - phase[index - Nx];
			
			// check if there is a jump larger than pi between this entry and the previous
			if (abs(diff) > PI)
			{
				// yes, determine the sign of the jump and adjust the offset
				if (diff < 0.0)
				{
					// we jumped DOWN, so we need to ADD 2 pi
					offset += 2 * PI;
				}
				else
				{
					// we jumped UP, so we need to SUBTRACT 2 pi
					offset -= 2 * PI;
				}
			}
			
			// add the current offset
			phase[index] += offset;
		}
	}

	// if there were no phase dislocations then the phase should be correctly unwrapped now
}


/// HIFN Performs a weighted average of two phase angles, taking the branch cut at 2 pi into account.
///      See Mathematica notebook 'PhaseAveraging.nb' for implementation details.
///      NOTE: the return value is always between 0 and 2 PI
double SLM_averagePhase(double t1, double t2, double a)
{
	// allocate variable to return
	double avgt;
	
	// make sure that t1 and t2 are between 0 and 2 PI
	t1 = fmod(t1, 2 * PI);
	t2 = fmod(t2, 2 * PI);
	
	// perform the averaging, by choosing the smallest angle between the two phases
	// (i.e. the average is taken either clockwise or anti-clockwise, whichever is the closest)
	if (t2 > t1)
	{
		if (t2 - t1 < PI)	
			avgt = t1 + a * (t2 - t1);
		else
			avgt = t1 - a * (2 * PI - (t2 - t1));
	}
	else
	{
		if (t1 - t2 < PI)
			avgt = t1 - a * (t1 - t2);
		else
			avgt = t1 + a * (2 * PI - (t1 - t2));
	}
	
	// return the averaged result, where we ensure it is valued between 0 and 2 PI
	return fmod(avgt, 2 * PI);
}



/// HIFN Calculates the total squared difference between two (2D) arrays of double
double calcSquaredDifference(double* f1, double *f2, int N)
{
	double fsq = 0.0;
	for (int k = 0; k < N; k++)
		fsq += pow(f1[k] - f2[k], 2);
	return fsq;
}


/// HIFN Calculates the RMS difference between the intensities represented by two arrays that represent amplitudes (so they will be squared),
///      and the difference is only taken within a certain window
double calcRMSdiff(fftw_complex* f, double* v, int Nx, int Ny, int SigXmin, int SigXmax, int SigYmin, int SigYmax)
{
	// first we compute the total intensity in the signal window
	double fsqtotal = 0.0, vsqtotal = 0.0;
	double vsqmax = 0.0;
	for (int k = SigXmin; k < SigXmax; k++)
	for (int l = SigYmin; l < SigYmax; l++)
	{
		// compute the intensities (square values)
		double fsq = pow(f[k + l * Nx][0], 2) + pow(f[k + l * Nx][1], 2); 	
		double vsq = v[k + l * Nx] * v[k + l * Nx];
		
		// update the maximum v^2
		if (vsq > vsqmax)
			vsqmax = vsq;
	
		// add current values to the total intensities
		vsqtotal += vsq;			
		fsqtotal += fsq;
	}
	
	// next, we compute the rms where we normalise each array to have the same total intensity in the window
	double rms = 0;
	for (int k = SigXmin; k < SigXmax; k++)
	for (int l = SigYmin; l < SigYmax; l++)
	{
		// compute the intensities (square values)
		double fsq = pow(f[k + l * Nx][0], 2) + pow(f[k + l * Nx][1], 2); 	
		double vsq = v[k + l * Nx] * v[k + l * Nx];
		
		// add the normalised squared difference
		rms += pow(fsq / fsqtotal - vsq / vsqtotal, 2);
	}
	
	// compute the rms
	rms = sqrt(rms / ((SigXmax - SigXmin) * (SigYmax - SigYmin))) / (vsqmax / vsqtotal);
	return rms;
}
