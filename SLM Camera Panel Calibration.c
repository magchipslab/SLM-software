//==============================================================================
//
// Title:       SLM Camera Panel Calibration
// Purpose:     Library of functions to let the computer know which pixel ends
//				up where in the camera plane.
//
// Created on:  29-07-2013 by Rick.
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

//==============================================================================
// Include files

#include <analysis.h>
#include <formatio.h>
#include <ansi_c.h>
#include <cvirte.h>     
#include <userint.h>
#include "SLM Control Panel.h"
#include "toolbox.h"
#include "SLM.h"
#include "Flycapture2/FlyCapture2_C.h"
#include "SLM Camera Panel.h"
#include "libmat/mat.h" 
#include "libmat/matrix.h"
#include "SLM_internal.h"
#include "SLM Control Panel Internal.h"


// indicates whether the user is currently marking the zeroth order spot
int gCalibrationZeroOrderMarking;

// the coordinates of a box containing the zeroth order spot (such that this region can be ignored while calibrating)
int gZeroOrderULX, gZeroOrderULY, gZeroOrderLRX, gZeroOrderLRY;

// the calibration points in SLM space, and the corresponding points in camera space
static double* CalibrationSLM_X;
static double* CalibrationSLM_Y;
static double* CalibrationCam_X;
static double* CalibrationCam_Y;
static int NumCalX;
static int NumCalY;
static int gCameraWindowOriginX, gCameraWindowOriginY;

// variable indicating whether the extensive calibration has been performed
int gCamExtendedCalibration;

// variable indicating whether the user is currently marking a calibration point on the camera canvas
int gReCalibratingPoint;


// Callback function for marking the zeroth order spot
int CVICALLBACK CamCalMarkZeroOrder_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			gCalibrationZeroOrderMarking = 1;

			break;
	}
	return 0;
}

int CVICALLBACK CamCalibration_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// get the number of calibration points
			GetCtrlVal(panel, PANEL_CamCalScanNx, &NumCalX);
			GetCtrlVal(panel, PANEL_CamCalScanNy, &NumCalY);
			
			// get the calibration area
			double width, height, offsetx, offsety;
			GetCtrlVal(panel, PANEL_CamCalScanHeight, &height);
			GetCtrlVal(panel, PANEL_CamCalScanWidth,  &width);
			GetCtrlVal(panel, PANEL_CamCalScanOffsetX,    &offsetx);
			GetCtrlVal(panel, PANEL_CamCalScanOffsetY,    &offsety);
			
			// convert from micron to meters
			width   *= 1.0e-6;
			height  *= 1.0e-6;
			offsetx *= 1.0e-6;
			offsety *= 1.0e-6;
			
			// allocate memory for the calibration points
			CalibrationSLM_X = (double*) realloc(CalibrationSLM_X, NumCalX * NumCalY * sizeof(double));
			CalibrationSLM_Y = (double*) realloc(CalibrationSLM_Y, NumCalX * NumCalY * sizeof(double));
			CalibrationCam_X = (double*) realloc(CalibrationCam_X, NumCalX * NumCalY * sizeof(double));
			CalibrationCam_Y = (double*) realloc(CalibrationCam_Y, NumCalX * NumCalY * sizeof(double));
			
			// should we clear the SLM?
			int clearSLM;
			GetCtrlVal(panel, PANEL_CamCalClearSLM, &clearSLM);
			if (clearSLM)
				SLM_Clear();
			
			// update with the 'new pattern'
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
			// store the current translation settings, as these define the current SLM origin
			double htrans = gHorizTrans;
			double vtrans = gVertTrans;
			
			// get the number of pixels surrounding the maximum to use for fitting gaussian spots
			int fitpix;
			GetCtrlVal(panel, PANEL_CamCalFitPix, &fitpix); 
			int numxfit = fitpix;
			int numyfit = fitpix;
			
			// arrays for holding the integrated intensities for fitting gaussian spots
			// and the associated camera coordinates
			double* tmpx = (double*) malloc((2 * numxfit + 1) * sizeof(double));
			double* tmpy = (double*) malloc((2 * numyfit + 1) * sizeof(double));
			double* xfit = (double*) malloc((2 * numxfit + 1) * sizeof(double));
			double* yfit = (double*) malloc((2 * numyfit + 1) * sizeof(double));
			
			// the initial guess for the spot parameters
			double sigmax0 = 0.33 * numxfit * gCamPixelSize;
			double sigmay0 = 0.33 * numxfit * gCamPixelSize;
			
			// fill the arrays of calibration points 
			for (int l = 0; l < NumCalY; l++)
			for (int k = 0; k < NumCalX; k++)
			{
				
				// update the progress indicator
				char tmp[100];
				sprintf(tmp, "%i / %i", k + l * NumCalX + 1, NumCalX * NumCalY);
				SetCtrlVal(pnlCamera, PANEL_CalibrationProgress, tmp);
				
				// set the calibration point in SLM space
				CalibrationSLM_X[k + l * NumCalX] = k * width  / ((double) (NumCalX - 1)) + offsetx;
				CalibrationSLM_Y[k + l * NumCalX] = l * height / ((double) (NumCalY - 1)) + offsety;
				
				// let the SLM place a spot at the calibration point
				SLM_setHorizTrans(htrans + CalibrationSLM_X[k + l * NumCalX]);
				SLM_setVertTrans( vtrans + CalibrationSLM_Y[k + l * NumCalX]);
				
				// update the SLM translation and lens phases (not the base pattern)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);
				
				// wait a short while to make sure the image is adjusted to the new situation
				// (let the liquid crystals of the SLM settle or whatever they seem to be doing when you set an image)
				Delay(0.5);
				
				// call the camera relevant part of the camera timer callback function manually, since the camera 
				// timer events are not processed when the feedback routine is running
				for (int c = 0; c < getCameraNumAveragingFrames(); c++)
				{
					// first wait one frame
					Delay(1.0 / getCameraFrameRate());
					
					// grab a frame
					grabCameraFrame();
				}
				
				// average the camera frames
				performCameraAveraging();
				
				// erase the zeroth order spot from the camera image
				for (int n = gZeroOrderULX; n < gZeroOrderLRX; n++)
				for (int m = gZeroOrderULY; m < gZeroOrderLRY; m++)
					gAvgFrame[n + gCamX * m] = 0;
				
				// find the maximum intensity
				double imax = 0.0;
				int xmax, ymax;
				for (int n = 0; n < gCamX; n++)
				for (int m = 0; m < gCamY; m++)
				{
					if (gAvgFrame[n + gCamX * m] > imax)
					{
						// store the new maximum and its location
						imax = gAvgFrame[n + gCamX * m];
						xmax = n;
						ymax = m;
					}
				}
				
				// extract and integrate the surrounding area
				memset(tmpx, 0, (2 * numxfit + 1) * sizeof(double));
				memset(tmpy, 0, (2 * numyfit + 1) * sizeof(double));
				for (int n = 0; n < 2 * numxfit + 1; n++)
				for (int m = 0; m < 2 * numyfit + 1; m++)					
				{
					// add the current intensity (check for overshoot on array boundaries first)
					int i = (n + xmax - numxfit);
					int j = (m + ymax - numyfit);
					if ((i > 0) && (j > 0) && (i < gCamX) && (j < gCamY))
					{
						tmpx[n] += gAvgFrame[i + j * gCamX];
						tmpy[m] += gAvgFrame[i + j * gCamX];						
					}
					
					// also compute the associated coordinates
					xfit[n] = (i + 0.5) * gCamPixelSize;
					yfit[m] = (j + 0.5) * gCamPixelSize;
				}
				
				// fit the integrated intensities with a 1D gaussian
				double x, y, sigmax, sigmay, amplx, amply, resx, resy;
				double c0x[3], c0y[3];
				c0x[0] = imax; c0x[1] = xmax * gCamPixelSize; c0x[2] = sigmax0;
				c0y[0] = imax; c0y[1] = ymax * gCamPixelSize; c0y[2] = sigmay0;
				GaussFit (xfit, tmpx, NULL, 2 * numxfit + 1, BISQUARE, 1e-4, c0x, NULL, &amplx, &x, &sigmax, &resx);
				GaussFit (yfit, tmpy, NULL, 2 * numyfit + 1, BISQUARE, 1e-4, c0y, NULL, &amply, &y, &sigmay, &resy);
				
				// store the new sigmas as a new initial guess for the next fit
				sigmax0 = sigmax;
				sigmay0 = sigmay;
				
				// get the center coordinates of the gauss, and store these in the camera calibration points array
				CalibrationCam_X[k + l * NumCalX] = x;
				CalibrationCam_Y[k + l * NumCalX] = y;
			}
			
			// finally, compute the pixel at the camera window origin
			double cwox, cwoy;
			transformCoordinateSLMToCamera(0.0, 0.0, &cwox, &cwoy);
			gCameraWindowOriginX = (int) (cwox / gCamPixelSize + 0.5);
			gCameraWindowOriginY = (int) (cwoy / gCamPixelSize + 0.5);
			
			// clean up the temporary arrays for the gaussian fit
			free(tmpx);
			free(tmpy);
			free(xfit);
			free(yfit);
			
			// indicate that the calibration has been performed
			gCamExtendedCalibration = 1;
			SetCtrlAttribute(panel, PANEL_CamCalToggle,  ATTR_CTRL_VAL, gCamExtendedCalibration);
			
			// reset the original translational settings
			SLM_setHorizTrans(htrans);
			SLM_setVertTrans(vtrans);
				
			// update the SLM translation and lens phases (not the base pattern), placing the spot at the origin
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 0);

			break;
	}
	return 0;
}


// get the coordinates of the pixel at the origin of the camera window
void getCameraWindowOriginPixel(int *x, int *y)
{
	*x = gCameraWindowOriginX;
	*y = gCameraWindowOriginY;
}


// get the SLM coordinates of a calibration point
void getCalibrationPointSLM(int p, double *x, double *y)
{
	*x = CalibrationSLM_X[p];
	*y = CalibrationSLM_Y[p];
}


// draws the camera calibration points on the canvas, using the current pen
void drawCameraCalibrationPoints(int panel, int canvas)
{
	for (int k = 0; k < NumCalX; k++)
	for (int l = 0; l < NumCalY; l++)
	{
		// compute the coordinates of the current calibration point on the bitmap
		double bx = (CalibrationCam_X[k + l * NumCalX] / gCamPixelSize);
		double by = (CalibrationCam_Y[k + l * NumCalX] / gCamPixelSize);
		
		int x = TransformBitmapDoubleToCamCanvasX(bx);
		int y = TransformBitmapDoubleToCamCanvasY(by);
		
		// draw the calibration point as a small cross
		CanvasDrawLine(panel, canvas, MakePoint(x, y - 2), MakePoint(x, y + 2));
		CanvasDrawLine(panel, canvas, MakePoint(x - 2, y), MakePoint(x + 2, y));
	}
}



int CVICALLBACK CamCalToggle_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			gCamExtendedCalibration = (gCamExtendedCalibration ? 0 : 1);
			SetCtrlAttribute(panel, PANEL_CamCalToggle,  ATTR_CTRL_VAL, gCamExtendedCalibration);

			break;
	}
	return 0;
}


// transforms a coordinate from SLM space to camera space
// NOTE: slow code, not intended for heavy realtime use (use for precalculation)
void transformCoordinateSLMToCamera(double xfrom, double yfrom, double* xto, double* yto)
{
	transformCoordinate(CalibrationSLM_X, CalibrationSLM_Y, CalibrationCam_X, CalibrationCam_Y, NumCalX, NumCalY, xfrom, yfrom, xto, yto);
}


// transforms a coordinate from camera space to SLM space
// NOTE: slow code, not intended for heavy realtime use (use for precalculation)
void transformCoordinateCameraToSLM(double xfrom, double yfrom, double* xto, double* yto)
{
	transformCoordinate(CalibrationCam_X, CalibrationCam_Y, CalibrationSLM_X, CalibrationSLM_Y, NumCalX, NumCalY, xfrom, yfrom, xto, yto);
}


// transforms a coordinate with linear interpolation based on two grids which are in correspondence to each other
// NOTE: slow code, not intended for heavy realtime use (use for precalculation)
void transformCoordinate(double* MapFromX, double* MapFromY, double* MapToX, double* MapToY, int Nx, int Ny, double xfrom, double yfrom, double* xto, double* yto)
{
	// find the indices of the nearest coordinate in the 'from' map by a simple linear search
	int ix, iy, i0;
	double distsq = 1.0e10;
	for (int k = 0; k < Nx; k++)
	for (int l = 0; l < Ny; l++)
	{
		double xf = MapFromX[k + l * Nx] - xfrom;
		double yf = MapFromY[k + l * Nx] - yfrom;
		double tmp = xf * xf + yf * yf;
		if (tmp < distsq)
		{
			distsq = tmp;
			ix = k; 
			iy = l;
			i0 = k + l * Nx;
		}
	}
	
	// get the actual nearest coordinate, we use this as the origin of a local coordinate system
	double x0 = MapFromX[i0];
	double y0 = MapFromY[i0];
	
	// get the position of our to-be-mapped coordinate with respect to this origin
	double x = xfrom - x0;
	double y = yfrom - y0;
	
	// get the four vectors connecting the adjacent points
	// for each vector we check whether the coordinate is contained in the map, else we extrapolate the opposite vector negatively
	double v1x, v1y, v2x, v2y, v3x, v3y, v4x, v4y;
	int i1, i2, i3, i4;
	
	// v1 (right)
	if (ix < Nx - 1)
	{
		i1 = (ix + 1) + iy * Nx;
		v1x = MapFromX[i1] - x0;
		v1y = MapFromY[i1] - y0;
	}
	else
	{
		i1 = (ix - 1) + iy * Nx;
		v1x = -(MapFromX[i1] - x0);
		v1y = -(MapFromY[i1] - y0);
	}
	
	// v2 (up)
	if (iy < Ny - 1)
	{
		i2 = ix + (iy + 1) * Nx;
		v2x = MapFromX[i2] - x0;
		v2y = MapFromY[i2] - y0;
	}
	else
	{
		i2 = ix + (iy - 1) * Nx;
		v2x = -(MapFromX[i2] - x0);
		v2y = -(MapFromY[i2] - y0);
	}
	
	// v3 (left)
	if (ix > 0)
	{
		i3 = (ix - 1) + iy * Nx;
		v3x = MapFromX[i3] - x0;
		v3y = MapFromY[i3] - y0;
	}
	else
	{
		i3 = (ix + 1) + iy * Nx;
		v3x = -(MapFromX[i3] - x0);
		v3y = -(MapFromY[i3] - y0);
	}
	
	// v4 (down)
	if (iy > 0)
	{
		i4 = ix + (iy - 1) * Nx;
		v4x = MapFromX[i4] - x0;
		v4y = MapFromY[i4] - y0;
	}
	else
	{
		i4 = ix + (iy + 1) * Nx;
		v4x = -(MapFromX[i4] - x0);
		v4y = -(MapFromY[i4] - y0);
	}
	
	// find out in which 'quadrant' we are by checking where the point lies with respect to the four vectors
	// connecting the adjacent points, by computing the cross products
	double q1, q2, q3, q4;
	q1 = v1x * y - v1y * x;
	q2 = v2x * y - v2y * x;
	q3 = v3x * y - v3y * x;
	q4 = v4x * y - v4y * x;
	
	
	
	// TODO: The quadrant search can be done also more efficiently, instead of first computing all the cross products,
	//       and then compare, we can also compute and compare immediately, and stop when we have found the sign change
	
	
	
	
	// we are in the quadrant where the cross product changes sign from plus to minus
	// (should work as long as the angles are less than 180 degrees?)
	// determine the four corner points of this quadrant
	double x1, y1, x2, y2, x3, y3, x4, y4;
	double u1, v1, u2, v2, u3, v3, u4, v4;
	if ((q1 >= 0) && (q2 < 0))
	{
		// first quadrant
		
		// set the three coordinates we already know
		x2 = v2x;
		y2 = v2y;
		x3 = 0;
		y3 = 0;
		x4 = v1x;
		y4 = v1y;
		
		// for the last one we have to check if it falls in the map
		if ((ix < Nx - 1) && (iy < Ny - 1))
		{
			// yes, copy from map
			x1 = MapFromX[ix + 1 + (iy + 1) * Nx] - x0;
			y1 = MapFromY[ix + 1 + (iy + 1) * Nx] - y0;
			
			// we can now also safely copy the 4 corresponding coordinates from the target space
			u1 = MapToX[ix + 1 + (iy + 1) * Nx];
			v1 = MapToY[ix + 1 + (iy + 1) * Nx];
			u2 = MapToX[ix + 0 + (iy + 1) * Nx];
			v2 = MapToY[ix + 0 + (iy + 1) * Nx];
			u3 = MapToX[ix + 0 + (iy + 0) * Nx];
			v3 = MapToY[ix + 0 + (iy + 0) * Nx];
			u4 = MapToX[ix + 1 + (iy + 0) * Nx];
			v4 = MapToY[ix + 1 + (iy + 0) * Nx];
		}
		else
		{
			// no, guess a bit by using the other coordinates
			x1 = x4;
			y1 = y2;
			
			// TODO: fix mapped points!!
			u3 = MapToX[ix + 0 + (iy + 0) * Nx];
			v3 = MapToY[ix + 0 + (iy + 0) * Nx];
			u1 = u3 + 1.0;
			v1 = v3 + 1.0;
			u2 = u3;
			v2 = v3 + 1.0;
			u4 = u3 + 1.0;
			v4 = v3;
		}
	}
	else if ((q2 >= 0) && (q3 < 0))
	{
		// second quadrant
		
		// set the coordinates we already know
		x1 = v2x;
		y1 = v2y;
		x3 = v3x;
		y3 = v3y;
		x4 = 0;
		y4 = 0;
		
		// for the last one we have to check if it falls in the map
		if ((ix > 0) && (iy < Ny - 1))
		{
			// yes, copy from map
			x2 = MapFromX[ix - 1 + (iy + 1) * Nx] - x0;
			y2 = MapFromY[ix - 1 + (iy + 1) * Nx] - y0;	
			
			// we can now also safely copy the 4 corresponding coordinates from the target space
			u1 = MapToX[ix + 0 + (iy + 1) * Nx];
			v1 = MapToY[ix + 0 + (iy + 1) * Nx];
			u2 = MapToX[ix - 1 + (iy + 1) * Nx];
			v2 = MapToY[ix - 1 + (iy + 1) * Nx];
			u3 = MapToX[ix - 1 + (iy + 0) * Nx];
			v3 = MapToY[ix - 1 + (iy + 0) * Nx];
			u4 = MapToX[ix + 0 + (iy + 0) * Nx];
			v4 = MapToY[ix + 0 + (iy + 0) * Nx];
		}
		else
		{
			// no, guess a bit by using the other coordinates
			x2 = x3;
			y2 = y1;
			
			// TODO: fix mapped points!!
			u4 = MapToX[ix + 0 + (iy + 0) * Nx];
			v4 = MapToY[ix + 0 + (iy + 0) * Nx];
			u1 = u4;
			v1 = v4 + 1.0;
			u2 = u4 - 1.0;
			v2 = v4 + 1.0;
			u3 = u4 - 1.0;
			v3 = v4;
		}
	}
	else if ((q3 >= 0) && (q4 < 0))
	{
		// third quadrant
		
		// set the coordinates we already know
		x1 = 0;
		y1 = 0;
		x2 = v3x;
		y2 = v3y;
		x4 = v4x;
		y4 = v4y;
		
		// for the last one we have to check if it falls in the map
		if ((ix > 0) && (iy > 0))
		{
			// yes, copy from map
			x3 = MapFromX[ix - 1 + (iy - 1) * Nx] - x0;
			y3 = MapFromY[ix - 1 + (iy - 1) * Nx] - y0;	
			
			// we can now also safely copy the 4 corresponding coordinates from the target space
			u1 = MapToX[ix + 0 + (iy + 0) * Nx];
			v1 = MapToY[ix + 0 + (iy + 0) * Nx];
			u2 = MapToX[ix - 1 + (iy + 0) * Nx];
			v2 = MapToY[ix - 1 + (iy + 0) * Nx];
			u3 = MapToX[ix - 1 + (iy - 1) * Nx];
			v3 = MapToY[ix - 1 + (iy - 1) * Nx];
			u4 = MapToX[ix + 0 + (iy - 1) * Nx];
			v4 = MapToY[ix + 0 + (iy - 1) * Nx];
		}
		else
		{
			// no, guess a bit by using the other coordinates
			x3 = x2;
			y3 = y4;
			
			// TODO: fix mapped points!!
			u1 = MapToX[ix + 0 + (iy + 0) * Nx];
			v1 = MapToY[ix + 0 + (iy + 0) * Nx];
			u2 = u1 - 1.0;
			v2 = v1;
			u3 = u1 - 1.0;
			v3 = v1 - 1.0;
			u4 = u1;
			v4 = v1 - 1.0;
		}
	}
	else if ((q4 >= 0) && (q1 < 0))
	{
		// fourth quadrant
		
		// set the coordinates we already know
		x1 = v1x;
		y1 = v1y;
		x2 = 0;
		y2 = 0;
		x3 = v4x;
		y3 = v4y;
		
		// for the last one we have to check if it falls in the map
		if ((ix < (Nx - 1)) && (iy > 0))
		{
			// yes, copy from map
			x4 = MapFromX[ix + 1 + (iy - 1) * Nx] - x0;
			y4 = MapFromY[ix + 1 + (iy - 1) * Nx] - y0;			
			
			// we can now also safely copy the 4 corresponding coordinates from the target space
			u1 = MapToX[ix + 1 + (iy + 0) * Nx];
			v1 = MapToY[ix + 1 + (iy + 0) * Nx];
			u2 = MapToX[ix + 0 + (iy + 0) * Nx];
			v2 = MapToY[ix + 0 + (iy + 0) * Nx];
			u3 = MapToX[ix + 0 + (iy - 1) * Nx];
			v3 = MapToY[ix + 0 + (iy - 1) * Nx];
			u4 = MapToX[ix + 1 + (iy - 1) * Nx];
			v4 = MapToY[ix + 1 + (iy - 1) * Nx];
		}
		else
		{
			// no, guess a bit by using the other coordinates
			x4 = x1;
			y4 = y3;
			
			// TODO: fix mapped points!!
			u2 = MapToX[ix + 0 + (iy + 0) * Nx];
			v2 = MapToY[ix + 0 + (iy + 0) * Nx];
			u1 = u2 + 1.0;
			v1 = v2;
			u3 = u2;
			v3 = v2 - 1.0;
			u4 = u2 + 1.0;
			v4 = v2 - 1.0;
		}
	}
	
	
	// compute the barycentric coordinates of the two triangles making up the quadrant
	// (the separation diagonal runs between points 1 and 3)
																										   
	// first the barycentric coordinates of the upper triangle consisting of vertices 1, 2, 3
	double detTu = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
	double lu1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detTu;
	double lu2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detTu;
	double lu3 = 1 - lu1 - lu2;
	
	// check if the point (x, y) lies within the upper triangle
	if (   (lu1 >= 0) && (lu1 <= 1.0)
		&& (lu2 >= 0) && (lu2 <= 1.0)
		&& (lu3 >= 0) && (lu3 <= 1.0) )
	{
		// yes, it is in the upper triangle
		
		// compute the mapped coordinates by changing to the barycentric coordinates of the mapped triangle
		*xto = lu1 * u1 + lu2 * u2 + lu3 * u3;
		*yto = lu1 * v1 + lu2 * v2 + lu3 * v3;
		
	}
	else
	{
		// no, we must be in the lower triangle (or outside altogether)
		
		// compute the barycentric coordinates in the lower triangle with vertices 1, 3, 4
		double detTd = (y3 - y4) * (x1 - x4) + (x4 - x3) * (y1 - y4);
		double ld1 = ((y3 - y4) * (x - x4) + (x4 - x3) * (y - y4)) / detTd;
		double ld2 = ((y4 - y1) * (x - x4) + (x1 - x4) * (y - y4)) / detTd;
		double ld3 = 1 - ld1 - ld2;
		
		// compute the mapped coordinates by changing to the barycentric coordinates of the mapped triangle
		*xto = ld1 * u1 + ld2 * u3 + ld3 * u4;
		*yto = ld1 * v1 + ld2 * v3 + ld3 * v4;
	}
}

int CVICALLBACK CamCalLoad_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			char matfilename[1024];
			int SelectionStatus;

			// let the user specify a filename for saving the calibration data
		    SelectionStatus = FileSelectPopup("", "*.mat", "*.mat", "Save Calibration Data", 
									VAL_OK_BUTTON, 0, 1, 1, 1, matfilename);

			// did the user select a valid, existing file?
			if (SelectionStatus == 1)
			{
				// create a Matlab file pointer for loading the calibration data
				MATFile *pmat;

				// open the .mat file
				pmat = matOpen(matfilename, "r");
				
				// if the calibration arrays were already allocated, free them first
				if (CalibrationSLM_X != NULL)
				{
					free(CalibrationSLM_X);
					free(CalibrationSLM_Y);
					free(CalibrationCam_X);
					free(CalibrationCam_Y);
				}
				
				// read the calibration arrays from the matfile
				CalibrationSLM_X = readMatDoubleArray(pmat, "CalibrationSLM_X", &NumCalX, &NumCalY);
				CalibrationSLM_Y = readMatDoubleArray(pmat, "CalibrationSLM_Y", &NumCalX, &NumCalY);
				CalibrationCam_X = readMatDoubleArray(pmat, "CalibrationCam_X", &NumCalX, &NumCalY);
				CalibrationCam_Y = readMatDoubleArray(pmat, "CalibrationCam_Y", &NumCalX, &NumCalY);
				
				// read the camera window origin  
				gCameraWindowOriginX = (int) readMatDoubleScalar(pmat, "CameraWindowOriginX");
				gCameraWindowOriginY = (int) readMatDoubleScalar(pmat, "CameraWindowOriginY");
				
				// also read the zero order spot data
				gZeroOrderULX = readMatDoubleScalar(pmat, "ZeroOrderULX");
				gZeroOrderULY = readMatDoubleScalar(pmat, "ZeroOrderULY");
				gZeroOrderLRX = readMatDoubleScalar(pmat, "ZeroOrderLRX");
				gZeroOrderLRY = readMatDoubleScalar(pmat, "ZeroOrderLRY");
				
				// close the mat file
				matClose(pmat);
			}

			break;
	}
	return 0;
}

int CVICALLBACK CamCalSave_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		char matfilename[1024];
		int SelectionStatus;

		case EVENT_COMMIT:

			// let the user specify a filename for saving the calibration data
		    SelectionStatus = FileSelectPopup("", "*.mat", "*.mat", "Save Calibration Data", 
									VAL_OK_BUTTON, 0, 1, 1, 1, matfilename);

			// did the user select a valid file?
			if (SelectionStatus > 0)
			{
				// create a Matlab file pointer for saving the data
				MATFile *pmat;

				// open the .mat file
				pmat = matOpen(matfilename, "wz");

				// save the calibration data
				writeMatDoubleArray(pmat, "CalibrationSLM_X", CalibrationSLM_X, NumCalX, NumCalY);
				writeMatDoubleArray(pmat, "CalibrationSLM_Y", CalibrationSLM_Y, NumCalX, NumCalY);
				writeMatDoubleArray(pmat, "CalibrationCam_X", CalibrationCam_X, NumCalX, NumCalY);
				writeMatDoubleArray(pmat, "CalibrationCam_Y", CalibrationCam_Y, NumCalX, NumCalY);
				writeMatDoubleScalar(pmat, "ZeroOrderULX", gZeroOrderULX);
				writeMatDoubleScalar(pmat, "ZeroOrderULY", gZeroOrderULY);
				writeMatDoubleScalar(pmat, "ZeroOrderLRX", gZeroOrderLRX);
				writeMatDoubleScalar(pmat, "ZeroOrderLRY", gZeroOrderLRY);
				writeMatDoubleScalar(pmat, "CameraWindowOriginX", (double) gCameraWindowOriginX);
				writeMatDoubleScalar(pmat, "CameraWindowOriginY", (double) gCameraWindowOriginY);
				
				// close the mat-file
				matClose(pmat);					
			}
			
			break;
	}
	return 0;
}

int CVICALLBACK CamCalDisplayPoint_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// clear the SLM
			SLM_Clear();
			
			// store the current translation settings, as these define the current SLM origin
			double htrans = gHorizTrans;
			double vtrans = gVertTrans;
			
			// get the desired point
			int p;
			GetCtrlVal(panel, PANEL_CamCalDisplayPoint, &p);
			
			if (p < NumCalX * NumCalY)
			{
				// let the SLM place a spot at the calibration point
				SLM_setHorizTrans(htrans + CalibrationSLM_X[p]);
				SLM_setVertTrans( vtrans + CalibrationSLM_Y[p]);
			}
			
			// update the SLM 
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
			// reset the original translational settings, but don't update the SLM
			// so the spot stays in place until a new command is issued
			SLM_setHorizTrans(htrans);
			SLM_setVertTrans(vtrans);
			
			
			
			break;
	}
	return 0;
}


int CVICALLBACK CamCalReCal_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// get the desired point
			int p;
			GetCtrlVal(panel, PANEL_CamCalDisplayPoint, &p);
			
			if (p < NumCalX * NumCalY)
			{
				
				
			}

			break;
	}
	return 0;
}






int CVICALLBACK CamCalReCalClick_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// get the desired point
			int p;
			GetCtrlVal(panel, PANEL_CamCalDisplayPoint, &p);
			
			if (p < NumCalX * NumCalY)
			{
				// indicate that we are calibrating a point (this is handled in CameraTimer_Callback)
				gReCalibratingPoint = p;
			}

			break;
	}
	return 0;
}


void updateCalibrationPoint(int p, double SLMx, double SLMy, double CamX, double CamY)
{
	CalibrationSLM_X[p] = SLMx;
	CalibrationSLM_Y[p] = SLMy;
	CalibrationCam_X[p] = CamX;
	CalibrationCam_Y[p] = CamY;
	
	// if we were re-calibrating by mouse, unset the status
	gReCalibratingPoint = -1;
	
}
