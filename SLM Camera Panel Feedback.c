//==============================================================================
//
// Title:       SLM Camera Panel Feedback
// Purpose:     Library of functions to couple camera output to SLM data and 
//				get feedback from the camera stream
//
// Created on:  22-02-2012 at 22:22:12 by Rick.
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

//==============================================================================
// Include files

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


// correspondence of simulation with camera feed: origin coordinates
int gOx = 0, gOy = 0;

// correspondence of simulation with camera feed: magnification factors (camera image size / expected image size)
// NOTE: these values are the ratios of the physical sizes, see gCamMPx and gCamMPy for the pixel ratios
double gCamMx = 0.0, gCamMy = 0.0;

// size of SLM focal unit in camera pixels
double gCamMPx = 1.0, gCamMPy = 1.0;

// indicator variable whether the camera has been calibrated
int gCamCalibrated = 0;



/// HIFN Generates a spot pattern on the SLM, and then tries to find this pattern back in the camera feed
int CVICALLBACK CalibrateCam_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		// the parameters of the spot pattern
		int Nx;
		int Ny;
		int spotxspacing;
		int spotyspacing;
		int sigxoffset;
		int sigyoffset;
		int spotxoffset;
		int spotyoffset;
		int phaseconstraint;
		double phasestep;
		
		case EVENT_COMMIT:
				
			// get the spot pattern parameters from the control panel
			GetCtrlVal(TabPage_1_2, TABPANEL_5_NumXSpots,    &Nx);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_NumYSpots,    &Ny);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXOffset,  &spotxoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYOffset,  &spotyoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SigXOffset,   &sigxoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SigYOffset,   &sigyoffset);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXSpacing, &spotxspacing);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYSpacing, &spotyspacing);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_PhaseConstraint, &phaseconstraint);
			GetCtrlVal(TabPage_1_2, TABPANEL_5_PhaseStep,    &phasestep);
			
			// get the camera feed
			double* image = (double*) malloc((gCamXmax - gCamXmin) * (gCamYmax - gCamYmin) * sizeof(double));
			for (int k = gCamXmin; k < gCamXmax; k++)
			for (int l = gCamYmin; l < gCamYmax; l++)
				image[(l - gCamYmin) * (gCamXmax - gCamXmin) + k - gCamXmin] = (double) gAvgFrame[l * gCamX + k];
			
			// find the spot pattern in the camera feed
			double ulxpos,  ulypos;
			double periodx, periody;
			SLM_findSpotPattern(image, gCamXmax - gCamXmin, gCamYmax - gCamYmin, Nx, Ny, &ulxpos, &ulypos, &periodx, &periody);
			
			// we now have the coordinates of the upper left spot in camera pixel space
			
			// compute the physical magnification factor of the image 
			// (how large is it on CCD / how large was it supposed to be)
			gCamMx = (gCamPixelSize * 1e6 * periodx) / ((spotxspacing + 1.0) * SLM_getFocalUnitX());
			gCamMy = (gCamPixelSize * 1e6 * periody) / ((spotyspacing + 1.0) * SLM_getFocalUnitY());
			
			// compute the pixel magnification factor
			// (how many pixels in the camera plane correspond to one pixel in the focal plane?)
			gCamMPx = periodx / (spotxspacing + 1.0);
			gCamMPy = periody / (spotyspacing + 1.0);
			
			// find the origin of the SLM focal plane in the camera pixel space
			gOx = ulxpos - (sigxoffset + spotxoffset) * gCamMPx;
			gOy = ulypos - (sigyoffset + spotyoffset) * gCamMPy;
			
			// indicate that the camera has been calibrated
			gCamCalibrated = 1;
			
			// enable the update button on the feedback panel 
			// (the updating can only be done from a calibrated camera, it needs to know the window coordinates)
			SetCtrlAttribute(TabPage_3, TABPANE_10_UpdateFeedbackWindow, ATTR_DIMMED, 0);

			break;
	}
	return 0;
}


/// HIFN Returns a window area in the camera plane, specified in focal plane coordinates (!)
double* getSignalWindowFromCamera(int ulx, int uly, int Wx, int Wy, int corrx, int corry)
{
	return getSignalWindowFromCameraResampled(NULL, ulx, uly, Wx, Wy, corrx, corry, Wx, Wy);
}



double* getSignalWindowFromCameraResampled(double* window, int ulx, int uly, int Wx, int Wy, int corrx, int corry, int resampleWx, int resampleWy)
{
	if (gCamCalibrated)
	{
		// calculate the coordinates in the camera plane, corresponding 
		// to the ulx, lrx, ... coordinates in the focal plane
		int ulxw = TransformFocalPlaneToBitmapX(ulx);
		int ulyw = TransformFocalPlaneToBitmapY(uly);
		int lrxw = TransformFocalPlaneToBitmapX(ulx + Wx);
		int lryw = TransformFocalPlaneToBitmapY(uly + Wy);
		int cxw = lrxw - ulxw + corrx;
		int cyw = lryw - ulyw + corry;
		
		// allocate memory for the portion of the camera plane (camera plane coordinates!)
		double* camwindow = (double*) calloc(cxw * cyw, sizeof(double));
	
		// copy the specified values from the camera frame to the window array
		for (int k = 0; k < cxw; k++)
		for (int l = 0; l < cyw; l++)
		{
			// check if we are within the camera bitmap, if so, copy the value
			if ((k + ulxw > 0) && (k + ulxw < gCamX) && (l + ulyw > 0) && (l + ulyw < gCamY))
				camwindow[k + l * cxw] = gAvgFrame[(k + ulxw) + (l + ulyw) * gCamX];	
		}
		
		// resize the camera window to have the requested dimensions
		if (window == NULL)
		{
			// no valid pointer was supplied, so we create a new one while resampling			
			window = SLM_resampleBitmap(camwindow, cxw, cyw, resampleWx, resampleWy);
		}
		else
		{
			// a valid pointer was supplied, so we resample into that memory
			SLM_resampleBitmapInPlace(camwindow, window, cxw, cyw, resampleWx, resampleWy);
		}
	
		// clean up
		free(camwindow);

		// return the pointer to the window array
		return window;
	
	}
	else
		return NULL;
}


/// HIFN get the camera space coordinate of a focal plane coordinate (X)
int TransformFocalPlaneToBitmapX(int fx)
{
	return (int) (gOx + gCamMPx * ((double) fx));
}


/// HIFN get the camera space coordinate of a focal plane coordinate (Y)
int TransformFocalPlaneToBitmapY(int fy)
{
	return (int) (gOy + gCamMPy * ((double) fy));
}


/// HIFN draws lines marking the spot pattern
void drawSpotLines(int Panel, int Canvas)
{
	// if the current SLM pattern is a spot pattern, we plot a grid
	// indicating the spot positions
	if (SLM_getCurrentPattern() == SLM_BEAMSPLIT_IFTA)
	{
		// get the canvas dimensions
		int CanvasX, CanvasY;
		GetCtrlAttribute(Panel, Canvas, ATTR_WIDTH,  &CanvasX);
		GetCtrlAttribute(Panel, Canvas, ATTR_HEIGHT, &CanvasY);
		
		// retrieve the current settings from the control panel
		int Nx, Ny, sigxoffset, sigyoffset, spotxoffset, spotyoffset, spotxspacing, spotyspacing;
		GetCtrlVal(TabPage_1_2, TABPANEL_5_NumXSpots,    &Nx);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_NumYSpots,    &Ny);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXOffset,  &spotxoffset);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYOffset,  &spotyoffset);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_SigXOffset,   &sigxoffset);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_SigYOffset,   &sigyoffset);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotXSpacing, &spotxspacing);
		GetCtrlVal(TabPage_1_2, TABPANEL_5_SpotYSpacing, &spotyspacing);
		
		// convenience variables for the lattice periods
		double periodx = spotxspacing + 1.0;
		double periody = spotyspacing + 1.0;
		
		// compute the upper left spot coordinates in camera space
		int ulxpos = gOx + gCamMPx * (spotxoffset + sigxoffset + 1);
		int ulypos = gOy + gCamMPy * (spotyoffset + sigyoffset + 1);
		
		// compute the lower right spot coordinates
		int lrxpos = gOx + gCamMPx * ((spotxoffset + sigxoffset + 1) + (Nx - 1) * (periodx));
		int lrypos = gOy + gCamMPy * ((spotyoffset + sigyoffset + 1) + (Ny - 1) * (periody));
		
		// set the pen color
		SetCtrlAttribute(Panel, Canvas, ATTR_PEN_COLOR, 0x00AAAAAA);
		
		// draw a cross at the origin
		int xo = TransformBitmapToCamCanvasX(gOx);
		int yo = TransformBitmapToCamCanvasY(gOy);
		CanvasDrawLine(Panel, Canvas, MakePoint(xo, 0), MakePoint(xo, CanvasY));
		CanvasDrawLine(Panel, Canvas, MakePoint(0, yo), MakePoint(CanvasX, yo));
		
		// draw vertical lines indicating the grid
		for (int k = 0; k <= Nx; k++)
		{
			// get the coordinates of the endpoints of the gridline
			int xc = TransformBitmapToCamCanvasX((int) (ulxpos + gCamMx * SLM_getFocalUnitX() * (-periodx / 2.0 + k * periodx) / (gCamPixelSize * 1e6)));
			int y1 = TransformBitmapToCamCanvasY((int) (ulypos - gCamMy * SLM_getFocalUnitY() * periody * 0.5 / (gCamPixelSize * 1e6)));
			int y2 = TransformBitmapToCamCanvasY((int) (lrypos + gCamMy * SLM_getFocalUnitY() * periody * 0.5 / (gCamPixelSize * 1e6)));

			// draw the line
			CanvasDrawLine(Panel, Canvas, MakePoint(xc, y1), MakePoint(xc, y2));
		}

		// draw horizontal lines indicating the grid
		for (int l = 0; l <= Ny; l++)
		{
			// get the coordinates of the endpoints of the gridline
			int yc = TransformBitmapToCamCanvasY((int) (ulypos + gCamMy * SLM_getFocalUnitY() * (-periody / 2.0 + l * periody) / (gCamPixelSize * 1e6)));
			int x1 = TransformBitmapToCamCanvasX((int) (ulxpos - gCamMx * SLM_getFocalUnitX() * periodx * 0.5 / (gCamPixelSize * 1e6)));
			int x2 = TransformBitmapToCamCanvasX((int) (lrxpos + gCamMx * SLM_getFocalUnitX() * periodx * 0.5 / (gCamPixelSize * 1e6)));

			// draw the line
			CanvasDrawLine(Panel, Canvas, MakePoint(x1, yc), MakePoint(x2, yc));
		}
	}	
}



/// HIFN gets a window from the camera data
void getCameraWindow(int x0, int y0, int width, int height, double* dest)
{
	double sum = 0.0;
	for (int k = 0; k < width; k++)
	for (int l = 0; l < height; l++)
	{
		double tmp = (double) gAvgFrame[(k + x0) + (l + y0) * gCamX];
		dest[k + l * width] = tmp;
		sum += tmp;
	}
}





