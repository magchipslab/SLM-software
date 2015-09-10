//==============================================================================
//
// Title:       SLM Controller Shack-Hartmann Aberration Correction
// Purpose:     Control Panel routines for correcting aberrations using Shack-Hartmann
//
// Created on:  09-12-2012 by Rick.
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
#include "SLM Control Panel Internal.h"
#include "toolbox.h"
#include "SLM.h"
#include "SLM_internal.h"
#include "SLM Camera Panel.h"
#include "SLM_Zernike.h"




// the physical dimensions of the SLM (in m)
const double HalfLxSLM = 0.5*1.5 * 0.01;
const double HalfLySLM = 0.5*0.8 * 0.01;




int CVICALLBACK ShackHartmann_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// read out the input parameters in um
			double spotx, spoty, spotdiameter, strayx, strayy;
			GetCtrlVal(panel, TABPANE_13_SHSpotSize, &spotdiameter);
			GetCtrlVal(panel, TABPANE_13_SHSpotX,    &spotx);
			GetCtrlVal(panel, TABPANE_13_SHSpotY,    &spoty);
			GetCtrlVal(panel, TABPANE_13_SHStrayX,   &strayx);
			GetCtrlVal(panel, TABPANE_13_SHStrayY,   &strayy);
			
			// convert to meter
			spotx *= 1.0e-6;
			spoty *= 1.0e-6;
			spotdiameter *= 1.0e-6;
			strayx *= 1.0e-6;
			strayy *= 1.0e-6;
			
			// set the Shack-Hartmann pattern
			SLM_SetShackHartmannPattern(spotx, spoty, spotdiameter, strayx, strayy);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}

			break;
	}
	return 0;
}


// Calibrate the Phase aberration using the spot projector 
int CVICALLBACK CalibratePhaseAberr_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// get the number of spots, their spacing, spot diameter, strayx and strayy from the SH control panel
			unsigned int Nx, Ny, Ntot;
			double xspacing, yspacing, randomampl, spotdiameter, strayx, strayy, spotx, spoty;
			GetCtrlVal(panel, TABPANE_13_SHSamplePointsX, &Nx);
			GetCtrlVal(panel, TABPANE_13_SHSamplePointsY, &Ny);
			GetCtrlVal(panel, TABPANE_13_SHSpotSize, &spotdiameter);
			GetCtrlVal(panel, TABPANE_13_SHStrayX,   &strayx);
			GetCtrlVal(panel, TABPANE_13_SHStrayY,   &strayy);
			Ntot = Nx * Ny;
			// Convert to meter and set spacing and spots
			spotdiameter *= 1.0e-6;
			strayx = 1.0e-6;
			strayy = 1.0e-6;
			xspacing = 2 * HalfLxSLM / Nx;
			yspacing = 2 * HalfLySLM / Ny;
			spotx = -HalfLxSLM + (0.5 * xspacing);
			spoty =  HalfLySLM - (0.5*yspacing);
			
			//FIRST retrieve (0,0) spot on the camera frame
			//Clear the SLM
			SLM_Clear();
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
			// find the maximum intensity inside the visible part of the camera canvas => zeroth order spot
			double imax = 0.0;
			int xmax, ymax;
			for (int n = gCamXmin; n < gCamXmax; n++)
			for (int m = gCamYmin; m < gCamYmax; m++)
			{
				if (gAvgFrame[n + gCamX * m] > imax)
				{
					// store the new maximum and its location
					imax = gAvgFrame[n + gCamX * m];
					xmax = n;
					ymax = m;
				}
			}
				
			/*
			// loop over the number of spots of the SH pattern
			for (int nx = 1; nx <= Nx; nx++)
			{
				for (int ny = 1; ny <= Ny; ny++)
				{
					//Produce the Shack Hartmann pattern and update the SLM
					SLM_SetShackHartmannPattern(spotx, spoty, spotdiameter, strayx, strayy);
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);

					//Wait for one second to let the pattern set on the camera
					//delay(1);
					//Retrieve signal from camera
					//Find spot position on the camera and save it in array
					//Move to the next spot position
					spoty -= yspacing;
					
				}
				// Move to the next spot position
				spotx += xspacing;
				spoty = HalfLySLM - 0.5 * yspacing;
				
			}
			*/	
					
				
			
			
			// get the signal area from the camera
			
			
			// get the number of Zernike polynomials we need for fitting
			
			
			break;
		}
	}
	return 0;
}

