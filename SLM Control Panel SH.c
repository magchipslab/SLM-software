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
#include <stdio.h>
#include <stdlib.h>
#include "stats.h"
#include <time.h>   //to see current date and time, probably not needed
#include <utility.h>
#include <string.h>


// panel, canvas, graphs and timer controls for displaying camera output
static int gPanel, gCanvas, gHistoX, gHistoY, gTimer;
// bitmap id
static int gBitmap = -1;
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
			double spotx, spoty, spotdiameter, strayx, strayy, xstep, ystep, ytop, xright;
			GetCtrlVal(panel, TABPANE_13_SHSpotSize, &spotdiameter);
			GetCtrlVal(panel, TABPANE_13_SHSpotX,    &spotx);
			GetCtrlVal(panel, TABPANE_13_SHSpotY,    &spoty);
			GetCtrlVal(panel, TABPANE_13_SHStrayX,   &strayx);
			GetCtrlVal(panel, TABPANE_13_SHStrayY,   &strayy);
			GetCtrlVal(panel, TABPANE_13_SHxStep,    &xstep); 
			GetCtrlVal(panel, TABPANE_13_SHyStep,    &ystep);
			GetCtrlVal(panel, TABPANE_13_SHyTop,     &ytop);
			GetCtrlVal(panel, TABPANE_13_SHxRight,   &xright);
			
			// convert to meter
			spotx *= 1.0e-6;
			spoty *= 1.0e-6;
			spotdiameter *= 1.0e-6;
			strayx *= 1.0e-6;
			strayy *= 1.0e-6;
			xstep*=1.0e-6;
			ystep*=1.0e-6;
			xright*=1.0e-6;
			ytop*=1.0e-6;
			
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
			double xspacing, yspacing, randomampl, spotdiameter, strayx, strayy, spotx, spoty, xstep, ystep, ytop, xright;
			//GetCtrlVal(panel, TABPANE_13_SHSamplePointsX, &Nx);
			//GetCtrlVal(panel, TABPANE_13_SHSamplePointsY, &Ny);
			GetCtrlVal(panel, TABPANE_13_SHSpotSize, &spotdiameter);
			GetCtrlVal(panel, TABPANE_13_SHStrayX, &strayx);
			GetCtrlVal(panel, TABPANE_13_SHStrayY, &strayy);
			GetCtrlVal(panel, TABPANE_13_SHxStep, &xstep); 
			GetCtrlVal(panel, TABPANE_13_SHyStep, &ystep);
			GetCtrlVal(panel, TABPANE_13_SHyTop,  &ytop);
			GetCtrlVal(panel, TABPANE_13_SHxRight, &xright);
			GetCtrlVal(panel, TABPANE_13_SHSpotX,  &spotx);
			GetCtrlVal(panel, TABPANE_13_SHSpotY,  &spoty);
			// Convert to meter and set spacing and spots
			spotx *=1.0e-6;
			spoty *=1.0e-6;
			spotdiameter *=1.0e-6;
			strayx *= 1.0e-6;
			strayy *= 1.0e-6;
			xstep *=1.0e-6;
			ystep *=1.0e-6;
			xright *=1.0e-6;
			ytop *=1.0e-6;
			
			Nx=abs((xright-spotx+1.0e-6)/xstep);		 // 1 added to avoid problem with the number of steps
			Ny=abs((spoty-ytop)/ystep);
			Ntot = Nx * Ny;
			
			//FIRST retrieve (0,0) spot on the camera frame
			//Clear the SLM
			SLM_Clear();
			SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			 /*
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
			*/	
			 //Getting current date and creating a folder  
			 int		   month, day, year, hour, minute,second;
		     GetSystemDate (&month, &day, &year);	 //getting the current date 
			 GetSystemTime (&hour, &minute, &second);
 				char bufferD[33];
				char bufferM[33];
				char bufferY[33];
				char bufferH[33];
				char bufferMi[33];
				Fmt (bufferD,"%s<%i", day); 
				Fmt (bufferM,"%s<%i", month);
				Fmt (bufferY,"%s<%i", year);
				Fmt (bufferH,"%s<%i", hour);
				Fmt (bufferMi,"%s<%i", minute);
				char NewFolderName[255]="C:\\Users\\SLM Heroes\\Desktop\\Data\\";   //constructor of the folder name, current date
				char FolderTime[255]="";
				strcat(NewFolderName,bufferD);   
				strcat(NewFolderName,".");  
				strcat(NewFolderName,bufferM);
				strcat(NewFolderName,".");
				strcat(NewFolderName,bufferY);
				MakeDir(NewFolderName);//makes a new directory in Desktop\\Data
				//creating the second folder
				strcat(FolderTime,NewFolderName);
				strcat(FolderTime,"\\"); 
				strcat(FolderTime,bufferH);
				strcat(FolderTime,".");
				strcat(FolderTime,bufferMi);
				MakeDir(FolderTime);//makes a new directory in Desktop\\Data\\Current Date  
				
		
			
			// loop over the number of spots of the SH pattern, this loop was changed by David 
					  
			for (int nx = 0; nx <= Nx; nx++)
			{
				for (int ny = 0; ny <= Ny; ny++)
				{
					//Produce the Shack Hartmann pattern and update the SLM
					SLM_SetShackHartmannPattern(spotx, spoty, spotdiameter, strayx, strayy);
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
					
					 //Grab the picture from the camera
					grabCameraFrame();
				
					//Wait for one second to let the pattern set on the camera
					Delay(1);// this should be a delay of 1 sec
							
					//this is needed for correct image saving, otherwise it will save the same image each time
					performCameraAveraging();

						
							 //Create File Name
							 char bufferX [33];
							 char bufferY [33];
							 int Xposition=RoundRealToNearestInteger(spotx*1.0e6);
							 int Yposition=RoundRealToNearestInteger(spoty*1.0e6);
							 Fmt (bufferX,"%s<%i", Xposition);
							 Fmt (bufferY,"%s<%i", Yposition);
							 char filename[129]="";
							 strcat(filename,FolderTime);
							 strcat(filename,"\\");
							 strcat(filename,bufferX);
							 strcat(filename,".");
							 strcat(filename, bufferY);
							 strcat(filename,".");
							 strcat(filename, "mat");
							  
							//save the data
							/*
							// generate a bitmap from the canvas
							
							int TmpBitmap = -1;
							GetCtrlDisplayBitmap (panel, gCanvas, 0, &TmpBitmap);
				
							// save it to file
							//SaveBitmapToPNGFile (TmpBitmap, filename);
								
							// discard the bitmap
							DiscardBitmap(TmpBitmap);
							
							// construct a filename from the .png screenshot filename,
							// where we assume that the file is always ending in .png
							char matfilename[1024];
							strcpy(matfilename, pngfilename);
							strcpy((matfilename + StringLength(pngfilename) - 3), "mat");
							 */
							 
							// save the raw data as a Matlab .mat file
							CameraSaveMatFile(filename, 1);
									 	 
					//Retrieve signal from camera
					//Find spot position on the camera and save it in array
					//Move to the next spot position
					spoty -= ystep;
					
				}
				// Move to the next spot position
				spotx += xstep;
				spoty = spoty+(Ny+1)*ystep;
				
			}
				
					
				
			
			
			// get the signal area from the camera
			
			
			// get the number of Zernike polynomials we need for fitting
			
			
			break;
		}
	}
	return 0;
}

