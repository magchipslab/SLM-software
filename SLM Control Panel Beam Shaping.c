//==============================================================================
//
// Title:       SLM Controller
// Purpose:     Callbacks associated with (geometric) beam shaping and 
// 			    speckle-free IFTAs
//
// Created on:  04-01-2012 by Rick.
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

// this file is only included for the character beam shaping modus (temporarily)
#include "SLM_internal.h"

//#include "Camera Controller.h"


/// HIFN Callback function for geometric beam shaping
int CVICALLBACK BeamShape_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		double signalwidth, signalheight, sigmax, sigmay;
		int btgaussian, btsquare;
		
		// number of signal intensity points
		int Ns;
			
		// create arrays for the input and signal intensities
		double* isx;
		double* isy;
			
		// total intensity in the signal
		double sumi;
		
		case EVENT_COMMIT:
		case EVENT_VAL_CHANGED:
			
			// read out beam parameters
			GetCtrlVal(panel, TABPANEL_7_SignalWidth, &signalwidth);
			GetCtrlVal(panel, TABPANEL_7_SignalHeight, &signalheight);
			GetCtrlVal(panel, TABPANEL_7_SigmaX, &sigmax);
			GetCtrlVal(panel, TABPANEL_7_SigmaY, &sigmay);
										 
			// get the type of beam desired
			GetCtrlVal(panel, TABPANEL_7_BeamTypeGaussian, &btgaussian);
			GetCtrlVal(panel, TABPANEL_7_BeamTypeSquare, &btsquare);
			
			// convert the input values from microns to meters
			signalwidth  = signalwidth  * 1.0e-6;
			signalheight = signalheight * 1.0e-6;
			sigmax = sigmax * 1.0e-6;
			sigmay = sigmay * 1.0e-6;
			
			// number of signal intensity points
			Ns = 2048;
			
			// create arrays for the input and signal intensities
			isx = (double*) malloc(Ns * sizeof(double));
			isy = (double*) malloc(Ns * sizeof(double));
			
			// total intensity in the signal
			sumi = 0.0;
			
			// fill the signal arrays for the x direction
			for (int k = 0; k < Ns; k++)
			{
				// calculate the relative coordinate
				double x = ((double) k) / ((double) Ns - 1.0) - 0.5;
				
				// compute the signal intensity
				double isignal;
				if (btgaussian == 1)
				{
					// take the signal intensity to be a Gaussian
					isignal = exp(-(signalwidth * signalwidth * x * x) / (sigmax * sigmax));
				}
				else
				{
					// take the signal intensity to be square
					double sigxsize = 0.5 * sigmax / signalwidth;
					isignal = (fabs(x) - sigxsize < 0.0 ? 1.0 : 0.0); 
				}
				
				// set the signal intensities
				isx[k] = isignal;
			}
			
			// fill the signal arrays for the y direction
			for (int l = 0; l < Ns; l++)
			{
				// calculate the relative coordinate
				double y = ((double) l) / ((double) Ns - 1.0) - 0.5;
				
				// compute the signal intensity
				double isignal;
				if (btgaussian == 1)
				{
					// take the signal intensity to be a Gaussian
					isignal = exp(-(signalheight * signalheight * y * y) / (sigmay * sigmay));
				}
				else
				{
					// take the signal intensity to be square
					double sigysize = 0.5 * sigmay / signalheight;
					isignal = (fabs(y) - sigysize < 0.0 ? 1.0 : 0.0); 
				}
				
				// set the signal intensities
				isy[l] = isignal;
			}
			
			//PlotY(panel, TABPANEL_7_GRAPH, isx, Ns, VAL_DOUBLE, VAL_THIN_LINE, VAL_NO_POINT, VAL_SOLID, 1, 0x00FF0000);
			
			// calculate the geometric beam shaping phase profile
			SLM_generateBeamShape(isx, isy, Ns, Ns, 0.0, signalwidth, 0.0, signalheight, 1);
			
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


/// HIFN Callback function for selecting the beam type in beam shaping
int CVICALLBACK BeamType_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// check which radio button was toggled
			switch (control)
			{
				case TABPANEL_7_BeamTypeSquare:
					
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_7_BeamTypeSquare, ATTR_CTRL_VAL, 1);
					SetCtrlAttribute(panel, TABPANEL_7_BeamTypeGaussian, ATTR_CTRL_VAL, 0);
					
					break;
					
				case TABPANEL_7_BeamTypeGaussian:
					
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_7_BeamTypeSquare, ATTR_CTRL_VAL, 0);
					SetCtrlAttribute(panel, TABPANEL_7_BeamTypeGaussian, ATTR_CTRL_VAL, 1);
					
					break;
			}
			
			break;
	}
	return 0;
}



/// HIFN Callback function for the beam shaping mode radio buttons
int CVICALLBACK rbShapingMode_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// check which radio button was toggled
			switch (control)
			{
				case TABPANEL_6_rbShapingModeArb:
					
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_6_rbShapingModeArb, ATTR_CTRL_VAL, 1);
					SetCtrlAttribute(panel, TABPANEL_6_rbShapingModeStd, ATTR_CTRL_VAL, 0);
					
					// bring the corresponding tabpage to the front
					SetCtrlAttribute(panel, TABPANEL_6_TAB, ATTR_CTRL_INDEX, 1);
					
					break;
					
				case TABPANEL_6_rbShapingModeStd:
					
					// make sure this button is switched on, the others off
					SetCtrlAttribute(panel, TABPANEL_6_rbShapingModeArb, ATTR_CTRL_VAL, 0);
					SetCtrlAttribute(panel, TABPANEL_6_rbShapingModeStd, ATTR_CTRL_VAL, 1);
					
					// bring the corresponding tabpage to the front
					SetCtrlAttribute(panel, TABPANEL_6_TAB, ATTR_CTRL_INDEX, 0);
					
					break;
			}
			
			break;
	}
	return 0;
}

void GeneratePhaseFromFile(char* FilePath, int panel, int refinephase)
{
	// open the file
	int FileHandle = OpenFile(FilePath, VAL_READ_ONLY, VAL_OPEN_AS_IS, VAL_BINARY);
	
	if (FileHandle != -1)
	{
		// display the picture in the picture control
		DisplayImageFile(panel, TABPANEL_8_Picture, FilePath);
	
		// get the bitmap from the picture control
		int BitmapID = -1;
		GetBitmapFromFile(FilePath, &BitmapID);
	
		// allocate variables to hold the bitmap data
		int BytesPerRow, PixelDepth, Width, Height;
		unsigned char* BitmapData;
		AllocBitmapData (BitmapID, NULL, &BitmapData, NULL);
	
		// get the bitmap pixel data
		GetBitmapData (BitmapID, &BytesPerRow, &PixelDepth, &Width, &Height, NULL, BitmapData, NULL);

		// allocate array to hold the signal intensity
		double* Signal = (double*) malloc(Width * Height * sizeof(double));
	
		// extract the GREEN channel pixel data, and convert it to a double between 0 and 1
		for (int k = 0; k < Width * Height; k++)
			Signal[k] = ((double) BitmapData[k * (PixelDepth >> 3) + 1]) / 256.0;
	
		// get the number of iterations
		int numPF, numAF;
		GetCtrlVal(panel, TABPANEL_8_IterationsPF, &numPF);
		GetCtrlVal(panel, TABPANEL_8_IterationsAF, &numAF);
	
		// get the signal amplitude
		double sigamp;
		GetCtrlVal(panel, TABPANEL_8_SigAmp, &sigamp);
	
		// get the pattern width and height
		double sxmin = 0.0;
		double symin = 0.0;
		double sxmax, symax;
		GetCtrlVal(panel, TABPANEL_8_PictureWidth, &sxmax);
		GetCtrlVal(panel, TABPANEL_8_PictureHeight, &symax);
	
		// convert units from micron to meter
		sxmax = sxmax * 1e-6;
		symax = symax * 1e-6;
		
		// compute the pixel coordinates in the focal plane of the signal window
		int wx = sxmin / SLM_getFocalUnitX();
		int wy = symin / SLM_getFocalUnitY();
		int wxsize = sxmax / SLM_getFocalUnitX() - wx;
		int wysize = symax / SLM_getFocalUnitY() - wy;
		
		// should we use soft operators?
		int softop;
		GetCtrlVal(panel, TABPANEL_8_SoftOp, &softop);
		
		// should we use the MRAF algorithm
		int MRAF;
		GetCtrlVal(panel, TABPANEL_8_MRAF, &MRAF);
		
		// generate the phase pattern for this image
		SLM_generatePhase(Signal, Width, Height, wx, wy, wxsize, wysize, numPF, numAF, sigamp, refinephase, softop, MRAF, 1);
	
		// close the file again
		CloseFile(FileHandle);
		
		// clean up
		free(BitmapData);
		free(Signal);
	}
}


/// HIFN Callback function for the "load file" button for beam shaping
int CVICALLBACK LoadFile_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		// buffer for the filename
		char FilePath[600];
		int SelectionStatus;
		
		case EVENT_COMMIT:
			
			// open a file dialog
			SelectionStatus = FileSelectPopup("", "*.jpg", "*.tif; *.pcx; *.bmp; *.dib; *.rle; *.ico; *.jpg; *.png; *.wmf; *.emf",
				"Open an image file", VAL_OK_BUTTON, 0, 0, 1, 0, FilePath);
		        
			if (SelectionStatus > 0)
			{
				// display the name of the selected file in the textbox
				SetCtrlVal(panel, TABPANEL_8_BeamShapeFilename, FilePath);
				
				// see if we want to keep the current phase pattern as an initial guess
				int refinephase;
				GetCtrlVal(panel, TABPANEL_8_BeamShapeRefinePhase, &refinephase);
				
				// start phase pattern generation
				GeneratePhaseFromFile(FilePath, panel, refinephase);	
				
				// check if we need to update
				if (eventData1 != SLM_NO_UPDATE)
				{
					// update the SLM panel and the simulation Panel (if toggled)
					SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
				}
			}
			
			break;
	}
	return 0;
}


/// HIFN Starts arbitrary pattern generation from the picture file specified in the filename textbox
int CVICALLBACK PicPattern_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		{
			// get the current filename from the text box
			char FilePath[600];
			GetCtrlVal(panel, TABPANEL_8_BeamShapeFilename, FilePath);
			
			// see if we want to keep the current phase pattern as an initial guess
			int refinephase;
			GetCtrlVal(panel, TABPANEL_8_BeamShapeRefinePhase, &refinephase);
			
			// start phase pattern generation
			GeneratePhaseFromFile(FilePath, panel, refinephase);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			
			break;
		}
	}
	return 0;
}


int CVICALLBACK CharacterNumber_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			char* characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.,;:?!@#$%^&*()+=/";
			
			int charnumber;
			GetCtrlVal(panel, TABPANEL_8_CharacterNumber, &charnumber);
			
			int bold;
			GetCtrlVal(panel, TABPANEL_8_CharacterBold, &bold);
			
			int fontsize;
			GetCtrlVal(panel, TABPANEL_8_CharacterFontsize, &fontsize);
			
			char fontname[400];
			GetCtrlVal(panel, TABPANEL_8_CharacterFont, fontname);
			
			char printstr[2];
			printstr[0] = characters[charnumber - 1];
			printstr[1] = '\0';
			
			double posx, posy;
			GetCtrlVal(panel, TABPANEL_8_CharacterPosX, &posx);
			GetCtrlVal(panel, TABPANEL_8_CharacterPosY, &posy);
			
			// first we clear the canvas
			SetCtrlAttribute(panel, TABPANEL_8_CharacterCanvas, ATTR_PEN_COLOR, VAL_BLACK);
			CanvasDrawRect(panel, TABPANEL_8_CharacterCanvas, MakeRect(0, 0, 190, 190), VAL_DRAW_FRAME_AND_INTERIOR);
			
			// set the pen properties for drawing characters			
			SetCtrlAttribute(panel, TABPANEL_8_CharacterCanvas, ATTR_PEN_COLOR, VAL_WHITE);
			SetCtrlAttribute(panel, TABPANEL_8_CharacterCanvas, ATTR_PEN_FILL_COLOR, VAL_BLACK);
			
			CreateMetaFont("RicksFont", fontname, fontsize, bold, 0, 0, 0);	
			CanvasDrawText(panel, TABPANEL_8_CharacterCanvas, printstr, "RicksFont", MakeRect(posy, posx, 2 * 190, 2 * 190), VAL_UPPER_LEFT);
			CanvasUpdate(panel, TABPANEL_8_CharacterCanvas, VAL_ENTIRE_OBJECT);
			
			gCurrentPattern = SLM_BEAMSHAPE_CHAR;

			break;
	}
	return 0;
}

int CVICALLBACK CharacterGo_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			
			// write a character on the canvas
			CharacterNumber_Callback (panel, TABPANEL_8_CharacterNumber, EVENT_COMMIT, NULL, 0, 0);
			
			// see if we want to keep the current phase pattern as an initial guess
			int refinephase;
			GetCtrlVal(panel, TABPANEL_8_BeamShapeRefinePhase, &refinephase);
			
			// get the canvas size
			int Width, Height;
			GetCtrlAttribute(panel, TABPANEL_8_CharacterCanvas, ATTR_WIDTH, &Width);
			GetCtrlAttribute(panel, TABPANEL_8_CharacterCanvas, ATTR_HEIGHT, &Height);
			
			// get the pixel values from the canvas
			double* Signal = (double*) malloc(Width * Height * sizeof(double));
			int PixelColor;
			for (int k = 0; k < Width; k++)
			for (int l = 0; l < Height; l++)
			{
				CanvasGetPixel(panel, TABPANEL_8_CharacterCanvas, MakePoint(k, l), &PixelColor);
				Signal[l * Width + k] = ((double)(PixelColor & VAL_GREEN)) / ((double) VAL_GREEN);
				//if (PixelColor > 0)
				//	printf("(%i, %i) [%.2f %.2f %.2f]\n", k, l, ((double)(PixelColor & VAL_RED)) / ((double) VAL_RED), ((double)(PixelColor & VAL_GREEN)) / ((double) VAL_GREEN), ((double)(PixelColor & VAL_BLUE)) / ((double) VAL_BLUE));
				
			}
	
			// get the number of iterations
			int numPF, numAF;
			GetCtrlVal(panel, TABPANEL_8_IterationsPF, &numPF);
			GetCtrlVal(panel, TABPANEL_8_IterationsAF, &numAF);
	
			// get the signal amplitude
			double sigamp;
			GetCtrlVal(panel, TABPANEL_8_SigAmp, &sigamp);
	
			// get the pattern width and height
			double sxmin = 0.0;
			double symin = 0.0;
			double sxmax, symax;
			GetCtrlVal(panel, TABPANEL_8_PictureWidth, &sxmax);
			GetCtrlVal(panel, TABPANEL_8_PictureHeight, &symax);
	
			// convert units from micron to meter
			sxmax = sxmax * 1e-6;
			symax = symax * 1e-6;
		
			// compute the pixel coordinates in the focal plane of the signal window
			int wx = sxmin / SLM_getFocalUnitX();
			int wy = symin / SLM_getFocalUnitY();
			int wxsize = sxmax / SLM_getFocalUnitX() - wx;
			int wysize = symax / SLM_getFocalUnitY() - wy;
		
			// should we use soft operators?
			int softop;
			GetCtrlVal(panel, TABPANEL_8_SoftOp, &softop);
		
			// should we use the MRAF algorithm
			int MRAF;
			GetCtrlVal(panel, TABPANEL_8_MRAF, &MRAF);
		
			// generate the phase pattern for this image
			SLM_generatePhase(Signal, Width, Height, wx, wy, wxsize, wysize, numPF, numAF, sigamp, refinephase, softop, MRAF, 1);
			
			// check if we need to update
			if (eventData1 != SLM_NO_UPDATE)
			{
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			}
			
			free(Signal);
			gCurrentPattern = SLM_BEAMSHAPE_CHAR;
	
			break;
	}
	return 0;
}
