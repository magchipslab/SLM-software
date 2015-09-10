//==============================================================================
//
// Title:       SLM Camera Panel
// Purpose:     Library of functions to display camera output on a panel
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


// FlyCapture gError variable to catch error messages
static fc2Error gError;

// bitmap id
static int gBitmap = -1;

// camera colormap
static int* gColorMap;

// camera context and guid (Point Grey data structures)
static fc2Context gContext;
static fc2PGRGuid gGuid;

// the camera image data structure (Point Grey)
static fc2Image gImage;

// indicator variable for indicating whether there is a camera present on this computer
static int gCamPresent = 0;

// panel, canvas, graphs and timer controls for displaying camera output
static int gPanel, gCanvas, gHistoX, gHistoY, gTimer;

// number of pixels of the canvas
static int gCanvasX, gCanvasY;

// number of pixels of canvas and camera feed
int gCamX, gCamY;

// the pixel position of the bitmap (camera feed) that should be centered in the canvas
int gCamCenterX = 0, gCamCenterY = 0;

// the pixel coordinates of the corners of the visible portion of the bitmap
// (convenience variables, as these can in principle be computed from center and zoom variables)
int gCamXmin, gCamXmax, gCamYmin, gCamYmax;

// zoom factor
double gCamZoomFactor = 1.0;

// convenience arrays with indices of pixels
static int* gPixX;
static int* gPixY;

// dummy data array for the case when there is no camera present
static unsigned char* gSignalChannel;

// histogram data arrays
static unsigned char* gHistoXdata;
static unsigned char* gHistoYdata;

// grid variables
static int gShowGrid = 1;
static int gGridSpacing = 100;
static int gShowCrosshairs = 1;

// frame buffer
static unsigned char* gFrameBuffer[MAXFRAMES];

// frame counter within gFrameBuffer array
static int gFrameCounter = 0;

// total frame number
int gFrameNumber = 0;

// the number of frames to average over
int gNumFrames = 1;

// the time-averaged camera frame
unsigned char* gAvgFrame;

// temporary frame (of larger data type) to hold the sum of all the frames in the buffer
static unsigned int* gTmpFrame;

// indicate whether to continue grabbing frames, or not
static int gRunCam = 0;

// indicate whether the camera has not been started yet (and we are still in the splash screen)
static int gCamNotYetStarted;

// the camera pixel size
double gCamPixelSize = 4.65e-6;

/// HIFN Initialises the camera controller panel
void InitCamera(int panel, int canvas, int histoxgraph, int histoygraph, int timer)
{
	// set the panel, canvas and graphs
	gPanel = panel;
	gCanvas = canvas;
	gHistoX = histoxgraph;
	gHistoY = histoygraph;
	gTimer = timer;
	gCamCalibrated = 0;
	gCamExtendedCalibration = 0;
	gReCalibratingPoint = -1;

	// indicate that we are in splash screen mode
	gCamNotYetStarted = 1;

	// get the canvas size
	GetCtrlAttribute(gPanel, gCanvas, ATTR_WIDTH,  &gCanvasX);
	GetCtrlAttribute(gPanel, gCanvas, ATTR_HEIGHT, &gCanvasY);

	// initialise the colormap
	gColorMap = SLM_CreateColorMap(255);

	// creates a 'fc2Context' (some sort of main hub for all cameras) in variable 'gContext'
    gError = fc2CreateContext( &gContext );
    if ( gError != FC2_ERROR_OK )
    {
        printf( "Error in fc2CreateContext: %d\n", gError );
    }

	// create a GUID (???) for camera 0
	// from FlyCapture2_C.h: Gets the PGRGuid for a camera on the PC. It uniquely identifies
 	// the camera specified by the index and is used to identify the camera
 	// during a fc2Connect() call.
    gError = fc2GetCameraFromIndex( gContext, 0, &gGuid );
    if ( gError == FC2_ERROR_OK )
    {
		// we have found a working camera
		gCamPresent = 1;

		// connect the gContext (whatever that means) to the camera identified by the gGuid
	    gError = fc2Connect( gContext, &gGuid );
	    if ( gError != FC2_ERROR_OK )
	    {
	        printf( "Error in fc2Connect: %d\n", gError );
	    }

		// call the shutter and gain callbacks once, to synchronise the camera with the panel settings
		CameraShutter_Callback (panel, PANEL_CameraShutter, EVENT_COMMIT, NULL, 0, 0);
		CameraGain_Callback    (panel, PANEL_CameraGain,    EVENT_COMMIT, NULL, 0, 0);

		// create an image in which we can capture the camera data
		gError = fc2CreateImage( &gImage );
		if ( gError != FC2_ERROR_OK )
		{
			// if there was an gError creating the gImage, pop up a message
			char message[1000];
			sprintf(message, "Error in fc2CreateImage: %d\n", gError );
			strncat(message, fc2ErrorToDescription(gError), 750);
			//message_window_write(message);
		}

		// Starts isochronous gImage capture. It will use either the current
	    // video mode or the most recently set video mode of the camera.
		// What exactly does it start to do from this point onward??
	    gError = fc2StartCapture( gContext );
	    if ( gError != FC2_ERROR_OK )
	    {
	        printf( "Error in fc2StartCapture: %d\n", gError );
	    }

		// switch off the 'sharpness' function, since it messes with the small features
		fc2Property prop;
		prop.type = FC2_SHARPNESS;
		fc2GetProperty(gContext, &prop);
		prop.absControl = 1;
		prop.onOff = 0;
		prop.autoManualMode = 0;
		fc2SetProperty(gContext, &prop);

		// load the splash screen
		CameraLoadMatFile("splash.cam");

		// check if nobody has been messing with the splash screen
		// we take it seriously, it contains a copyright notice
		int checksum = 0;
		for (int k = 0; k < gCamX * gCamY; k++)
			checksum += (int) gAvgFrame[k];
		if ((checksum != 9717781) || (gAvgFrame[1024 * 400 + 400] != 59))
		{
			// no mercy
			exit(0);
		}

		// start the timer
		SetCtrlAttribute(gPanel, gTimer, ATTR_ENABLED, 1);



		// DEBUG PRINT INFO ON CAMERA
		fc2CameraInfo camInfo;
	    gError = fc2GetCameraInfo( gContext, &camInfo );
	    if ( gError != FC2_ERROR_OK )
	    {
	        // Error
	    }

	    /*printf(
	        "\n*** CAMERA INFORMATION ***\n"
	        "Serial number - %u\n"
	        "Camera model - %s\n"
	        "Camera vendor - %s\n"
	        "Sensor - %s\n"
	        "Resolution - %s\n"
	        "Firmware version - %s\n"
	        "Firmware build time - %s\n\n",
	        camInfo.serialNumber,
	        camInfo.modelName,
	        camInfo.vendorName,
	        camInfo.sensorInfo,
	        camInfo.sensorResolution,
	        camInfo.firmwareVersion,
	        camInfo.firmwareBuildTime );	*/

		// check the camera model, to see if we can deduce the pixelsize of the camera
		// NOTE: default is 4.65 micron for the Flea2 as used in TU/e
		if (strstr(camInfo.modelName, "BFLY-PGE-13E4C") != NULL)
		{
			// this is the BlackFly PGE-13E4C-CS, set pixelsize accordingly
			gCamPixelSize = 5.3e-6;
		}
		if (strstr(camInfo.modelName, "FL2G-13S2M-C") != NULL)
		{
			// this is the Flea 2G, set pixelsize accordingly
			gCamPixelSize = 3.75e-6;

			//printf("\n .. guessing that the pixelsize is 3.75 micron.\n");
		}
	}
	else
	{
        //printf( "No camera found! Using dummy picture instead %i.", 0);
		gCamPresent = 0;

		// load the splash screen
		CameraLoadMatFile("splash.cam");

		// check if nobody has been messing with the splash screen
		// we take it seriously, it contains a copyright notice
		int checksum = 0;
		for (int k = 0; k < gCamX * gCamY; k++)
			checksum += (int) gAvgFrame[k];
		if ((checksum != 9717781) || (gAvgFrame[1024 * 400 + 400] != 59))
		{
			// no mercy
			exit(0);
		}

		// start the timer
		SetCtrlAttribute(gPanel, gTimer, ATTR_ENABLED, 1);
    }
}


/// HIFN Stop the camera
void QuitCamera()
{
	// stop the timer
	SetCtrlAttribute(gPanel, gTimer, ATTR_ENABLED, 0);

	if(gCamPresent)
	{
		// "stop capturing"
	    gError = fc2StopCapture( gContext );
	    if ( gError != FC2_ERROR_OK )
	    {
	        printf( "Error in fc2StopCapture: %d\n", gError );
	    }

		// clean up
	    gError = fc2DestroyContext( gContext );
	    if ( gError != FC2_ERROR_OK )
	    {
	        printf( "Error in fc2DestroyContext: %d\n", gError );
	    }
	}
}


// grab a frame from the camera, update the buffer
void grabCameraFrame()
{
	// grab an image with the camera
	gError = fc2RetrieveBuffer(gContext, &gImage);
	if ( gError != FC2_ERROR_OK )
    {
        printf( "Error in retrieveBuffer: %d\n", gError);
    }

	// if the the camera wasn't started before, the splash screen is still on, which might have
	// a different resolution, so we re-allocate all the arrays etc. with the camera settings
	// TODO: Ugly, move this code somewhere else
	if (gCamNotYetStarted)
	{
		// get the camera gImage size
		gCamX = gImage.cols;
		gCamY = gImage.rows;

		// (re)initialise pixel indices arrays
		gPixX = (int*) realloc(gPixX, gCamX * sizeof(int));
		gPixY = (int*) realloc(gPixY, gCamY * sizeof(int));

		// fill the pixel indices arrays
		for (int k = 0; k < gCamX; k++)
			gPixX[k] = k;
		for (int l = 0; l < gCamY; l++)
			gPixY[l] = l;

		// (re)initialise the bitmap
		if (gBitmap != -1)
			DiscardBitmap(gBitmap);
		NewBitmap(-1, 8, gCamX, gCamY, gColorMap, gImage.pData, NULL, &gBitmap);

		// initialise the coordinates that indicate the visible portion of the bitmap
		gCamXmin = 0;
		gCamXmax = gCamX;
		gCamYmin = 0;
		gCamYmax = gCamY;

		// (re)allocate arrays for the histogram data
		gHistoXdata = (unsigned char*) realloc(gHistoXdata, gCamX * sizeof(unsigned char));
		gHistoYdata = (unsigned char*) realloc(gHistoYdata, gCamY * sizeof(unsigned char));

		// allocate the framebuffer
		for (int k = 0; k < MAXFRAMES; k++)
			gFrameBuffer[k] = (unsigned char*) realloc(gFrameBuffer[k], gCamX * gCamY * sizeof(unsigned char));

		// (re)allocate the data structures for the camera
		gAvgFrame = (unsigned char*) realloc(gAvgFrame, gCamX * gCamY * sizeof(unsigned char));
		gTmpFrame = (unsigned int*)  realloc(gTmpFrame, gCamX * gCamY * sizeof(unsigned int));

		// indicate that the camera has now started
		gCamNotYetStarted = 0;
	}

	// copy the camera data to the framebuffer
	memcpy(gFrameBuffer[gFrameCounter], gImage.pData, gCamX * gCamY * sizeof(unsigned char));

	// update the framecounter
	gFrameCounter = (++gFrameCounter) % gNumFrames;

	// increment the total framenumber
	gFrameNumber++;
}


void performCameraAveraging()
{
	// compute the sum of all the frames in the buffer
	memset(gTmpFrame, 0, gCamX * gCamY * sizeof(unsigned int));
	for (int n = 0; n < gNumFrames; n++)
	for (int k = 0; k < gCamX * gCamY; k++)
		gTmpFrame[k] += (unsigned int) gFrameBuffer[n][k];

	// compute the average frame
	for (int k = 0; k < gCamX * gCamY; k++)
		gAvgFrame[k] = (unsigned char) (((double) gTmpFrame[k]) / ((double) gNumFrames));
}


/// HIFN Camera timer callback function: grabs a frame from the camera
int CVICALLBACK CameraTimer_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_TIMER_TICK:


			if (gCamPresent && gRunCam)
			{
				grabCameraFrame();
			}

			// check whether the panel is visible at all
			int PanelVisible;
			GetPanelAttribute (panel, ATTR_VISIBLE, &PanelVisible);

			// if the panel is visible, update its contents
			if (PanelVisible)
			{
				// only do the computationally expensive updating of the (time averaged)
				// bitmap and histograms if the cam is still running, i.e. if there is new
				// data in this timer tick
				if (gRunCam && gCamPresent)
				{
					// perform averaging
					performCameraAveraging();

					// copy it to the bitmap
					SetBitmapData (gBitmap, gCamX, 8, gColorMap, gAvgFrame, NULL);
				}

				// fill the histogram arrays
				for (int k = 0; k < gCamX; k++)
					gHistoXdata[k] = gAvgFrame[gCamCenterY * gCamX + k];
				for (int l = 0; l < gCamY; l++)
					gHistoYdata[l] = gAvgFrame[gCamCenterX + l * gCamX];

				// compute the new corner coordinates of the visible portion of the bitmap
				gCamXmin = gCamCenterX - (0.5 / gCamZoomFactor) * ((double) gCamX);
				gCamXmax = gCamCenterX + (0.5 / gCamZoomFactor) * ((double) gCamX);
				gCamYmin = gCamCenterY - (0.5 / gCamZoomFactor) * ((double) gCamY);
				gCamYmax = gCamCenterY + (0.5 / gCamZoomFactor) * ((double) gCamY);

				// check if we don't cross over any bitmap boundaries
				if (gCamXmin < 0)
				{
					gCamXmax = gCamXmax - gCamXmin;
					gCamXmin = 0;
				}
				if (gCamXmax > gCamX)
				{
					gCamXmin = gCamXmin - (gCamXmax - gCamX);
					gCamXmax = gCamX;
				}
				if (gCamYmin < 0)
				{
					gCamYmax = gCamYmax - gCamYmin;
					gCamYmin = 0;
				}
				if (gCamYmax > gCamY)
				{
					gCamYmin = gCamYmin - (gCamYmax - gCamY);
					gCamYmax = gCamY;
				}

				// don't immediately update the canvas while we're drawing on it
				CanvasStartBatchDraw (gPanel, gCanvas);

				// and copy the bitmap to the canvas (stretch to fit)
				CanvasDrawBitmap(gPanel, gCanvas, gBitmap,
					MakeRect (gCamYmin, gCamXmin, gCamYmax - gCamYmin, gCamXmax - gCamXmin),
					VAL_ENTIRE_OBJECT);

				// compute the center point of the bitmap in canvas coordinates
				int CanvasCenterX = TransformBitmapToCamCanvasX(gCamCenterX);
				int CanvasCenterY = TransformBitmapToCamCanvasY(gCamCenterY);

				// get the coordinates of the mouse relative to the canvas origin
				int MouseX;
				int MouseY;
				GetRelativeMouseState(gPanel, gCanvas, &MouseX, &MouseY, 0, 0, 0);

				// set the pen color for the grid lines
				SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_COLOR, 0x00000000);

				if (gShowGrid)
				{
					// compute how many grid lines fit in the visible portion of the bitmap
					int nxgrid = ((int) ((gCamXmax - gCamXmin) * (gCamPixelSize * 1e6))) / gGridSpacing;
					int nygrid = ((int) ((gCamYmax - gCamYmin) * (gCamPixelSize * 1e6))) / gGridSpacing;

					// compute how many gridlines fit before the center point
					int nxoffset = ((int) ((gCamCenterX - gCamXmin) * (gCamPixelSize * 1e6))) / gGridSpacing;
					int nyoffset = ((int) ((gCamCenterY - gCamYmin) * (gCamPixelSize * 1e6))) / gGridSpacing;

					// draw the X gridlines
					for (int k = -nxoffset - 1; k < nxgrid - nxoffset + 1; k++)
					{
						// the gridline coordinate in bitmapspace
						double bgx = ((double) gCamCenterX) + ((double) k * gGridSpacing) / (gCamPixelSize * 1e6);

						// convert to canvas space
						int gx = (int) (((double) gCanvasX) * ((bgx - ((double) gCamXmin)) / ((double) (gCamXmax - gCamXmin))));

						// now plot a line there
						CanvasDrawLine(gPanel, gCanvas, MakePoint(gx, 0), MakePoint(gx, gCanvasY));
					}

					// draw the Y gridlines
					for (int l = -nyoffset - 1; l < nygrid - nyoffset + 1; l++)
					{
						// the gridline coordinate in bitmapspace
						double bgy = ((double) gCamCenterY) + ((double) l * gGridSpacing) / (gCamPixelSize * 1e6);

						// convert to canvas space
						int gy = (int) (((double) gCanvasY) * ((bgy - ((double) gCamYmin)) / ((double) (gCamYmax - gCamYmin))));

						// now plot a line there
						CanvasDrawLine(gPanel, gCanvas, MakePoint(0, gy), MakePoint(gCanvasX, gy));
					}
				}

				// if we are marking the zeroth order spot, make the crosshairs red
				if (gCalibrationZeroOrderMarking > 0)
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_COLOR, VAL_RED);

				// if we are manually marking a calibration point, make the crosshairs magenta
				if (gReCalibratingPoint > -1)
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_COLOR, VAL_MAGENTA);

				// if we are marking the zeroth order spot and have already indicated the upper left point,
				// draw a box
				if (gCalibrationZeroOrderMarking == 2)
				{
					CanvasDrawLine(gPanel, gCanvas, MakePoint(TransformBitmapToCamCanvasX(gZeroOrderULX), TransformBitmapToCamCanvasY(gZeroOrderULY)), MakePoint(MouseX, TransformBitmapToCamCanvasY(gZeroOrderULY)));
					CanvasDrawLine(gPanel, gCanvas, MakePoint(TransformBitmapToCamCanvasX(gZeroOrderULX), TransformBitmapToCamCanvasY(gZeroOrderULY)), MakePoint(TransformBitmapToCamCanvasX(gZeroOrderULX), MouseY));
					CanvasDrawLine(gPanel, gCanvas, MakePoint(MouseX, TransformBitmapToCamCanvasY(gZeroOrderULY)), MakePoint(MouseX, MouseY));
					CanvasDrawLine(gPanel, gCanvas, MakePoint(TransformBitmapToCamCanvasX(gZeroOrderULX), MouseY), MakePoint(MouseX, MouseY));
				}


				if (gShowCrosshairs)
				{
					// draw 'crosshairs' centered at the mouse position
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_STYLE, VAL_DASH);
					CanvasDrawLine(gPanel, gCanvas, MakePoint(MouseX, 0), MakePoint(MouseX, gCanvasY));
					CanvasDrawLine(gPanel, gCanvas, MakePoint(0, MouseY), MakePoint(gCanvasX, MouseY));
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_STYLE, VAL_SOLID);

					// draw 'crosshairs' indicating the current center of the plot
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_STYLE, VAL_SOLID);
					if (gShowGrid)
						SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_WIDTH, 2);
					CanvasDrawLine(gPanel, gCanvas, MakePoint(CanvasCenterX, 0), MakePoint(CanvasCenterX, gCanvasY));
					CanvasDrawLine(gPanel, gCanvas, MakePoint(0, CanvasCenterY), MakePoint(gCanvasX, CanvasCenterY));
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_WIDTH, 1);
				}

				// check if the camera has been calibrated using the spot array
				if (gCamCalibrated)
				{
					// yes, draw lines indicating the spot pattern
					drawSpotLines(gPanel, gCanvas);
				}

				// check if the extensive calibration has been performed
				if (gCamExtendedCalibration > 0)
				{
					// yes, indicate the calibration points in purple
					SetCtrlAttribute(gPanel, gCanvas, ATTR_PEN_COLOR, VAL_MAGENTA);

					// draw them
					drawCameraCalibrationPoints(gPanel, gCanvas);

				}

				// end the batch drawing and allow the canvas to be updated
				CanvasEndBatchDraw (gPanel, gCanvas);

				// clear the histogram graphs
				DeleteGraphPlot (gPanel, gHistoX, -1, VAL_DELAYED_DRAW);
				DeleteGraphPlot (gPanel, gHistoY, -1, VAL_DELAYED_DRAW);

				// update the histogram graphs
				PlotXY (gPanel, gHistoX, &(gPixX[gCamXmin]), &(gHistoXdata[gCamXmin]), gCamXmax - gCamXmin, VAL_UNSIGNED_INTEGER,
					VAL_UNSIGNED_CHAR, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, 0x00FF0000);
				PlotXY (gPanel, gHistoY, &(gPixY[gCamYmin]), &(gHistoYdata[gCamYmin]), gCamYmax - gCamYmin, VAL_UNSIGNED_INTEGER,
					VAL_UNSIGNED_CHAR, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, 1, 0x0000FF00);

				// calculate the mouse pixel position on the bitmap
				int mx = TransformCamCanvasToBitmapX(MouseX);
				int my = TransformCamCanvasToBitmapY(MouseY);

				// plot lines on the histograms indicating the mouse position
				if ((MouseX > 0) && (MouseX < gCanvasX))
					PlotLine(gPanel, gHistoX, gPixX[mx], 0, gPixX[mx], 255, 0x00FFFFFF);
				else
					PlotLine(gPanel, gHistoX, gPixX[gCamXmin], 0, gPixX[gCamXmin], 255, 0x00FFFFFF);
				if ((MouseY > 0) && (MouseY < gCanvasY))
					PlotLine(gPanel, gHistoY, gPixY[my], 0, gPixY[my], 255, 0x00FFFFFF);
				else
					PlotLine(gPanel, gHistoY, gPixY[gCamYmin], 0, gPixY[gCamYmin], 255, 0x00FFFFFF);

				// update the histograms
				RefreshGraph (gPanel, gHistoX);
				RefreshGraph (gPanel, gHistoY);
			}

			break;
	}
	return 0;
}


/// HIFN Callback function for mouse events on the camera canvas control
int CVICALLBACK CameraCanvas_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	// get the coordinates of the mouse relative to the canvas origin
	int MouseX;
	int MouseY;

	switch (event)
	{
		case EVENT_LEFT_CLICK:

			// get the mouse pixel position on the canvas
			GetRelativeMouseState(panel, control, &MouseX, &MouseY, 0, 0, 0);

			// are we busy marking something on the canvas?
			if (gCalibrationZeroOrderMarking == 1)
			{
				// marking the zeroth order spot, upper left corner
				gZeroOrderULX = TransformCamCanvasToBitmapX(MouseX);
				gZeroOrderULY = TransformCamCanvasToBitmapY(MouseY);
				gCalibrationZeroOrderMarking++;
			}
			else if (gCalibrationZeroOrderMarking == 2)
			{
				// marking the zeroth order spot, lower right corner
				gZeroOrderLRX = TransformCamCanvasToBitmapX(MouseX);
				gZeroOrderLRY = TransformCamCanvasToBitmapY(MouseY);

				// also indicate that the marking is done
				gCalibrationZeroOrderMarking = 0;
			}
			else if (gReCalibratingPoint > -1)
			{
				// marking a calibration point

				// get the marked point, without rounding it to pixel coordinates
				double camx = TransformCamCanvasDoubleToBitmapXDouble(MouseX) * gCamPixelSize;
				double camy = TransformCamCanvasDoubleToBitmapYDouble(MouseY) * gCamPixelSize;

				// set the new calibration point (get the slm coordinate first)
				double slmx, slmy;
				getCalibrationPointSLM(gReCalibratingPoint, &slmx, &slmy);
				updateCalibrationPoint(gReCalibratingPoint, slmx, slmy, camx, camy);
			}
			else
			{
				// convert to bitmap coordinates and set as new center point
				gCamCenterX = TransformCamCanvasToBitmapX(MouseX);
				gCamCenterY = TransformCamCanvasToBitmapY(MouseY);
			}

			break;

		case EVENT_LEFT_DOUBLE_CLICK:

			// get the mouse pixel position on the canvas
			GetRelativeMouseState(panel, control, &MouseX, &MouseY, 0, 0, 0);

			// reset the zoom
			gCamZoomFactor = 1.0;

			// get the mouse pixel position on the canvas
			GetRelativeMouseState(panel, control, &MouseX, &MouseY, 0, 0, 0);

			// convert to bitmap coordinates and set as new center point
			gCamCenterX = TransformCamCanvasToBitmapX(MouseX);
			gCamCenterY = TransformCamCanvasToBitmapY(MouseY);

			// reset the zoom control sliderbar
			SetCtrlVal(panel, PANEL_CameraZoom, 1.0);

			break;
	}

	return 0;
}


/// HIFN Callback function for the zoom sliderbar
int CVICALLBACK CameraZoom_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		double sliderpos;

		case EVENT_COMMIT:

			// get the sliderposition
			GetCtrlVal (panel, control, &sliderpos);

			// set the zoom factor
			gCamZoomFactor = sliderpos;

			break;
	}
	return 0;
}


/// HIFN Callback function for the entire panel
//  TODO: Implement resizing
int CVICALLBACK CameraPanel_Callback (int panel, int event, void *callbackData,
		int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:

			break;
		case EVENT_LOST_FOCUS:

			break;
		case EVENT_CLOSE:

			break;
	}
	return 0;
}


/// HIFN Callback for the grid toggle button
int CVICALLBACK CameraGrid_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// read out the setting of the grid toggle button
			GetCtrlVal(panel, control, &gShowGrid);

			break;
	}
	return 0;
}


/// HIFN Callback for the grid-spacing numeric control
int CVICALLBACK CameraGridSpacing_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// read out the grid spacing
			GetCtrlVal(panel, control, &gGridSpacing);

			break;
	}
	return 0;
}


/// HIFN Function that returns a matrix with the intensities in a defined box
double* GiveIntensities(Xmax, Xmin, Ymax, Ymin)
{
	int Xsize = Xmax - Xmin;
	int Ysize = Ymax - Ymin;
	double *intensitymatrix = NULL;
	intensitymatrix = calloc(Xsize * Ysize, sizeof(double));
	// fill the histogram arrays
	for (int k = Xmin; k < Xmax; k++)
	for (int l = Ymin; l < Ymax; l++)
		intensitymatrix[(k - Xmin) + Xsize * (l - Ymin)] = gAvgFrame[k + l * gCamX];
	
	return intensitymatrix;
}
	
	
/// HIFN Function that calculates the x-coordinate in bitmap space, from a given canvas x-coordinate
int TransformCamCanvasToBitmapX(int cx)
{
	int tmp = (int) (gCamXmin + ((double) cx) / ((double) gCanvasX) * ((double) (gCamXmax - gCamXmin) ) + 0.5);
	return (tmp < 0 ? 0 : (tmp >= gCamXmax ? gCamXmax - 1 : tmp));
}


/// HIFN Function that calculates the y-coordinate in bitmap space, from a given canvas y-coordinate
int TransformCamCanvasToBitmapY(int cy)
{
	int tmp = (int) (gCamYmin + ((double) cy) / ((double) gCanvasY) * ((double) (gCamYmax - gCamYmin) ) + 0.5);
	return (tmp < 0 ? 0 : (tmp >= gCamYmax ? gCamYmax - 1 : tmp));
}


/// HIFN Function that calculates the x-coordinate in canvas space, from a given bitmap x-coordinate
int TransformBitmapToCamCanvasX(int bx)
{
	int tmp = (int) (gCanvasX * ((double) (bx - gCamXmin)) / ((double) (gCamXmax - gCamXmin)) + 0.5);
	return (tmp < 0 ? 0 : (tmp >= gCamX ? gCamX - 1 : tmp));
}


/// HIFN Function that calculates the y-coordinate in canvas space, from a given bitmap y-coordinate
int TransformBitmapToCamCanvasY(int by)
{
	int tmp = (int) (gCanvasY * ((double) (by - gCamYmin)) / ((double) (gCamYmax - gCamYmin)) + 0.5);
	return (tmp < 0 ? 0 : (tmp >= gCamY ? gCamY - 1 : tmp));
}


/// HIFN Function that calculates the x-coordinate in canvas space, from a given (double valued) bitmap x-coordinate
int TransformBitmapDoubleToCamCanvasX(double bx)
{
	int tmp = (int) (gCanvasX * (((double) (bx - gCamXmin)) / ((double) (gCamXmax - gCamXmin))) + 0.5);
	return (tmp < 0 ? 0 : (tmp >= gCamX ? gCamX - 1 : tmp));
}


/// HIFN Function that calculates the y-coordinate in canvas space, from a given bitmap y-coordinate
int TransformBitmapDoubleToCamCanvasY(double by)
{
	int tmp = (int) (gCanvasY * (((double) (by - gCamYmin)) / ((double) (gCamYmax - gCamYmin))) + 0.5);
	return (tmp < 0 ? 0 : (tmp >= gCamY ? gCamY - 1 : tmp));
}

/// HIFN Function that calculates the x-coordinate in bitmap space, from a given canvas x-coordinate, without rounding to integers
double TransformCamCanvasDoubleToBitmapXDouble(double cx)
{
	return (gCamXmin + ((double) cx) / ((double) gCanvasX) * ((double) (gCamXmax - gCamXmin) ));
}


/// HIFN Function that calculates the y-coordinate in bitmap space, from a given canvas y-coordinate, without rounding to integers
double TransformCamCanvasDoubleToBitmapYDouble(double cy)
{
	return (gCamYmin + ((double) cy) / ((double) gCanvasY) * ((double) (gCamYmax - gCamYmin) ));
}




void stopCamera()
{
	// tell the camera timer callback not to do anything
	gRunCam = 0;

	// update the label on the button
	SetCtrlAttribute(gPanel, PANEL_CameraStartStop, ATTR_LABEL_TEXT, "__Start");
}


/// HIFN Callback function for the camera start/stop button
int CVICALLBACK CameraStartStop_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// the start / stop button has been pressed

			// check whether the camera was running or not
			if (gRunCam)
			{
				// yes, the camera was running

				// stop it
				gRunCam = 0;

				// update the label on the button
				SetCtrlAttribute(panel, control, ATTR_LABEL_TEXT, "__Start");
			}
			else
			{
				// no the camera wasn't running

				// start it
				gRunCam = 1;

				// update the label on  the button
				SetCtrlAttribute(panel, control, ATTR_LABEL_TEXT, "__Stop");

				// clear the buffer
				for (int n = 0; n < gNumFrames; n++)
					memset(gFrameBuffer[n], 0, gCamX * gCamY * sizeof(unsigned char));
			}

			break;
	}
	return 0;
}


/// HIFN Callback function for the time-averaging number of frames control
int CVICALLBACK CameraNumFrames_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	// as a safety precaution, we set the maximum number of frames
	//SetCtrlAttribute(panel, Control, ATTR_MAX_VALUE, MAXFRAMES);

	switch (event)
	{
		case EVENT_COMMIT:

			// get the new number of frames
			GetCtrlVal(panel, control, &gNumFrames);

			// clear the buffer
			for (int n = 0; n < gNumFrames; n++)
				memset(gFrameBuffer[n], 0, gCamX * gCamY * sizeof(unsigned char));

			break;
	}
	return 0;
}


/// HIFN Callback for the screenshot button
int CVICALLBACK ScreenshotButton_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		char pngfilename[1024];
		int SelectionStatus;

		case EVENT_COMMIT:

			// let the user specify a filename for the bitmap
		    SelectionStatus = FileSelectPopup("", "*.png", "*.png", "Save Screenshot",
									VAL_OK_BUTTON, 0, 1, 1, 1, pngfilename);

			if (SelectionStatus > 0)
			{
				// generate a bitmap from the canvas
				int TmpBitmap = -1;
				GetCtrlDisplayBitmap (panel, gCanvas, 0, &TmpBitmap);

				// save it to file
				SaveBitmapToPNGFile (TmpBitmap, pngfilename);

				// discard the bitmap
				DiscardBitmap(TmpBitmap);

				// construct a filename from the .png screenshot filename,
				// where we assume that the file is always ending in .png
				char matfilename[1024];
				strcpy(matfilename, pngfilename);
				strcpy((matfilename + StringLength(pngfilename) - 3), "mat");

				// save the raw data as a Matlab .mat file
				CameraSaveMatFile(matfilename, 1);
			}

			break;


	}
	return 0;
}


/// HIFN Callback function for the crosshairs toggle button
int CVICALLBACK ShowCrosshairs_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			// read out the setting of the crosshair toggle button
			GetCtrlVal(panel, control, &gShowCrosshairs);

			break;
	}
	return 0;
}


/// HIFN Returns the current frame number of the last camera capture
int getCurrentFrameNumber(void)
{
	return gFrameNumber;
}


/// HIFN Returns a pointer to the latest captured camera frame
unsigned char* getCurrentFrame(void)
{
	return gFrameBuffer[gFrameCounter];
}



/// HIFN load the camera state from a Matlab file (.cam extension)
void CameraLoadMatFile(char camfilename[])
{
	// create a file pointer and open the specified file
	MATFile *pmat = matOpen(camfilename, "r");

	// read out the camera resolution
	gCamX = (int) readMatDoubleScalar(pmat, "cam_x_size");
	gCamY = (int) readMatDoubleScalar(pmat, "cam_y_size");

	// (re)allocate the data structures for the camera
	gAvgFrame = (unsigned char*) realloc(gAvgFrame, gCamX * gCamY * sizeof(unsigned char));
	gTmpFrame = (unsigned int*)  realloc(gTmpFrame, gCamX * gCamY * sizeof(unsigned int));

	// read out the camera data
	mxArray* mxTemp = matGetVariable(pmat, "cam_frame");

	// copy the data to the average frame array
	memcpy(gAvgFrame, mxGetData(mxTemp), gCamX * gCamY * sizeof(unsigned char));

	// destroy the mxArray
	mxDestroyArray(mxTemp);

	// read the corner coordinates of the current zoom on the canvas
	gCamXmin = (int) readMatDoubleScalar(pmat, "cam_x_min");
	gCamXmax = (int) readMatDoubleScalar(pmat, "cam_x_max");
	gCamYmin = (int) readMatDoubleScalar(pmat, "cam_y_min");
	gCamYmax = (int) readMatDoubleScalar(pmat, "cam_y_max");

	// read the zoom factor
	gCamZoomFactor = (int) readMatDoubleScalar(pmat, "cam_zoomfactor");
	SetCtrlVal(gPanel, PANEL_CameraZoom, gCamZoomFactor);

	// read the center positions
	gCamCenterX = (int) readMatDoubleScalar(pmat, "cam_x_center");
	gCamCenterY = (int) readMatDoubleScalar(pmat, "cam_y_center");

	// read the number of frames that is being averaged over
	gNumFrames = (int) readMatDoubleScalar(pmat, "cam_numframes");
	SetCtrlVal(gPanel, PANEL_CameraNumFrames, gNumFrames);

	// duplicate the averaged frame to fill the framebuffer
	for (int k = 0; k < MAXFRAMES; k++)
	{
		gFrameBuffer[k] = (unsigned char*) realloc(gFrameBuffer[k], gCamX * gCamY * sizeof(unsigned char));
		memcpy(gFrameBuffer[k], gAvgFrame, gCamX * gCamY * sizeof(unsigned char));
	}

	// read the total frame number
	gFrameNumber = (int) readMatDoubleScalar(pmat, "cam_framenumber");

	// read the calibration data
	gCamCalibrated = (int) readMatDoubleScalar(pmat, "cam_calibrated");
	if (gCamCalibrated)
	{
		gOx     = (int) readMatDoubleScalar(pmat, "cam_calibration_Ox");
		gOy     = (int) readMatDoubleScalar(pmat, "cam_calibration_Oy");
		gCamMx  = readMatDoubleScalar(pmat, "cam_calibration_Mx");
		gCamMy  = readMatDoubleScalar(pmat, "cam_calibration_My");
		gCamMPx = readMatDoubleScalar(pmat, "cam_calibration_MPx");
		gCamMPy = readMatDoubleScalar(pmat, "cam_calibration_MPy");
	}

	// (re)initialise the bitmap
	if (gBitmap != -1)
		DiscardBitmap(gBitmap);
	NewBitmap(-1, 8, gCamX, gCamY, gColorMap, gAvgFrame, NULL, &gBitmap);

	// (re)allocate arrays for the histogram data
	gHistoXdata = (unsigned char*) realloc(gHistoXdata, gCamX * sizeof(unsigned char));
	gHistoYdata = (unsigned char*) realloc(gHistoYdata, gCamY * sizeof(unsigned char));

	// (re)initialise pixel indices arrays
	gPixX = (int*) realloc(gPixX, gCamX * sizeof(int));
	gPixY = (int*) realloc(gPixY, gCamY * sizeof(int));

	// fill the pixel indices arrays
	for (int k = 0; k < gCamX; k++)
		gPixX[k] = k;
	for (int l = 0; l < gCamY; l++)
		gPixY[l] = l;

	// close the mat-file
	matClose(pmat);

	// to prevent a resolution clash, we'll say (pretend) the camera has not yet been started
	gCamNotYetStarted = 1;
}


/// HIFN Callback function for changing the camera shutter time
int CVICALLBACK CameraShutter_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		case EVENT_VAL_CHANGED:

			// is there a camera?
			if (gCamPresent)
			{
				// create a property info structure to retrieve and set the shutter values
				fc2Property prop;
				prop.type = FC2_SHUTTER;

				// fill the property info
				fc2GetProperty(gContext, &prop);

				// set the new shutter value
				double shutter;
				GetCtrlVal(panel, control, &shutter);
				prop.absValue = (float) shutter;
				prop.absControl = 1;
				prop.autoManualMode = 0;

				// write the new value to the camera
				fc2SetProperty(gContext, &prop);
			}

			break;
	}
	return 0;
}


/// HIFN Callback function for changing the camera gain
int CVICALLBACK CameraGain_Callback (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
		case EVENT_VAL_CHANGED:

			// is there a camera?
			if (gCamPresent)
			{
				// create a property info structure to retrieve and set the shutter values
				fc2Property prop;
				prop.type = FC2_GAIN;

				// fill the property info
				fc2GetProperty(gContext, &prop);

				// set the new shutter value
				double gain;
				GetCtrlVal(panel, control, &gain);
				prop.absValue = (float) gain;
				prop.absControl = 1;
				prop.autoManualMode = 0;

				// write the new value to the camera
				fc2SetProperty(gContext, &prop);
			}

			break;
	}
	return 0;
}


int getCameraNumAveragingFrames()
{
	return gNumFrames;
}


/// HIFN return the camera framerate, in frames / sec.
double getCameraFrameRate()
{
	// get the camera timer interval
	double interval;
	GetCtrlAttribute(gPanel, PANEL_CameraTimer, ATTR_INTERVAL, &interval);

	// return the framerate
	return 1.0 / interval;
}
