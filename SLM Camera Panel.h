//==============================================================================
//
// Title:       Camera Controller.h
// Purpose:     Header file for "Camera Controller.c"
//
// Created on:  13-12-2011 at 10:45:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================

#ifndef MAXFRAMES
#define MAXFRAMES 15
#endif

// the pixel size in micron
extern double gCamPixelSize;

// number of pixels of canvas and camera feed
extern int gCamX, gCamY;

// the pixel position of the bitmap (camera feed) that should be centered in the canvas
extern int gCamCenterX, gCamCenterY;

// the pixel coordinates of the corners of the visible portion of the bitmap
// (convenience variables, as these can in principle be computed from center and zoom variables)
extern int gCamXmin, gCamXmax, gCamYmin, gCamYmax;

// the camera zoom factor
extern double gCamZoomFactor;

// the number of frames to average over
extern int gNumFrames;

// the time-averaged camera frame
extern unsigned char* gAvgFrame;

// total frame number
extern int gFrameNumber;

// correspondence of simulation with camera feed: origin coordinates
extern int gOx, gOy;

// correspondence of simulation with camera feed: magnification factors (camera image size / expected image size)
// NOTE: these values are the ratios of the physical sizes, see gMPx and gMPy for the pixel ratios
extern double gCamMx, gCamMy;

// size of SLM focal unit in camera pixels
extern double gCamMPx, gCamMPy;

// indicator variable whether the camera has been calibrated
extern int gCamCalibrated;

// variable indicating whether the extensive calibration has been performed
extern int gCamExtendedCalibration;

// indicates whether the user is currently marking the zeroth order spot
extern int gCalibrationZeroOrderMarking;

// indicates whether the user is currently marking a calibrationpoint
extern int gReCalibratingPoint;

// the coordinates of a box containing the zeroth order spot (such that this region can be ignored while calibrating)
extern int gZeroOrderULX, gZeroOrderULY, gZeroOrderLRX, gZeroOrderLRY;

// get the SLM coordinates of a calibration point   
void getCalibrationPointSLM(int p, double *x, double *y);

// change a calibration point
void updateCalibrationPoint(int p, double SLMx, double SLMy, double CamX, double CamY);

void InitCamera(int panel, int canvas, int histox, int histoy, int timer);
void QuitCamera(void);

// convert coordinates between the camera feed (bitmap) and the canvas to display it on
int TransformCamCanvasToBitmapX(int cx);
int TransformCamCanvasToBitmapY(int cy);
int TransformBitmapToCamCanvasX(int bx);
int TransformBitmapToCamCanvasY(int by);
int TransformBitmapDoubleToCamCanvasX(double bx);
int TransformBitmapDoubleToCamCanvasY(double by);
double TransformCamCanvasDoubleToBitmapXDouble(double cx);
double TransformCamCanvasDoubleToBitmapYDouble(double cy);


// get the camera space coordinate of a focal plane coordinate
int TransformFocalPlaneToBitmapX(int fx);
int TransformFocalPlaneToBitmapY(int fy);

int getCurrentFrameNumber(void);
unsigned char* getCurrentFrame(void);

void CameraSaveMatFile(char matfilename[], int saveslmdata);

void CameraLoadMatFile(char matfilename[]);

double* getSignalWindowFromCamera(int ulx, int uly, int Wx, int Wy, int corrx, int corry);
double* getSignalWindowFromCameraResampled(double* dest, int ulx, int uly, int Wx, int Wy, int corrx, int corry, int Rx, int Ry);

void drawSpotLines(int Panel, int Canvas);

void stopCamera(void);
void grabCameraFrame(void);
void performCameraAveraging(void);

/// HIFN return the number of frames the camera uses for averaging
int getCameraNumAveragingFrames(void);

/// HIFN return a matrix with intensities from the camera in a specified box
double* GiveIntensities(int Xmax, int Xmin, int Ymax, int Ymin);

/// HIFN return the camera framerate, in frames / sec.
double getCameraFrameRate(void);

/// HIFN draws the camera calibration points on the canvas, using the current pen
void drawCameraCalibrationPoints(int panel, int canvas);

/// HIFN transform a coordinate according to two sets of coordinates which are in correspondence to each other
void transformCoordinate(double* MapFromX, double* MapFromY, double* MapToX, double* MapToY, int Nx, int Ny, double xfrom, double yfrom, double* xto, double* yto);

/// HIFN transform a coordinate from SLM space to camera space
// NOTE: slow code, not intended for heavy realtime use (use for precalculation)
void transformCoordinateSLMToCamera(double xfrom, double yfrom, double* xto, double* yto);

/// HIFN transform a coordinate from camera space to SLM space
// NOTE: slow code, not intended for heavy realtime use (use for precalculation)
void transformCoordinateCameraToSLM(double xfrom, double yfrom, double* xto, double* yto);

/// HIFN get the coordinates of the pixel at the origin of the camera window
void getCameraWindowOriginPixel(int *x, int *y);

/// HIFN gets a 1:1 window from the camera data
void getCameraWindow(int x0, int y0, int width, int height, double* dest);
