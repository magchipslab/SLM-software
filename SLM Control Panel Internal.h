//==============================================================================
//
// Title:       SLM Control Panel Internal.h
// Purpose:     Header file for "SLM Pixel Panel.c", "SLM Simulation Panel.c" 
// 				and all "SLM Control Panel *.c" source files
//				All definitions in this file would ordinarily be placed in 
//				"SLM Control Panel.h", but LabWindows does not allow / support
//				additions to that file.
//
// Created on:  04-01-2012 at 11:25:02 by Rick van Bijnen
// Copyright:   Technische Universiteit Eindhoven. All Rights Reserved.
//
//==============================================================================


// includes needed for the MATfile type
#include "libmat/mat.h" 
#include "libmat/matrix.h"


// Flag for not updating the SLM pixels and simulation, such that we can call 
// the control callbacks ourselves
#define SLM_NO_UPDATE -1

// non-callback (convenience) functions 
void setPhaseGratingFromParameters(void);
void setSpotPatternFromParameters(void);
void setDefaultValuesFromPanel(void);
void UpdateSLMinfo(void);
void UpdateFeedbackCanvas(int Panel, int Canvas);
void WritePatternSettings(MATFile *pmat);
void LoadSLMSettings(char slmfilename[]);

// panel IDs
extern int pnlControl, pnlSLMpixels, pnlSimPanel, pnlCamera, pnlDebug;

// tabpanel IDs
extern int TabPage_0, TabPage_1, TabPage_2, TabPage_3, TabPage_4, TabPage_5, TabPage_6, TabPage_7, TabPage_1_0, TabPage_1_1, TabPage_1_2;
extern int TabPage_2, TabPage_2_0, TabPage_2_1;

	
