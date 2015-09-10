//==============================================================================
//
// Title:       SLM Pixel Panel
// Purpose:     Callback routines for the panel displaying the SLM pixels 
//              (i.e., the phase pattern).
//
// Created on:  23-02-2012 by Rick.
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

/// HIFN Callback for the SLM pixels panel 
int CVICALLBACK SLMpixels_Callback (int panel, int event, void *callbackData,
		int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_PANEL_SIZE:
		{	
			// get the new panel size
			int xsize, ysize;
			GetPanelAttribute(panel, ATTR_WIDTH,  &xsize);
			GetPanelAttribute(panel, ATTR_HEIGHT, &ysize);
			
			// check if we are not dealing with an event fired after a programmatical resize
			if ((xsize != SLM_getXres()) || (ysize != SLM_getYres()))
			{
			
				// manually check for a maximize event, which sets the y size of the panel to 1061
				// so we use that value to detect a maximization (NOTE: This code is highly specific to 
				// the UCP lab computer and/or the Windows version running on that)
				if (ysize == 1061)
				{
					// yes, apparently we tried to maximize the panel
				
					// however, 1061 is not what we want, so we
					// set the size manually again to the desired SLM resolution
					SetPanelSize (panel, 1080, 1920);
				
					// set the upper left corner to the upper left corner of the screen,
					// and hope we have banned the titlebar
					SetPanelPos(panel, 0, 1920);
				}
			
				// adjust the canvas size accordingly
				SetCtrlAttribute(panel, SLMpixels_SLMcanvas, ATTR_WIDTH,  xsize);
				SetCtrlAttribute(panel, SLMpixels_SLMcanvas, ATTR_HEIGHT, ysize);
			
				// get the size of the simulation canvas
				int SimX, SimY;
				GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_WIDTH,  &SimX);
				GetCtrlAttribute(pnlSimPanel, SimPanel_CANVAS, ATTR_HEIGHT, &SimY);
				
				// get the subsampling factor
				int subsample;
				GetCtrlVal(TabPage_0, TABPANEL_SubSample, &subsample);

				// all data structures of the SLM need to be resized to acommodate the new 
				// number of pixels, so we just reintialise the SLM in its entirety			
				SLM_initialise(xsize, ysize, SimX, SimY, subsample);
			
				// update the SLM panel and the simulation Panel (if toggled)
				SLM_update(pnlSLMpixels, SLMpixels_SLMcanvas, pnlSimPanel, SimPanel_CANVAS, 1);
			
				// update the displayed information about the SLM 
				UpdateSLMinfo();
			}
				
			break;
		}
	}
	return 0;
}
