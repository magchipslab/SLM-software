/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/* Copyright (c) National Instruments 2011. All Rights Reserved.          */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  PANEL                            1       /* callback function: panelCB */
#define  PANEL_CLIENTMODECOMMAND          2       /* callback function: clientmodecallback */
#define  PANEL_EXAMPLE_NUMERIC            3       /* callback function: changecallback */
#define  PANEL_SCANCONTR_NUMERIC          4       /* callback function: changecallback */
#define  PANEL_EXP_MON                    5
#define  PANEL_MATLAB_SCRIPT              6       /* callback function: matlab_script_callback */
#define  PANEL_MATLAB_FIT_VALUE           7
#define  PANEL_ANALYSE_IN_MATLAB          8       /* callback function: analyze_in_matlab_callback */
#define  PANEL_RUN_MATLAB_SCRIPT          9       /* callback function: run_matlab_script_callback */


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK analyze_in_matlab_callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK changecallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK clientmodecallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK matlab_script_callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK panelCB(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK run_matlab_script_callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
