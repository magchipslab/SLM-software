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

#define  MESSAGE                          1       /* callback function: message_window_callback */
#define  MESSAGE_BUTTON_HIDE              2       /* callback function: Button_Hide */
#define  MESSAGE_BUTTON_CLOSE             3       /* callback function: Button_Close */
#define  MESSAGE_BUTTON_CLEAR             4       /* callback function: Button_Clear */
#define  MESSAGE_TEXTBOX                  5


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK Button_Clear(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Button_Close(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Button_Hide(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK message_window_callback(int panel, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
