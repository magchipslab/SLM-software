/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  PANEL                           1
#define  PANEL_DATA                      2
#define  PANEL_ACCEPT_DATA               3       /* callback function: Accept */
#define  PANEL_CONFIDENCE                4
#define  PANEL_SAMPLES                   5
#define  PANEL_DONE                      6       /* callback function: Done */
#define  PANEL_OUTPUT                    7       /* callback function: OutputToSTDIO */
#define  PANEL_DEVIATION                 8
#define  PANEL_FREEDOM                   9
#define  PANEL_MEAN                      10
#define  PANEL_HELP                      11      /* callback function: HelpCallback */
#define  PANEL_QUIT                      12      /* callback function: Quit */
#define  PANEL_DECORATION                13
#define  PANEL_DECORATION_2              14
#define  PANEL_TEXTMSG                   15
#define  PANEL_DECORATION_3              16


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */ 

int  CVICALLBACK Accept(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Done(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK HelpCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OutputToSTDIO(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Quit(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
