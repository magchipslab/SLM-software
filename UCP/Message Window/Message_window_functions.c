/*------------------------------------------------------------------------------------------------------------------------
 *	File		Message_window_functions.c
 *	Author		Wouter Engelen
 *	Version		January 2010
 *		
 *	Details		UI for displaying messages
 *
------------------------------------------------------------------------------------------------------------------------*/
#include <utility.h>
#include <ansi_c.h>
#include "Message_window.h"

#define max_number_of_message_lines 500

static int message_window_Handle = 0;

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_open
	Arguments		-
	Return value	-
					
	Description		Opens the Message window UI
------------------------------------------------------------------------------------------------------------------------*/
void message_window_open(void)
{
	int program_title_length, panel_left_position, panel_top_position, panel_height, panel_width;
	char window_title[215] = "Message Window "; //200+15
	
	GetPanelAttribute(1, ATTR_TITLE_LENGTH, &program_title_length); //Main panel has always a panelhandle of 1
	if (program_title_length < 200)
		GetPanelAttribute(1, ATTR_TITLE, window_title+15);
	
	//#define message_window_path_uir_file "S:\\bunches\\cvi\\General\\Message_window.uir"
	//message_window_Handle = LoadPanel(0, message_window_path_uir_file, MESSAGE);
	message_window_Handle = LoadPanelEx(0, "Message_window.uir", MESSAGE, __CVIUserHInst); //When compiling as a dll
	SetPanelAttribute(message_window_Handle, ATTR_TITLE, window_title); //Set title of panel
	
	//Set message window position so its vertical center equals the program vertical center, 
	//and the right coordinate of the message window is 100 pixels more than the right coordinate of the program
	//so you can always see the message window, even if there is a message generated before the program gui is displayed
	GetPanelAttribute(1, ATTR_LEFT,   &panel_left_position);
	GetPanelAttribute(1, ATTR_TOP,    &panel_top_position);
	GetPanelAttribute(1, ATTR_WIDTH,  &panel_width);
	GetPanelAttribute(1, ATTR_HEIGHT, &panel_height);
	SetPanelAttribute(message_window_Handle, ATTR_LEFT, panel_left_position+panel_width+100-500); //message window with = 500
	SetPanelAttribute(message_window_Handle, ATTR_TOP,  panel_top_position+panel_height/2-85); //message window height = 190, so subtract 85 to overlap the center of the message window with the center of the program
	
    DisplayPanel(message_window_Handle);
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_write_main_thread
	Arguments		-
	Return value	-
					
	Description		Write a message to the Message window UI (in main thread)
------------------------------------------------------------------------------------------------------------------------*/
void CVICALLBACK message_window_write_main_thread(void *callbackData)
{
	int number_of_message_lines;
	char *message;
	
	message = callbackData;

	if (message_window_Handle == 0) 
		message_window_open();
	SetPanelAttribute(message_window_Handle, ATTR_VISIBLE, 1); //Make sure the UI is visible and on top
	
	GetNumTextBoxLines(message_window_Handle, MESSAGE_TEXTBOX, &number_of_message_lines);
	if (number_of_message_lines != 0) //Begin new line
		SetCtrlVal(message_window_Handle, MESSAGE_TEXTBOX, "\n");
	
	SetCtrlVal(message_window_Handle, MESSAGE_TEXTBOX, message); //Set Message

	// Clean up the message box by deleting the first lines when the maximum number of lines is reached
	if (number_of_message_lines > max_number_of_message_lines)
	{
		DeleteTextBoxLines(message_window_Handle, MESSAGE_TEXTBOX, 0, number_of_message_lines - max_number_of_message_lines);
		ReplaceTextBoxLine(message_window_Handle, MESSAGE_TEXTBOX, 0, "<...older messages have been deleted...>");
	}
	free(message); //clear memory of string
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_write_main_thread_dummy
	Arguments		-
	Return value	-
					
	Description		Dummy function needed to call message_window_write_main_thread as a main thread function
------------------------------------------------------------------------------------------------------------------------*/
//void message_window_write_main_thread_dummy(void)
//{
//}

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_write_main_thread
	Arguments		-
	Return value	-
					
	Description		Write a message to the Message window UI (in main thread)
------------------------------------------------------------------------------------------------------------------------*/
/*void CVICALLBACK message_window_write_main_thread(int poolHandle, int functionID, unsigned int event, int value, void *callbackData)
{
	int number_of_message_lines;
	char *message_with_time;
	char *message;
	
	message = callbackData;
	message_with_time = calloc(strlen(message)+13, sizeof(char)); //set size of string: message + 6 for time + 3 for spaces + 1 for terminating null + 1 for possible newline character for the first line in the text box

	if (message_window_Handle == 0) 
		message_window_open();
	SetPanelAttribute(message_window_Handle, ATTR_VISIBLE, 1); //Make sure the UI is visible and on top
	
	GetNumTextBoxLines(message_window_Handle, MESSAGE_TEXTBOX, &number_of_message_lines);
	
	//Add the current time to the message
	if (number_of_message_lines != 0)
		strcat(message_with_time, "\n");
	strcat(message_with_time, TimeStr());
	strcat(message_with_time, "   ");
	strcat(message_with_time, message);

	SetCtrlVal(message_window_Handle, MESSAGE_TEXTBOX, message_with_time); //Set Message

	// Clean up the message box by deleting the first lines when the maximum number of lines is reached
	if (number_of_message_lines > max_number_of_message_lines)
	{
		DeleteTextBoxLines(message_window_Handle, MESSAGE_TEXTBOX, 0, number_of_message_lines - max_number_of_message_lines);
		ReplaceTextBoxLine(message_window_Handle, MESSAGE_TEXTBOX, 0, "<...older messages have been deleted...>");
	}
	free(message_with_time); //clear memory of string
}*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_write
	Arguments		message		message to display
	Return value	-
					
	Description		Write a message to the Message window UI
					If this function is called from a different thread than the main thread, the panel 
					will not update. Therefore the message_window_write_main_thread() function, which
					will write the message, is always called from the main thread.
------------------------------------------------------------------------------------------------------------------------*/
void message_window_write(const char *message)
{
	/*if ( CmtGetMainThreadID() == CmtGetCurrentThreadID() ) //message_window_write() is called from the main thread
		message_window_write_main_thread(0, 0, 0, 0, message); //Call message_window_write_main_thread() directly
	else //Call message_window_write_main_thread() from the main thread
	{
		int Thread_ID;
		CmtScheduleThreadPoolFunctionAdv(DEFAULT_THREAD_POOL_HANDLE, (ThreadFunctionPtr) message_window_write_main_thread_dummy,
			NULL, THREAD_PRIORITY_NORMAL, (ThreadFunctionCallbackPtr) message_window_write_main_thread, EVENT_TP_THREAD_FUNCTION_BEGIN, message, CmtGetMainThreadID(), &Thread_ID );	
		CmtWaitForThreadPoolFunctionCompletion (DEFAULT_THREAD_POOL_HANDLE, Thread_ID, OPT_TP_PROCESS_EVENTS_WHILE_WAITING);
		CmtReleaseThreadPoolFunctionID(DEFAULT_THREAD_POOL_HANDLE, Thread_ID);
		//Wrong code: the function above waits until message_window_write_main_thread_dummy is finished, not message_window_write_main_thread 
		//then the message could already be freed
		
		PostDeferredCall(DeferredCallbackFunction, message);
		ProcessSystemEvents();
	}*/ 
	
	char *message_with_time;
	
	//Add the current time to the message
	message_with_time = calloc(strlen(message)+13, sizeof(char)); //set size of string: message + 6 for time + 3 for spaces + 1 for terminating null + 1 for possible newline character for the first line in the text box
	strcat(message_with_time, TimeStr());
	strcat(message_with_time, "   ");
	strcat(message_with_time, message);
	
	if ( CmtGetMainThreadID() == CmtGetCurrentThreadID() ) //message_window_write() is called from the main thread
		message_window_write_main_thread(message_with_time); //Call message_window_write_main_thread() directly
	else //Call message_window_write_main_thread() from the main thread
	{
		PostDeferredCall(message_window_write_main_thread, message_with_time);
		ProcessSystemEvents();
	}
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		Button_Clear
					
	Description		Button Clear pressed
					All message lines are cleared
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK Button_Clear (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			DeleteTextBoxLines(message_window_Handle, MESSAGE_TEXTBOX, 0, -1);
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		Button_Hide
					
	Description		Button Hide messages pressed
					Hide the UI
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK Button_Hide (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			SetPanelAttribute(message_window_Handle, ATTR_VISIBLE, 0);
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		Button_Close
					
	Description		Button close pressed
					UI will close
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK Button_Close (int panel, int control, int event,
		void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			DiscardPanel(message_window_Handle);
			message_window_Handle = 0;
			break;
	}
	return 0;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_callback
					
	Description		User closes or resized the GUI
------------------------------------------------------------------------------------------------------------------------*/
int CVICALLBACK message_window_callback (int panel, int event, void *callbackData,
        int eventData1, int eventData2)
{
    if (event == EVENT_CLOSE)
	{
        DiscardPanel(message_window_Handle);
		message_window_Handle = 0;
	}
	else if (event == EVENT_PANEL_SIZE)
	{
		int size;
		GetPanelAttribute(message_window_Handle, ATTR_WIDTH, &size);
		SetCtrlAttribute(message_window_Handle, MESSAGE_TEXTBOX, ATTR_WIDTH, size-20);
		GetPanelAttribute(message_window_Handle, ATTR_HEIGHT, &size);
		SetCtrlAttribute(message_window_Handle, MESSAGE_TEXTBOX, ATTR_HEIGHT, (size/13)*13-50); //Control can only be multiple of 13 pixel in height
		SetCtrlAttribute(message_window_Handle, MESSAGE_BUTTON_CLEAR, ATTR_TOP, size-30);
		SetCtrlAttribute(message_window_Handle, MESSAGE_BUTTON_HIDE,  ATTR_TOP, size-30);
		SetCtrlAttribute(message_window_Handle, MESSAGE_BUTTON_CLOSE, ATTR_TOP, size-30);
	}
    return 0;
}
