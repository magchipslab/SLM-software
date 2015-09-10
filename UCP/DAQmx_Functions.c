/*------------------------------------------------------------------------------------------------------------------------
 *  File 		DAQmx_Functions.c
 *  Author		Wouter Engelen
 *	Version		September 2009
 *		
 *  Details		Control I/O on NI USB-6259 with DAQmx
 *
------------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
 	Channel Names: 	dev1/ + 
					ao0-3			Coaxial Analog Out 0-3
					ai0-31			Coaxial Analog In 0-7 16-23, rest unknown
					ctr0-1			Counter. Counter Out is connected to P2.4-5 (For other counter connections see Chapter 7-30 DAQmx manual)
					port0/line0-31  Right Panel Digital and Timing IO (In/Out), P0.0-31
					port1/line0-7	Coaxial Digital and Timing IO (In/Out), P1.0-7
					port2/line0-7	Right Panel Digital and Timing IO (In/Out), P2.0-7
					(list not complete)
 ------------------------------------------------------------------------------------------------------------------------*/
#include <utility.h>
#include <ansi_c.h>
#include <NIDAQmx.h>
#include <DAQmxIOctrl.h>

void error_occurred(void); // Announce function for DAQmxErrChk
#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) error_occurred()
#define number_of_ai_samples 10

int error;

/*------------------------------------------------------------------------------------------------------------------------
	Function		create_task_AO
	Arguments		*taskHandle	Pointer to taskhandle
					Channel		AO channel
					min_value	Minimum voltage (V) that can be set
					max_value	Maximum voltage (V) that can be set
	Return value	-
					
	Description		Create a task that can set a voltage on an AO (Analog Output)
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_create_task_AO(TaskHandle *taskHandle, char *channel, double min_value, double max_value)
{
	DAQmxErrChk (DAQmxCreateTask("", taskHandle));
	DAQmxErrChk (DAQmxCreateAOVoltageChan(*taskHandle, channel, "", min_value, max_value, DAQmx_Val_Volts, ""));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		run_task_AO
	Arguments		taskHandle	taskhandle
					voltage		Voltage (V)
	Return value	-
					
	Description		Set a voltage on an AO (Analog Output)
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_run_task_AO(TaskHandle taskHandle, float64 voltage)
{
	DAQmxErrChk (DAQmxStartTask(taskHandle));
	DAQmxErrChk (DAQmxWriteAnalogF64(taskHandle, 1, 1, 10.0, DAQmx_Val_GroupByChannel, &voltage, NULL, NULL));
	DAQmxErrChk (DAQmxStopTask(taskHandle));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		create_task_AI
	Arguments		*taskHandle	Pointer to taskhandle
					Channel		AI channel
					min_value	Minimum voltage (V) that can be get
					max_value	Maximum voltage (V) that can be get
	Return value	-
					
	Description		Create a task that can get the voltage on an AI (Analog Input)
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_create_task_AI(TaskHandle *taskHandle, char *channel, double min_value, double max_value)
{
	DAQmxErrChk (DAQmxCreateTask("", taskHandle));
	DAQmxErrChk (DAQmxCreateAIVoltageChan(*taskHandle, channel, "", DAQmx_Val_Cfg_Default, min_value, max_value, DAQmx_Val_Volts, NULL));
	DAQmxErrChk (DAQmxCfgSampClkTiming(*taskHandle, "", 1000, DAQmx_Val_Rising, DAQmx_Val_FiniteSamps, number_of_ai_samples));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		run_task_AI
	Arguments		taskHandle	taskhandle
	Return value	Voltage (V)
					
	Description		Get the voltage on an AI (Analog Input)
------------------------------------------------------------------------------------------------------------------------*/
double DAQmx_run_task_AI(TaskHandle taskHandle)
{
	double sum = 0;
	float64 *data = NULL;
	int i;
	if( (data=malloc(number_of_ai_samples*sizeof(float64))) == NULL ) MessagePopup("Error","Not enough memory");
	
	DAQmxErrChk (DAQmxStartTask(taskHandle));
	DAQmxErrChk (DAQmxReadAnalogF64(taskHandle,number_of_ai_samples,10.0,DAQmx_Val_GroupByChannel,data,number_of_ai_samples,NULL,NULL));
	DAQmxErrChk (DAQmxStopTask(taskHandle));

    for (i = 0; i < number_of_ai_samples; ++i) sum += data[i];
	free(data);
	return sum/number_of_ai_samples;
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		create_task_DO
	Arguments		*taskHandle	Pointer to taskhandle
					Channel		DO channel
	Return value	-
					
	Description		Create a task that can set 0 / 5 V on an DO (Digital Output)
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_create_task_DO(TaskHandle *taskHandle, char *channel)
{
	DAQmxErrChk (DAQmxCreateTask("",taskHandle));
	DAQmxErrChk (DAQmxCreateDOChan(*taskHandle, channel, "", DAQmx_Val_ChanPerLine));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		run_task_DO
	Arguments		taskHandle	taskhandle
					channel_on	0	Channel off, 0 V
								1	Channel on,  5 V
	Return value	-
					
	Description		Set 0 / 5 V on a DO (Digital Output)
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_run_task_DO(TaskHandle taskHandle, uInt8 channel_on)
{
	DAQmxErrChk (DAQmxStartTask(taskHandle));
	DAQmxErrChk (DAQmxWriteDigitalLines(taskHandle, 1, FALSE, 1.0, DAQmx_Val_GroupByChannel, &channel_on, NULL, NULL));
	DAQmxErrChk (DAQmxStopTask(taskHandle));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		create_task_Pulse_Train
	Arguments		*taskHandle	Pointer to taskhandle
					Channel		Channel
					frequency   The frequency at which to generate pulses 
	Return value	-
					
	Description		Create a task that can generate digital pulses
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_create_task_Pulse_Train(TaskHandle *taskHandle, char *channel, float64 frequency)
{
	float64 Initialdelay = 0.0; //The amount of time in seconds to wait before generating the first pulse
	float64 dutycycle = 0.5;	//The width of the pulse divided by the pulse period, i.e. percentage of high time with respect to duration of one pulse
	
	DAQmxErrChk (DAQmxCreateTask("",taskHandle)); 
	DAQmxErrChk (DAQmxCreateCOPulseChanFreq(*taskHandle,channel,"", DAQmx_Val_Hz, DAQmx_Val_Low, Initialdelay, frequency, dutycycle));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		run_task_Pulse_Train
	Arguments		taskHandle			Taskhandle
					number_of_pulses	Number of pulses in the pulse train
	Return value	-
					
	Description		Make a pulse train
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_run_task_Pulse_Train(TaskHandle taskHandle, uInt64 number_of_pulses)
{
	DAQmxErrChk (DAQmxCfgImplicitTiming(taskHandle, DAQmx_Val_FiniteSamps, number_of_pulses));
	DAQmxErrChk (DAQmxStartTask(taskHandle));
	DAQmxErrChk (DAQmxWaitUntilTaskDone(taskHandle,-1));
	DAQmxErrChk (DAQmxStopTask(taskHandle));
}

 /*------------------------------------------------------------------------------------------------------------------------
	Function		DAQmx_create_task_trigger
	Arguments		*taskHandle	Pointer to taskhandle
					Channel		DO channel
	Return value	-
					
	Description		Create a task that gives a TTL trigger pulse on a DO (Digital Output) channel
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_create_task_trigger(TaskHandle *taskHandle, char *channel)
{
	DAQmxErrChk (DAQmxCreateTask("",taskHandle));
	DAQmxErrChk (DAQmxCreateDOChan(*taskHandle, channel, "", DAQmx_Val_ChanPerLine));
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		DAQmx_run_task_trigger
	Arguments		taskHandle			taskhandle
					trigger_duration	Duration of trigger. Minimal 4 ms. 
	Return value	-
					
	Description		Give a TTL trigger pulse. Output pulse duration jitter +- 1 ms
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_run_task_trigger(TaskHandle taskHandle, double trigger_duration)
{
	uInt8 channel_state;
		
	DAQmxErrChk (DAQmxStartTask(taskHandle));
	channel_state = 1;
	DAQmxErrChk (DAQmxWriteDigitalLines(taskHandle, 1, FALSE, 1.0, DAQmx_Val_GroupByChannel, &channel_state, NULL, NULL)); //Set channel to high output
	Delay(trigger_duration-2e-3); //Delay duration of trigger, minus 2 ms offset
	channel_state = 0;
	DAQmxErrChk (DAQmxWriteDigitalLines(taskHandle, 1, FALSE, 1.0, DAQmx_Val_GroupByChannel, &channel_state, NULL, NULL)); //Set channel to low output  
	DAQmxErrChk (DAQmxStopTask(taskHandle));
}
	
/*------------------------------------------------------------------------------------------------------------------------
	Function		delete_task
	Arguments		taskHandle	taskhandle
	Return value	-
					
	Description		Delete a DAQmx task
------------------------------------------------------------------------------------------------------------------------*/
void DAQmx_delete_task(TaskHandle taskHandle)
{
	DAQmxErrChk (DAQmxClearTask(taskHandle));	
}

/*------------------------------------------------------------------------------------------------------------------------
	Function		error_occurred
	Arguments		-
	Return value	-
					
	Description		Displays information about an occurred error
------------------------------------------------------------------------------------------------------------------------*/
void error_occurred(void)
{
	char errBuff[2048]={'\0'};
	DAQmxGetExtendedErrorInfo(errBuff,2048);
	MessagePopup("DAQmx Error",errBuff);
}
