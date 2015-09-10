typedef void* TaskHandle;
typedef double float64;
typedef unsigned char uInt8;
typedef unsigned __int64 uInt64;

void DAQmx_create_task_AO(TaskHandle *taskHandle, char *channel, double min_value, double max_value);
void DAQmx_run_task_AO(TaskHandle taskHandle, float64 voltage);
void DAQmx_create_task_AI(TaskHandle *taskHandle, char *channel, double min_value, double max_value);
double DAQmx_run_task_AI(TaskHandle taskHandle);
void DAQmx_create_task_DO(TaskHandle *taskHandle, char *channel);
void DAQmx_run_task_DO(TaskHandle taskHandle, uInt8 channel_on);
void DAQmx_create_task_Pulse_Train(TaskHandle *taskHandle, char *channel, float64 frequency);
void DAQmx_run_task_Pulse_Train(TaskHandle taskHandle, uInt64 number_of_pulses);
void DAQmx_create_task_trigger(TaskHandle *taskHandle, char *channel);
void DAQmx_run_task_trigger(TaskHandle taskHandle, double trigger_duration);
void DAQmx_delete_task(TaskHandle taskHandle);
