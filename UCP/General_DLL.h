/*------------------------------------------------------------------------------------------------------------------------
 *
 *			ScanController
 *
-----------------------------------------------------------------------------------------------------------------------*/
#define ScanController_strlen        100                               

// client types
#define ScanController_clienttype_hw   1	//hardware client
#define ScanController_clienttype_acq  2	//acquisition client
#define ScanController_clienttype_mon  3	//monitor client
#define ScanController_clienttype_stat 4	//for a program to remotely check the status of all clients

// communication structure for adding a variable from a client to the ScanController
typedef struct
{
	char     name[ScanController_strlen+1]; 	// name
	double   initvalue;				// initial value
	double   initmin;				// min value
	double   initmax;				// max value
	char     unit[ScanController_strlen+1];	// unit
	unsigned clienthandle;			// handle of client belonging to the variable
	int      clienttype;			// clienttype_ hw / acq / mon
	int		 vartype;				// int: VAL_INTEGER  / 1, double: VAL_DOUBLE / 4. For commumication between client and ScanController a double variable is used. In the ScanController and client this is typecasted back to an integer
} ScanController_Struct_VarInfo;

// communication structure from the ScanController to a client with info about saving data during an acquisition callback
typedef struct
{
	char     path[2*ScanController_strlen+1];		  //Main path\series path, e.g. S:\Bunches\Data\2010-01-01 Wavelength Scan\20100101_083000_480nm
	char     filenameprefix[ScanController_strlen+1]; //timestamp_seriesname_currentstep, e.g. 20100101_083000_480nm_S001
} ScanController_Struct_AcqInfo;

typedef void (*tvar_cli_cb_start)       (void);                 	//pointer to callback function for "start" event of clients
typedef void (*tvar_cli_cb_startacq)    (ScanController_Struct_AcqInfo); //pointer to callback function for "startacq" event of clients
typedef void (*tvar_cli_cb_updatevar)   (unsigned, char*, double);  //pointer to callback function for "client updatevar" event
typedef void (*tvar_cli_cb_disconnect)  (void);                     //pointer to callback function for "client disconnect" event

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_Connect
	Arguments		clienttype		ScanController_clienttype_hw	1	hardware client
									ScanController_clienttype_acq	2	acquisition client
									ScanController_clienttype_mon   3	monitor client
					clientname		Display name in scancontroller
					cb_start		Callback function for start event, which is called when a variable value is changed. This event can be for
									example updating hardware. Acquisition clients can both have a start AND acquisition callback function
					cb_startacq		If the client is an acquisition client, this callback function is called when data acquisition should begin
					cb_disconnect	Callback function when the ScanController disconnects. A message window automatically pops up 
									in the client to display when this happens
					setup_code		1 for UCP setup, 2 for UCIS setup  
	Return value	1	client connected successfully
					0	error in connecting
					
	Description		Connect client to ScanController
------------------------------------------------------------------------------------------------------------------------*/
char ScanController_Client_Connect(int clienttype, char *clientname, tvar_cli_cb_start cb_start, tvar_cli_cb_startacq cb_startacq, tvar_cli_cb_disconnect cb_disconnect, char setup_code);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_Disconnect
	Arguments		-
	Return value	-
					
	Description		Disconnect client from ScanController
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Client_Disconnect(void);   

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_AddVar
	Arguments		varinfo			min, max, unit, name info of var
					variable		Variable. When the variable value is changed in the ScanController, 
									its new value is automatically updated in the client 
									and the 'updated' char is set to 1 if it is used
					variable_size   size of variable, pass sizeof(variable) to function
					updated			pointer to char that indicates if the variable is updated (optional, if not needed pass NULL)
	Return value	-
					
	Description		Add variable to ScanController
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Client_AddVar(ScanController_Struct_VarInfo varinfo, void *variable, int variable_size, char *updated);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_ClearAllVars
	Arguments		-
	Return value	-
					
	Description		Clear all variables from this client on ScanController
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Client_ClearAllVars(void);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_ClientReady
	Arguments		-
	Return value	-
					
	Description		Signal ScanController that client is ready after updating
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Client_ClientReady(void);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_ClientError
	Arguments		-
	Return value	-
					
	Description		Signal ScanController that client has an error
------------------------------------------------------------------------------------------------------------------------*/
void ScanController_Client_ClientError(void);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ScanController_Client_DelVar
	Arguments		varname		variable name
	Return value	-
					
	Description		Delete variable from ScanController
------------------------------------------------------------------------------------------------------------------------*/
//void  ScanController_Client_DelVar(char *varname); //currently not used



/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Experiment Monitor
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		ExpMonitor_Client_Add_Var
	Arguments		Client_Name		Name of the client. Only alpha-numerical characters are allowed
					Variable_Name	Variable name
			   Out  Var_Handle 		Handle for the variable, used internally. The value of this handle is the number of
			   						the connected variable, i.e. this is 1 for first connected variable, 2 for second...
									Add Var1, handle = 1, add Var2, handle = 2. Del Var1, handle Var2 is changed to 1
					Panel_Handle	PanelHandle of the control belonging to the variable
					Control_ID		Control ID of the control belonging to the variable
					setup_code		1 for UCP setup, 2 for UCIS setup
	Return value	-
					
	Description		Add a variable to the Experiment Monitor
------------------------------------------------------------------------------------------------------------------------*/
void ExpMonitor_Client_Add_Var(char *Client_Name, char *Variable_Name, int *Var_Handle, int Panel_Handle, int Control_ID, char setup_code);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ExpMonitor_Client_Del_Var
	Arguments		Var_Handle	Handle for the variable
	Return value	-
					
	Description		Delete a variable from the Experiment Monitor
------------------------------------------------------------------------------------------------------------------------*/
void ExpMonitor_Client_Del_Var(int VarHandle);

/*------------------------------------------------------------------------------------------------------------------------
	Function		ExpMonitor_Client_Del_All_Vars
	Arguments		-
	Return value	-
					
	Description		Deletes all variables of this client from the Experiment Monitor 
------------------------------------------------------------------------------------------------------------------------*/
void ExpMonitor_Client_Del_All_Vars(void);



/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Message window
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		message_window_write
	Arguments		message		message to display
	Return value	-
					
	Description		Write a message to the Message window UI
------------------------------------------------------------------------------------------------------------------------*/
void message_window_write(const char *message);



/*------------------------------------------------------------------------------------------------------------------------
 *
 *			NLLS Fitting
 *
-----------------------------------------------------------------------------------------------------------------------*/
/*
The Levenberg-Marquardt (LM) algorithm is an iterative technique that
locates the minimum of a function that is expressed as the sum of squares
of nonlinear functions. It has become a standard technique for nonlinear
least-squares problems and can be thought of as a combination of steepest
descent and the Gauss-Newton method.
*/

#define LM_OPTS_SZ    	 5 /* max(4, 5) */
#define LM_INIT_MU    	 1E-03
#define LM_STOP_THRESH	 1E-17

/* 
 * This function seeks the parameter vector p that best describes the measurements vector x.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized.
 *
 * This function requires an analytic Jacobian. In case the latter is unavailable,
 * use LEVMAR_DIF() bellow
 *
 * Returns the number of iterations (>=0) if successful, LM_ERROR if failed
 *
 * For more details, see K. Madsen, H.B. Nielsen and O. Tingleff's lecture notes on 
 * non-linear least squares at http://www.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 *
int LEVMAR_DER(
  void (*func)(LM_REAL *p, LM_REAL *hx, int m, int n, void *adata), // functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n
  void (*jacf)(LM_REAL *p, LM_REAL *j, int m, int n, void *adata),  // function to evaluate the Jacobian \part x / \part p  
  LM_REAL *p,         /* I/O: initial parameter estimates. On output has the estimated solution 
  LM_REAL *x,         /* I: measurement vector. NULL implies a zero vector 
  int m,              /* I: parameter vector dimension (i.e. #unknowns) 
  int n,              /* I: measurement vector dimension 
  int itmax,          /* I: maximum number of iterations 
  LM_REAL opts[4],    /* I: minim. options [\mu, \epsilon1, \epsilon2, \epsilon3]. Respectively the scale factor for initial \mu,
                       * stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2. Set to NULL for defaults to be used
                    
  LM_REAL info[LM_INFO_SZ],
					  /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]= ||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small Dp
                      *                                 3 - stopped by itmax
                      *                                 4 - singular matrix. Restart from current p with increased mu 
                      *                                 5 - no further error reduction is possible. Restart with increased mu
                      *                                 6 - stopped by small ||e||_2
                      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                      * info[7]= # function evaluations
                      * info[8]= # Jacobian evaluations
                      * info[9]= # linear systems solved, i.e. # attempts for reducing error
                      
  LM_REAL *work,     /* working memory at least LM_DER_WORKSZ() reals large, allocated if NULL 
  LM_REAL *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed.
  void *adata)       /* pointer to possibly additional data, passed uninterpreted to func & jacf.
                      * Set to NULL if not needed
                      */

extern int dlevmar_der(
      void (*func)(double *p, double *hx, int m, int n, void *adata),
      void (*jacf)(double *p, double *j, int m, int n, void *adata),
      double *p, double *x, int m, int n, int itmax, double *opts,
      double *info, double *work, double *covar, void *adata);	  

/* Secant version of the LEVMAR_DER() function above: the Jacobian is approximated with 
 * the aid of finite differences (forward or central, see the comment for the opts argument)

int LEVMAR_DIF(
  void (*func)(LM_REAL *p, LM_REAL *hx, int m, int n, void *adata), / functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n 
  LM_REAL *p,          I/O: initial parameter estimates. On output has the estimated solution 
  LM_REAL *x,         /* I: measurement vector. NULL implies a zero vector 
  int m,              /* I: parameter vector dimension (i.e. #unknowns) 
  int n,              /* I: measurement vector dimension 
  int itmax,          /* I: maximum number of iterations 
  LM_REAL opts[5],    /* I: opts[0-4] = minim. options [\mu, \epsilon1, \epsilon2, \epsilon3, \delta]. Respectively the
                       * scale factor for initial \mu, stopping thresholds for ||J^T e||_inf, ||Dp||_2 and ||e||_2 and
                       * the step used in difference approximation to the Jacobian. Set to NULL for defaults to be used.
                       * If \delta<0, the Jacobian is approximated with central differences which are more accurate
                       * (but slower!) compared to the forward differences employed by default. 
                       *
  LM_REAL info[LM_INFO_SZ],
					           /* O: information regarding the minimization. Set to NULL if don't care
                      * info[0]= ||e||_2 at initial p.
                      * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                      * info[5]= # iterations,
                      * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                      *                                 2 - stopped by small Dp
                      *                                 3 - stopped by itmax
                      *                                 4 - singular matrix. Restart from current p with increased mu 
                      *                                 5 - no further error reduction is possible. Restart with increased mu
                      *                                 6 - stopped by small ||e||_2
                      *                                 7 - stopped by invalid (i.e. NaN or Inf) "func" values. This is a user error
                      * info[7]= # function evaluations
                      * info[8]= # Jacobian evaluations
                      * info[9]= # linear systems solved, i.e. # attempts for reducing error
                      *
  LM_REAL *work,     /* working memory at least LM_DIF_WORKSZ() reals large, allocated if NULL
  LM_REAL *covar,    /* O: Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed.
  void *adata)       /* pointer to possibly additional data, passed uninterpreted to func.
                      * Set to NULL if not needed
                      */
/*
extern int dlevmar_dif(
      void (*func)(double *p, double *hx, int m, int n, void *adata),
      double *p, double *x, int m, int n, int itmax, double *opts,
      double *info, double *work, double *covar, void *adata);
*/




/*  standard deviation (standard error) of the best-fit parameter i.
 *  You can calculate 95% confidence intervals from the standard errors.
 *  Hereby we assume that the difference beteen the data and fit is normally distributed.
 *  The integral of a Gaussian distribution between - 1.96 sigma and + 1.96 sigma
 *  gives 0.950004 , so 1.96*standard_error give the 95% confidence intervals
 *
 *  covar is the mxm covariance matrix of the best-fit parameters (see also LEVMAR_COVAR()).
 *  The standard deviation is computed as \sigma_{i} = \sqrt{Covar_{ii}} 
 */
//extern double dlevmar_stddev( double *covar, int m, int i);


/* coefficient of determination. R squared
 * see  http://en.wikipedia.org/wiki/Coefficient_of_determination
 
  void (*func)(LM_REAL *p, LM_REAL *hx, int m, int n, void *adata), // functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n
  LM_REAL *p,         /* I: parameter estimates
  LM_REAL *x,         /* I: measurement vector. NULL implies a zero vector 
  int m,              /* I: parameter vector dimension (i.e. #unknowns) 
  int n,              /* I: measurement vector dimension 
  void *adata)        /* pointer to possibly additional data, passed uninterpreted to func & jacf.
                       * Set to NULL if not needed
                       */
//extern double dlevmar_R2(void (*func)(double *p, double *hx, int m, int n, void *adata), double *p, double *x, int m, int n, void *adata);


/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Matlab Automation Server
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		StartStopMatlabServer
	Arguments		start_matlab			1	Start matlab
											0	Stop matlab
					panelHandle				Panel Handle
					ControlIDScriptRing		Control ID of the Ring to set the script names in
					ControlIDResultsNumeric	Control ID of the Numeric control to show the script results in
	Return value	0	Success
					-1	Error
	
	Description		Start / Stop the Matlab server
					In the panel you need
					- A ring control (initially hidden) where the available matlab scripts are shown
					- A numeric control (initially hidden) to show the results of the script
					  There should be space in the panel to copy the control MATLAB_DISPLAY_RESULTS_ROWS-1
					  times below the current location, and MATLAB_DISPLAY_RESULTS_COLUMS-1 times right 
					  of the current location (+ x times below that)
------------------------------------------------------------------------------------------------------------------------*/
char StartStopMatlabServer(char start_matlab, int panelHandle, int ControlIDScriptRing, int ControlIDResultsNumeric);

/*------------------------------------------------------------------------------------------------------------------------
	Function		MatlabRunScript
	Arguments		arguments	Arguments to pass to script
	Return value	-
					
	Description		Call a script with (optional) arguments and show results.
------------------------------------------------------------------------------------------------------------------------*/
void MatlabRunScript(char *arguments);

/*------------------------------------------------------------------------------------------------------------------------
	Function		SelectMatlabScript
	Arguments		-
	Return value	-
	
	Description		Select the active MATLAB Script in the ring index
------------------------------------------------------------------------------------------------------------------------*/
void SelectMatlabScript(void);


/*------------------------------------------------------------------------------------------------------------------------
 *
 *			Misc General Functions
 *
-----------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------
	Function		check_for_duplicate_program_running
	Arguments		0	Program is not running
					-1	Another instance of this program is already running
	Return value	-
					
	Description		Check if onother copy of the program is already running 
------------------------------------------------------------------------------------------------------------------------*/
char check_for_duplicate_program_running(void);

/*------------------------------------------------------------------------------------------------------------------------
	Function		install_main_callback
	Arguments		-
	Return value	-
					
	Description		Install a Main Callback function
					When the user closes the program by right-clicking -> close in the taskbar, this callback makes
					sure that the functions defined after RunUserInterface() are executed, which is not the case 
					when no Main Callback is installed
------------------------------------------------------------------------------------------------------------------------*/
void install_main_callback(void);

/*------------------------------------------------------------------------------------------------------------------------
	Function		recall_panel_state
	Arguments		-
	Return value	-
					
	Description		Recall a previous panel state, saved in the file 'panelstate.dat'
------------------------------------------------------------------------------------------------------------------------*/
void recall_panel_state(void);
	
/*------------------------------------------------------------------------------------------------------------------------
	Function		save_panel_state
	Arguments		-
	Return value	-
					
	Description		Save the current panel state to the file 'panelstate.dat'
------------------------------------------------------------------------------------------------------------------------*/
void save_panel_state(void);

/*------------------------------------------------------------------------------------------------------------------------
	Function		save_panel_state
	Arguments		path_name			Path name of the directory to make
	Return value	-
					
	Description		Make a directory. The standard MakeDir function can only make one directory at a time
					E.g. if D:\ is empty it cannot make D:\Dir1\Subdir1 directly
------------------------------------------------------------------------------------------------------------------------*/
void make_directory(char *path_name);
