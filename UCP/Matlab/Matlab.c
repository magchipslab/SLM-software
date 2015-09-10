#include "Matlab.h"

typedef interface tagMLApp_IMLApp_Interface MLApp_IMLApp_Interface;

typedef struct tagMLApp_IMLApp_VTable
{
	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( MLApp_IMLApp_Interface __RPC_FAR * This, 
	                                                         REFIID riid, 
	                                                         void __RPC_FAR *__RPC_FAR *ppvObject);

	ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( MLApp_IMLApp_Interface __RPC_FAR * This);

	ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( MLApp_IMLApp_Interface __RPC_FAR * This);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetFullMatrix_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                        BSTR name, 
	                                                        BSTR workspace, 
	                                                        SAFEARRAY **pr, 
	                                                        SAFEARRAY **pi);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PutFullMatrix_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                        BSTR name, 
	                                                        BSTR workspace, 
	                                                        SAFEARRAY *pr, 
	                                                        SAFEARRAY *pi);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Execute_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                  BSTR name);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *MinimizeCommandWindow_) (MLApp_IMLApp_Interface __RPC_FAR *This);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *MaximizeCommandWindow_) (MLApp_IMLApp_Interface __RPC_FAR *This);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Quit_) (MLApp_IMLApp_Interface __RPC_FAR *This);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetCharArray_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                       BSTR name, 
	                                                       BSTR workspace, 
	                                                       BSTR *mlString);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PutCharArray_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                       BSTR name, 
	                                                       BSTR workspace, 
	                                                       BSTR charArray);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetVisible_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                     long *pVal);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *SetVisible_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                     long pVal);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetWorkspaceData_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                           BSTR name, 
	                                                           BSTR workspace, 
	                                                           VARIANT *pdata);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PutWorkspaceData_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                           BSTR name, 
	                                                           BSTR workspace, 
	                                                           VARIANT data);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Feval_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                BSTR bstrName, 
	                                                long nargout, 
	                                                VARIANT *pvarArgOut, 
	                                                VARIANT arg1, 
	                                                VARIANT arg2, 
	                                                VARIANT arg3, 
	                                                VARIANT arg4, 
	                                                VARIANT arg5, 
	                                                VARIANT arg6, 
	                                                VARIANT arg7, 
	                                                VARIANT arg8, 
	                                                VARIANT arg9, 
	                                                VARIANT arg10, 
	                                                VARIANT arg11, 
	                                                VARIANT arg12, 
	                                                VARIANT arg13, 
	                                                VARIANT arg14, 
	                                                VARIANT arg15, 
	                                                VARIANT arg16, 
	                                                VARIANT arg17, 
	                                                VARIANT arg18, 
	                                                VARIANT arg19, 
	                                                VARIANT arg20, 
	                                                VARIANT arg21, 
	                                                VARIANT arg22, 
	                                                VARIANT arg23, 
	                                                VARIANT arg24, 
	                                                VARIANT arg25, 
	                                                VARIANT arg26, 
	                                                VARIANT arg27, 
	                                                VARIANT arg28, 
	                                                VARIANT arg29, 
	                                                VARIANT arg30, 
	                                                VARIANT arg31, 
	                                                VARIANT arg32);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetVariable_) (MLApp_IMLApp_Interface __RPC_FAR *This, 
	                                                      BSTR name, 
	                                                      BSTR workspace, 
	                                                      VARIANT *pdata);

} MLApp_IMLApp_VTable;

typedef interface tagMLApp_IMLApp_Interface
{
	CONST_VTBL MLApp_IMLApp_VTable __RPC_FAR *lpVtbl;
} MLApp_IMLApp_Interface;

typedef interface tagMLApp_IMLEval_Interface MLApp_IMLEval_Interface;

typedef struct tagMLApp_IMLEval_VTable
{
	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( MLApp_IMLEval_Interface __RPC_FAR * This, 
	                                                         REFIID riid, 
	                                                         void __RPC_FAR *__RPC_FAR *ppvObject);

	ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( MLApp_IMLEval_Interface __RPC_FAR * This);

	ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( MLApp_IMLEval_Interface __RPC_FAR * This);

	HRESULT ( STDMETHODCALLTYPE __RPC_FAR *XLEval_) (MLApp_IMLEval_Interface __RPC_FAR *This, 
	                                                 BSTR bstrName, 
	                                                 long nargout, 
	                                                 VARIANT *pvarArgOut, 
	                                                 long nargin, 
	                                                 VARIANT varArgIn);

} MLApp_IMLEval_VTable;

typedef interface tagMLApp_IMLEval_Interface
{
	CONST_VTBL MLApp_IMLEval_VTable __RPC_FAR *lpVtbl;
} MLApp_IMLEval_Interface;

const IID MLApp_IID_DIMLApp =
	{
		0x669CEC93, 0x6E22, 0x11CF, 0xA4, 0xD6, 0x0, 0xA0, 0x24, 0x58, 0x3C, 0x19
	};

const IID MLApp_IID_IMLApp =
	{
		0x669CEC92, 0x6E22, 0x11CF, 0xA4, 0xD6, 0x0, 0xA0, 0x24, 0x58, 0x3C, 0x19
	};

const IID MLApp_IID_IMLEval =
	{
		0x27BEA9CE, 0xA20C, 0x430F, 0x87, 0xC2, 0x1C, 0xC8, 0xBA, 0x31, 0xC3, 0xA8
	};

const IID MLApp_IID_DIMLEval =
	{
		0x6E813920, 0x23FE, 0x4D6D, 0x91, 0xF8, 0x56, 0xFA, 0xB0, 0x6C, 0x5D, 0x13
	};

HRESULT CVIFUNC MLApp_NewDIMLApp (const char *server, int supportMultithreading,
                                  LCID locale, int reserved,
                                  CAObjHandle *objectHandle)
{
	HRESULT __result = S_OK;
	GUID clsid = {0x3A579EB9, 0x8A43, 0x4022, 0xA1, 0xAF, 0x4C, 0x13, 0xDD,
	              0xEC, 0x40, 0x79};

	__result = CA_CreateObjectByClassIdEx (&clsid, server, &IID_IDispatch,
	                                       supportMultithreading, locale,
	                                       reserved, objectHandle);

	return __result;
}

HRESULT CVIFUNC MLApp_OpenDIMLApp (const char *fileName, const char *server,
                                   int supportMultithreading, LCID locale,
                                   int reserved, CAObjHandle *objectHandle)
{
	HRESULT __result = S_OK;
	GUID clsid = {0x3A579EB9, 0x8A43, 0x4022, 0xA1, 0xAF, 0x4C, 0x13, 0xDD,
	              0xEC, 0x40, 0x79};

	__result = CA_LoadObjectFromFileByClassIdEx (fileName, &clsid, server,
	                                             &IID_IDispatch,
	                                             supportMultithreading, locale,
	                                             reserved, objectHandle);

	return __result;
}

HRESULT CVIFUNC MLApp_ActiveDIMLApp (const char *server,
                                     int supportMultithreading, LCID locale,
                                     int reserved, CAObjHandle *objectHandle)
{
	HRESULT __result = S_OK;
	GUID clsid = {0x3A579EB9, 0x8A43, 0x4022, 0xA1, 0xAF, 0x4C, 0x13, 0xDD,
	              0xEC, 0x40, 0x79};

	__result = CA_GetActiveObjectByClassIdEx (&clsid, server,
	                                          &IID_IDispatch,
	                                          supportMultithreading, locale,
	                                          reserved, objectHandle);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppGetFullMatrix (CAObjHandle objectHandle,
                                            ERRORINFO *errorInfo,
                                            const char *name,
                                            const char *workspace,
                                            SAFEARRAY **pr, SAFEARRAY **pi)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING,
	                               CAVT_DOUBLE | CAVT_ARRAY | CAVT_BYREFIO,
	                               CAVT_DOUBLE | CAVT_ARRAY | CAVT_BYREFIO};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010000, CAVT_EMPTY, NULL, 4,
	                              __paramTypes, name, workspace, pr, pi);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppPutFullMatrix (CAObjHandle objectHandle,
                                            ERRORINFO *errorInfo,
                                            const char *name,
                                            const char *workspace, SAFEARRAY *pr,
                                            SAFEARRAY *pi)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING,
	                               CAVT_DOUBLE | CAVT_ARRAY,
	                               CAVT_DOUBLE | CAVT_ARRAY};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010001, CAVT_EMPTY, NULL, 4,
	                              __paramTypes, name, workspace, pr, pi);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppExecute (CAObjHandle objectHandle,
                                      ERRORINFO *errorInfo, const char *name,
                                      char **returnValue)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010002, CAVT_CSTRING, returnValue, 1,
	                              __paramTypes, name);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppMinimizeCommandWindow (CAObjHandle objectHandle,
                                                    ERRORINFO *errorInfo)
{
	HRESULT __result = S_OK;

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010003, CAVT_EMPTY, NULL, 0, NULL);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppMaximizeCommandWindow (CAObjHandle objectHandle,
                                                    ERRORINFO *errorInfo)
{
	HRESULT __result = S_OK;

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010004, CAVT_EMPTY, NULL, 0, NULL);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppQuit (CAObjHandle objectHandle,
                                   ERRORINFO *errorInfo)
{
	HRESULT __result = S_OK;

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010005, CAVT_EMPTY, NULL, 0, NULL);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppGetCharArray (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace,
                                           char **returnValue)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010006, CAVT_CSTRING, returnValue, 2,
	                              __paramTypes, name, workspace);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppPutCharArray (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace,
                                           const char *charArray)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING, CAVT_CSTRING};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010007, CAVT_EMPTY, NULL, 3,
	                              __paramTypes, name, workspace, charArray);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppGetVisible (CAObjHandle objectHandle,
                                         ERRORINFO *errorInfo, long *returnValue)
{
	HRESULT __result = S_OK;

	__result = CA_InvokeHelperEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010008, DISPATCH_PROPERTYGET, CAVT_LONG,
	                              returnValue, 0, NULL);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppSetVisible (CAObjHandle objectHandle,
                                         ERRORINFO *errorInfo, long newValue)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_LONG};

	__result = CA_InvokeHelperEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x60010008, DISPATCH_PROPERTYPUT, CAVT_EMPTY,
	                              NULL, 1, __paramTypes, newValue);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppGetWorkspaceData (CAObjHandle objectHandle,
                                               ERRORINFO *errorInfo,
                                               const char *name,
                                               const char *workspace,
                                               VARIANT *pdata)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING,
	                               CAVT_VARIANT | CAVT_BYREFO};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x6001000A, CAVT_EMPTY, NULL, 3,
	                              __paramTypes, name, workspace, pdata);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppPutWorkspaceData (CAObjHandle objectHandle,
                                               ERRORINFO *errorInfo,
                                               const char *name,
                                               const char *workspace,
                                               VARIANT data)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING, CAVT_VARIANT};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x6001000B, CAVT_EMPTY, NULL, 3,
	                              __paramTypes, name, workspace, data);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppFeval (CAObjHandle objectHandle,
                                    ERRORINFO *errorInfo, const char *bstrName,
                                    long nargout, VARIANT *pvarArgOut,
                                    VARIANT arg1, VARIANT arg2, VARIANT arg3,
                                    VARIANT arg4, VARIANT arg5, VARIANT arg6,
                                    VARIANT arg7, VARIANT arg8, VARIANT arg9,
                                    VARIANT arg10, VARIANT arg11, VARIANT arg12,
                                    VARIANT arg13, VARIANT arg14, VARIANT arg15,
                                    VARIANT arg16, VARIANT arg17, VARIANT arg18,
                                    VARIANT arg19, VARIANT arg20, VARIANT arg21,
                                    VARIANT arg22, VARIANT arg23, VARIANT arg24,
                                    VARIANT arg25, VARIANT arg26, VARIANT arg27,
                                    VARIANT arg28, VARIANT arg29, VARIANT arg30,
                                    VARIANT arg31, VARIANT arg32)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_LONG,
	                               CAVT_VARIANT | CAVT_BYREFO, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT, CAVT_VARIANT, CAVT_VARIANT,
	                               CAVT_VARIANT};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x6001000C, CAVT_EMPTY, NULL, 35,
	                              __paramTypes, bstrName, nargout, pvarArgOut,
	                              arg1, arg2, arg3, arg4, arg5, arg6, arg7,
	                              arg8, arg9, arg10, arg11, arg12, arg13,
	                              arg14, arg15, arg16, arg17, arg18, arg19,
	                              arg20, arg21, arg22, arg23, arg24, arg25,
	                              arg26, arg27, arg28, arg29, arg30, arg31,
	                              arg32);

	return __result;
}

HRESULT CVIFUNC MLApp_DIMLAppGetVariable (CAObjHandle objectHandle,
                                          ERRORINFO *errorInfo, const char *name,
                                          const char *workspace,
                                          VARIANT *returnValue)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_CSTRING};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &IID_IDispatch,
	                              0x6001000D, CAVT_VARIANT, returnValue, 2,
	                              __paramTypes, name, workspace);

	return __result;
}

HRESULT CVIFUNC MLApp_NewIMLApp (const char *server, int supportMultithreading,
                                 LCID locale, int reserved,
                                 CAObjHandle *objectHandle)
{
	HRESULT __result = S_OK;
	GUID clsid = {0x3A579EB9, 0x8A43, 0x4022, 0xA1, 0xAF, 0x4C, 0x13, 0xDD,
	              0xEC, 0x40, 0x79};

	__result = CA_CreateObjectByClassIdEx (&clsid, server, &MLApp_IID_IMLApp,
	                                       supportMultithreading, locale,
	                                       reserved, objectHandle);

	return __result;
}

HRESULT CVIFUNC MLApp_OpenIMLApp (const char *fileName, const char *server,
                                  int supportMultithreading, LCID locale,
                                  int reserved, CAObjHandle *objectHandle)
{
	HRESULT __result = S_OK;
	GUID clsid = {0x3A579EB9, 0x8A43, 0x4022, 0xA1, 0xAF, 0x4C, 0x13, 0xDD,
	              0xEC, 0x40, 0x79};

	__result = CA_LoadObjectFromFileByClassIdEx (fileName, &clsid, server,
	                                             &MLApp_IID_IMLApp,
	                                             supportMultithreading, locale,
	                                             reserved, objectHandle);

	return __result;
}

HRESULT CVIFUNC MLApp_ActiveIMLApp (const char *server,
                                    int supportMultithreading, LCID locale,
                                    int reserved, CAObjHandle *objectHandle)
{
	HRESULT __result = S_OK;
	GUID clsid = {0x3A579EB9, 0x8A43, 0x4022, 0xA1, 0xAF, 0x4C, 0x13, 0xDD,
	              0xEC, 0x40, 0x79};

	__result = CA_GetActiveObjectByClassIdEx (&clsid, server,
	                                          &MLApp_IID_IMLApp,
	                                          supportMultithreading, locale,
	                                          reserved, objectHandle);

	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppGetFullMatrix (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace, SAFEARRAY **pr,
                                           SAFEARRAY **pi)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->GetFullMatrix_ (__vtblIFacePtr,
	                                                    name__AutoType,
	                                                    workspace__AutoType,
	                                                    pr, pi));

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	if (FAILED(__result))
		{
		if (pr)
			{
			CA_SafeArrayDestroy (*pr);
			*pr = 0;
			}
		if (pi)
			{
			CA_SafeArrayDestroy (*pi);
			*pi = 0;
			}
		}
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppPutFullMatrix (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace, SAFEARRAY *pr,
                                           SAFEARRAY *pi)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->PutFullMatrix_ (__vtblIFacePtr,
	                                                    name__AutoType,
	                                                    workspace__AutoType,
	                                                    pr, pi));

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppExecute (CAObjHandle objectHandle,
                                     ERRORINFO *errorInfo, const char *name)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->Execute_ (__vtblIFacePtr,
	                                              name__AutoType));

Error:
	CA_FreeBSTR (name__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppMinimizeCommandWindow (CAObjHandle objectHandle,
                                                   ERRORINFO *errorInfo)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->MinimizeCommandWindow_ (__vtblIFacePtr));

Error:
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppMaximizeCommandWindow (CAObjHandle objectHandle,
                                                   ERRORINFO *errorInfo)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->MaximizeCommandWindow_ (__vtblIFacePtr));

Error:
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppQuit (CAObjHandle objectHandle, ERRORINFO *errorInfo)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->Quit_ (__vtblIFacePtr));

Error:
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppGetCharArray (CAObjHandle objectHandle,
                                          ERRORINFO *errorInfo, const char *name,
                                          const char *workspace, char **mlString)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;
	BSTR mlString__AutoType = 0;

	if (mlString)
		*mlString = 0;

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->GetCharArray_ (__vtblIFacePtr,
	                                                   name__AutoType,
	                                                   workspace__AutoType,
	                                                   &mlString__AutoType));

	if (mlString)
		__caErrChk (CA_BSTRGetCString (mlString__AutoType, mlString));

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	CA_FreeBSTR (mlString__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	if (FAILED(__result))
		{
		if (mlString)
			{
			CA_FreeMemory (*mlString);
			*mlString = 0;
			}
		}
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppPutCharArray (CAObjHandle objectHandle,
                                          ERRORINFO *errorInfo, const char *name,
                                          const char *workspace,
                                          const char *charArray)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;
	BSTR charArray__AutoType = 0;

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));
	__caErrChk (CA_CStringToBSTR (charArray, &charArray__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->PutCharArray_ (__vtblIFacePtr,
	                                                   name__AutoType,
	                                                   workspace__AutoType,
	                                                   charArray__AutoType));

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	CA_FreeBSTR (charArray__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppGetVisible (CAObjHandle objectHandle,
                                        ERRORINFO *errorInfo, long *pVal)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	long pVal__Temp;

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->GetVisible_ (__vtblIFacePtr,
	                                                 &pVal__Temp));

	if (pVal)
		{
		*pVal = pVal__Temp;
		}

Error:
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppSetVisible (CAObjHandle objectHandle,
                                        ERRORINFO *errorInfo, long pVal)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->SetVisible_ (__vtblIFacePtr, pVal));

Error:
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppGetWorkspaceData (CAObjHandle objectHandle,
                                              ERRORINFO *errorInfo,
                                              const char *name,
                                              const char *workspace,
                                              VARIANT *pdata)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;
	VARIANT pdata__Temp;

	if (pdata)
		CA_VariantSetEmpty (pdata);
	CA_VariantSetEmpty (&pdata__Temp);

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->GetWorkspaceData_ (__vtblIFacePtr,
	                                                       name__AutoType,
	                                                       workspace__AutoType,
	                                                       &pdata__Temp));

	if (pdata)
		{
		*pdata = pdata__Temp;
		CA_VariantSetEmpty (&pdata__Temp);
		}

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	CA_VariantClear (&pdata__Temp);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	if (FAILED(__result))
		{
		if (pdata)
			CA_VariantClear (pdata);
		}
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppPutWorkspaceData (CAObjHandle objectHandle,
                                              ERRORINFO *errorInfo,
                                              const char *name,
                                              const char *workspace,
                                              VARIANT data)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->PutWorkspaceData_ (__vtblIFacePtr,
	                                                       name__AutoType,
	                                                       workspace__AutoType,
	                                                       data));

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppFeval (CAObjHandle objectHandle,
                                   ERRORINFO *errorInfo, const char *bstrName,
                                   long nargout, VARIANT *pvarArgOut,
                                   VARIANT arg1, VARIANT arg2, VARIANT arg3,
                                   VARIANT arg4, VARIANT arg5, VARIANT arg6,
                                   VARIANT arg7, VARIANT arg8, VARIANT arg9,
                                   VARIANT arg10, VARIANT arg11, VARIANT arg12,
                                   VARIANT arg13, VARIANT arg14, VARIANT arg15,
                                   VARIANT arg16, VARIANT arg17, VARIANT arg18,
                                   VARIANT arg19, VARIANT arg20, VARIANT arg21,
                                   VARIANT arg22, VARIANT arg23, VARIANT arg24,
                                   VARIANT arg25, VARIANT arg26, VARIANT arg27,
                                   VARIANT arg28, VARIANT arg29, VARIANT arg30,
                                   VARIANT arg31, VARIANT arg32)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR bstrName__AutoType = 0;
	VARIANT pvarArgOut__Temp;

	if (pvarArgOut)
		CA_VariantSetEmpty (pvarArgOut);
	CA_VariantSetEmpty (&pvarArgOut__Temp);

	__caErrChk (CA_CStringToBSTR (bstrName, &bstrName__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->Feval_ (__vtblIFacePtr,
	                                            bstrName__AutoType, nargout,
	                                            &pvarArgOut__Temp, arg1, arg2,
	                                            arg3, arg4, arg5, arg6, arg7,
	                                            arg8, arg9, arg10, arg11,
	                                            arg12, arg13, arg14, arg15,
	                                            arg16, arg17, arg18, arg19,
	                                            arg20, arg21, arg22, arg23,
	                                            arg24, arg25, arg26, arg27,
	                                            arg28, arg29, arg30, arg31,
	                                            arg32));

	if (pvarArgOut)
		{
		*pvarArgOut = pvarArgOut__Temp;
		CA_VariantSetEmpty (&pvarArgOut__Temp);
		}

Error:
	CA_FreeBSTR (bstrName__AutoType);
	CA_VariantClear (&pvarArgOut__Temp);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	if (FAILED(__result))
		{
		if (pvarArgOut)
			CA_VariantClear (pvarArgOut);
		}
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLAppGetVariable (CAObjHandle objectHandle,
                                         ERRORINFO *errorInfo, const char *name,
                                         const char *workspace, VARIANT *pdata)
{
	HRESULT __result = S_OK;
	MLApp_IMLApp_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR name__AutoType = 0;
	BSTR workspace__AutoType = 0;
	VARIANT pdata__Temp;

	if (pdata)
		CA_VariantSetEmpty (pdata);
	CA_VariantSetEmpty (&pdata__Temp);

	__caErrChk (CA_CStringToBSTR (name, &name__AutoType));
	__caErrChk (CA_CStringToBSTR (workspace, &workspace__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLApp,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->GetVariable_ (__vtblIFacePtr,
	                                                  name__AutoType,
	                                                  workspace__AutoType,
	                                                  &pdata__Temp));

	if (pdata)
		{
		*pdata = pdata__Temp;
		CA_VariantSetEmpty (&pdata__Temp);
		}

Error:
	CA_FreeBSTR (name__AutoType);
	CA_FreeBSTR (workspace__AutoType);
	CA_VariantClear (&pdata__Temp);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	if (FAILED(__result))
		{
		if (pdata)
			CA_VariantClear (pdata);
		}
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLApp, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_IMLEvalXLEval (CAObjHandle objectHandle,
                                     ERRORINFO *errorInfo, const char *bstrName,
                                     long nargout, VARIANT *pvarArgOut,
                                     long nargin, VARIANT varArgIn)
{
	HRESULT __result = S_OK;
	MLApp_IMLEval_Interface * __vtblIFacePtr = 0;
	int __didAddRef;
	int __errorInfoPresent = 0;
	BSTR bstrName__AutoType = 0;

	__caErrChk (CA_CStringToBSTR (bstrName, &bstrName__AutoType));

	__caErrChk (CA_GetInterfaceFromObjHandle (objectHandle, &MLApp_IID_IMLEval,
	                                          0, &__vtblIFacePtr, &__didAddRef));
	__caErrChk (__vtblIFacePtr->lpVtbl->XLEval_ (__vtblIFacePtr,
	                                             bstrName__AutoType, nargout,
	                                             pvarArgOut, nargin, varArgIn));

Error:
	CA_FreeBSTR (bstrName__AutoType);
	if (__vtblIFacePtr && __didAddRef)
		__vtblIFacePtr->lpVtbl->Release (__vtblIFacePtr);
	if (FAILED(__result))
		{
		if (pvarArgOut)
			CA_VariantClear (pvarArgOut);
		}
	CA_FillErrorInfoEx (objectHandle, &MLApp_IID_IMLEval, __result, errorInfo,
	                    &__errorInfoPresent);
	if (__errorInfoPresent)
		__result = DISP_E_EXCEPTION;
	return __result;
}

HRESULT CVIFUNC MLApp_DIMLEvalXLEval (CAObjHandle objectHandle,
                                      ERRORINFO *errorInfo, const char *bstrName,
                                      long nargout, VARIANT *pvarArgOut,
                                      long nargin, VARIANT varArgIn)
{
	HRESULT __result = S_OK;
	unsigned int __paramTypes[] = {CAVT_CSTRING, CAVT_LONG,
	                               CAVT_VARIANT | CAVT_BYREFIO, CAVT_LONG,
	                               CAVT_VARIANT};

	__result = CA_MethodInvokeEx (objectHandle, errorInfo, &MLApp_IID_DIMLEval,
	                              0x60010000, CAVT_EMPTY, NULL, 5,
	                              __paramTypes, bstrName, nargout, pvarArgOut,
	                              nargin, varArgIn);

	return __result;
}
