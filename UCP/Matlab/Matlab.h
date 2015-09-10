#ifndef _MLAPP_H
#define _MLAPP_H

#if defined(INCLUDE_AFTER_WINDOWS_H) && !defined(_INC_WINDOWS)
#error  This header must be included before utility.h and formatio.h
#error  because it includes cviauto.h which includes Windows SDK headers.
#endif /* INCLUDE_AFTER_WINDOWS_H */

#include <cviauto.h>

#ifdef __cplusplus
    extern "C" {
#endif
extern const IID MLApp_IID_DIMLApp;

extern const IID MLApp_IID_IMLApp;
extern const IID MLApp_IID_IMLEval;
extern const IID MLApp_IID_DIMLEval;

HRESULT CVIFUNC MLApp_NewDIMLApp (const char *server, int supportMultithreading,
                                  LCID locale, int reserved,
                                  CAObjHandle *objectHandle);

HRESULT CVIFUNC MLApp_OpenDIMLApp (const char *fileName, const char *server,
                                   int supportMultithreading, LCID locale,
                                   int reserved, CAObjHandle *objectHandle);

HRESULT CVIFUNC MLApp_ActiveDIMLApp (const char *server,
                                     int supportMultithreading, LCID locale,
                                     int reserved, CAObjHandle *objectHandle);

HRESULT CVIFUNC MLApp_DIMLAppGetFullMatrix (CAObjHandle objectHandle,
                                            ERRORINFO *errorInfo,
                                            const char *name,
                                            const char *workspace,
                                            SAFEARRAY **pr, SAFEARRAY **pi);

HRESULT CVIFUNC MLApp_DIMLAppPutFullMatrix (CAObjHandle objectHandle,
                                            ERRORINFO *errorInfo,
                                            const char *name,
                                            const char *workspace, SAFEARRAY *pr,
                                            SAFEARRAY *pi);

HRESULT CVIFUNC MLApp_DIMLAppExecute (CAObjHandle objectHandle,
                                      ERRORINFO *errorInfo, const char *name,
                                      char **returnValue);

HRESULT CVIFUNC MLApp_DIMLAppMinimizeCommandWindow (CAObjHandle objectHandle,
                                                    ERRORINFO *errorInfo);

HRESULT CVIFUNC MLApp_DIMLAppMaximizeCommandWindow (CAObjHandle objectHandle,
                                                    ERRORINFO *errorInfo);

HRESULT CVIFUNC MLApp_DIMLAppQuit (CAObjHandle objectHandle,
                                   ERRORINFO *errorInfo);

HRESULT CVIFUNC MLApp_DIMLAppGetCharArray (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace,
                                           char **returnValue);

HRESULT CVIFUNC MLApp_DIMLAppPutCharArray (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace,
                                           const char *charArray);

HRESULT CVIFUNC MLApp_DIMLAppGetVisible (CAObjHandle objectHandle,
                                         ERRORINFO *errorInfo, long *returnValue);

HRESULT CVIFUNC MLApp_DIMLAppSetVisible (CAObjHandle objectHandle,
                                         ERRORINFO *errorInfo, long newValue);

HRESULT CVIFUNC MLApp_DIMLAppGetWorkspaceData (CAObjHandle objectHandle,
                                               ERRORINFO *errorInfo,
                                               const char *name,
                                               const char *workspace,
                                               VARIANT *pdata);

HRESULT CVIFUNC MLApp_DIMLAppPutWorkspaceData (CAObjHandle objectHandle,
                                               ERRORINFO *errorInfo,
                                               const char *name,
                                               const char *workspace,
                                               VARIANT data);

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
                                    VARIANT arg31, VARIANT arg32);

HRESULT CVIFUNC MLApp_DIMLAppGetVariable (CAObjHandle objectHandle,
                                          ERRORINFO *errorInfo, const char *name,
                                          const char *workspace,
                                          VARIANT *returnValue);

HRESULT CVIFUNC MLApp_NewIMLApp (const char *server, int supportMultithreading,
                                 LCID locale, int reserved,
                                 CAObjHandle *objectHandle);

HRESULT CVIFUNC MLApp_OpenIMLApp (const char *fileName, const char *server,
                                  int supportMultithreading, LCID locale,
                                  int reserved, CAObjHandle *objectHandle);

HRESULT CVIFUNC MLApp_ActiveIMLApp (const char *server,
                                    int supportMultithreading, LCID locale,
                                    int reserved, CAObjHandle *objectHandle);

HRESULT CVIFUNC MLApp_IMLAppGetFullMatrix (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace, SAFEARRAY **pr,
                                           SAFEARRAY **pi);

HRESULT CVIFUNC MLApp_IMLAppPutFullMatrix (CAObjHandle objectHandle,
                                           ERRORINFO *errorInfo,
                                           const char *name,
                                           const char *workspace, SAFEARRAY *pr,
                                           SAFEARRAY *pi);

HRESULT CVIFUNC MLApp_IMLAppExecute (CAObjHandle objectHandle,
                                     ERRORINFO *errorInfo, const char *name);

HRESULT CVIFUNC MLApp_IMLAppMinimizeCommandWindow (CAObjHandle objectHandle,
                                                   ERRORINFO *errorInfo);

HRESULT CVIFUNC MLApp_IMLAppMaximizeCommandWindow (CAObjHandle objectHandle,
                                                   ERRORINFO *errorInfo);

HRESULT CVIFUNC MLApp_IMLAppQuit (CAObjHandle objectHandle, ERRORINFO *errorInfo);

HRESULT CVIFUNC MLApp_IMLAppGetCharArray (CAObjHandle objectHandle,
                                          ERRORINFO *errorInfo, const char *name,
                                          const char *workspace, char **mlString);

HRESULT CVIFUNC MLApp_IMLAppPutCharArray (CAObjHandle objectHandle,
                                          ERRORINFO *errorInfo, const char *name,
                                          const char *workspace,
                                          const char *charArray);

HRESULT CVIFUNC MLApp_IMLAppGetVisible (CAObjHandle objectHandle,
                                        ERRORINFO *errorInfo, long *pVal);

HRESULT CVIFUNC MLApp_IMLAppSetVisible (CAObjHandle objectHandle,
                                        ERRORINFO *errorInfo, long pVal);

HRESULT CVIFUNC MLApp_IMLAppGetWorkspaceData (CAObjHandle objectHandle,
                                              ERRORINFO *errorInfo,
                                              const char *name,
                                              const char *workspace,
                                              VARIANT *pdata);

HRESULT CVIFUNC MLApp_IMLAppPutWorkspaceData (CAObjHandle objectHandle,
                                              ERRORINFO *errorInfo,
                                              const char *name,
                                              const char *workspace,
                                              VARIANT data);

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
                                   VARIANT arg31, VARIANT arg32);

HRESULT CVIFUNC MLApp_IMLAppGetVariable (CAObjHandle objectHandle,
                                         ERRORINFO *errorInfo, const char *name,
                                         const char *workspace, VARIANT *pdata);

HRESULT CVIFUNC MLApp_IMLEvalXLEval (CAObjHandle objectHandle,
                                     ERRORINFO *errorInfo, const char *bstrName,
                                     long nargout, VARIANT *pvarArgOut,
                                     long nargin, VARIANT varArgIn);

HRESULT CVIFUNC MLApp_DIMLEvalXLEval (CAObjHandle objectHandle,
                                      ERRORINFO *errorInfo, const char *bstrName,
                                      long nargout, VARIANT *pvarArgOut,
                                      long nargin, VARIANT varArgIn);
#ifdef __cplusplus
    }
#endif
#endif /* _MLAPP_H */
