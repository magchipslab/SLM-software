/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/* Copyright (c) National Instruments 2015. All Rights Reserved.          */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  Control                          1       /* callback function: Control_Callback */
#define  Control_TAB                      2       /* control type: tab, callback function: (none) */

#define  DebugPanel                       2
#define  DebugPanel_CANVAS                2       /* control type: canvas, callback function: (none) */
#define  DebugPanel_CANVAS_2              3       /* control type: canvas, callback function: (none) */
#define  DebugPanel_CANVAS_3              4       /* control type: canvas, callback function: (none) */
#define  DebugPanel_CANVAS_4              5       /* control type: canvas, callback function: (none) */
#define  DebugPanel_CANVAS_5              6       /* control type: canvas, callback function: (none) */
#define  DebugPanel_CANVAS_6              7       /* control type: canvas, callback function: (none) */
#define  DebugPanel_DBGGRAPH3             8       /* control type: graph, callback function: (none) */
#define  DebugPanel_DBGGRAPH2             9       /* control type: graph, callback function: (none) */
#define  DebugPanel_DBGGRAPH              10      /* control type: graph, callback function: (none) */

#define  PANEL                            3       /* callback function: CameraPanel_Callback */
#define  PANEL_CameraCanvas               2       /* control type: canvas, callback function: CameraCanvas_Callback */
#define  PANEL_CameraGain                 3       /* control type: scale, callback function: CameraGain_Callback */
#define  PANEL_CameraShutter              4       /* control type: scale, callback function: CameraShutter_Callback */
#define  PANEL_CameraZoom                 5       /* control type: scale, callback function: CameraZoom_Callback */
#define  PANEL_HistoY                     6       /* control type: graph, callback function: (none) */
#define  PANEL_HistoX                     7       /* control type: graph, callback function: (none) */
#define  PANEL_CameraGridSpacing          8       /* control type: numeric, callback function: CameraGridSpacing_Callback */
#define  PANEL_CameraGrid                 9       /* control type: radioButton, callback function: CameraGrid_Callback */
#define  PANEL_SpotSequence               10      /* control type: command, callback function: SpotSequence_Callback */
#define  PANEL_CameraStartStop            11      /* control type: command, callback function: CameraStartStop_Callback */
#define  PANEL_CameraNumFrames            12      /* control type: numeric, callback function: CameraNumFrames_Callback */
#define  PANEL_ScreenshotButton           13      /* control type: command, callback function: ScreenshotButton_Callback */
#define  PANEL_ShowCrosshairs             14      /* control type: radioButton, callback function: ShowCrosshairs_Callback */
#define  PANEL_SpotSeqFilePrefix          15      /* control type: string, callback function: (none) */
#define  PANEL_CameraTimer                16      /* control type: timer, callback function: CameraTimer_Callback */
#define  PANEL_SpotSeqTimer               17      /* control type: timer, callback function: SpotSeqTimer_Callback */
#define  PANEL_CalibrateCam               18      /* control type: command, callback function: CalibrateCam_Callback */
#define  PANEL_CamCalFitPix               19      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalScanNy               20      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalScanOffsetY          21      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalScanOffsetX          22      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalDisplayPoint         23      /* control type: numeric, callback function: CamCalDisplayPoint_Callback */
#define  PANEL_CamCalScanNx               24      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalScanHeight           25      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalScanWidth            26      /* control type: numeric, callback function: (none) */
#define  PANEL_CamCalibration             27      /* control type: command, callback function: CamCalibration_Callback */
#define  PANEL_CamCalReCalClick           28      /* control type: command, callback function: CamCalReCalClick_Callback */
#define  PANEL_CamCalReCal                29      /* control type: command, callback function: CamCalReCal_Callback */
#define  PANEL_CamCalSave                 30      /* control type: command, callback function: CamCalSave_Callback */
#define  PANEL_CamCalLoad                 31      /* control type: command, callback function: CamCalLoad_Callback */
#define  PANEL_CamCalMarkZeroOrder        32      /* control type: command, callback function: CamCalMarkZeroOrder_Callback */
#define  PANEL_CamCalClearSLM             33      /* control type: radioButton, callback function: (none) */
#define  PANEL_CamCalToggle               34      /* control type: radioButton, callback function: CamCalToggle_Callback */
#define  PANEL_DECORATION                 35      /* control type: deco, callback function: (none) */
#define  PANEL_CalibrationLabel           36      /* control type: textMsg, callback function: (none) */
#define  PANEL_CalibrationProgress        37      /* control type: textMsg, callback function: (none) */
#define  PANEL_CurrSpotSeq                38      /* control type: textMsg, callback function: (none) */

#define  SimPanel                         4       /* callback function: SimPanel_Callback */
#define  SimPanel_CANVAS                  2       /* control type: canvas, callback function: (none) */
#define  SimPanel_SimZoom                 3       /* control type: scale, callback function: SimZoom_Callback */
#define  SimPanel_SimZPos                 4       /* control type: numeric, callback function: SimZPos_Callback */
#define  SimPanel_SimGridSpacing          5       /* control type: numeric, callback function: SimGridSpacing_Callback */
#define  SimPanel_SimGrid                 6       /* control type: radioButton, callback function: SimGrid_Callback */
#define  SimPanel_SimHelmholtz            7       /* control type: radioButton, callback function: SimHelmholtz_Callback */
#define  SimPanel_SimPhaseToggle          8       /* control type: radioButton, callback function: SimToggle_Callback */
#define  SimPanel_SimToggle               9       /* control type: radioButton, callback function: SimToggle_Callback */
#define  SimPanel_SimSaturation           10      /* control type: scale, callback function: SimSaturation_Callback */
#define  SimPanel_SaveSim                 11      /* control type: command, callback function: SaveSim_Callback */
#define  SimPanel_AmplMod                 12      /* control type: scale, callback function: AmplMod_Callback */
#define  SimPanel_SimSuperSample          13      /* control type: radioButton, callback function: SimSuperSample_Callback */
#define  SimPanel_SLMpixeltoggle          14      /* control type: radioButton, callback function: SimToggle_Callback */

#define  SLMpixels                        5       /* callback function: SLMpixels_Callback */
#define  SLMpixels_SLMcanvas              2       /* control type: canvas, callback function: (none) */

     /* tab page panel controls */
#define  TABPANE_10_FeedbackCanvas        2       /* control type: canvas, callback function: (none) */
#define  TABPANE_10_rbFeedbackCorrection  3       /* control type: radioButton, callback function: rbFeedbackToggle_Callback */
#define  TABPANE_10_rbFeedbackAdjusted    4       /* control type: radioButton, callback function: rbFeedbackToggle_Callback */
#define  TABPANE_10_rbFeedbackDifference  5       /* control type: radioButton, callback function: rbFeedbackToggle_Callback */
#define  TABPANE_10_rbFeedbackCamera      6       /* control type: radioButton, callback function: rbFeedbackToggle_Callback */
#define  TABPANE_10_rbFeedbackSignal      7       /* control type: radioButton, callback function: rbFeedbackToggle_Callback */
#define  TABPANE_10_UpdateFeedbackWindow  8       /* control type: command, callback function: UpdateFeedbackWindow_Callback */
#define  TABPANE_10_RefinePhase           9       /* control type: command, callback function: RefinePhase_Callback */
#define  TABPANE_10_WinSizeY              10      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_WinSizeX              11      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_EdgeTreshold          12      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_EdgeDetection         13      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_EdgeBorder            14      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_WinBorder             15      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_SignalWeight          16      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_AdjustmentFactor      17      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_IterationsAF          18      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_Repeat                19      /* control type: numeric, callback function: (none) */
#define  TABPANE_10_ClearCorr             20      /* control type: command, callback function: ClearCorr_Callback */
#define  TABPANE_10_FixSignal             21      /* control type: command, callback function: FixSignal_Callback */
#define  TABPANE_10_txtSignalQuality      22      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_txtSignalQualityVal   23      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_txtWinSize            24      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_rbPhaseUpdate         25      /* control type: radioButton, callback function: PhaseOption_Callback */
#define  TABPANE_10_rbPhaseSignal         26      /* control type: radioButton, callback function: PhaseOption_Callback */
#define  TABPANE_10_rbPhaseReset          27      /* control type: radioButton, callback function: PhaseOption_Callback */
#define  TABPANE_10_txtWinSizeVal         28      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_txtPhaseOption        29      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_MRAF                  30      /* control type: radioButton, callback function: (none) */
#define  TABPANE_10_FeedbackFilename      31      /* control type: string, callback function: (none) */
#define  TABPANE_10_txtIterationUpdate    32      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_itadjmax              33      /* control type: numeric, callback function: itcalc */
#define  TABPANE_10_nadj                  34      /* control type: numeric, callback function: itcalc */
#define  TABPANE_10_adjmax                35      /* control type: numeric, callback function: itcalc */
#define  TABPANE_10_adjmin                36      /* control type: numeric, callback function: itcalc */
#define  TABPANE_10_RunBatch              37      /* control type: command, callback function: RunBatchFeedback */
#define  TABPANE_10_totalit               38      /* control type: textMsg, callback function: (none) */
#define  TABPANE_10_DECORATION            39      /* control type: deco, callback function: (none) */

     /* tab page panel controls */
#define  TABPANE_11_InitPhaseError        2       /* control type: command, callback function: InitPhaseError_Callback */
#define  TABPANE_11_PhaseErrorClear       3       /* control type: command, callback function: PhaseErrorClearCallback */
#define  TABPANE_11_PhaseErrorBatchEval   4       /* control type: command, callback function: PhaseErrorBatchEvalCallback */
#define  TABPANE_11_PhaseErrorBatchApply  5       /* control type: command, callback function: PhaseErrorBatchApplyCallback */
#define  TABPANE_11_PhaseErrorBatchLoad   6       /* control type: command, callback function: PhaseErrorBatchLoadCallback */
#define  TABPANE_11_PhaseErrorBatchSave   7       /* control type: command, callback function: PhaseErrorBatchSaveCallback */
#define  TABPANE_11_PhaseErrorApply       8       /* control type: command, callback function: PhaseErrorApplyCallback */
#define  TABPANE_11_PhaseErrorBatchStop   9       /* control type: command, callback function: PhaseBatchStopCallback */
#define  TABPANE_11_PhaseErrorBatchGo     10      /* control type: command, callback function: PhaseBatchGoCallback */
#define  TABPANE_11_PhaseErrorGo          11      /* control type: command, callback function: PhaseErrorGoCallback */
#define  TABPANE_11_SSRDenominator        12      /* control type: numeric, callback function: SSRUpdate */
#define  TABPANE_11_SSRNumerator          13      /* control type: numeric, callback function: SSRUpdate */
#define  TABPANE_11_Slash                 14      /* control type: textMsg, callback function: (none) */
#define  TABPANE_11_SamplingPoints        15      /* control type: textMsg, callback function: (none) */
#define  TABPANE_11_SamplingPointsValue   16      /* control type: textMsg, callback function: (none) */
#define  TABPANE_11_LoadCJ                17      /* control type: command, callback function: LoadCJ_Callback */
#define  TABPANE_11_ClearCJ               18      /* control type: command, callback function: ClearCJ_Callback */
#define  TABPANE_11_RandomizeCJ           19      /* control type: command, callback function: RandomizeCJ_Callback */
#define  TABPANE_11_PhaseErrorNoiseCutoff 20      /* control type: numeric, callback function: (none) */
#define  TABPANE_11_PhaseErrorBeta        21      /* control type: numeric, callback function: (none) */
#define  TABPANE_11_PhaseErrorBatchSelect 22      /* control type: numeric, callback function: PhaseErrorBatchSelect_Callback */
#define  TABPANE_11_PhaseErrorBatchNumIt  23      /* control type: numeric, callback function: (none) */
#define  TABPANE_11_PhaseErrorPrjInterval 24      /* control type: numeric, callback function: (none) */
#define  TABPANE_11_PhaseErrorNumZernike  25      /* control type: numeric, callback function: (none) */
#define  TABPANE_11_PhaseErrorIterations  26      /* control type: numeric, callback function: (none) */
#define  TABPANE_11_PreallocateZernike    27      /* control type: radioButton, callback function: (none) */
#define  TABPANE_11_SamplingPointsValCam  28      /* control type: textMsg, callback function: (none) */
#define  TABPANE_11_BatchScore            29      /* control type: textMsg, callback function: (none) */
#define  TABPANE_11_BatchProgress         30      /* control type: textMsg, callback function: (none) */
#define  TABPANE_11_ManualSlide           31      /* control type: scale, callback function: ManualSlide_Callback */
#define  TABPANE_11_ManualCj              32      /* control type: numeric, callback function: ManualCj_Callback */
#define  TABPANE_11_PhaseErrorBatchTimer  33      /* control type: timer, callback function: PhaseErrorBatchTimer_Callback */
#define  TABPANE_11_AbbCorrMirrorY        34      /* control type: radioButton, callback function: (none) */
#define  TABPANE_11_AbbCorrMirrorX        35      /* control type: radioButton, callback function: (none) */
#define  TABPANE_11_ZerCorrFactor         36      /* control type: numeric, callback function: (none) */

     /* tab page panel controls */
#define  TABPANE_13_SHSpotSize            2       /* control type: numeric, callback function: ShackHartmann_Callback */
#define  TABPANE_13_SHSamplePointsY       3       /* control type: numeric, callback function: (none) */
#define  TABPANE_13_SHSamplePointsX       4       /* control type: numeric, callback function: (none) */
#define  TABPANE_13_NumZernikePolynomials 5       /* control type: numeric, callback function: (none) */
#define  TABPANE_13_SHStrayY              6       /* control type: scale, callback function: ShackHartmann_Callback */
#define  TABPANE_13_CalibratePhaseAberr   7       /* control type: command, callback function: CalibratePhaseAberr_Callback */
#define  TABPANE_13_SHStrayX              8       /* control type: scale, callback function: ShackHartmann_Callback */
#define  TABPANE_13_SHSpotY               9       /* control type: numeric, callback function: ShackHartmann_Callback */
#define  TABPANE_13_DECORATION            10      /* control type: deco, callback function: (none) */
#define  TABPANE_13_SHSpotX               11      /* control type: numeric, callback function: ShackHartmann_Callback */
#define  TABPANE_13_txtSHUnits            12      /* control type: textMsg, callback function: (none) */
#define  TABPANE_13_txtSHUnits_2          13      /* control type: textMsg, callback function: (none) */
#define  TABPANE_13_txtSHUnits_3          14      /* control type: textMsg, callback function: (none) */

     /* tab page panel controls */
#define  TABPANE_14_About1                2       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About3_3              3       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About3_2              4       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About3_4              5       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About3_5              6       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About3                7       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About2                8       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About3_6              9       /* control type: textMsg, callback function: (none) */
#define  TABPANE_14_About2_3              10      /* control type: textMsg, callback function: (none) */

     /* tab page panel controls */
#define  TABPANEL_HorizTrans              2       /* control type: scale, callback function: HorizTrans_Callback */
#define  TABPANEL_VertTrans               3       /* control type: scale, callback function: VertTrans_Callback */
#define  TABPANEL_LensXphase              4       /* control type: scale, callback function: Lens_Callback */
#define  TABPANEL_LensYphase              5       /* control type: scale, callback function: Lens_Callback */
#define  TABPANEL_LinkLensPhases          6       /* control type: radioButton, callback function: Lens_Callback */
#define  TABPANEL_Xpos                    7       /* control type: numeric, callback function: InputIntensity_Callback */
#define  TABPANEL_Ypos                    8       /* control type: numeric, callback function: InputIntensity_Callback */
#define  TABPANEL_SigmaX                  9       /* control type: numeric, callback function: InputIntensity_Callback */
#define  TABPANEL_SigmaY                  10      /* control type: numeric, callback function: InputIntensity_Callback */
#define  TABPANEL_LensFocalLength         11      /* control type: numeric, callback function: LensFocalLength_Callback */
#define  TABPANEL_Wavelength              12      /* control type: numeric, callback function: Wavelength_Callback */
#define  TABPANEL_Bias                    13      /* control type: numeric, callback function: Bias_Callback */
#define  TABPANEL_SubSample               14      /* control type: numeric, callback function: SubSample_Callback */
#define  TABPANEL_ClearSLM                15      /* control type: command, callback function: ClearSLM_Callback */
#define  TABPANEL_ClientMode              16      /* control type: command, callback function: ClientMode_Callback */
#define  TABPANEL_DECORATION              17      /* control type: deco, callback function: (none) */
#define  TABPANEL_DECORATION_2            18      /* control type: deco, callback function: (none) */
#define  TABPANEL_FocalUnitLabel          19      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_SLMResolutionLabel      20      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_SLMResolution           21      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_FocalUnit               22      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_SubSampleInv            23      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS                   24      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_3                 25      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_2                 26      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_4                 27      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_6                 28      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_7                 29      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_8                 30      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_10                31      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_UNITS_9                 32      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_rbInputLoadFromFile     33      /* control type: radioButton, callback function: rbInput_Callback */
#define  TABPANEL_rbInputAnnular          34      /* control type: radioButton, callback function: rbInput_Callback */
#define  TABPANEL_rbInputGaussian         35      /* control type: radioButton, callback function: rbInput_Callback */
#define  TABPANEL_txtInputLoadFromFile    36      /* control type: string, callback function: (none) */
#define  TABPANEL_btnInputLoadFromFile    37      /* control type: command, callback function: InputLoadFromFile_Callback */
#define  TABPANEL_txtInputIntensity       38      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_DECORATION_3            39      /* control type: deco, callback function: (none) */
#define  TABPANEL_UNITS_5                 40      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_Magnification           41      /* control type: numeric, callback function: Magnification_Callback */
#define  TABPANEL_Focallengthlabel        42      /* control type: textMsg, callback function: (none) */
#define  TABPANEL_EffectiveFocalLength    43      /* control type: textMsg, callback function: (none) */

     /* tab page panel controls */
#define  TABPANEL_2_TAB                   2       /* control type: tab, callback function: (none) */
#define  TABPANEL_2_rbPhaseRetrieval      3       /* control type: radioButton, callback function: rbPattern_Callback */
#define  TABPANEL_2_rbPhaseGrating        4       /* control type: radioButton, callback function: rbPattern_Callback */
#define  TABPANEL_2_rbSpotProjector       5       /* control type: radioButton, callback function: rbPattern_Callback */
#define  TABPANEL_2_SLMPatternMode        6       /* control type: textMsg, callback function: (none) */

     /* tab page panel controls */
#define  TABPANEL_3_inputNy               2       /* control type: numeric, callback function: SpotControls_Callback */
#define  TABPANEL_3_inputNx               3       /* control type: numeric, callback function: SpotControls_Callback */
#define  TABPANEL_3_Yspacing              4       /* control type: scale, callback function: SpotControls_Callback */
#define  TABPANEL_3_RandomAmpl            5       /* control type: scale, callback function: SpotControls_Callback */
#define  TABPANEL_3_Spotlens              6       /* control type: scale, callback function: SpotControls_Callback */
#define  TABPANEL_3_Xspacing              7       /* control type: scale, callback function: SpotControls_Callback */

     /* tab page panel controls */
#define  TABPANEL_4_Angle                 2       /* control type: scale, callback function: PhaseGrating_Callback */
#define  TABPANEL_4_Period2               3       /* control type: scale, callback function: PhaseGrating_Callback */
#define  TABPANEL_4_Amplitude2            4       /* control type: scale, callback function: PhaseGrating_Callback */
#define  TABPANEL_4_Period1               5       /* control type: scale, callback function: PhaseGrating_Callback */
#define  TABPANEL_4_Amplitude1            6       /* control type: scale, callback function: PhaseGrating_Callback */
#define  TABPANEL_4_GaussianPeaks         7       /* control type: radioButton, callback function: GaussianPeaks_Callback */

     /* tab page panel controls */
#define  TABPANEL_5_PhaseRetrievalGo      2       /* control type: command, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_Output                3       /* control type: textBox, callback function: (none) */
#define  TABPANEL_5_SpotYSpacing          4       /* control type: scale, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_SpotXSpacing          5       /* control type: scale, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_NumYSpots             6       /* control type: numeric, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_NumXSpots             7       /* control type: numeric, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_SpotXOffset           8       /* control type: numeric, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_SpotYOffset           9       /* control type: numeric, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_SigYOffset            10      /* control type: numeric, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_SigXOffset            11      /* control type: numeric, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_PhaseConstraint       12      /* control type: radioButton, callback function: PhaseRetrieval_Callback */
#define  TABPANEL_5_PhaseStep             13      /* control type: numeric, callback function: PhaseRetrieval_Callback */

     /* tab page panel controls */
#define  TABPANEL_6_TAB                   2       /* control type: tab, callback function: (none) */
#define  TABPANEL_6_rbShapingModeArb      3       /* control type: radioButton, callback function: rbShapingMode_Callback */
#define  TABPANEL_6_rbShapingModeStd      4       /* control type: radioButton, callback function: rbShapingMode_Callback */
#define  TABPANEL_6_ShapingMode           5       /* control type: textMsg, callback function: (none) */

     /* tab page panel controls */
#define  TABPANEL_7_SignalHeight          2       /* control type: numeric, callback function: (none) */
#define  TABPANEL_7_SignalWidth           3       /* control type: numeric, callback function: (none) */
#define  TABPANEL_7_SigmaY                4       /* control type: numeric, callback function: (none) */
#define  TABPANEL_7_SigmaX                5       /* control type: numeric, callback function: (none) */
#define  TABPANEL_7_BeamShapeOK           6       /* control type: command, callback function: BeamShape_Callback */
#define  TABPANEL_7_GRAPH                 7       /* control type: graph, callback function: (none) */
#define  TABPANEL_7_BeamTypeSquare        8       /* control type: radioButton, callback function: BeamType_Callback */
#define  TABPANEL_7_BeamTypeGaussian      9       /* control type: radioButton, callback function: BeamType_Callback */
#define  TABPANEL_7_DECORATION            10      /* control type: deco, callback function: (none) */
#define  TABPANEL_7_TEXTMSG               11      /* control type: textMsg, callback function: (none) */

     /* tab page panel controls */
#define  TABPANEL_8_BeamShapeFilename     2       /* control type: string, callback function: (none) */
#define  TABPANEL_8_LoadFile              3       /* control type: command, callback function: LoadFile_Callback */
#define  TABPANEL_8_PicPatternGo          4       /* control type: command, callback function: PicPattern_Callback */
#define  TABPANEL_8_PictureHeight         5       /* control type: numeric, callback function: (none) */
#define  TABPANEL_8_PictureWidth          6       /* control type: numeric, callback function: (none) */
#define  TABPANEL_8_IterationsPF          7       /* control type: numeric, callback function: (none) */
#define  TABPANEL_8_IterationsAF          8       /* control type: numeric, callback function: (none) */
#define  TABPANEL_8_TEXTBOX               9       /* control type: textBox, callback function: (none) */
#define  TABPANEL_8_SigAmp                10      /* control type: numeric, callback function: (none) */
#define  TABPANEL_8_BeamShapeRefinePhase  11      /* control type: radioButton, callback function: (none) */
#define  TABPANEL_8_Picture               12      /* control type: picture, callback function: (none) */
#define  TABPANEL_8_MRAF                  13      /* control type: radioButton, callback function: (none) */
#define  TABPANEL_8_SoftOp                14      /* control type: radioButton, callback function: (none) */
#define  TABPANEL_8_CharacterGo           15      /* control type: command, callback function: CharacterGo_Callback */
#define  TABPANEL_8_CharacterNumber       16      /* control type: numeric, callback function: CharacterNumber_Callback */
#define  TABPANEL_8_CharacterCanvas       17      /* control type: canvas, callback function: (none) */
#define  TABPANEL_8_CharacterBold         18      /* control type: radioButton, callback function: CharacterNumber_Callback */
#define  TABPANEL_8_CharacterFont         19      /* control type: string, callback function: CharacterNumber_Callback */
#define  TABPANEL_8_CharacterFontsize     20      /* control type: numeric, callback function: CharacterNumber_Callback */
#define  TABPANEL_8_CharacterPosY         21      /* control type: numeric, callback function: CharacterNumber_Callback */
#define  TABPANEL_8_CharacterPosX         22      /* control type: numeric, callback function: CharacterNumber_Callback */
#define  TABPANEL_8_DECORATION            23      /* control type: deco, callback function: (none) */

     /* tab page panel controls */
#define  TABPANEL_9_SaveSLMPattern        2       /* control type: command, callback function: SaveSLMPattern_Callback */
#define  TABPANEL_9_LoadSLMPattern        3       /* control type: command, callback function: LoadSLMPattern_Callback */
#define  TABPANEL_9_LoadFilename          4       /* control type: string, callback function: (none) */
#define  TABPANEL_9_LoadKeepInpInt        5       /* control type: radioButton, callback function: (none) */
#define  TABPANEL_9_LoadKeepSettings      6       /* control type: radioButton, callback function: (none) */
#define  TABPANEL_9_LoadSaveCameraState   7       /* control type: radioButton, callback function: (none) */
#define  TABPANEL_9_FactoryCorrection     8       /* control type: string, callback function: (none) */
#define  TABPANEL_9_LoadFactoryCorrection 9       /* control type: command, callback function: LoadFactoryCorrection_Callback */
#define  TABPANEL_9_FactorySwitch         10      /* control type: radioButton, callback function: FactorySwitch_Callback */
#define  TABPANEL_9_SHcorrection          11      /* control type: string, callback function: (none) */
#define  TABPANEL_9__LoadSHcorrection     12      /* control type: command, callback function: LoadSHcorrection_Callback */
#define  TABPANEL_9_SHswitch              13      /* control type: radioButton, callback function: SHswitch_Callback */
#define  TABPANEL_9_Xmin                  14      /* control type: numeric, callback function: (none) */
#define  TABPANEL_9_Xmax                  15      /* control type: numeric, callback function: (none) */
#define  TABPANEL_9_Ymin                  16      /* control type: numeric, callback function: (none) */
#define  TABPANEL_9_Ymax                  17      /* control type: numeric, callback function: (none) */
#define  TABPANEL_9_IntensityMatrix       18      /* control type: command, callback function: IntensityMatrix_Callback */


     /* Control Arrays: */

#define  CTRLARRAY                        1

     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK AmplMod_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK BeamShape_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK BeamType_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Bias_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CalibrateCam_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CalibratePhaseAberr_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalDisplayPoint_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalibration_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalLoad_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalMarkZeroOrder_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalReCal_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalReCalClick_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalSave_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CamCalToggle_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraCanvas_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraGain_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraGrid_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraGridSpacing_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraNumFrames_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraPanel_Callback(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraShutter_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraStartStop_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraTimer_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CameraZoom_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CharacterGo_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK CharacterNumber_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ClearCJ_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ClearCorr_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ClearSLM_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ClientMode_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Control_Callback(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK FactorySwitch_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK FixSignal_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK GaussianPeaks_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK HorizTrans_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK InitPhaseError_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK InputIntensity_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK InputLoadFromFile_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK IntensityMatrix_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK itcalc(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Lens_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK LensFocalLength_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK LoadCJ_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK LoadFactoryCorrection_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK LoadFile_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK LoadSHcorrection_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK LoadSLMPattern_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Magnification_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ManualCj_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ManualSlide_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseBatchGoCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseBatchStopCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorApplyCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorBatchApplyCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorBatchEvalCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorBatchLoadCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorBatchSaveCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorBatchSelect_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorBatchTimer_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorClearCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseErrorGoCallback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseGrating_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseOption_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PhaseRetrieval_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK PicPattern_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK RandomizeCJ_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK rbFeedbackToggle_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK rbInput_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK rbPattern_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK rbShapingMode_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK RefinePhase_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK RunBatchFeedback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveSim_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveSLMPattern_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ScreenshotButton_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ShackHartmann_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK ShowCrosshairs_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SHswitch_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimGrid_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimGridSpacing_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimHelmholtz_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimPanel_Callback(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimSaturation_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimSuperSample_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimToggle_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimZoom_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SimZPos_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SLMpixels_Callback(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SpotControls_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SpotSeqTimer_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SpotSequence_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SSRUpdate(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SubSample_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK UpdateFeedbackWindow_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK VertTrans_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK Wavelength_Callback(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif
