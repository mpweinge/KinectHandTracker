/****************************************************************************
 *                                                                           *
 *  OpenNI 1.x Alpha                                                         *
 *  Copyright (C) 2011 PrimeSense Ltd.                                       *
 *                                                                           *
 *  This file is part of OpenNI.                                             *
 *                                                                           *
 *  OpenNI is free software: you can redistribute it and/or modify           *
 *  it under the terms of the GNU Lesser General Public License as published *
 *  by the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  OpenNI is distributed in the hope that it will be useful,                *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU Lesser General Public License *
 *  along with OpenNI. If not, see <http://www.gnu.org/licenses/>.           *
 *                                                                           *
 ****************************************************************************/
//---------------------------------------------------------------------------
// Includes
//---------------------------------------------------------------------------
#include "NiHandViewer.h"
#include "XnCppWrapper.h"

#if (XN_PLATFORM == XN_PLATFORM_LINUX_X86 || XN_PLATFORM == XN_PLATFORM_LINUX_ARM)
#define UNIX
#define GLX_GLXEXT_LEGACY
#endif

#if (XN_PLATFORM == XN_PLATFORM_MACOSX)
#define MACOS
#endif

//#define GLH_EXT_SINGLE_FILE
//#pragma warning(push, 3)
//#include "glh/glh_obs.h"
//#include "glh/glh_glut2.h"
//#pragma warning(pop)

#if (XN_PLATFORM == XN_PLATFORM_MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

//using namespace glh;

#include "XnLog.h"
#include "Device.h"
#include "Capture.h"
#include "Draw.h"
#include "Keyboard.h"
#include "Menu.h"
#include "Statistics.h"
#include "MouseInput.h"
#include "NiHandTracker.h"


#if (XN_PLATFORM == XN_PLATFORM_WIN32)
#include <conio.h>
#include <direct.h>
#elif (XN_PLATFORM == XN_PLATFORM_LINUX_X86 || XN_PLATFORM == XN_PLATFORM_LINUX_ARM || XN_PLATFORM == XN_PLATFORM_MACOSX)
#define _getch() getchar()
#endif

//---------------------------------------------------------------------------
// Defines
//---------------------------------------------------------------------------
#define LENGTHOF(arr)			(sizeof(arr)/sizeof(arr[0]))
#define FOR_ALL(arr, action)	{for(int i = 0; i < LENGTHOF(arr); ++i){action(arr[i])}}

#define ADD_GESTURE(name)		{if(gGestureGenerator.AddGesture(name, NULL) != XN_STATUS_OK){printf("Unable to add gesture"); exit(1);}}
#define REMOVE_GESTURE(name)	{if(gGestureGenerator.RemoveGesture(name) != XN_STATUS_OK){printf("Unable to remove gesture"); exit(1);}}

#define ADD_ALL_GESTURES		FOR_ALL(cGestures, ADD_GESTURE)
#define REMOVE_ALL_GESTURES		FOR_ALL(cGestures, REMOVE_GESTURE)


//---------------------------------------------------------------------------
// Consts
//---------------------------------------------------------------------------
// Gestures to track
static const char			cClickStr[] = "Click";
static const char			cWaveStr[] = "Wave";
static const char           cSteadyStr[] = "Steady";
static const char           cRaiseHand[] = "RaiseHand";
static const char           cMovingHand[] = "MovingHand";
static const char* const	cGestures[] =
{
	cClickStr,
	cWaveStr/*,
    cRaiseHand*/
};

//---------------------------------------------------------------------------
// Defines
//---------------------------------------------------------------------------
#define SAMPLE_XML_PATH "/Users/Mikey/KinectHandTracker_CommandLine/KinectHandTracker_CommandLine/SamplesConfig.xml"

// --------------------------------
// Types
// --------------------------------
enum {
	ERR_OK,
	ERR_USAGE,
	ERR_DEVICE,
	ERR_UNKNOWN
};

//---------------------------------------------------------------------------
// Globals
//---------------------------------------------------------------------------
xn::ScriptNode	g_scriptNodeMain;
xn::DepthGenerator g_DepthGenerator;
xn::UserGenerator g_UserGenerator;
HandTracker*	m_HandTracker;
HandsGenerator gHandsGenerator;
GestureGenerator gGestureGenerator;

XnBool g_bNeedPose = FALSE;
XnChar g_strPose[20] = "";

#define CHECK_RC(nRetVal, what)					    \
if (nRetVal != XN_STATUS_OK)				    \
{								    \
printf("%s failed: %s\n", what, xnGetStatusString(nRetVal));    \
return nRetVal;						    \
}

/* When true, frames will not be read from the device. */
bool g_bPause = false;
/* When true, only a single frame will be read, and then reading will be paused. */
bool g_bStep = false;

//glut_perspective_reshaper reshaper;
//glut_callbacks cb;

IntPair mouseLocation;
IntPair windowSize;

// --------------------------------
// Utilities
// --------------------------------
void MotionCallback(int x, int y)
{
	mouseInputMotion(int((double)x/windowSize.X*WIN_SIZE_X), int((double)y/windowSize.Y*WIN_SIZE_Y));
}

void MouseCallback(int button, int state, int x, int y)
{
	mouseInputButton(button, state, int((double)x/windowSize.X*WIN_SIZE_X), int((double)y/windowSize.Y*WIN_SIZE_Y));
}

void KeyboardCallback(unsigned char key, int /*x*/, int /*y*/)
{
	if (isCapturing())
	{
		captureStop(0);
	}
	else
	{
		handleKey(key);
	}
    
	//glutPostRedisplay();
}

void ReshapeCallback(int width, int height)
{
	windowSize.X = width;
	windowSize.Y = height;
	windowReshaped(width, height);
}

void IdleCallback()
{
	XnStatus nRetVal = XN_STATUS_OK;
	if (g_bPause != TRUE)
	{
		// read a frame
		readFrame();
        
		// capture if needed
		nRetVal = captureFrame();
		if (nRetVal != XN_STATUS_OK)
		{
			displayMessage("Error capturing frame: '%s'", xnGetStatusString(nRetVal));
		}
        
		// add to statistics
		statisticsAddFrame();
	}
    
	if (g_bStep == TRUE)
	{
		g_bStep = FALSE;
		g_bPause = TRUE;
	}
    
	glutPostRedisplay();
}

void seek(int nDiff)
{
	if (!isPlayerOn())
	{
		displayMessage("Seeking is only supported in playback mode!");
		return;
	}
    
	seekFrame(nDiff);
    
	// now step the last one (that way, if seek is not supported, as in sensor, at least one frame
	// will be read).
	g_bPause = false;
	g_bStep = true;
}

void init_opengl()
{
	glClearStencil(128);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_NORMALIZE);
    
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	GLfloat ambient[] = {0.5, 0.5, 0.5, 1};
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightf (GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.1f);
    glEnableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
}

void closeSample(int errCode)
{
	captureStop(0);
	closeDevice();
    
	if (errCode != ERR_OK)
	{
		printf("Press any key to continue . . .\n");
		_getch();
	}
    
	exit(errCode);
}

void togglePause(int)
{
	g_bPause = !g_bPause;
}

void resetCropping(MapGenerator* pGenerator)
{
	XnCropping crop;
	crop.bEnabled = FALSE;
	crop.nXOffset = crop.nYOffset = 0;
	crop.nXSize = crop.nYSize = 0;
	setStreamCropping(pGenerator, &crop);
}

void resetDepthCropping(int)
{
	resetCropping(getDepthGenerator());
}

void resetImageCropping(int)
{
	resetCropping(getImageGenerator());
}

void resetIRCropping(int)
{
	resetCropping(getIRGenerator());
}

void resetAllCropping(int)
{
	if (getDepthGenerator() != NULL)
		resetDepthCropping(0);
    
	if (getImageGenerator() != NULL)
		resetImageCropping(0);
    
	if (getIRGenerator() != NULL)
		resetIRCropping(0);
}

void step(int)
{
	g_bStep = true;
	g_bPause = false;
}

void startCapture(int delay)
{
	if (g_bPause)
	{
		displayMessage("Cannot record when paused!");
	}
	else
	{
		captureStart(delay);
	}
}

void createKeyboardMap()
{
	startKeyboardMap();
	{
		startKeyboardGroup(KEYBOARD_GROUP_PRESETS);
		{
			registerKey('`', getPresetName(0), setPreset, 0);
			registerKey('1', getPresetName(1), setPreset, 1);
			registerKey('2', getPresetName(2), setPreset, 2);
			registerKey('3', getPresetName(3), setPreset, 3);
			registerKey('4', getPresetName(4), setPreset, 4);
			registerKey('5', getPresetName(5), setPreset, 5);
			registerKey('6', getPresetName(6), setPreset, 6);
			registerKey('7', getPresetName(7), setPreset, 7);
			registerKey('8', getPresetName(8), setPreset, 8);
			registerKey('9', getPresetName(9), setPreset, 9);
			registerKey('0', getPresetName(10), setPreset, 10);
			registerKey('-', getPresetName(11), setPreset, 11);
			registerKey('=', getPresetName(12), setPreset, 12);
		}
		endKeyboardGroup();
		startKeyboardGroup(KEYBOARD_GROUP_DEVICE);
		{
			registerKey('m', "Mirror on/off", toggleMirror, 0);
			registerKey('/', "Reset all croppings", resetAllCropping, 0);
		}
		endKeyboardGroup();
		startKeyboardGroup(KEYBOARD_GROUP_CAPTURE);
		{
			registerKey('s', "Start", startCapture, 0);
			registerKey('d', "Start (5 sec delay)", startCapture, 5);
			registerKey('x', "Stop", captureStop, 0);
			registerKey('c', "Capture current frame only", captureSingleFrame, 0);
		}
		endKeyboardGroup();
		startKeyboardGroup(KEYBOARD_GROUP_DISPLAY);
		{
			registerKey('p', "Pointer Mode On/Off", togglePointerMode, 0);
			registerKey('f', "Full Screen On/Off", toggleFullScreen, 0);
			registerKey('?', "Show/Hide Help screen", toggleHelpScreen, 0);
		}
		endKeyboardGroup();
		startKeyboardGroup(KEYBOARD_GROUP_GENERAL);
		{
			registerKey('z', "Start/Stop Collecting Statistics", toggleStatistics, 0);
			registerKey('?', "Show/Hide help screen", toggleHelpScreen, 0);
			registerKey(27, "Exit", closeSample, ERR_OK);
		}
		endKeyboardGroup();
		startKeyboardGroup(KEYBOARD_GROUP_PLAYER);
		{
			registerKey('o', "Pause/Resume", togglePause, 0);
			registerKey('l', "Jump 1 frame forward", seek, 1);
			registerKey('L', "Jump 10 frames forward", seek, 10);
			registerKey('k', "Jump 1 frame backwards", seek, -1);
			registerKey('K', "Jump 10 frames backwards", seek, -10);
			registerKey(';', "Read one frame", step, 0);
			registerKey('[', "Decrease playback speed", setPlaybackSpeed, -1);
			registerKey(']', "Increase playback speed", setPlaybackSpeed, 1);
		}
		endKeyboardGroup();
	}
	endKeyboardMap();
}

void createMenu()
{
	startMenu();
	{
		startSubMenu("View");
		{
			startSubMenu("Presets");
			{
				for (int i = 0; i < PRESET_COUNT; ++i)
				{
					createMenuEntry(getPresetName(i), setPreset, i);
				}
			}
			endSubMenu();
			startSubMenu("Screen Layout");
			{
				createMenuEntry("Side by Side", setScreenLayout, (int)SIDE_BY_SIDE);
				createMenuEntry("Overlay", setScreenLayout, (int)OVERLAY);
			}
			endSubMenu();
			startSubMenu("Depth");
			{
				for (int i = 0; i < NUM_OF_DEPTH_TYPES; ++i)
				{
					createMenuEntry(g_DepthColoring[i], setDepthDrawing, i);
				}
			}
			endSubMenu();
			startSubMenu("Image");
			{
				for (int i = 0; i < NUM_OF_IMAGE_TYPES; ++i)
				{
					createMenuEntry(g_ImageColoring[i], setImageDrawing, i);
				}
			}
			endSubMenu();
			createMenuEntry("Pointer Mode On/Off", togglePointerMode, 0);
			createMenuEntry("Show/Hide Background", toggleBackground, 0);
			createMenuEntry("Show/Hide Help Screen", toggleHelpScreen, 0);
		}
		endSubMenu();
		startSubMenu("Device");
		{
			startSubMenu("Streams");
			{
				startSubMenu("Depth");
				{
					createMenuEntry("On/Off", toggleDepthState, 0);
					startSubMenu("Registration");
					{
						for (int i = 0; i < g_Registration.nValuesCount; ++i)
						{
							unsigned int nValue = g_Registration.pValues[i];
							createMenuEntry(g_Registration.pValueToName[nValue], changeRegistration, nValue);
						}
					}
					endSubMenu();
					startSubMenu("Resolution");
					{
						createMenuEntry("QVGA", setDepthResolution, XN_RES_QVGA);
						createMenuEntry("VGA", setDepthResolution, XN_RES_VGA);
					}
					endSubMenu();
					startSubMenu("FPS");
					{
						createMenuEntry("25", setDepthFPS, 25);
						createMenuEntry("30", setDepthFPS, 30);
						createMenuEntry("60", setDepthFPS, 60);
					}
					endSubMenu();
					createMenuEntry("Reset Cropping", resetDepthCropping, 0);
				}
				endSubMenu();
				startSubMenu("Image");
				{
					createMenuEntry("On/Off", toggleImageState, 0);
					startSubMenu("Resolution");
					{
						createMenuEntry("QVGA", setImageResolution, XN_RES_QVGA);
						createMenuEntry("VGA", setImageResolution, XN_RES_VGA);
						createMenuEntry("SXGA", setImageResolution, XN_RES_SXGA);
						createMenuEntry("UXGA", setImageResolution, XN_RES_UXGA);
					}
					endSubMenu();
					startSubMenu("FPS");
					{
						createMenuEntry("25", setImageFPS, 25);
						createMenuEntry("30", setImageFPS, 30);
						createMenuEntry("60", setImageFPS, 60);
					}
					endSubMenu();
					createMenuEntry("Reset Cropping", resetImageCropping, 0);
				}
				endSubMenu();
				startSubMenu("IR");
				{
					createMenuEntry("On/Off", toggleIRState, 0);
					startSubMenu("Resolution");
					{
						createMenuEntry("QVGA", setIRResolution, XN_RES_QVGA);
						createMenuEntry("VGA", setIRResolution, XN_RES_VGA);
						createMenuEntry("SXGA", setIRResolution, XN_RES_SXGA);
					}
					endSubMenu();
					startSubMenu("FPS");
					{
						createMenuEntry("25", setIRFPS, 25);
						createMenuEntry("30", setIRFPS, 30);
						createMenuEntry("60", setIRFPS, 60);
					}
					endSubMenu();
					createMenuEntry("Reset Cropping", resetIRCropping, 0);
				}
				endSubMenu();
				startSubMenu("Primary Stream");
				{
					for (int i = 0; i < g_PrimaryStream.nValuesCount; ++i)
					{
						createMenuEntry(g_PrimaryStream.pValues[i], changePrimaryStream, i);
					}
				}
				endSubMenu();
			}
			endSubMenu();
            
			createMenuEntry("Mirror", toggleMirror, 0);
		}
		endSubMenu();
		startSubMenu("Capture");
		{
			startSubMenu("Depth Format");
			{
				for (int i = 0; i < g_DepthFormat.nValuesCount; ++i)
				{
					unsigned int nValue = g_DepthFormat.pValues[i];
					createMenuEntry(g_DepthFormat.pIndexToName[i], captureSetDepthFormat, nValue);
				}
			}
			endSubMenu();
			startSubMenu("Image Format");
			{
				for (int i = 0; i < g_ImageFormat.nValuesCount; ++i)
				{
					unsigned int nValue = g_ImageFormat.pValues[i];
					createMenuEntry(g_ImageFormat.pIndexToName[i], captureSetImageFormat, nValue);
				}
			}
			endSubMenu();
			startSubMenu("IR Format");
			{
				for (int i = 0; i < g_IRFormat.nValuesCount; ++i)
				{
					unsigned int nValue = g_IRFormat.pValues[i];
					createMenuEntry(g_IRFormat.pIndexToName[i], captureSetIRFormat, nValue);
				}
			}
			endSubMenu();
			createMenuEntry("Browse", captureBrowse, 0);
			createMenuEntry("Start", captureStart, 0);
			createMenuEntry("Start (5 sec delay)", captureStart, 5);
			createMenuEntry("Restart", captureRestart, 0);
			createMenuEntry("Stop", captureStop, 0);
		}
		endSubMenu();
		startSubMenu("Statistics");
		{
			createMenuEntry("Start Collecting", statisticsStart, 0);
			createMenuEntry("Stop Collecting", statisticsStop, 0);
			createMenuEntry("Clear", statisticsClear, 0);
		}
		endSubMenu();
		startSubMenu("Player");
		{
			createMenuEntry("Pause/Resume", togglePause, 0);
			createMenuEntry("Skip 1 frame forward", seek, 1);
			createMenuEntry("Skip 10 frame forward", seek, 10);
			createMenuEntry("Skip 1 frame backwards", seek, -1);
			createMenuEntry("Skip 10 frame backwards", seek, -10);
		}
		endSubMenu();
		createMenuEntry("Quit", closeSample, ERR_OK);
	}
	endMenu();
}

void onExit()
{
	captureStop(0);
}

int changeDirectory(char* arg0)
{
	// get dir name
	XnChar strDirName[XN_FILE_MAX_PATH];
	XnStatus nRetVal = xnOSGetDirName(arg0, strDirName, XN_FILE_MAX_PATH);
	XN_IS_STATUS_OK(nRetVal);
    
	// now set current directory to it
	nRetVal = xnOSSetCurrentDir(strDirName);
	XN_IS_STATUS_OK(nRetVal);
    
	return 0;
}

void XN_CALLBACK_TYPE UserCalibration_CalibrationComplete(xn::SkeletonCapability& /*capability*/, XnUserID nId, XnCalibrationStatus eStatus, void* /*pCookie*/)
{
    XnUInt32 epochTime = 0;
    xnOSGetEpochTime(&epochTime);
    if (eStatus == XN_CALIBRATION_STATUS_OK)
    {
        // Calibration succeeded
        printf("%d Calibration complete, start tracking user %d\n", epochTime, nId);
        g_UserGenerator.GetSkeletonCap().StartTracking(nId);
    }
    else
    {
        // Calibration failed
        printf("%d Calibration failed for user %d\n", epochTime, nId);
        if(eStatus==XN_CALIBRATION_STATUS_MANUAL_ABORT)
        {
            printf("Manual abort occured, stop attempting to calibrate!");
            return;
        }
        if (g_bNeedPose)
        {
            g_UserGenerator.GetPoseDetectionCap().StartPoseDetection(g_strPose, nId);
        }
        else
        {
            g_UserGenerator.GetSkeletonCap().RequestCalibration(nId, TRUE);
        }
    }
}

// Callback: New user was detected
void XN_CALLBACK_TYPE User_NewUser(xn::UserGenerator& /*generator*/, XnUserID nId, void* /*pCookie*/)
{
    XnUInt32 epochTime = 0;
    xnOSGetEpochTime(&epochTime);
    printf("%d New User %d\n", epochTime, nId);
    // New user found
    if (g_bNeedPose)
    {
        g_UserGenerator.GetPoseDetectionCap().StartPoseDetection(g_strPose, nId);
    }
    else
    {
        g_UserGenerator.GetSkeletonCap().RequestCalibration(nId, TRUE);
    }
}
// Callback: An existing user was lost
void XN_CALLBACK_TYPE User_LostUser(xn::UserGenerator& /*generator*/, XnUserID nId, void* /*pCookie*/)
{
    XnUInt32 epochTime = 0;
    xnOSGetEpochTime(&epochTime);
    printf("%d Lost user %d\n", epochTime, nId);
}
// Callback: Detected a pose
void XN_CALLBACK_TYPE UserPose_PoseDetected(xn::PoseDetectionCapability& /*capability*/, const XnChar* strPose, XnUserID nId, void* /*pCookie*/)
{
    XnUInt32 epochTime = 0;
    xnOSGetEpochTime(&epochTime);
    printf("%d Pose %s detected for user %d\n", epochTime, strPose, nId);
    g_UserGenerator.GetPoseDetectionCap().StopPoseDetection(nId);
    g_UserGenerator.GetSkeletonCap().RequestCalibration(nId, TRUE);
}
// Callback: Started calibration
void XN_CALLBACK_TYPE UserCalibration_CalibrationStart(xn::SkeletonCapability& /*capability*/, XnUserID nId, void* /*pCookie*/)
{
    XnUInt32 epochTime = 0;
    xnOSGetEpochTime(&epochTime);
    printf("%d Calibration started for user %d\n", epochTime, nId);
}

void XN_CALLBACK_TYPE Hand_Create(	xn::HandsGenerator& /*generator*/,
                                               XnUserID			nId,
                                               const XnPoint3D*	pPosition,
                                               XnFloat				/*fTime*/,
                                               void*				pCookie)
{
	printf("New Hand: %d @ (%f,%f,%f)\n", nId, pPosition->X, pPosition->Y, pPosition->Z);
    
	/*HandTracker*	pThis = static_cast<HandTracker*>(pCookie);
	if(sm_Instances.Find(pThis) == sm_Instances.End())
	{
		printf("Dead HandTracker: skipped!\n");
		return;
	}
    
	pThis->m_History[nId].Push(*pPosition);*/
}

void XN_CALLBACK_TYPE Hand_Update(	xn::HandsGenerator& /*generator*/,
                                               XnUserID			nId,
                                               const XnPoint3D*	pPosition,
                                               XnFloat				/*fTime*/,
                                               void*				pCookie)
{
	//HandTracker*	pThis = static_cast<HandTracker*>(pCookie);
	/*if(sm_Instances.Find(pThis) == sm_Instances.End())
	{
		printf("Dead HandTracker: skipped!\n");
		return;
	}
    
	// Add to this user's hands history
	TrailHistory::Iterator it = pThis->m_History.Find(nId);
	if (it == pThis->m_History.End())
	{
		printf("Dead hand update: skipped!\n");
		return;
	}
    
	it->Value().Push(*pPosition);*/
    printf("New Position: %d @ (%f,%f,%f)\n", nId, pPosition->X, pPosition->Y, pPosition->Z);
}

void XN_CALLBACK_TYPE Hand_Destroy(	xn::HandsGenerator& /*generator*/,
                                                XnUserID			nId,
                                                XnFloat				/*fTime*/,
                                                void*				pCookie)
{
	printf("Lost Hand: %d\n", nId);
    
}

void XN_CALLBACK_TYPE Gesture_Recognized(	xn::GestureGenerator&	/*generator*/,
                                                      const XnChar*			strGesture,
                                                      const XnPoint3D*		pIDPosition,
                                                      const XnPoint3D*		pEndPosition,
                                                      void*					pCookie)
{
	printf("Gesture recognized: %s\n", strGesture);
    //REMOVE_ALL_GESTURES;
     gHandsGenerator.StartTracking(*pEndPosition);
}

static void XN_CALLBACK_TYPE Gesture_Process(	xn::GestureGenerator&	/*generator*/,
                                             const XnChar*			/*strGesture*/,
                                             const XnPoint3D*		/*pPosition*/,
                                             XnFloat					/*fProgress*/,
                                             void*					/*pCookie*/)	{}

void drawFunctionMain()
{
    drawFrame();
    //DisplayPostDraw(m_HandTracker);
}

int main(int argc, char* argv[])
{
    
    XnBool bChooseDevice = false;
	const char* csRecordingName = NULL;
    
	if (argc > 1)
	{
		if (strcmp(argv[1], "-devices") == 0)
		{
			bChooseDevice = TRUE;
		}
		else
		{
			csRecordingName = argv[1];
		}
	}
    
	if (csRecordingName != NULL)
	{
		// check if running from a different directory. If so, we need to change directory
		// to the real one, so that path to INI file will be OK (for log initialization, for example)
		if (0 != changeDirectory(argv[0]))
		{
			return(ERR_DEVICE);
		}
	}
    
	// Xiron Init
	XnStatus rc = XN_STATUS_OK;
	EnumerationErrors errors;
    
	if (csRecordingName != NULL)
	{
		xnLogInitFromXmlFile(SAMPLE_XML_PATH);
		rc = openDeviceFile(argv[1]);
	}
	else if (bChooseDevice)
	{
		rc = openDeviceFromXmlWithChoice(SAMPLE_XML_PATH, errors);
	}
	else
	{
		rc = openDeviceFromXml(SAMPLE_XML_PATH, errors);
	}
    
	if (rc == XN_STATUS_NO_NODE_PRESENT)
	{
		XnChar strError[1024];
		errors.ToString(strError, 1024);
		printf("%s\n", strError);
		closeSample(ERR_DEVICE);
		return (rc);
	}
	else if (rc != XN_STATUS_OK)
	{
		printf("Open failed: %s\n", xnGetStatusString(rc));
		closeSample(ERR_DEVICE);
	}
    
	captureInit();
	statisticsInit();
    
	//reshaper.zNear = 1;
	//reshaper.zFar = 100;
	//glut_add_interactor(&reshaper);
    
	//cb.mouse_function = MouseCallback;
	//cb.motion_function = MotionCallback;
	//cb.passive_motion_function = MotionCallback;
	//cb.keyboard_function = KeyboardCallback;
	//cb.reshape_function = ReshapeCallback;
	//glut_add_interactor(&cb);
    
    glutInit(&argc, argv);
	glutInitDisplayString("stencil double rgb");
	glutInitWindowSize(WIN_SIZE_X, WIN_SIZE_Y);
	glutCreateWindow("OpenNI Viewer");
	glutFullScreen();
	glutSetCursor(GLUT_CURSOR_NONE);
    
    init_opengl();
    
    glutIdleFunc(IdleCallback);
	glutDisplayFunc(drawFunctionMain);
    
	//createKeyboardMap();
	//createMenu();

    atexit(onExit);
    HandTracker mainHandTracker(g_Context);
    m_HandTracker = &mainHandTracker;
    
    drawInit(m_HandTracker);
    
    mainHandTracker.Init();
    mainHandTracker.Run();
    xn::ImageGenerator test;
    g_Context.FindExistingNode(XN_NODE_TYPE_IMAGE, test);
    xn::DepthGenerator depth;
    g_Context.FindExistingNode(XN_NODE_TYPE_DEPTH, depth);
    
    depth.GetAlternativeViewPointCap().SetViewPoint(test);
    /*rc = g_Context.FindExistingNode(XN_NODE_TYPE_HANDS,gHandsGenerator);
    bool bHandsGenerator = (rc != XN_STATUS_OK);
    if (bHandsGenerator)
    {
        rc = gHandsGenerator.Create(g_Context);
        //CHECK_RC(rc, bHandsGenerator, "Find hands generator");
        XnCallbackHandle hHandsCallbacks;
        gHandsGenerator.RegisterHandCallbacks(Hand_Create, Hand_Update, Hand_Destroy, 0, hHandsCallbacks);
        /*XnCallbackHandle h;
        if (gHandsGenerator.IsCapabilitySupported(XN_CAPABILITY_HAND_TOUCHING_FOV_EDGE))
        {
            gHandsGenerator.GetHandTouchingFOVEdgeCap().RegisterToHandTouchingFOVEdge(TouchingCallback, NULL, h);
        }*/
    //}
    // Create generators
	/*rc = gGestureGenerator.Create(g_Context);
	if (rc != XN_STATUS_OK)
	{
		printf("Unable to create GestureGenerator.");
		return rc;
	}
    
    XnCallbackHandle hHandsCallbacks;
    ADD_ALL_GESTURES
    rc = gGestureGenerator.RegisterGestureCallbacks(Gesture_Recognized, Gesture_Process, 0, hHandsCallbacks);
    if (rc != XN_STATUS_OK)
    {
        printf("Unable to register gesture callbacks.");
        return rc;
    }*/
    
    //m_HandTracker->Init(startingLeftHand, startingRightHand);
    //m_HandTracker->Run();
    
    //Take depth snapshot
    //Look for localized changes in depth
    //If you see one, use that as a starting point for hand tracking
    
    const DepthMetaData* pDepthMD = getDepthMetaData();
    const XnDepthPixel* pDepth = pDepthMD->Data();
    //Diff this against our current depth
    
    glutMainLoop();
    
	/*SimpleViewer& viewer = HandViewer::CreateInstance(g_Context);
    XnSkeletonJointTransformation leftHand;
        XnSkeletonJointTransformation rightHand;
    viewer.setHandStartingLocations(leftHand, rightHand);
	rc = viewer.Init(argc, argv);
	if (rc != XN_STATUS_OK)
	{
		printf("Viewer init failed: %s\n", xnGetStatusString(rc));
		return 1;
	}
    
	rc = viewer.Run();
	if (rc != XN_STATUS_OK)
	{
		printf("Viewer run failed: %s\n", xnGetStatusString(rc));
		return 1;
	}*/
    
	return 0;
    
}

