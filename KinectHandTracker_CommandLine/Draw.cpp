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
// --------------------------------
// Includes
// --------------------------------
#include "Draw.h"
#include "Device.h"
#include "Keyboard.h"
#include "Capture.h"
#include "NiHandTracker.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#if (XN_PLATFORM == XN_PLATFORM_MACOSX)
	#include <GLUT/glut.h>
	#include <OpenGL/gl.h>
#else
	#include <GL/gl.h>
	#include <GL/glut.h>
#endif
#include "Statistics.h"
#include "MouseInput.h"

#if (XN_PLATFORM == XN_PLATFORM_WIN32)
	#ifdef __INTEL_COMPILER
		#include <ia32intrin.h>
	#else
		#include <intrin.h>
	#endif
#endif

// --------------------------------
// Defines
// --------------------------------
#define YUV422_U  0
#define YUV422_Y1 1
#define YUV422_V  2
#define YUV422_Y2 3
#define YUV422_BPP 4
#define YUV_RED   0
#define YUV_GREEN 1
#define YUV_BLUE  2
#define YUV_ALPHA  3
#define YUV_RGBA_BPP 4
#define LENGTHOF(arr)			(sizeof(arr)/sizeof(arr[0]))


#define XDepthGrid 64
#define YDepthGrid 48
#define PI 3.141593

// --------------------------------
// Types
// --------------------------------
typedef struct  
{
	StreamsDrawConfig Streams; 
	bool bShowPointer;
	bool bShowMessage;
	bool bHelp;
	const XnChar* strErrorState;
	IntRect DepthLocation;
	IntRect ImageLocation;
} DrawConfig;

typedef struct
{
	const char* csName;
	StreamsDrawConfig Config;
} DrawConfigPreset;

typedef struct XnTextureMap
{
	IntPair Size;
	IntPair OrigSize;
	unsigned char* pMap;
	unsigned int nBytesPerPixel;
	GLuint nID;
	GLenum nFormat;
	bool bInitialized;
	IntPair CurSize;
} XnTextureMap;

// --------------------------------
// Global Variables
// --------------------------------
DrawConfig g_DrawConfig;

XnUInt8 PalletIntsR [256] = {0};
XnUInt8 PalletIntsG [256] = {0};
XnUInt8 PalletIntsB [256] = {0};

/* Linear Depth Histogram */
float* g_pDepthHist;

HandTracker* m_handTracker;

const char* g_DepthColoring[NUM_OF_DEPTH_TYPES];
const char* g_ImageColoring[NUM_OF_IMAGE_TYPES];

typedef struct DrawUserInput
{
	SelectionState State;
	IntRect SoapRect;
    IntRect TowelRect;
    int CurrentSelect;
	IntPair Cursor;
} DrawUserInput;

DrawUserInput g_DrawUserInput;

int g_nMaxDepth = 0;
bool bFirstFrame = true;
int pFirstFrame [2] = {0};
float pDepthMap [XDepthGrid][YDepthGrid] = { 0.0f };
int pDepthDiffMap [XDepthGrid][YDepthGrid] = { 0 };
std::vector<int> DepthOnSoap;
std::vector<int> DepthOnTowel;

DrawConfigPreset g_Presets[PRESET_COUNT] = 
{
	// NAME,								BACKGRD, { Depth_Type, Transparency}, { Image_Type  }			Arrangement }}
	{ "Standard Deviation",					{ false, { STANDARD_DEVIATION,	1 }, { IMAGE_OFF },				OVERLAY } },
	{ "Depth Histogram",					{ false, { LINEAR_HISTOGRAM,	1 }, { IMAGE_OFF },				OVERLAY } },
	{ "Psychedelic Depth [Centimeters]",	{ false, { PSYCHEDELIC,			1 }, { IMAGE_OFF },				OVERLAY } },
	{ "Psychedelic Depth [Millimeters]",	{ false, { PSYCHEDELIC_SHADES,	1 }, { IMAGE_OFF },				OVERLAY } },
	{ "Rainbow Depth",						{ false, { CYCLIC_RAINBOW_HISTOGRAM,1 },{ IMAGE_OFF },			OVERLAY } },
	{ "Depth masked Image",					{ false, { DEPTH_OFF,			1 }, { DEPTH_MASKED_IMAGE },	OVERLAY } },
	{ "Background Removal",					{ true,	 { DEPTH_OFF,			1 }, { DEPTH_MASKED_IMAGE },	OVERLAY } },
	{ "Side by Side",						{ false, { LINEAR_HISTOGRAM,	1 }, { IMAGE_NORMAL },			SIDE_BY_SIDE } },
	{ "Depth on Image",						{ false, { LINEAR_HISTOGRAM,	1 }, { IMAGE_NORMAL },			OVERLAY } },
	{ "Transparent Depth on Image",			{ false, { LINEAR_HISTOGRAM,  0.6 }, { IMAGE_NORMAL },			OVERLAY } },
	{ "Rainbow Depth on Image",				{ false, { RAINBOW,			  0.6 }, { IMAGE_NORMAL },			OVERLAY } },
	{ "Cyclic Rainbow Depth on Image",		{ false, { CYCLIC_RAINBOW,	  0.6 }, { IMAGE_NORMAL },			OVERLAY } },
	{ "Image Only",							{ false, { DEPTH_OFF,			1 }, { IMAGE_NORMAL },			OVERLAY } },
};

/* Texture maps for depth and image */
XnTextureMap g_texDepth = {0};
XnTextureMap g_texImage = {0};
XnTextureMap g_texBackground = {0};

/* A user message to be displayed. */
char g_csUserMessage[256];

bool g_bFullScreen = true;
bool g_bFirstTimeNonFull = true;
IntPair g_NonFullWinSize = { WIN_SIZE_X, WIN_SIZE_Y };

// --------------------------------
// Textures
// --------------------------------
int GetPowerOfTwo(int num)
{
	int result = 1;

	while (result < num)
		result <<= 1;

	return result;
}

void TextureMapInit(XnTextureMap* pTex, int nSizeX, int nSizeY, unsigned int nBytesPerPixel, int nCurX, int nCurY)
{
	// check if something changed
	if (pTex->bInitialized && pTex->OrigSize.X == nSizeX && pTex->OrigSize.Y == nSizeY)
	{
		if (pTex->CurSize.X != nCurX || pTex->CurSize.Y != nCurY)
		{
			// clear map
			xnOSMemSet(pTex->pMap, 0, pTex->Size.X * pTex->Size.Y * pTex->nBytesPerPixel);

			// update
			pTex->CurSize.X = nCurX;
			pTex->CurSize.Y = nCurY;
			return;
		}
	}

	// free memory if it was allocated
	if (pTex->pMap != NULL)
		delete pTex->pMap;

	// update it all
	pTex->OrigSize.X = nSizeX;
	pTex->OrigSize.Y = nSizeY;
	pTex->Size.X = GetPowerOfTwo(nSizeX);
	pTex->Size.Y = GetPowerOfTwo(nSizeY);
	pTex->nBytesPerPixel = nBytesPerPixel;
	pTex->CurSize.X = nCurX;
	pTex->CurSize.Y = nCurY;
	pTex->pMap = new unsigned char[pTex->Size.X * pTex->Size.Y * nBytesPerPixel];
	xnOSMemSet(pTex->pMap, 0, pTex->Size.X * pTex->Size.Y * nBytesPerPixel);
	
	if (!pTex->bInitialized)
	{
		glGenTextures(1, &(pTex->nID));
		glBindTexture(GL_TEXTURE_2D, pTex->nID);

		switch (pTex->nBytesPerPixel)
		{
		case 3:
			pTex->nFormat = GL_RGB;
			break;
		case 4:
			pTex->nFormat = GL_RGBA;
			break;
		}

		glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		pTex->bInitialized = TRUE;
	}
}

inline unsigned char* TextureMapGetLine(XnTextureMap* pTex, unsigned int nLine)
{
	return &pTex->pMap[nLine * pTex->Size.X * pTex->nBytesPerPixel];
}

void TextureMapSetPixel(XnTextureMap* pTex, int x, int y, int red, int green, int blue)
{
	if (x < 0 || y < 0 || x >= (int)pTex->OrigSize.X || y >= (int)pTex->OrigSize.Y)
		return;

	unsigned char* pPixel = TextureMapGetLine(pTex, y) + x * pTex->nBytesPerPixel;
	pPixel[0] = red;
	pPixel[1] = green;
	pPixel[2] = blue;

	if (pTex->nBytesPerPixel > 3)
		pPixel[3] = 255;
}

void TextureMapDrawCursor(XnTextureMap* pTex, IntPair cursor)
{
	// marked pixel
	TextureMapSetPixel(pTex, cursor.X, cursor.Y, 255, 0, 0);

	// top left marker
	TextureMapSetPixel(pTex, cursor.X-2, cursor.Y-2, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X-2, cursor.Y-1, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X-1, cursor.Y-2, 255, 0, 0);

	// top right marker
	TextureMapSetPixel(pTex, cursor.X+2, cursor.Y-2, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X+2, cursor.Y-1, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X+1, cursor.Y-2, 255, 0, 0);

	// bottom left marker
	TextureMapSetPixel(pTex, cursor.X-2, cursor.Y+2, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X-2, cursor.Y+1, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X-1, cursor.Y+2, 255, 0, 0);

	// bottom right marker
	TextureMapSetPixel(pTex, cursor.X+2, cursor.Y+2, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X+2, cursor.Y+1, 255, 0, 0);
	TextureMapSetPixel(pTex, cursor.X+1, cursor.Y+2, 255, 0, 0);
}

void TextureMapUpdate(XnTextureMap* pTex)
{
	// set current texture object
	glBindTexture(GL_TEXTURE_2D, pTex->nID);

	// set the current image to the texture
	glTexImage2D(GL_TEXTURE_2D, 0, pTex->nFormat, pTex->Size.X, pTex->Size.Y, 0, pTex->nFormat, GL_UNSIGNED_BYTE, pTex->pMap);
}

void TextureMapDraw(XnTextureMap* pTex, IntRect* pLocation)
{
	// set current texture object
	glBindTexture(GL_TEXTURE_2D, pTex->nID);

	// turn on texture mapping
	glEnable(GL_TEXTURE_2D);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// set drawing mode to rectangles
	glBegin(GL_QUADS);

	// set the color of the polygon
	glColor4f(1, 1, 1, 1);

	// upper left
	glTexCoord2f(0, 0);
	glVertex2f(pLocation->uLeft, pLocation->uBottom);
	// upper right
	glTexCoord2f((float)pTex->OrigSize.X/(float)pTex->Size.X, 0);
	glVertex2f(pLocation->uRight, pLocation->uBottom);
	// bottom right
	glTexCoord2f((float)pTex->OrigSize.X/(float)pTex->Size.X, (float)pTex->OrigSize.Y/(float)pTex->Size.Y);
	glVertex2f(pLocation->uRight, pLocation->uTop);
	// bottom left
	glTexCoord2f(0, (float)pTex->OrigSize.Y/(float)pTex->Size.Y);
	glVertex2f(pLocation->uLeft, pLocation->uTop);

	glEnd();

	// turn off texture mapping
	glDisable(GL_TEXTURE_2D);

	glDisable(GL_BLEND);
}

// --------------------------------
// Code
// --------------------------------

//Iterate through area of interest and apply a gaussian kernel at each point to figure out the strength of
//The points around it. Using this to filter out noise
void GaussianFilter(int pDepthDiffMap[XDepthGrid][YDepthGrid], float pOutput[XDepthGrid][YDepthGrid], int kernelSize)
{
    //Construct our kernel
    float mean = kernelSize / 2.0f;
    float sigmaSQ = kernelSize / 2.0f; // This should work
    float sigma = sqrt(sigmaSQ);
    
    float frontCoeff = 1 / (sigma * sqrt(2 * PI));
    
    std::vector< std::vector<double> > kernel;
    
    // Set up sizes. (kernelSize x kernelSize)
    kernel.resize(kernelSize);
    
    for (int i = 0; i < kernelSize; ++i)
        kernel[i].resize(kernelSize);
    
    //Generate our 2D gaussian kernel
    for (int i = 0; i < kernelSize; i++)
    {
        for (int j = 0; j < kernelSize; j++)
        {
            float dist = pow((i - mean), 2) + pow((j - mean), 2);
            
            //Calculate the strength at each point (from the center of the kernel)
            float gaussDist = frontCoeff * exp( - dist / (2 * sigmaSQ));
            
            kernel[i][j] = gaussDist;
        }
    }
    
    //Iterate through our map and apply this kernel. Return the sum of strengths at each point
    for (int i = 0; i < XDepthGrid - kernelSize; i++)
    {
        for (int j = 0; j < YDepthGrid - kernelSize; j++)
        {
            //Ok multiply our kernel by the value at each point and sum to calculate energy output
            float energyOutput = 0.0f;
            for (int kernX = 0; kernX < kernelSize; kernX++)
            {
                for (int kernY = 0; kernY < kernelSize; kernY++)
                {
                    energyOutput += kernel[kernX][kernY] * pDepthDiffMap[i + kernX][j + kernY];
                }
            }
            pOutput[i + kernelSize / 2][j + kernelSize / 2] = energyOutput;
        }
    }
}

//Calculate a linear regression based off a series of input points
void CalculateLinearRegression(float pInput[XDepthGrid][YDepthGrid], float & xSlope, float & ySlope, float & xIntercept, float & yIntercept)
{
    //Input has a 1 where there is a point in the 2D grid
    float xMean = 0;
    float yMean = 0;
    
    int numPoints = 0;
    
    //Based off of http://en.wikipedia.org/wiki/Simple_linear_regression#Fitting_the_regression_line
    for (int i = 0; i < XDepthGrid; i++)
    {
        for (int j = 0; j < YDepthGrid; j++)
        {
            if (pInput[i][j] == 1)
            {
                numPoints++;
                xMean += i;
                yMean += j;
            }
        }
    }
    
    xMean /= numPoints;
    yMean /= numPoints;
    
    float xSQSum = 0;
    float multSum = 0;
    float ySQSum = 0;
    
    for (int i = 0; i < XDepthGrid; i++)
    {
        for (int j = 0; j < YDepthGrid; j++)
        {
            if (pInput[i][j] == 1)
            {
                xSQSum += pow((i - xMean), 2);
                ySQSum += pow((j - yMean), 2);
                multSum += (i - xMean) * (j - yMean);
            }
        }
    }
    
    xSlope = multSum / xSQSum;
    ySlope = multSum / ySQSum;
    
    ySlope = 1 / ySlope;
    
    xIntercept = yMean - xSlope * xMean;
    
    yIntercept = yMean - ySlope * xMean;
    
    if (ySQSum > xSQSum)
    {
        xSlope = ySlope;
        xIntercept = yIntercept;
    }
    else
    {
        ySlope = xSlope;
        yIntercept = xIntercept;
    }
}

void CreateRainbowPallet()
{
	unsigned char r, g, b;
	for (int i=1; i<255; i++)
	{
		if (i<=29)
		{
			r = (unsigned char)(129.36-i*4.36);
			g = 0;
			b = (unsigned char)255;
		}
		else if (i<=86)
		{
			r = 0;
			g = (unsigned char)(-133.54+i*4.52);
			b = (unsigned char)255;
		}
		else if (i<=141)
		{
			r = 0;
			g = (unsigned char)255;
			b = (unsigned char)(665.83-i*4.72);
		}
		else if (i<=199)
		{
			r = (unsigned char)(-635.26+i*4.47);
			g = (unsigned char)255;
			b = 0;
		}
		else
		{
			r = (unsigned char)255;
			g = (unsigned char)(1166.81-i*4.57);
			b = 0;
		}

		PalletIntsR[i] = r;
		PalletIntsG[i] = g;
		PalletIntsB[i] = b;
	}
}

void glPrintString(void *font, const char *str)
{
	int i,l = (int)strlen(str);

	for(i=0; i<l; i++)
	{
		glutBitmapCharacter(font,*str++);
	}
}

void drawConfigChanged()
{
	// recalculate registration
	bool bRegistration = 
		(g_DrawConfig.Streams.ScreenArrangement == OVERLAY) && 
		(g_DrawConfig.Streams.Image.Coloring != IMAGE_OFF) &&
		(g_DrawConfig.Streams.Depth.Coloring != DEPTH_OFF || g_DrawConfig.Streams.Image.Coloring == DEPTH_MASKED_IMAGE);

	changeRegistration(bRegistration);
}

void setPreset(int preset)
{
	g_DrawConfig.Streams = g_Presets[preset].Config;
	drawConfigChanged();
}

const char* getPresetName(int preset)
{
	return g_Presets[preset].csName;
}

void setScreenLayout(int layout)
{
	g_DrawConfig.Streams.ScreenArrangement = (ScreenArrangementType)layout;
	drawConfigChanged();
}

void windowReshaped(int width, int height)
{
	g_NonFullWinSize.X = width;
	g_NonFullWinSize.Y = height;
}

void toggleFullScreen(int)
{
	if (g_bFullScreen)
	{
		if (g_bFirstTimeNonFull)
		{
			g_NonFullWinSize.X = WIN_SIZE_X/2;
			g_NonFullWinSize.Y = WIN_SIZE_Y/2;
			g_bFirstTimeNonFull = false;
		}

		glutReshapeWindow(g_NonFullWinSize.X, g_NonFullWinSize.Y);
		g_bFullScreen = false;
	}
	else
	{
		glutFullScreen();
		g_bFullScreen = true;
	}
}

void displayMessage(const char* csFormat, ...)
{
	g_DrawConfig.bShowMessage = true;
	va_list args;
	va_start(args, csFormat);
	vsprintf(g_csUserMessage, csFormat, args);
	va_end(args);
}

void setErrorState(const char* strMessage)
{
	g_DrawConfig.strErrorState = strMessage;
}

void drawCropStream(MapGenerator* pGenerator, IntRect location, IntRect selection, int dividedBy)
{
	if (!pGenerator->IsCapabilitySupported(XN_CAPABILITY_CROPPING))
	{
		return;
	}

	XnMapOutputMode Mode;
	pGenerator->GetMapOutputMode(Mode);

	// check if entire selection is in location
	if (selection.uLeft >= location.uLeft &&
		selection.uRight <= location.uRight &&
		selection.uBottom >= location.uBottom &&
		selection.uTop <= location.uTop)
	{
		IntRect cropRect;
		cropRect.uBottom = Mode.nYRes * (selection.uBottom - location.uBottom) / (location.uTop - location.uBottom);
		cropRect.uTop = Mode.nYRes * (selection.uTop - location.uBottom) / (location.uTop - location.uBottom);
		cropRect.uLeft = Mode.nXRes * (selection.uLeft - location.uLeft) / (location.uRight - location.uLeft);
		cropRect.uRight = Mode.nXRes * (selection.uRight - location.uLeft) / (location.uRight - location.uLeft);

		XnCropping cropping;
		cropping.bEnabled = TRUE;
		cropping.nXOffset = cropRect.uLeft;
		cropping.nYOffset = cropRect.uBottom;
		cropping.nXSize = cropRect.uRight - cropRect.uLeft;
		cropping.nYSize = cropRect.uTop	- cropRect.uBottom;

		if ((cropping.nXOffset % dividedBy) != 0)
			cropping.nXOffset -= (cropping.nXOffset % dividedBy);
		if ((cropping.nXSize % dividedBy) != 0)
			cropping.nXSize += dividedBy - (cropping.nXSize % dividedBy);

		setStreamCropping(pGenerator, &cropping);
	}
}

void drawSelectionChanged(SelectionState state, IntRect SoapSelection, IntRect TowelSelection, int iActiveSelection)
{
	g_DrawUserInput.State = state;
	g_DrawUserInput.SoapRect = SoapSelection;
    g_DrawUserInput.TowelRect = TowelSelection;
    g_DrawUserInput.CurrentSelect = iActiveSelection;
    
    bFirstFrame = true;

	/*if (state == SELECTION_DONE)
	{
		// Crop depth
		if (getDepthGenerator() != NULL && g_DrawConfig.Streams.Depth.Coloring != DEPTH_OFF)
		{
			drawCropStream(getDepthGenerator(), g_DrawConfig.DepthLocation, selection, 2);
		}

		// Crop image
		if (getImageGenerator() != NULL && g_DrawConfig.Streams.Image.Coloring != IMAGE_OFF)
		{
			drawCropStream(getImageGenerator(), g_DrawConfig.ImageLocation, selection, 4);
		}

		// Crop IR
		if (getIRGenerator() != NULL && g_DrawConfig.Streams.Image.Coloring != IMAGE_OFF)
		{
			drawCropStream(getIRGenerator(), g_DrawConfig.ImageLocation, selection, 4);
		}
	}*/
}

void drawCursorMoved(IntPair location)
{
	g_DrawUserInput.Cursor = location;
}

void drawInit(HandTracker* handTracker)
{
	g_DepthColoring[DEPTH_OFF] = "Off";
	g_DepthColoring[LINEAR_HISTOGRAM] = "Linear Histogram";
	g_DepthColoring[PSYCHEDELIC] = "Psychedelic";
	g_DepthColoring[PSYCHEDELIC_SHADES] = "Psychedelic (Millimeters)";
	g_DepthColoring[RAINBOW] = "Rainbow";
	g_DepthColoring[CYCLIC_RAINBOW] = "Cyclic Rainbow";
	g_DepthColoring[CYCLIC_RAINBOW_HISTOGRAM] = "Cyclic Rainbow Histogram";
	g_DepthColoring[STANDARD_DEVIATION] = "Standard Deviation";

	g_ImageColoring[IMAGE_OFF] = "Off";
	g_ImageColoring[IMAGE_NORMAL] = "Normal";
	g_ImageColoring[DEPTH_MASKED_IMAGE] = "Depth Masked Image";

	CreateRainbowPallet();

	setPreset(7);

	TextureMapInit(&g_texBackground, 1024, 1024, 3, 1024, 1024);

	// load background image
	xnOSLoadFile("..\\..\\..\\..\\Data\\RGBViewer\\back.raw", TextureMapGetLine(&g_texBackground, 0), 1024*1024*3);

	TextureMapUpdate(&g_texBackground);

	mouseInputRegisterForSelectionRectangle(drawSelectionChanged);
	mouseInputRegisterForCursorMovement(drawCursorMoved);
    
    m_handTracker = handTracker;
}

void togglePointerMode(int)
{
	g_DrawConfig.bShowPointer = !g_DrawConfig.bShowPointer;
}

void toggleHelpScreen(int)
{
	g_DrawConfig.bHelp = !g_DrawConfig.bHelp;
}

void toggleBackground(int)
{
	g_DrawConfig.Streams.bBackground = !g_DrawConfig.Streams.bBackground;
}

void calculateHistogram()
{
	DepthGenerator* pDepthGen = getDepthGenerator();

	if (pDepthGen == NULL)
		return;

	XnUInt32 nZRes = pDepthGen->GetDeviceMaxDepth() + 1;
	if (g_pDepthHist == NULL)
	{
		g_pDepthHist = new float[nZRes];
	}

	xnOSMemSet(g_pDepthHist, 0, nZRes*sizeof(float));
	int nNumberOfPoints = 0;

	XnDepthPixel nValue;

	const XnDepthPixel* pDepth = pDepthGen->GetDepthMap();
	const XnDepthPixel* pDepthEnd = pDepth + (pDepthGen->GetDataSize() / sizeof(XnDepthPixel));

	while (pDepth != pDepthEnd)
	{
		nValue = *pDepth;

		XN_ASSERT(nValue <= nZRes);

		if (nValue != 0)
		{
			g_pDepthHist[nValue]++;
			nNumberOfPoints++;
		}

		pDepth++;
	}

	XnUInt32 nIndex;
	for (nIndex = 1; nIndex < nZRes; nIndex++)
	{
		g_pDepthHist[nIndex] += g_pDepthHist[nIndex-1];
	}
	for (nIndex = 1; nIndex < nZRes; nIndex++)
	{
		if (g_pDepthHist[nIndex] != 0)
		{
			g_pDepthHist[nIndex] = (nNumberOfPoints-g_pDepthHist[nIndex]) / nNumberOfPoints;
		}
	}
}

// --------------------------------
// Drawing
// --------------------------------
#if (XN_PLATFORM == XN_PLATFORM_WIN32)

void YUV422ToRGB888(const XnUInt8* pYUVImage, XnUInt8* pRGBAImage, XnUInt32 nYUVSize, XnUInt32 nRGBSize)
{
	const XnUInt8* pYUVLast = pYUVImage + nYUVSize - 8;
	XnUInt8* pRGBLast = pRGBAImage + nRGBSize - 16;

	const __m128 minus128 = _mm_set_ps1(-128);
	const __m128 plus113983 = _mm_set_ps1(1.13983F);
	const __m128 minus039466 = _mm_set_ps1(-0.39466F);
	const __m128 minus058060 = _mm_set_ps1(-0.58060F);
	const __m128 plus203211 = _mm_set_ps1(2.03211F);
	const __m128 zero = _mm_set_ps1(0);
	const __m128 plus255 = _mm_set_ps1(255);

	// define YUV floats
	__m128 y;
	__m128 u;
	__m128 v;

	__m128 temp;

	// define RGB floats
	__m128 r;
	__m128 g;
	__m128 b;

	// define RGB integers
	__m128i iR;
	__m128i iG;
	__m128i iB;

	XnUInt32* piR = (XnUInt32*)&iR;
	XnUInt32* piG = (XnUInt32*)&iG;
	XnUInt32* piB = (XnUInt32*)&iB;

	while (pYUVImage <= pYUVLast && pRGBAImage <= pRGBLast)
	{
		// process 4 pixels at once (values should be ordered backwards)
		y = _mm_set_ps(pYUVImage[YUV422_Y2 + YUV422_BPP], pYUVImage[YUV422_Y1 + YUV422_BPP], pYUVImage[YUV422_Y2], pYUVImage[YUV422_Y1]);
		u = _mm_set_ps(pYUVImage[YUV422_U + YUV422_BPP],  pYUVImage[YUV422_U + YUV422_BPP],  pYUVImage[YUV422_U],  pYUVImage[YUV422_U]);
		v = _mm_set_ps(pYUVImage[YUV422_V + YUV422_BPP],  pYUVImage[YUV422_V + YUV422_BPP],  pYUVImage[YUV422_V],  pYUVImage[YUV422_V]);

		u = _mm_add_ps(u, minus128); // u -= 128
		v = _mm_add_ps(v, minus128); // v -= 128

		/*

		http://en.wikipedia.org/wiki/YUV

		From YUV to RGB:
		R =     Y + 1.13983 V
		G =     Y - 0.39466 U - 0.58060 V
		B =     Y + 2.03211 U

		*/ 

		temp = _mm_mul_ps(plus113983, v);
		r = _mm_add_ps(y, temp);

		temp = _mm_mul_ps(minus039466, u);
		g = _mm_add_ps(y, temp);
		temp = _mm_mul_ps(minus058060, v);
		g = _mm_add_ps(g, temp);

		temp = _mm_mul_ps(plus203211, u);
		b = _mm_add_ps(y, temp);

		// make sure no value is smaller than 0
		r = _mm_max_ps(r, zero);
		g = _mm_max_ps(g, zero);
		b = _mm_max_ps(b, zero);

		// make sure no value is bigger than 255
		r = _mm_min_ps(r, plus255);
		g = _mm_min_ps(g, plus255);
		b = _mm_min_ps(b, plus255);

		// convert floats to int16 (there is no conversion to uint8, just to int8).
		iR = _mm_cvtps_epi32(r);
		iG = _mm_cvtps_epi32(g);
		iB = _mm_cvtps_epi32(b);

		// extract the 4 pixels RGB values.
		// because we made sure values are between 0 and 255, we can just take the lower byte
		// of each INT16
		pRGBAImage[0] = piR[0];
		pRGBAImage[1] = piG[0];
		pRGBAImage[2] = piB[0];
		pRGBAImage[3] = 255;

		pRGBAImage[4] = piR[1];
		pRGBAImage[5] = piG[1];
		pRGBAImage[6] = piB[1];
		pRGBAImage[7] = 255;

		pRGBAImage[8] = piR[2];
		pRGBAImage[9] = piG[2];
		pRGBAImage[10] = piB[2];
		pRGBAImage[11] = 255;

		pRGBAImage[12] = piR[3];
		pRGBAImage[13] = piG[3];
		pRGBAImage[14] = piB[3];
		pRGBAImage[15] = 255;

		// advance the streams
		pYUVImage += 8;
		pRGBAImage += 16;
	}
}

#else // not Win32

void YUV444ToRGBA(XnUInt8 cY, XnUInt8 cU, XnUInt8 cV,
					XnUInt8& cR, XnUInt8& cG, XnUInt8& cB, XnUInt8& cA)
{
	XnInt32 nC = cY - 16;
	XnInt16 nD = cU - 128;
	XnInt16 nE = cV - 128;

	nC = nC * 298 + 128;

	cR = XN_MIN(XN_MAX((nC            + 409 * nE) >> 8, 0), 255);
	cG = XN_MIN(XN_MAX((nC - 100 * nD - 208 * nE) >> 8, 0), 255);
	cB = XN_MIN(XN_MAX((nC + 516 * nD           ) >> 8, 0), 255);
	cA = 255;
}

void YUV422ToRGB888(const XnUInt8* pYUVImage, XnUInt8* pRGBImage, XnUInt32 nYUVSize, XnUInt32 nRGBSize)
{
	const XnUInt8* pCurrYUV = pYUVImage;
	XnUInt8* pCurrRGB = pRGBImage;
	const XnUInt8* pLastYUV = pYUVImage + nYUVSize - YUV422_BPP;
	XnUInt8* pLastRGB = pRGBImage + nRGBSize - YUV_RGBA_BPP;

	while (pCurrYUV <= pLastYUV && pCurrRGB <= pLastRGB)
	{
		YUV444ToRGBA(pCurrYUV[YUV422_Y1], pCurrYUV[YUV422_U], pCurrYUV[YUV422_V],
						pCurrRGB[YUV_RED], pCurrRGB[YUV_GREEN], pCurrRGB[YUV_BLUE], pCurrRGB[YUV_ALPHA]);
		pCurrRGB += YUV_RGBA_BPP;
		YUV444ToRGBA(pCurrYUV[YUV422_Y2], pCurrYUV[YUV422_U], pCurrYUV[YUV422_V],
						pCurrRGB[YUV_RED], pCurrRGB[YUV_GREEN], pCurrRGB[YUV_BLUE], pCurrRGB[YUV_ALPHA]);
		pCurrRGB += YUV_RGBA_BPP;
		pCurrYUV += YUV422_BPP;
	}
}

#endif

void drawClosedStream(IntRect* pLocation, const char* csStreamName)
{
	char csMessage[512];
	sprintf(csMessage, "%s stream is OFF", csStreamName);
	void* pFont = GLUT_BITMAP_TIMES_ROMAN_24;

	int nWidth = glutBitmapLength(pFont, (const unsigned char*)csMessage);
	int nXLocation = (pLocation->uRight + pLocation->uLeft - nWidth) / 2;
	int nYLocation = (pLocation->uTop + pLocation->uBottom) / 2;

	glColor3f(1.0, 0, 0);
	glRasterPos2i(nXLocation, nYLocation);
	glPrintString(pFont, csMessage);
}

void drawColorImage(IntRect* pLocation, IntPair* pPointer)
{
	if (g_DrawConfig.Streams.bBackground)
		TextureMapDraw(&g_texBackground, pLocation);

	if (g_DrawConfig.Streams.Image.Coloring == IMAGE_OFF)
		return;

	if (!isImageOn() && !isIROn())
	{
		drawClosedStream(pLocation, "Image");
		return;
	}

	const MapMetaData* pImageMD;
	const XnUInt8* pImage = NULL;

	if (isImageOn())
	{
		pImageMD = getImageMetaData();
		pImage = getImageMetaData()->Data();
	}
	else if (isIROn())
	{
		pImageMD = getIRMetaData();
		pImage = (const XnUInt8*)getIRMetaData()->Data();
	}
	else
		return;

	if (pImageMD->FrameID() == 0)
	{
		return;
	}

	const DepthMetaData* pDepthMetaData = getDepthMetaData();

	for (XnUInt16 nY = pImageMD->YOffset(); nY < pImageMD->YRes() + pImageMD->YOffset(); nY++)
	{
		XnUInt8* pTexture = TextureMapGetLine(&g_texImage, nY) + pImageMD->XOffset()*4;

		if (pImageMD->PixelFormat() == XN_PIXEL_FORMAT_YUV422)
		{
			YUV422ToRGB888(pImage, pTexture, pImageMD->XRes()*2, g_texImage.Size.X*g_texImage.nBytesPerPixel);
			pImage += pImageMD->XRes()*2;
		}
		else
		{
			for (XnUInt16 nX = 0; nX < pImageMD->XRes(); nX++, pTexture+=4)
			{
				XnInt32 nDepthIndex = 0;

				if (pDepthMetaData != NULL)
				{
					XnDouble dRealX = (nX + pImageMD->XOffset()) / (XnDouble)pImageMD->FullXRes();
					XnDouble dRealY = nY / (XnDouble)pImageMD->FullYRes();

					XnUInt32 nDepthX = dRealX * pDepthMetaData->FullXRes() - pDepthMetaData->XOffset();
					XnUInt32 nDepthY = dRealY * pDepthMetaData->FullYRes() - pDepthMetaData->YOffset();

					if (nDepthX >= pDepthMetaData->XRes() || nDepthY >= pDepthMetaData->YRes())
					{
						nDepthIndex = -1;
					}
					else
					{
						nDepthIndex = nDepthY*pDepthMetaData->XRes() + nDepthX;
					}
				}

				switch (pImageMD->PixelFormat())
				{
				case XN_PIXEL_FORMAT_RGB24:
					pTexture[0] = pImage[0];
					pTexture[1] = pImage[1];
					pTexture[2] = pImage[2];
					pImage+=3; 
					break;
				case XN_PIXEL_FORMAT_GRAYSCALE_8_BIT:
					pTexture[0] = pTexture[1] = pTexture[2] = *pImage;
					pImage+=1; 
					break;
				case XN_PIXEL_FORMAT_GRAYSCALE_16_BIT:
					XnUInt16* p16 = (XnUInt16*)pImage;
					pTexture[0] = pTexture[1] = pTexture[2] = (*p16) >> 2;
					pImage+=2; 
					break;
				}

				// decide if pixel should be lit or not
				if (g_DrawConfig.Streams.Image.Coloring == DEPTH_MASKED_IMAGE &&
					(pDepthMetaData == NULL || nDepthIndex == -1 || pDepthMetaData->Data()[nDepthIndex] == 0))
				{
					pTexture[3] = 0;
				}
				else
				{
					pTexture[3] = 255;
				}
			}
		}
	}

	if (pPointer != NULL)
	{
		TextureMapDrawCursor(&g_texImage, *pPointer);
	}

	TextureMapUpdate(&g_texImage);
	TextureMapDraw(&g_texImage, pLocation);
}


//Most of the drawing, filter code is in here
//This is called whenever a depth frame is processed
void drawDepth(IntRect* pLocation, IntPair* pPointer)
{
	if (g_DrawConfig.Streams.Depth.Coloring != DEPTH_OFF)
	{
		if (!isDepthOn())
		{
			drawClosedStream(pLocation, "Depth");
			return;
		}

		const DepthMetaData* pDepthMD = getDepthMetaData();
		const XnDepthPixel* pDepth = pDepthMD->Data();
		XN_ASSERT(pDepth);
        
        std::vector<int> diffDepthPoints;
        
		if (pDepthMD->FrameID() == 0)
		{
			return;
		}

		if (g_DrawConfig.Streams.Depth.Coloring == STANDARD_DEVIATION)
		{
			XnPixelStatistics* pStatistics = g_PixelStatistics;

			//for (XnUInt16 nY = pDepthMD->YOffset(); nY < pDepthMD->YRes() + pDepthMD->YOffset(); nY++)
            for (XnUInt16 nY = 0; nY < pDepthMD->YRes() ; nY++)
			{
				XnUInt8* pTexture = TextureMapGetLine(&g_texDepth, nY) + pDepthMD->XOffset()*4;
				for (XnUInt16 nX = 0; nX < pDepthMD->XRes(); nX++, pTexture+=4, pStatistics++)
				{
					pTexture[0] = pTexture[1] = XN_MIN((int)pStatistics->dStdDev, 255);
					pTexture[2] = 0;
					pTexture[3] = g_DrawConfig.Streams.Depth.fTransparency*255;
				}
			}
		}
		else
		{
			// copy depth into texture-map
			for (XnUInt16 nY = pDepthMD->YOffset(); nY < pDepthMD->YRes() + pDepthMD->YOffset(); nY++)
			{
				XnUInt8* pTexture = TextureMapGetLine(&g_texDepth, nY) + pDepthMD->XOffset()*4;
				for (XnUInt16 nX = 0; nX < pDepthMD->XRes(); nX++, pDepth++, pTexture+=4)
				{
					XnUInt8 nRed = 0;
					XnUInt8 nGreen = 0;
					XnUInt8 nBlue = 0;
					XnUInt8 nAlpha = g_DrawConfig.Streams.Depth.fTransparency*255;

					XnUInt16 nColIndex;
                    //g_DrawConfig.Streams.Depth.Coloring = CYCLIC_RAINBOW_HISTOGRAM;
					switch (g_DrawConfig.Streams.Depth.Coloring)
					{
					case LINEAR_HISTOGRAM:
						nRed = nGreen = g_pDepthHist[*pDepth]*255;
						break;
					case PSYCHEDELIC_SHADES:
						nAlpha *= (((XnFloat)(*pDepth % 10) / 20) + 0.5);
					case PSYCHEDELIC:

						switch ((*pDepth/10) % 10)
						{
						case 0:
							nRed = 255;
							break;
						case 1:
							nGreen = 255;
							break;
						case 2:
							nBlue = 255;
							break;
						case 3:
							nRed = 255;
							nGreen = 255;
							break;
						case 4:
							nGreen = 255;
							nBlue = 255;
							break;
						case 5:
							nRed = 255;
							nBlue = 255;
							break;
						case 6:
							nRed = 255;
							nGreen = 255;
							nBlue = 255;
							break;
						case 7:
							nRed = 127;
							nBlue = 255;
							break;
						case 8:
							nRed = 255;
							nBlue = 127;
							break;
						case 9:
							nRed = 127;
							nGreen = 255;
							break;
						}
						break;
					case RAINBOW:
						nColIndex = (XnUInt16)((*pDepth / (g_nMaxDepth / 256.)));
						nRed = PalletIntsR[nColIndex];
						nGreen = PalletIntsG[nColIndex];
						nBlue = PalletIntsB[nColIndex];
						break;
					case CYCLIC_RAINBOW:
						nColIndex = (*pDepth % 256);
						nRed = PalletIntsR[nColIndex];
						nGreen = PalletIntsG[nColIndex];
						nBlue = PalletIntsB[nColIndex];
						break;
					case CYCLIC_RAINBOW_HISTOGRAM:
						float fHist = g_pDepthHist[*pDepth];
						nColIndex = (*pDepth % 256);
						nRed = PalletIntsR[nColIndex]   * fHist;
						nGreen = PalletIntsG[nColIndex] * fHist;
						nBlue = PalletIntsB[nColIndex]  * fHist;
						break;
					}

					pTexture[0] = nRed;
					pTexture[1] = nGreen;
					pTexture[2] = nBlue;

					if (*pDepth == 0)
						pTexture[3] = 0;
					else
						pTexture[3] = nAlpha;
				}
			}
		} // not STANDRARD_DEVIATION

		if (pPointer != NULL)
		{
			TextureMapDrawCursor(&g_texDepth, *pPointer);
		}

		TextureMapUpdate(&g_texDepth);
		TextureMapDraw(&g_texDepth, pLocation);
        
        /*glBegin(GL_LINES);
        glLineWidth(3.0);
        glColor3f( 1.0, 1.0, 1.0);
        for (int i = 0; i < diffDepthPoints.size(); i++)
        {
            int xCoord = (double)diffDepthPoints[i] / pDepthMD->YRes();
            int yCoord = diffDepthPoints[i] % pDepthMD->YRes();
            glVertex2f( xCoord, yCoord );
        }
        glEnd();*/
        //If we don't have an initial depth reading, copy it in
         //else
        //
         {
            //They have selected a region of interest
            if (g_DrawUserInput.State == SELECTION_DONE)
            {
                if (bFirstFrame)
                {
                    //Take a depth snapshot
                    //see if our depth has changed?
                    float averageDepth = 0;
                    int numPointsReadIn = 0;
                    
                    DepthGenerator* pDepthGen = getDepthGenerator();
                    
                    if (pDepthGen == NULL)
                        return;
                    
                    const XnDepthPixel* pDepthCopy = pDepthGen->GetDepthMap();
                    
                    //Do a depth snapshot of what we just read in (soap or towel)
                    numPointsReadIn = 0;
                    averageDepth = 0;
                    IntRect currRect;
                    //What rectangle are we working with right now
                    if (g_DrawUserInput.CurrentSelect == SOAP)
                        currRect = g_DrawUserInput.SoapRect;
                    else
                        currRect = g_DrawUserInput.TowelRect;
                    
                    //iterate over the selected rectangle and store the average depth over the rectangle
                    for (int i = currRect.uLeft; i< currRect.uRight; i++)
                    {
                        for (int j = currRect.uBottom; j < currRect.uTop; j++)
                        {
                            numPointsReadIn++;
                            int newDepth = *(pDepthCopy + pDepthMD->XRes() * j + i);
                            
                            averageDepth += (float)newDepth;
                        }
                    }
                    
                    averageDepth = averageDepth / numPointsReadIn;
                    pFirstFrame[g_DrawUserInput.CurrentSelect] = averageDepth;
                    bFirstFrame = false;
                    
                    //XDepthGrid stores the number of x points we are collecting for depth (breaking down
                    //the depth image into a xdepthgrid x ydepthgrid granulated region for efficiency)
                    
                    //640x480 is hard coded resolution of the depth image coming from the kinect
                    int width = 640 / XDepthGrid;
                    int height = 480 / YDepthGrid;
                    
                    for (int i = 0; i < XDepthGrid; i++)
                    {
                        for (int j=0 ; j < YDepthGrid; j++)
                        {
                            numPointsReadIn = 0;
                            averageDepth = 0;
                            //Store average depth across xdepthgrid x ydepthgrid
                            for (int k = i * width; k < (i+1) * width; k++)
                            {
                                for (int l = j * height; l < (j+1) * height; l++)
                                {
                                    numPointsReadIn++;
                                    int newDepth = *(pDepthCopy + pDepthMD->XRes() * l + k);
                                    
                                    averageDepth += (float)newDepth;
                                }
                            }
                            
                            averageDepth /= numPointsReadIn;
                            pDepthMap[i][j] = averageDepth;
                        }
                    }
                    
                }
                else //not the first frame (where we collect a snapshot)
                {
                    //see if our depth has changed?
                    float averageDepth = 0;
                    int numPointsReadIn = 0;
                    
                    DepthGenerator* pDepthGen = getDepthGenerator();
                    
                    if (pDepthGen == NULL)
                        return;
                    
                    const XnDepthPixel* pDepthCopy = pDepthGen->GetDepthMap();
                    
                    numPointsReadIn = 0;
                    averageDepth = 0;
                    
                    IntRect currRect;
                    for (int select = 0; select < 2; select++)
                    {
                        numPointsReadIn = 0;
                        averageDepth = 0;
                        
                        if (select == 0)
                            currRect = g_DrawUserInput.SoapRect;
                        else
                            currRect = g_DrawUserInput.TowelRect;
                        
                        //Collect the average depth across the rectangle
                        for (int i = currRect.uLeft; i< currRect.uRight; i++)
                        {
                            for (int j = currRect.uBottom; j < currRect.uTop; j++)
                            {
                                numPointsReadIn++;
                                int newDepth = *(pDepthCopy + pDepthMD->XRes() * j + i);
                                
                                averageDepth += (float)newDepth;
                            }
                        }
                        averageDepth = averageDepth / numPointsReadIn;
                        
                        int diff = averageDepth - pFirstFrame[select];
                        
                        //Has the depth changed?
                        if (select == 0)
                        {
                            DepthOnSoap.push_back(diff);
                            
                            if (DepthOnSoap.size() > 25)
                                DepthOnSoap.erase(DepthOnSoap.begin());
                        }
                        else
                        {
                            DepthOnTowel.push_back(diff);
                            
                            if (DepthOnTowel.size() > 25)
                                DepthOnTowel.erase(DepthOnTowel.begin());
                        }
                    }
                }
                //Put "hand above soap" message 
                /*char buf[512] = "";
                if (bSelectedRegionDepthChange)
                {
                    sprintf(buf, "Hand at selected region? YES");
                }
                else
                {
                    sprintf(buf, "Hand at selected region? NO");
                }
                int nYLocation = WIN_SIZE_Y - 60;
                int nXLocation =  10;
                glColor3f(1,0,0);
                glRasterPos2i(nXLocation,nYLocation);
                glPrintString(GLUT_BITMAP_HELVETICA_12, buf);*/
                
                //Plot the soap depth over time as a line graph (to show whether the soap rectangle is changing)
                char buf[512] = "";
                int nYLocation = WIN_SIZE_Y - 60;
                int nXLocation =  10;
                sprintf(buf, "Soap Depth Over Time: ");
                //nYLocation += 20;
                nXLocation =  WIN_SIZE_X - 350;
                glColor3f(1,0,0);
                glRasterPos2i(nXLocation,nYLocation);
                glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                
                glBegin(GL_LINE_LOOP);
                for (int i = 0; i < DepthOnSoap.size(); i++)
                {
                    glVertex2f(WIN_SIZE_X - 200 + i*5, nYLocation + DepthOnSoap[i]);
                }
                glEnd();
                if (DepthOnSoap.size() > 0)
                {
                    //The depth change is significant? Assume hand is over soap
                    if (DepthOnSoap[DepthOnSoap.size() - 1] < -5)
                    {
                        //Hand overtop of soap
                        sprintf(buf, "Hand over soap");
                        nYLocation = WIN_SIZE_Y - 80;
                        glColor3f(1,0,0);
                        glRasterPos2i(nXLocation,nYLocation);
                        glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                    }
                    else if (DepthOnSoap[DepthOnSoap.size() -1] > 5)
                    {
                        //Depth decreased, so must have dispensed some soap
                        sprintf(buf, "Hand depressed soap");
                        nYLocation = WIN_SIZE_Y - 80;
                        glColor3f(1,0,0);
                        glRasterPos2i(nXLocation,nYLocation);
                        glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                    }
                }
                
                sprintf(buf, "Towel Depth Over Time: ");
                nYLocation = WIN_SIZE_Y - 60;
                nXLocation =  WIN_SIZE_X - 850;
                glColor3f(1,0,0);
                glRasterPos2i(nXLocation,nYLocation);
                glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                
                glBegin(GL_LINE_LOOP);
                for (int i = 0; i < DepthOnTowel.size(); i++)
                {
                    glVertex2f(WIN_SIZE_X - 700 + i*5, nYLocation + DepthOnTowel[i]);
                }
                glEnd();
                
                if (DepthOnTowel.size() > 0)
                {
                    if (DepthOnTowel[DepthOnSoap.size() - 1] < -5)
                    {
                        //Hand overtop of soap
                        sprintf(buf, "Hand over towel");
                        nYLocation = WIN_SIZE_Y - 80;
                        glColor3f(1,0,0);
                        glRasterPos2i(nXLocation,nYLocation);
                        glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                    }
                    else if (DepthOnTowel[DepthOnSoap.size() -1] > 5)
                    {
                        //Depth decreased, so must have dispensed some soap
                        sprintf(buf, "Towel moved");
                        nYLocation = WIN_SIZE_Y - 80;
                        glColor3f(1,0,0);
                        glRasterPos2i(nXLocation,nYLocation);
                        glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                    }
                }
                
                //Not the first frame
                if (!bFirstFrame)
                {
                    int numPointsReadIn = 0;
                    float averageDepth = 0;
                    
                    DepthGenerator* pDepthGen = getDepthGenerator();
                    
                    if (pDepthGen == NULL)
                        return;
                    
                    const XnDepthPixel* pDepthCopy = pDepthGen->GetDepthMap();
                    
                    int minLeftWallDist = 64;
                    int minRightWallDist = 64;
                    int minTopWallDist = 64;
                    int minBottomWallDist = 64;
                    
                    int width = 640 / XDepthGrid;
                    int height = 480 / YDepthGrid;
                    
                    for (int i = 0; i < XDepthGrid; i++)
                    {
                        for (int j=0 ; j < YDepthGrid; j++)
                        {
                            pDepthDiffMap[i][j] = 0;
                            
                            numPointsReadIn = 0;
                            averageDepth = 0;
                            for (int k = i * width; k < (i+1) * width; k++)
                            {
                                for (int l = j * height; l < (j+1) * height; l++)
                                {
                                    numPointsReadIn++;
                                    int newDepth = *(pDepthCopy + pDepthMD->XRes() * l + k);
                                    
                                    averageDepth += (float)newDepth;
                                }
                            }
                            averageDepth /= numPointsReadIn;
                            //1000 is noise (too far for kinect to reliably pick up, so don't draw it)
                            //See if the depth has changed? If so, likely a body part (ie. arm)
                            if ((abs(averageDepth - pDepthMap[i][j]) > 50) && (averageDepth < 1000))
                            {
								glBegin(GL_LINE_LOOP);
                                {
                                    glVertex2i(i * width, j * height);
                                    glVertex2i(i * width, (j+1) * height);
                                    glVertex2i( (i+1) * width, (j+1) * height);
                                    glVertex2i( (i+1) * width, j * height);
                                    glVertex2i( i * width, j * height);
                                }
                                glEnd();
                                
                                pDepthDiffMap[i][j] = 1;
                                
                                /*char buf[512] = "";
                                int nYLocation = (j+1) * height - 20;
                                int nXLocation =  i * width + 10;
                                sprintf(buf, "D: %3.0f", averageDepth);
                                glColor3f(1,0,0);
                                glRasterPos2i(nXLocation,nYLocation);
                                glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
                                
                                nYLocation += 20;
                                sprintf(buf, "A: %3.0f", pDepthMap[i][j]);
                                glColor3f(1,0,0);
                                glRasterPos2i(nXLocation,nYLocation);
                                glPrintString(GLUT_BITMAP_HELVETICA_12, buf);*/
                            }
                        }
                    }
                    
                    //Now we have this pDepthDiffMap holding all the indices where
                    //There was a difference in Depth. We want to filter this to identify the hand center
                    //Use a window filter with a variable sized kernel
                    
                    int kernelSize = 8;
                    float pIntensityMap[XDepthGrid][YDepthGrid] = {0};
                    GaussianFilter(pDepthDiffMap, pIntensityMap, kernelSize);
                    
					float pHighIntensityPts[XDepthGrid][YDepthGrid] = { 0 };
					std::vector< IntPair > highIntensityVector;
                    
                    //We need to find out what side of the rectangular grid we are coming from. Calculate min distance from each wall and proceed that way.
                    
                    for (int i = 0; i < XDepthGrid; i++)
                    {
                        for (int j = 0; j < YDepthGrid; j++)
                        {
                            glEnable(GL_BLEND);
                            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                            
                            //Draw in shades of blue the depth values that have changed, and
                            //change the alpha value depending on the value of the difference
                            glBegin(GL_QUADS);
                                glColor4f(0, pIntensityMap[i][j] / 5, pIntensityMap[i][j] / 5, pIntensityMap[i][j] / 5);
                                glVertex2i( i * width , j * height);
                                glVertex2i( i * width , (j+1) * height);
                                glVertex2i( (i+1) * width , (j+1) * height);
                                glVertex2i( (i+1) * width , j * height);
                            glEnd();
                            
                            glDisable(GL_BLEND);
                            
                            if ((pIntensityMap[i][j] / 5) > 0.5)
                            {
                                pHighIntensityPts[i][j] = 1;
								IntPair currPoint;
								currPoint.X = i;
								currPoint.Y = j;
								highIntensityVector.push_back(currPoint);
                            }
                        }
                    }
					
					//We want to perform a round of k-means clustering to identify arms
					//We want a starting estimate of the center of the two. To do this, pick the two points that 
					//Are furthest away and set as center
					double maxDist = 0;
					int maxIndex1 = 0;
					int maxIndex2 = 0;
					for (int i = 0; i < highIntensityVector.size(); i++)
					{
						for (int j = i+1; j < highIntensityVector.size(); j++)
						{
							double dist = pow((highIntensityVector[i].X -highIntensityVector[j].X), 2)
							+ pow((highIntensityVector[i].Y -highIntensityVector[j].Y), 2);
							if (dist > maxDist)
							{
								maxDist = dist;
								maxIndex1 = i;
								maxIndex2 = j;
							}
						}
					}
                    
                    //Calculate linear regression from data points
                    float xSlope, ySlope;
                    float xIntercept, yIntercept;
                    CalculateLinearRegression(pHighIntensityPts, xSlope, ySlope, xIntercept, yIntercept);
                    
                    glPointSize(10);
					glBegin(GL_POINTS);
					glColor3f(1, 0, 1);
					
					//These are our clusters
					if (highIntensityVector.size() > 0)
					{
						glVertex2f(highIntensityVector[maxIndex1].X * width, highIntensityVector[maxIndex1].Y * height);
						glVertex2f(highIntensityVector[maxIndex2].X * width, highIntensityVector[maxIndex2].Y * height);
					}
                    glEnd();
                    glPointSize(1);
					
					//Run through an iteration of k-means clustering
					int * colourAffiliations = new int[highIntensityVector.size()];
					for (int i = 0; i < highIntensityVector.size(); i++)
					{
						double dist1 = pow((highIntensityVector[i].X -highIntensityVector[maxIndex1].X), 2)
						+ pow((highIntensityVector[i].Y -highIntensityVector[maxIndex1].Y), 2);

						double dist2 = pow((highIntensityVector[i].X -highIntensityVector[maxIndex2].X), 2)
						+ pow((highIntensityVector[i].Y -highIntensityVector[maxIndex2].Y), 2);
						
						if (dist1 > dist2)
						{
							colourAffiliations[i] = 0;
						}
						else {
							colourAffiliations[i] = 1;
						}
					}
					
					//Recalculate our centers
					double xTotal[2];
					double yTotal[2];
					int numPoints[2];
					xTotal[0] = xTotal[1] = yTotal[0] = yTotal[1] = numPoints[0] = numPoints[1]= 0;

					for (int i = 0; i < highIntensityVector.size(); i++)
					{
						if (colourAffiliations[i] == 0)
						{
							xTotal[0] += highIntensityVector[i].X;
							yTotal[0] += highIntensityVector[i].Y;
							numPoints[0]++;
						}
						else {
							xTotal[1] += highIntensityVector[i].X;
							yTotal[1] += highIntensityVector[i].Y;
							numPoints[1]++;
						}

					}
					
					double xCentroid[2];
					double yCentroid[2];
					xCentroid[0] = xTotal[0] / numPoints[0];
					xCentroid[1] = xTotal[1] / numPoints[1];
					yCentroid[0] = yTotal[0] / numPoints[0];
					yCentroid[1] = yTotal[1] / numPoints[1];
					
					glPointSize(10);
					glBegin(GL_POINTS);
					glColor3f(1, 0.5, 1);
					
					//These are our clusters
					if (highIntensityVector.size() > 0)
					{
						glVertex2f(xCentroid[0] * width, yCentroid[0] * height);
						glVertex2f(xCentroid[1] * width, yCentroid[1] * height);
					}
                    glEnd();
                    glPointSize(1);
					
                    //Classify based off euclidian distance
					for (int i = 0; i < highIntensityVector.size(); i++)
					{
						double dist1 = pow((highIntensityVector[i].X -xCentroid[0]), 2)
						+ pow((highIntensityVector[i].Y - yCentroid[0]), 2);
						
						double dist2 = pow((highIntensityVector[i].X -xCentroid[1]), 2)
						+ pow((highIntensityVector[i].Y - yCentroid[1]), 2);
						
						if (dist1 > dist2)
						{
							colourAffiliations[i] = 0;
						}
						else {
							colourAffiliations[i] = 1;
						}
					}
					
					glPointSize(10);
					glBegin(GL_POINTS);
					glColor3f(1, 1, 0);
					for (int i = 0; i < highIntensityVector.size(); i++)
					{
						if (colourAffiliations[i] == 0)
						{
							glVertex2f(highIntensityVector[i].X * width, highIntensityVector[i].Y * height);
						}
					}
					glEnd(); 
                    
                    //We know which wall we are coming from, ie. where the hand is reaching from, find the min distance from each boundary
					//IntPair xMax1 = 0, xMin1 = 64;
					IntPair xMax1, xMin1; 
					xMax1.X = 0;
					xMin1.X = XDepthGrid;
                    //IntPair yMax1 = 0, yMin1 = 64;
					IntPair yMax1, yMin1;
					yMax1.Y = 0;
					yMin1.Y = YDepthGrid;
					
					//IntPair xMax2 = 0, xMin2 = 64;
					IntPair xMax2, xMin2; 
					xMax2.X = 0;
					xMin2.X = XDepthGrid;
                    //IntPair yMax2 = 0, yMin2 = 64;
					IntPair yMax2, yMin2;
					yMax2.Y = 0;
					yMin2.Y = YDepthGrid;
					
					for (int i = 0; i < highIntensityVector.size(); i++)
					{
						if  (colourAffiliations[i] == 0)
						{
							if (highIntensityVector[i].X > xMax1.X)
							{
								xMax1 = highIntensityVector[i];
							}
							if (highIntensityVector[i].X < xMin1.X)
							{
								xMin1 = highIntensityVector[i];
							}
							
							if (highIntensityVector[i].Y > yMax1.Y)
							{
								yMax1 = highIntensityVector[i];
							}
							if (highIntensityVector[i].Y < yMin1.Y)
							{
								yMin1 = highIntensityVector[i];
							}
						}
						else 
						{
							if (highIntensityVector[i].X > xMax2.X)
							{
								xMax2 = highIntensityVector[i];
							}
							if (highIntensityVector[i].X < xMin2.X)
							{
								xMin2 = highIntensityVector[i];
							}
							
							if (highIntensityVector[i].Y > yMax2.Y)
							{
								yMax2 = highIntensityVector[i];
							}
							if (highIntensityVector[i].Y < yMin2.Y)
							{
								yMin2 = highIntensityVector[i];
							}
						}
						
						//TODO: Change to accumulated distances instead for better noise detection!
						if (highIntensityVector[i].X < minLeftWallDist)
							minLeftWallDist = highIntensityVector[i].X;
						
						if ( (XDepthGrid - highIntensityVector[i].X) < minRightWallDist )
							minRightWallDist = (XDepthGrid -highIntensityVector[i].X);
						
						if (highIntensityVector[i].Y < minTopWallDist)
							minTopWallDist = highIntensityVector[i].Y;
						
						if ((YDepthGrid - highIntensityVector[i].Y) < minBottomWallDist)
							minBottomWallDist = (YDepthGrid - highIntensityVector[i].Y);
						
					}
                    
					/*glBegin(GL_LINE_LOOP);
                    glLineWidth(2.0f);
                    glColor3f(0.5f, 0.35f, 0.05f);
					//Draw line of regression slope
                    for (int i = 0; i < XDepthGrid; i++)
                    {
                        for (int j = 0; j < YDepthGrid; j++)
                        {
                            if (pHighIntensityPts[i][j] == 1)
                            {
                                float yCoord = i * ySlope + yIntercept;
                                if ((yCoord < yMax) && (yCoord > yMin))
                                    glVertex2f( (i+0.5) * width, (yCoord + 0.5) * height);
                                
                                float xCoord = (j - yIntercept) / ySlope;
                                if ((xCoord < xMax) && (xCoord > xMin))
                                    glVertex2f( (xCoord+0.5) * width, (j + 0.5) * height);
                            }
                        }
                    }
                    glLineWidth(1.0f);
                    glEnd();*/
					
					IntPair handTip1;
					IntPair handTip2;
					
					bool bHandTipOverSoap = false;
                    
                    glPointSize(10);
                    glBegin(GL_POINTS);
                    glColor3f(0,0,0);
                    
                    //Are we reaching from the top of the grid
                    if ((minTopWallDist <= minRightWallDist)
						&& (minTopWallDist <= minLeftWallDist))
					{
						//int yCoord = yMax ;
                        //float xCoord = (yCoord - yIntercept) / ySlope;
                        //glVertex2f( (xCoord + 0.5f) * width, (yCoord + 0.5) * height);
						handTip1 = yMax1;
						handTip2 = yMax2;
						glVertex2f( (yMax1.X + 0.5f) * width, (yMax1.Y + 0.5) * height);
						glVertex2f( (yMax2.X + 0.5f) * width, (yMax2.Y + 0.5) * height);
						
						int bottomSoap = g_DrawUserInput.SoapRect.uBottom;
						int topSoap = g_DrawUserInput.SoapRect.uTop;
						int leftSoap = g_DrawUserInput.SoapRect.uLeft;
						int rightSoap = g_DrawUserInput.SoapRect.uRight;
						
                        //Draw a square at the furthest point from the top of the grid (for each cluster aka each hand)
						for (int i = 0; i < highIntensityVector.size(); i++)
						{
							if (colourAffiliations[i] == 0)
							{
								if (highIntensityVector[i].Y == yMax1.Y)
								{
                                    //Plot the point for each cluster at the extremity
									glVertex2f( (highIntensityVector[i].X + 0.5f) * width, (highIntensityVector[i].Y + 0.5) * height);
									
									int scaledY = highIntensityVector[i].Y * 480 / YDepthGrid;
									int scaledX = highIntensityVector[i].X * 640 / XDepthGrid;
									if ( ( scaledY < topSoap) && (scaledY > bottomSoap) 
										&& (scaledX > leftSoap) && (scaledX < rightSoap) )
									{
										bHandTipOverSoap = true;
									}
								}
							}
							else {
								if (highIntensityVector[i].Y == yMax2.Y)
								{
									glVertex2f( (highIntensityVector[i].X + 0.5f) * width, (highIntensityVector[i].Y + 0.5) * height);
									
									int scaledY = highIntensityVector[i].Y * WIN_SIZE_Y / YDepthGrid;
									int scaledX = highIntensityVector[i].X * WIN_SIZE_X / XDepthGrid;
									if ( ( scaledY < topSoap) && (scaledY > bottomSoap) 
										&& (scaledX > leftSoap) && (scaledX < rightSoap) )
									{
										bHandTipOverSoap = true;
									}
								}
							}

						}
					}
					else if ((minBottomWallDist < minRightWallDist) && (minBottomWallDist < minTopWallDist)
							 && (minBottomWallDist < minLeftWallDist))
                    {
                        //TODO: Mirror the work from the top wall for each wall
                        //int yCoord = yMin ;
                        //float xCoord = (yCoord - yIntercept) / ySlope;
                        //glVertex2f( (xCoord + 0.5f) * width, (yCoord + 0.5) * height);
						handTip1 = yMin1;
						handTip2 = yMin2;
						glVertex2f( (yMin1.X + 0.5f) * width, (yMin1.Y + 0.5) * height);
						glVertex2f( (yMin2.X + 0.5f) * width, (yMin2.Y + 0.5) * height);
                    }
					else if ((minLeftWallDist < minBottomWallDist) && ( minLeftWallDist < minTopWallDist)
                        && (minLeftWallDist < minRightWallDist))
                    {
                        //TODO: Mirror the work from the top wall for each wall
                        //int xCoord = xMax;
                        //float yCoord = xMax * ySlope + yIntercept;
                        //glVertex2f( (xCoord + 0.5f) * width, (yCoord + 0.5) * height);
						
						handTip1 = xMax1;
						handTip2 = xMax2;
						
						glVertex2f( (xMax1.X + 0.5f) * width, (xMax1.Y + 0.5) * height);
						glVertex2f( (xMax2.X + 0.5f) * width, (xMax2.Y + 0.5) * height);
                    }
                    else if ((minRightWallDist < minBottomWallDist) && (minRightWallDist < minTopWallDist)
							 && (minRightWallDist < minLeftWallDist))
                    {
                        //TODO: Mirror the work from the top wall for each wall
                        //int xCoord = xMin ;
                        //float yCoord = xMax * ySlope + yIntercept;
                        //glVertex2f( (xCoord + 0.5f) * width, (yCoord + 0.5) * height);
						
						handTip1 = xMin1;
						handTip2 = xMin2;
						glVertex2f( (xMin1.X + 0.5f) * width, (xMin1.Y + 0.5) * height);
						glVertex2f( (xMin2.X + 0.5f) * width, (xMin2.Y + 0.5) * height);
                    }
                    else
                    {
                        //TODO: Mirror the work from the top wall for each wall
                        //int yCoord = yMax ;
                        //float xCoord = (yCoord - yIntercept) / ySlope;
                        //glVertex2f( (xCoord + 0.5f) * width, (yCoord + 0.5) * height);
						handTip1 = yMax1;
						handTip2 = yMax2;
						glVertex2f( (yMax1.X + 0.5f) * width, (yMax1.Y + 0.5) * height);
						glVertex2f( (yMax2.X + 0.5f) * width, (yMax2.Y + 0.5) * height);
                    }
					glEnd();
					
					//Hand overtop of soap
					if (bHandTipOverSoap)
					{
						nXLocation =  WIN_SIZE_X - 350;
						sprintf(buf, "Hand Tip Over Soap");
                        nYLocation = WIN_SIZE_Y - 140;
                        glColor3f(1,0,0);
                        glRasterPos2i(nXLocation,nYLocation);
                        glPrintString(GLUT_BITMAP_HELVETICA_12, buf);
					}
                }
            }
        }
	}
}

void drawPointerMode(IntPair* pPointer)
{
	char buf[512] = "";
	int nCharWidth = glutBitmapWidth(GLUT_BITMAP_HELVETICA_18, '0');
	int nPointerValue = 0;

	XnDouble dTimestampDivider = 1E6;

	const DepthMetaData* pDepthMD = getDepthMetaData();

	if (pDepthMD != NULL)
	{
		// Print the scale black background
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glBegin(GL_QUADS);
		glColor4f(0, 0, 0, 0.7);
		glVertex2i(0, WIN_SIZE_Y); // lower left
		glVertex2i(WIN_SIZE_X, WIN_SIZE_Y);
		glVertex2i(WIN_SIZE_X, WIN_SIZE_Y - 135);
		glVertex2i(0, WIN_SIZE_Y - 135);
		glEnd();

		glDisable(GL_BLEND);

		// set a large point size (for the scale)
		glPointSize(15);

		// Print the scale data
		glBegin(GL_POINTS);
		for (int i = 0; i < pDepthMD->ZRes(); i+=1)
		{
			float fNewColor = g_pDepthHist[i];
			if ((fNewColor > 0.004) && (fNewColor < 0.996))
			{
				glColor3f(fNewColor, fNewColor, 0);
				glVertex3f(((i/10)*2), WIN_SIZE_Y - 23, 1);
			}
		}
		glEnd();

		// Print the pointer scale data
		if (pPointer != NULL)
		{
			// make sure pointer in on a depth pixel (take in mind cropping might be in place)
			IntPair pointerInDepth = *pPointer;
			pointerInDepth.X -= pDepthMD->XOffset();
			pointerInDepth.Y -= pDepthMD->YOffset();

			if (pointerInDepth.X < (int)pDepthMD->XRes() && pointerInDepth.Y < (int)pDepthMD->YRes())
			{
				nPointerValue = (*pDepthMD)(pointerInDepth.X, pointerInDepth.Y);

				glBegin(GL_POINTS);
				glColor3f(1,0,0);
				glVertex3f(10 + ((nPointerValue/10)*2), WIN_SIZE_Y - 70, 1);
				glEnd();
			}
		}

		// Print the scale texts
		for (int i = 0; i < pDepthMD->ZRes()/10; i+=25)
		{
			int xPos = i*2 + 10;

			// draw a small line in this position
			glBegin(GL_LINES);
			glColor3f(0, 1, 0);
			glVertex2i(xPos, WIN_SIZE_Y - 54);
			glVertex2i(xPos, WIN_SIZE_Y - 62);
			glEnd();

			// place a label under, and in the middle of, that line.
			int chars = sprintf(buf, "%d", i);
			glColor3f(1,0,0);
			glRasterPos2i(xPos - chars*nCharWidth/2, WIN_SIZE_Y - 40);
			glPrintString(GLUT_BITMAP_HELVETICA_18,buf);
		}

		sprintf(buf, "%s - Frame %4u, Timestamp %.3f", getDepthGenerator()->GetInfo().GetInstanceName(), pDepthMD->FrameID(), (double)pDepthMD->Timestamp()/dTimestampDivider);
	}

	const ImageMetaData* pImageMD = getImageMetaData();
	if (pImageMD != NULL)
	{
		if (buf[0] != '\0')
			sprintf(buf + strlen(buf), " | ");

		sprintf(buf + strlen(buf), "%s - Frame %4u, Timestamp %.3f", getImageGenerator()->GetInfo().GetInstanceName(), pImageMD->FrameID(), (double)pImageMD->Timestamp()/dTimestampDivider);
	}

	const IRMetaData* pIRMD = getIRMetaData();
	if (pIRMD != NULL)
	{
		if (buf[0] != '\0')
			sprintf(buf + strlen(buf), " | ");

		sprintf(buf + strlen(buf), "%s - Frame %4u, Timestamp %.3f", getIRGenerator()->GetInfo().GetInstanceName(), pIRMD->FrameID(), (double)pIRMD->Timestamp()/dTimestampDivider);
	}

	const AudioMetaData* pAudioMD = getAudioMetaData();
	if (pAudioMD != NULL)
	{
		if (buf[0] != '\0')
			sprintf(buf + strlen(buf), " | ");

		sprintf(buf + strlen(buf), "%s - Timestamp %.3f", getAudioGenerator()->GetInfo().GetInstanceName(), (double)pAudioMD->Timestamp()/dTimestampDivider);
	}

	int nYLocation = WIN_SIZE_Y - 88;
	glColor3f(1,0,0);
	glRasterPos2i(10,nYLocation);
	glPrintString(GLUT_BITMAP_HELVETICA_18, buf);
	nYLocation -= 26;

	if (pPointer != NULL && isStatisticsActive())
	{
		XnPixelStatistics* pStatistics = &g_PixelStatistics[pPointer->Y * pDepthMD->XRes() + pPointer->X];
		sprintf(buf, "Collected: %3u, Min: %4u Max: %4u Avg: %6.2f StdDev: %6.2f", 
			pStatistics->nCount, pStatistics->nMin, pStatistics->nMax, pStatistics->dAverage, pStatistics->dStdDev);
		glRasterPos2i(10,nYLocation);
		glPrintString(GLUT_BITMAP_HELVETICA_18, buf);
		nYLocation -= 26;
	}

	if (pPointer != NULL)
	{
		// Print the pointer text
		XnUInt64 nCutOffMin = 0;
		XnUInt64 nCutOffMax = (pDepthMD != NULL) ? g_nMaxDepth : 0;

		XnChar sPointerValue[100];
		if (nPointerValue != g_nMaxDepth)
		{
			sprintf(sPointerValue, "%.1f", (float)nPointerValue/10);
		}
		else
		{
			sprintf(sPointerValue, "-");
		}

		sprintf(buf, "Pointer Value: %s (X:%d Y:%d) Cutoff: %llu-%llu.", 
			sPointerValue, pPointer->X, pPointer->Y, nCutOffMin, nCutOffMax);

		glRasterPos2i(10,nYLocation);
		glPrintString(GLUT_BITMAP_HELVETICA_18, buf);
		nYLocation -= 26;
	}
}

void drawCenteredMessage(void* font, int y, const char* message, float fRed, float fGreen, float fBlue)
{
	const XnUInt32 nMaxLines = 5;
	XnChar buf[512];
	XnChar* aLines[nMaxLines];
	XnUInt32 anLinesWidths[nMaxLines];
	XnUInt32 nLine = 0;
	XnUInt32 nLineLengthChars = 0;
	XnUInt32 nLineLengthPixels = 0;
	XnUInt32 nMaxLineLength = 0;
	
	aLines[0] = buf;
	
	// parse message to lines
	const char* pChar = message;
	for (;;)
	{
		if (*pChar == '\n' || *pChar == '\0')
		{
			if (nLineLengthChars > 0)
			{
				aLines[nLine][nLineLengthChars++] = '\0';
				aLines[nLine+1] = &aLines[nLine][nLineLengthChars];
				anLinesWidths[nLine] = nLineLengthPixels;
				nLine++;
				if (nLineLengthPixels > nMaxLineLength)
				{
					nMaxLineLength = nLineLengthPixels;
				}
				nLineLengthPixels = 0;
				nLineLengthChars = 0;
			}

			if (nLine >= nMaxLines || *pChar == '\0')
			{
				break;
			}
		}
		else
		{
			aLines[nLine][nLineLengthChars++] = *pChar;
			nLineLengthPixels += glutBitmapWidth(font, *pChar);
		}
		pChar++;
	}
	
	XnUInt32 nHeight = 26;
	int nXLocation = XN_MAX(0, (WIN_SIZE_X - nMaxLineLength) / 2);
	int nYLocation = y;

	// Draw black background
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glBegin(GL_QUADS);
	glColor4f(0, 0, 0, 0.6);
	glVertex2i(nXLocation - 5, nYLocation - nHeight - 5);
	glVertex2i(nXLocation + nMaxLineLength + 5, nYLocation - nHeight - 5);
	glVertex2i(nXLocation + nMaxLineLength + 5, nYLocation + nHeight * nLine + 5);
	glVertex2i(nXLocation - 5, nYLocation + nHeight * nLine + 5);
	glEnd();

	glDisable(GL_BLEND);

	// show message
	glColor3f(fRed, fGreen, fBlue);
	for (XnUInt32 i = 0; i < nLine; ++i)
	{
		glRasterPos2i(nXLocation + (nMaxLineLength - anLinesWidths[i])/2, nYLocation + i * nHeight);
		glPrintString(font, aLines[i]);
	}
}

void drawUserMessage()
{
	static XnUInt64 nStartShowMessage = 0;
	if (g_DrawConfig.bShowMessage)
	{
		g_DrawConfig.bShowMessage = false;
		xnOSGetTimeStamp(&nStartShowMessage);
	}
	
	XnUInt64 nNow;
	xnOSGetTimeStamp(&nNow);

	if (nNow - nStartShowMessage < 3000)
	{
		drawCenteredMessage(GLUT_BITMAP_TIMES_ROMAN_24, WIN_SIZE_Y * 4 / 5, g_csUserMessage, 0, 1, 0);
	}
}

void printRecordingInfo()
{
	char csMessage[256];
	getCaptureMessage(csMessage);

	if (csMessage[0] != 0)
		drawCenteredMessage(GLUT_BITMAP_TIMES_ROMAN_24, 30, csMessage, 1, 0, 0);

	sprintf(csMessage, "Capture Formats - Depth: %s | Image: %s | IR: %s | Audio: %s",
		captureGetDepthFormatName(), 
		captureGetImageFormatName(), 
		captureGetIRFormatName(), 
		captureGetAudioFormatName());

	drawCenteredMessage(GLUT_BITMAP_HELVETICA_12, WIN_SIZE_Y - 3, csMessage, 0, 1, 0);
}

void printStatisticsInfo()
{
	char csMessage[256];
	getStatisticsMessage(csMessage);

	if (csMessage[0] != 0)
		drawCenteredMessage(GLUT_BITMAP_HELVETICA_18, 20, csMessage, 0, 1, 0);
}

void printHelpGroup(int nXLocation, int* pnYLocation, const char* csGroup)
{
	int nYLocation = *pnYLocation;

	unsigned char aKeys[20];
	const char* aDescs[20];
	int nCount;

	getGroupItems(csGroup, aKeys, aDescs, &nCount);

	glColor3f(0, 1, 0);
	glRasterPos2i(nXLocation, nYLocation);
	glPrintString(GLUT_BITMAP_TIMES_ROMAN_24, csGroup);
	nYLocation += 30;

	for (int i = 0; i < nCount; ++i, nYLocation += 22)
	{
		char buf[256];
		switch (aKeys[i])
		{
		case 27:
			sprintf(buf, "Esc");
			break;
		default:
			sprintf(buf, "%c", aKeys[i]);
			break;
		}

		glColor3f(1, 0, 0);
		glRasterPos2i(nXLocation, nYLocation);
		glPrintString(GLUT_BITMAP_HELVETICA_18, buf);

		glRasterPos2i(nXLocation + 40, nYLocation);
		glPrintString(GLUT_BITMAP_HELVETICA_18, aDescs[i]);
	}

	*pnYLocation = nYLocation + 20;
}

void drawErrorState()
{
	if (g_DrawConfig.strErrorState == NULL)
		return;

	// place a black rect on entire screen
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		
	glBegin(GL_QUADS);
	glColor4f(0, 0, 0, 0.8);
	glVertex2i(0, 0);
	glVertex2i(WIN_SIZE_X, 0);
	glVertex2i(WIN_SIZE_X, WIN_SIZE_Y);
	glVertex2i(0, WIN_SIZE_Y);
	glEnd();
	glDisable(GL_BLEND);

	int nYLocation = WIN_SIZE_Y/2 - 30;

	drawCenteredMessage(GLUT_BITMAP_TIMES_ROMAN_24, nYLocation, "ERROR!", 1, 0, 0);
	nYLocation += 40;
	drawCenteredMessage(GLUT_BITMAP_TIMES_ROMAN_24, nYLocation, g_DrawConfig.strErrorState, 1, 0, 0);
}

void drawHelpScreen()
{
	int nXStartLocation = WIN_SIZE_X/8;
	int nYStartLocation = WIN_SIZE_Y/5;
	int nXEndLocation = WIN_SIZE_X*7/8;
	int nYEndLocation = WIN_SIZE_Y*4/5;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		

	glBegin(GL_QUADS);
	glColor4f(0, 0, 0, 0.8);
	glVertex2i(nXStartLocation, nYStartLocation);
	glVertex2i(nXStartLocation, nYEndLocation);
	glVertex2i(nXEndLocation, nYEndLocation);
	glVertex2i(nXEndLocation, nYStartLocation);
	glEnd();

	glDisable(GL_BLEND);

	// set color to red
	glColor3f(1, 0, 0);

	// leave some margins
	nYStartLocation += 30;
	nXStartLocation += 30;

	// print left pane
	int nXLocation = nXStartLocation;
	int nYLocation = nYStartLocation;
	printHelpGroup(nXLocation, &nYLocation, KEYBOARD_GROUP_PRESETS);
	printHelpGroup(nXLocation, &nYLocation, KEYBOARD_GROUP_DISPLAY);
	printHelpGroup(nXLocation, &nYLocation, KEYBOARD_GROUP_DEVICE);

	// print right pane
	nXLocation = WIN_SIZE_X/2;
	nYLocation = nYStartLocation;
	printHelpGroup(nXLocation, &nYLocation, KEYBOARD_GROUP_PLAYER);
	printHelpGroup(nXLocation, &nYLocation, KEYBOARD_GROUP_CAPTURE);
	printHelpGroup(nXLocation, &nYLocation, KEYBOARD_GROUP_GENERAL);
}

void drawUserInput(bool bCursor)
{
	if (bCursor)
	{
		// draw cursor
		IntPair cursor = g_DrawUserInput.Cursor;
		glPointSize(1);
		glBegin(GL_POINTS);
		glColor3f(1,0,0);
		glVertex2i(cursor.X, cursor.Y);

		// upper left marker
		glVertex2i(cursor.X - 2, cursor.Y - 2);
		glVertex2i(cursor.X - 2, cursor.Y - 1);
		glVertex2i(cursor.X - 1, cursor.Y - 2);

		// bottom left marker
		glVertex2i(cursor.X - 2, cursor.Y + 2);
		glVertex2i(cursor.X - 2, cursor.Y + 1);
		glVertex2i(cursor.X - 1, cursor.Y + 2);

		// upper right marker
		glVertex2i(cursor.X + 2, cursor.Y - 2);
		glVertex2i(cursor.X + 2, cursor.Y - 1);
		glVertex2i(cursor.X + 1, cursor.Y - 2);

		// lower right marker
		glVertex2i(cursor.X + 2, cursor.Y + 2);
		glVertex2i(cursor.X + 2, cursor.Y + 1);
		glVertex2i(cursor.X + 1, cursor.Y + 2);

		glEnd();
	}

	// draw selection frame
	if (g_DrawUserInput.State == SELECTION_ACTIVE)
	{
        
        IntRect completedRect, drawingRect;
        if (g_DrawUserInput.CurrentSelect == 0)
        {
            completedRect = g_DrawUserInput.TowelRect;
            drawingRect = g_DrawUserInput.SoapRect;
        }
        else
        {
            completedRect = g_DrawUserInput.SoapRect;
            drawingRect = g_DrawUserInput.TowelRect;
        }
        
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		

		glBegin(GL_QUADS);
		glColor4f(0, 1, 1, 0.5);
        
		glVertex2i(drawingRect.uLeft, drawingRect.uTop); // Upper left
		glVertex2i(drawingRect.uRight, drawingRect.uTop); // upper right
		glVertex2i(drawingRect.uRight, drawingRect.uBottom); // lower right
		glVertex2i(drawingRect.uLeft, drawingRect.uBottom); // lower left
		glEnd();

		glDisable(GL_BLEND);
        
        glBegin(GL_LINES);
            glVertex2f(completedRect.uLeft, completedRect.uBottom);
            glVertex2f(completedRect.uLeft, completedRect.uTop);
            
            glVertex2f(completedRect.uLeft, completedRect.uBottom);
            glVertex2f(completedRect.uRight, completedRect.uBottom);
            
            glVertex2f(completedRect.uRight, completedRect.uTop);
            glVertex2f(completedRect.uRight, completedRect.uBottom);
            
            glVertex2f(completedRect.uRight, completedRect.uTop);
            glVertex2f(completedRect.uLeft, completedRect.uTop);
        glEnd();
	}
    else if (g_DrawUserInput.State == SELECTION_DONE)
    {
        //Draw a box over what has been selected
        
        glLineWidth(2.0f);
        glBegin(GL_LINES);
        {
            glColor3f(0, 1, 1);
            IntRect currRect;
            for (int i = 0; i < 2; i++)
            {
                if (i == 0)
                    currRect = g_DrawUserInput.SoapRect;
                else
                    currRect = g_DrawUserInput.TowelRect;
                
                glVertex2f(currRect.uLeft, currRect.uBottom);
                glVertex2f(currRect.uLeft, currRect.uTop);
                
                glVertex2f(currRect.uLeft, currRect.uBottom);
                glVertex2f(currRect.uRight, currRect.uBottom);
                
                glVertex2f(currRect.uRight, currRect.uTop);
                glVertex2f(currRect.uRight, currRect.uBottom);
                
                glVertex2f(currRect.uRight, currRect.uTop);
                glVertex2f(currRect.uLeft, currRect.uTop);
            }
        }
        glEnd();
        
        char buf[50] = "";
        int nYLocation = g_DrawUserInput.SoapRect.uBottom - 5;
        sprintf(buf, "S");
        //nYLocation += 20;
        int nXLocation =  g_DrawUserInput.SoapRect.uLeft - 5;
        glColor3f(0,1,1);
        glRasterPos2i(nXLocation,nYLocation);
        glPrintString(GLUT_BITMAP_HELVETICA_18, buf);
        
        nYLocation = g_DrawUserInput.TowelRect.uBottom - 5;
        sprintf(buf, "T");
        //nYLocation += 20;
        nXLocation =  g_DrawUserInput.TowelRect.uLeft - 5;
        glColor3f(0,1,1);
        glRasterPos2i(nXLocation,nYLocation);
        glPrintString(GLUT_BITMAP_HELVETICA_18, buf);
    }
}

void fixLocation(IntRect* pLocation, int xRes, int yRes)
{
	double resRatio = (double)xRes / yRes;

	double locationRatio = (pLocation->uRight - pLocation->uLeft) / (pLocation->uTop - pLocation->uBottom);

	if (locationRatio > resRatio) 
	{
		// location is wider. use height as reference.
		double width = (pLocation->uTop - pLocation->uBottom) * resRatio;
		pLocation->uRight = (pLocation->uLeft + width);
	}
	else if (locationRatio < resRatio)
	{
		// res is wider. use width as reference.
		double height = (pLocation->uRight - pLocation->uLeft) / resRatio;
		pLocation->uTop = (pLocation->uBottom + height);
	}
}

bool isPointInRect(IntPair point, IntRect* pRect)
{
	return (point.X >= pRect->uLeft && point.X <= pRect->uRight &&
		point.Y >= pRect->uBottom && point.Y <= pRect->uTop);
}

void drawPlaybackSpeed()
{
	XnDouble dSpeed = getPlaybackSpeed();
	if (dSpeed != 1.0)
	{
		XnChar strSpeed[30];
		int len = sprintf(strSpeed, "x%g", dSpeed);
		int width = 0;
		for (int i = 0; i < len; ++i)
			width += glutBitmapWidth(GLUT_BITMAP_TIMES_ROMAN_24, strSpeed[i]);

		glColor3f(0, 1, 0);
		glRasterPos2i(WIN_SIZE_X - width - 3, 30);
		glPrintString(GLUT_BITMAP_TIMES_ROMAN_24, strSpeed);
	}
}

void ScalePoint(XnPoint3D& point)
{
	point.X *= WIN_SIZE_X;
	point.X /= getDepthMetaData()->XRes();
    
	point.Y *= WIN_SIZE_Y;
	point.Y /= getDepthMetaData()->YRes();
}

void DisplayPostDraw(HandTracker& m_HandTracker)
{
	typedef TrailHistory			History;
	typedef History::ConstIterator	HistoryIterator;
	typedef Trail::ConstIterator	TrailIterator;
    
	static const float colours[][3] =
	{
		{ 0.5f, 0.5f, 0.5f},
		{ 0.0f, 1.0f, 0.0f},
		{ 0.0f, 0.5f, 1.0f},
		{ 1.0f, 1.0f, 0.0f},
		{ 1.0f, 0.5f, 0.0f},
		{ 1.0f, 0.0f, 1.0f}
	};
	const TrailHistory&	history = m_HandTracker.GetHistory();
    
	// History points coordinates buffer
	XnFloat	coordinatesOne[3 * MAX_HAND_TRAIL_LENGTH];
    XnFloat coordinatesTwo[3 * MAX_HAND_TRAIL_LENGTH];
    
	const HistoryIterator	hend = history.End();
	for(HistoryIterator		hit = history.Begin(); hit != hend; ++hit)
	{
        
		// Dump the history to local buffer
		int				numpoints = 0;
		const Trail&	trail = hit->Value();
        
		const TrailIterator	tend = trail.End();
		for(TrailIterator	tit = trail.Begin(); tit != tend; ++tit)
		{
			XnPoint3D	point = *tit;
			DepthGenerator* pDepthGen = getDepthGenerator();
            pDepthGen->ConvertRealWorldToProjective(1, &point, &point);
			ScalePoint(point);
			coordinatesOne[numpoints * 3] = point.X / 2;
			coordinatesOne[numpoints * 3 + 1] = point.Y / 2;
			coordinatesOne[numpoints * 3 + 2] = 0;
            
            //ImageGenerator* pImageGen = getImageGenerator();
            //pImageGen->ConvertRealWorldToProjective(1, &point, &point);
			//ScalePoint(point);
            coordinatesTwo[numpoints * 3] = point.X / 2 + WIN_SIZE_X/2;
			coordinatesTwo[numpoints * 3 + 1] = point.Y / 2;
			coordinatesTwo[numpoints * 3 + 2] = 0;
            
			++numpoints;
		}
		assert(numpoints <= MAX_HAND_TRAIL_LENGTH);
        
		// Draw the hand trail history
		XnUInt32 nColor = hit->Key() % LENGTHOF(colours);
		glColor4f(colours[nColor][0],
                  colours[nColor][1],
                  colours[nColor][2],
                  1.0f);
		glPointSize(2);
		glVertexPointer(3, GL_FLOAT, 0, coordinatesOne);
		glDrawArrays(GL_LINE_STRIP, 0, numpoints);
        glPointSize(8);
		glDrawArrays(GL_POINTS, 0, 1);
        glVertexPointer(3, GL_FLOAT, 0, coordinatesTwo);
		glDrawArrays(GL_LINE_STRIP, 0, numpoints);
		// Current point as a larger dot
		glPointSize(8);
		glDrawArrays(GL_POINTS, 0, 1);
		glFlush();
	}
}

void setDepthDrawing(int nColoring)
{
	g_DrawConfig.Streams.Depth.Coloring	= (DepthColoringType)nColoring;
}

void setImageDrawing(int nColoring)
{
	g_DrawConfig.Streams.Image.Coloring	= (ImageColoringType)nColoring;
}

void drawFrame()
{
	// calculate locations
	g_DrawConfig.DepthLocation.uBottom = 0;
	g_DrawConfig.DepthLocation.uTop = WIN_SIZE_Y - 1;
	g_DrawConfig.DepthLocation.uLeft = 0;
	g_DrawConfig.DepthLocation.uRight = WIN_SIZE_X - 1;
    
	g_DrawConfig.ImageLocation.uBottom = 0;
	g_DrawConfig.ImageLocation.uTop = WIN_SIZE_Y - 1;
	g_DrawConfig.ImageLocation.uLeft = 0;
	g_DrawConfig.ImageLocation.uRight = WIN_SIZE_X - 1;
    
	if (g_DrawConfig.Streams.ScreenArrangement == SIDE_BY_SIDE)
	{
		g_DrawConfig.DepthLocation.uTop = WIN_SIZE_Y / 2 - 1;
		g_DrawConfig.DepthLocation.uRight = WIN_SIZE_X / 2 - 1;
		g_DrawConfig.ImageLocation.uTop = WIN_SIZE_Y / 2 - 1;
		g_DrawConfig.ImageLocation.uLeft = WIN_SIZE_X / 2;
	}
    
	// Texture map init
	const DepthMetaData* pDepthMD = getDepthMetaData();
	if (isDepthOn())
	{
		g_nMaxDepth = getDepthGenerator()->GetDeviceMaxDepth();
		TextureMapInit(&g_texDepth, pDepthMD->FullXRes(), pDepthMD->FullYRes(), 4, pDepthMD->XRes(), pDepthMD->YRes());
		fixLocation(&g_DrawConfig.DepthLocation, pDepthMD->FullXRes(), pDepthMD->FullYRes());
	}
    
	const MapMetaData* pImageMD = NULL;
    
	if (isImageOn())
	{
		pImageMD = getImageMetaData();
	}
	else if (isIROn())
	{
		pImageMD = getIRMetaData();
	}
    
	if (pImageMD != NULL)
	{
		TextureMapInit(&g_texImage, pImageMD->FullXRes(), pImageMD->FullYRes(), 4, pImageMD->XRes(), pImageMD->YRes());
		fixLocation(&g_DrawConfig.ImageLocation, pImageMD->FullXRes(), pImageMD->FullYRes());
	}
    
	// check if pointer is over a map
	bool bOverDepth = (pDepthMD != NULL) && isPointInRect(g_DrawUserInput.Cursor, &g_DrawConfig.DepthLocation);
	bool bOverImage = (pImageMD != NULL) && isPointInRect(g_DrawUserInput.Cursor, &g_DrawConfig.ImageLocation);
    
	IntPair pointerInDepth;
	IntPair pointerInImage;
    
	if (bOverDepth)
	{
		pointerInDepth.X = (double)(g_DrawUserInput.Cursor.X - g_DrawConfig.DepthLocation.uLeft) / (g_DrawConfig.DepthLocation.uRight - g_DrawConfig.DepthLocation.uLeft + 1) * pDepthMD->FullXRes();
		pointerInDepth.Y = (double)(g_DrawUserInput.Cursor.Y - g_DrawConfig.DepthLocation.uBottom) / (g_DrawConfig.DepthLocation.uTop - g_DrawConfig.DepthLocation.uBottom + 1) * pDepthMD->FullYRes();
	}
    
	if (bOverImage)
	{
		pointerInImage.X = (double)(g_DrawUserInput.Cursor.X - g_DrawConfig.ImageLocation.uLeft) / (g_DrawConfig.ImageLocation.uRight - g_DrawConfig.ImageLocation.uLeft + 1) * pImageMD->FullXRes();
		pointerInImage.Y = (double)(g_DrawUserInput.Cursor.Y - g_DrawConfig.ImageLocation.uBottom) / (g_DrawConfig.ImageLocation.uTop - g_DrawConfig.ImageLocation.uBottom + 1) * pImageMD->FullYRes();
	}
    
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    
	// Setup the opengl env for fixed location view
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,WIN_SIZE_X,WIN_SIZE_Y,0,-1.0,1.0);
	glDisable(GL_DEPTH_TEST);
    
	if (g_DrawConfig.Streams.Depth.Coloring == CYCLIC_RAINBOW_HISTOGRAM || g_DrawConfig.Streams.Depth.Coloring == LINEAR_HISTOGRAM || g_DrawConfig.bShowPointer)
		calculateHistogram();
    
	drawColorImage(&g_DrawConfig.ImageLocation, bOverImage ? &pointerInImage : NULL);
    
	drawDepth(&g_DrawConfig.DepthLocation, bOverDepth ? &pointerInDepth : NULL);
    
	printStatisticsInfo();
	printRecordingInfo();
    
	if (g_DrawConfig.bShowPointer)
		drawPointerMode(bOverDepth ? &pointerInDepth : NULL);
    
	drawUserInput(!bOverDepth && !bOverImage);
    
	drawUserMessage();
	drawPlaybackSpeed();
    
	if (g_DrawConfig.strErrorState != NULL)
		drawErrorState();
    
	if (g_DrawConfig.bHelp)
		drawHelpScreen();
    
    //DisplayPostDraw(*m_handTracker);
    
    glutSwapBuffers();
}
