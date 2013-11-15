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
#include "XnPlatform.h"
#include "MouseInput.h"
#if (XN_PLATFORM == XN_PLATFORM_MACOSX)
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif

// --------------------------------
// Types
// --------------------------------
typedef struct MouseInput
{
	SelectionState nSelectionState;
	IntPair StartSelection;
	IntPair LastLocation;
	SelectionRectangleChangedPtr pSelectionCallback;
	CursorMovedPtr pCursorCallback;
} MouseInput;

MouseInput g_MouseInput = {SELECTION_NONE, {0,0}, {0,0}, NULL, NULL};

int iActiveSelection = SOAP;

IntRect Selections[2];

void mouseInputCallSelection()
{
	if (g_MouseInput.pSelectionCallback != NULL)
	{
		Selections[iActiveSelection].uBottom = XN_MIN(g_MouseInput.StartSelection.Y, g_MouseInput.LastLocation.Y);
		Selections[iActiveSelection].uTop = XN_MAX(g_MouseInput.StartSelection.Y, g_MouseInput.LastLocation.Y);
		Selections[iActiveSelection].uLeft = XN_MIN(g_MouseInput.StartSelection.X, g_MouseInput.LastLocation.X);
		Selections[iActiveSelection].uRight = XN_MAX(g_MouseInput.StartSelection.X, g_MouseInput.LastLocation.X);

		g_MouseInput.pSelectionCallback(g_MouseInput.nSelectionState, Selections[0], Selections[1], iActiveSelection);
	}
}

void mouseInputMotion(int x, int y)
{
	g_MouseInput.LastLocation.X = x;
	g_MouseInput.LastLocation.Y	= y;

	if (g_MouseInput.pCursorCallback != NULL)
		g_MouseInput.pCursorCallback(g_MouseInput.LastLocation);

	if (g_MouseInput.nSelectionState == SELECTION_ACTIVE)
		mouseInputCallSelection();
}

void mouseInputButton(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			g_MouseInput.nSelectionState = SELECTION_ACTIVE;
			g_MouseInput.StartSelection.X = x;
			g_MouseInput.StartSelection.Y = y;
		}
		else if (state == GLUT_UP && g_MouseInput.nSelectionState == SELECTION_ACTIVE)
		{
			// this is only a selection if mouse has moved
			if (x != g_MouseInput.StartSelection.X && y != g_MouseInput.StartSelection.Y)
			{
				g_MouseInput.nSelectionState = SELECTION_DONE;
            
				mouseInputCallSelection();
                
                iActiveSelection++;
                iActiveSelection = iActiveSelection % 2;
			}

			g_MouseInput.nSelectionState = SELECTION_NONE;
		}
	}
}

void mouseInputRegisterForSelectionRectangle(SelectionRectangleChangedPtr pFunc)
{
	g_MouseInput.pSelectionCallback = pFunc;
}

void mouseInputRegisterForCursorMovement(CursorMovedPtr pFunc)
{
	g_MouseInput.pCursorCallback = pFunc;
}

