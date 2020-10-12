#pragma once


/*!
\file   graphics.h
\brief  This file contains common data types used by various 2D and 3D graphic algorithms
\author JAH
\date   11-12-2007
*/

/*
 * Copyright ?2009 BA Logix, Inc. All rights reserved. 
 * Use of this software is governed by the terms and conditions of the end user license agreement (EULA)
 * that accompanies this software. EXCEPT AS WARRANTED IN THE LICENSE AGREEMENT, BA LOGIX HEREBY 
 * DISCLAIMS ALL WARRANTIES AND CONDITIONS WITH REGARD TO THE SOFTWARE, INCLUDING ALL IMPLIED 
 * WARRANTIES AND CONDITIONS OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND 
 * NON-INFRINGEMENT.
 */

#include "BALogix.h"


/*!
\brief  A Structure to hold the pixel coordinate location
\author JAH
\date   09-27-2007
*/
typedef struct _PixelCoordinate {
	INT32 row;  //!< The horizontal row of the Pixel location
	INT32 col;  //!< The vertical column of the Pixel location
} PixelCoordinate;


/*!
\brief  A Structure to hold a two dimensional coordinate
\author JAH
\date   09-27-2007
*/
typedef struct _Coordinate2D {
	INT32 x;   //!< The x-axis coordinate
	INT32 y;   //!< The y-axis coordinate
} Coordinate2D;


/*!
\brief  A Structure to hold a three dimensional coordinate
\author JAH
\date   09-27-2007
*/
typedef struct _Coordinate3D {
	INT32 x;  //!< The x-axis coordinate
	INT32 y;  //!< The y-axis coordinate
	INT32 z;  //!< The z-axis coordinate
} Coordinate3D;


