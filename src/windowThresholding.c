/*!
\file   windowThresholding.c
\brief  This file contains the implementation of the threshold algorithms to convert likelihood maps to edge maps
\author JAH
\date   11-13-2007
*/

/*
 * Copyright ?2009 BA Logix, Inc. All rights reserved. 
 * Use of this software is governed by the terms and conditions of the end user license agreement (EULA)
 * that accompanies this software. EXCEPT AS WARRANTED IN THE LICENSE AGREEMENT, BA LOGIX HEREBY 
 * DISCLAIMS ALL WARRANTIES AND CONDITIONS WITH REGARD TO THE SOFTWARE, INCLUDING ALL IMPLIED 
 * WARRANTIES AND CONDITIONS OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND 
 * NON-INFRINGEMENT.
 */

#include "windowThresholding.h"
#include <malloc.h>
#include <string.h>
#include <stdio.h>



/*!
 \brief  Converts a likelihood map into an edge map
 \par    Detailed Description:
         This function takes a likelihood map and generates a binary
         edge map.
 \param  pData   - the data to be thresholded
 \param  beta    - the beta factor to the threshold
 \param  pResult -  The resultant thresholded data
 \return true if successful, false otherwise
 \remark This function requires that the input image and output images are of the same
         dimensions and data type (i.e. same number of bytes per pixel)
 \author JAH
 \date   2009-02-10
*/
BOOL WindowThresholdBeta(BAL_sImage *pData, double beta, BAL_sImage *pResult) {
	double   min,max, minmax;
	double   value;
	double   num, den;
	double   threshold = 0.0;
	UINT32   row,col;
	INT32    row2,col2;

	// Validate the input structure
	if (pData == NULL) {
		fprintf(stderr, "Input image structure is NULL\n");
		return false;
	}

	// Validate the ouptut structure
	if (pResult == NULL) {
		fprintf(stderr, "Output image structure is NULL\n");
		return false;
	}

	// Validate the input buffer
	if (pData->scan0 == NULL) {
		fprintf(stderr, "Input buffer is NULL\n");
		return false;
	}

	// Validate the output buffer
	if (pData->scan0 == NULL) {
		fprintf(stderr, "Output buffer is NULL\n");
		return false;
	}

	// Do we have data points to be thresholded?
	if ((pData->height == 0) || (pData->width == 0)) {
		fprintf(stderr, "Buffer dimensions must be greater than 0\n");
		return false;
	}

	// Make sure the input and output buffers are the same dimensions
	if ((pData->height != pResult->height) || (pData->width != pResult->width)) {
		fprintf(stderr, "The input and output buffers must have the same dimensions\n");
		return false;
	}


	UINT32* pfscan0 = (UINT32*)pData->scan0;
	threshold = (beta >= 0.0) ? beta : 0.0;
	for (col = 1; col < pData->width-1; col++) {
		for (row = 1; row < pData->height-1; row++) {
			min = MAXDOUBLE;
			max = MINDOUBLE;
			for (col2 = -1; col2 <= 1; col2++) {
				for (row2 = -1; row2 <= 1; row2++) {
					if ((row2 == 0) && (col2 == 0))
						continue;

					value = (double)(pfscan0[(col + col2) * pData->height + row + row2]);
					
					// Find the max neighbor
					if (value > max) max = value;

					// Find the min neighbor
					if (value < min) min = value;
				}
			}
			minmax = min + max;
			num = 2.0 * pfscan0[col * pData->height + row] - minmax;
			den = 2.0 * pfscan0[col * pData->height + row] + minmax;

			if ((den != 0.0) && ((num / den) > threshold)) {
				pResult->scan0[col * pData->height + row] = 1;
			}
		}
	}

	return true;
} // end WindowThresholdBeta

