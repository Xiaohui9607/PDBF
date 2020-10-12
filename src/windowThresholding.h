#pragma once


/*!
 \file   windowThresholding.h
 \brief  This file contains the prototype of the threshold algorithms to convert likelihood maps to edge maps
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

#include "BALogix.h"


BOOL WindowThresholdBeta(BAL_sImage *pData, double beta, BAL_sImage *pResult);
