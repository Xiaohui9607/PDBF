#pragma once

/*!
 \file   bwmorph.h
 \brief  This file contains code that perform morphological operations on binary images
 \author JAH
 \date   11-19-2007
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
#include "graphics.h"

BOOL BW_Thin(BAL_sImage *pData, BAL_sImage *pResult);



BAL_sImage* Bw_Dilate(BAL_sImage* pData);