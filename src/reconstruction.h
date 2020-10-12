#pragma once

#include "bwmorph.h"
#include <math.h>
#include <malloc.h>
#include "BALogix.h"

BOOL Walk2Edge(BAL_sImage* recon, BAL_PixelCoordinate* start, BAL_PixelCoordinate* end, UINT32 rowIncr, UINT32 colIncr);

BAL_sImage* Reconstruction(BAL_sImage* img, BAL_sImage* edgeMap);