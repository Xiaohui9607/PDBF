#pragma once

#include "BALogix.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

BOOL guassian(UINT32 kernel_size, double sigma, BAL_sImage* in, BAL_sImage* out);