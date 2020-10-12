//
// Created by golf on 20-5-16.
//

#include <lzma.h>
#include "BALogix.h"
#include "reconstruction.h"

void edgedetect(const void * indatav, UINT32 width, UINT32 height, UINT32 nbitplanes, float beta, UINT32 winsize,
        float sigma, UINT32 kernelsize, UINT32 usegaussian, void * outdatav, UINT32 p_code, UINT32 n_code, UINT32 decomp) {

    //void cfun(const double * indata, int rowcount, int colcount, double * outdata) {
    const BYTE *indata = (BYTE *) indatav;
    BYTE *outdata = (BYTE *) outdatav;
    //BYTE *outrecon = (BYTE *) outreconv;
    BAL_sImage *pData, *pEdge, *pRecon;
    UINT32      bitplaneMask;
    BOOL        enableMorph;
//    UINT32      height, width, bpp;
    BYTE        grayValue;
    UINT32      sqrDim;
//    UINT32      startRow, endRow, startCol, endCol;
    UINT32      row, col;
//    UINT32      nOnes;
//    UINT32		p_code = 0;
//    UINT32		n_code = 8;
//    double		beta = 0.0;
//    UINT32		winsize = 2;
//    double		sigma = 1;
//    UINT32		kernelsize = 3;
//    UINT32		usegausian = 0;
    BAL_DECOMPOSITION_METHOD decompositionMethod;
    if (decomp==1)
        decompositionMethod = BAL_FIBONACCI_DECOMPOSITION;
    else
        decompositionMethod = BAL_BITPLANE_DECOMPOSITION;


    pData = pEdge = pRecon = NULL;
    enableMorph = TRUE;
    pData = BAL_NewImage(height, width, 1);
    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            pData->scan0[col *pData->height + row] = (BYTE)indata[col*pData->height + row] ; //grayValue;
        } // end for col
    }

    char* msg = NULL;

    pEdge = BAL_EdgeDetection(pData, nbitplanes, decompositionMethod, enableMorph, p_code, n_code, beta, winsize, sigma,kernelsize, usegaussian, &msg);


    for (row = 0; row < height; row++) {
        for (col = 0; col < width; col++) {
            outdata[col *pData->height + row] = (BYTE)pEdge->scan0[col*pData->height + row] ; //grayValue;
        } // end for col
    }

 

    BAL_DeleteImage(pData);
    BAL_DeleteImage(pEdge);
}
