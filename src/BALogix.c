//
// Created by golf on 20-5-16.
//

#include "BALogix.h"
#include "bwmorph.h"
#include "windowThresholding.h"
#include <malloc.h>
#include <stdio.h>
#include "fibonaccilut.h"
#include <math.h>
#include "gaussian.h"
#include "reconstruction.h"
#include <time.h>
#include <string.h>

static const UINT32 BAL_MAJOR_VERSION  = 1;
static const UINT32 BAL_MINOR_VERSION  = 0;
static const UINT32 BAL_RELEASE_SERIAL = 0;
static const char*  BAL_VERSION_STR    = "1.0.0";


static const BYTE lookupTable2x2[16][2][2] = {
        {{ 0, 0}, { 0, 0}},
        {{ 0, 1}, { 1, 0}},
        {{ 1, 0}, { 0, 1}},
        {{ 1, 1}, { 0, 0}},
        {{ 1, 0}, { 0, 1}},
        {{ 1, 0}, { 1, 0}},
        {{ 1, 0}, { 0, 1}},
        {{ 1, 0}, { 0, 0}},
        {{ 0, 1}, { 1, 0}},
        {{ 0, 1}, { 1, 0}},
        {{ 0, 1}, { 0, 1}},
        {{ 0, 1}, { 0, 0}},
        {{ 0, 0}, { 1, 1}},
        {{ 0, 0}, { 1, 0}},
        {{ 0, 0}, { 0, 1}},
        {{ 0, 0}, { 0, 0}}
};


BOOL BAL_PartialDerivativesBinary(BAL_sImage *pPlane, BAL_sImage *pEdge, UINT32 winsize);
BOOL BAL_PartialDerivativesBinary2x2(BAL_sImage* pPlane, BAL_sImage* pEdge);
BOOL BAL_BitplaneDecomposition(BAL_sImage *pData, UINT32 nPlane, BAL_sImage *pPlane);
BOOL BAL_FibonacciDecomposition(BAL_sImage* pData, UINT32 nPlane, BAL_sImage* pPlane, UINT32 lutIndex);
BOOL BAL_BitplaneFusion(BAL_sImage *pLikelihoodMap, BAL_sImage *pFusion, UINT32 nPlane);

/*!
 \brief  Returns the version string of the BA Logix DLL
 \return the version string
 \author JAH
 \date   2009-03-15
 */
const char* BAL_GetVersion() {
    return BAL_VERSION_STR;
} // end BAL_GetVersion




/*!
 \brief  Returns the major version of the BA Logix DLL
 \return the major version
 \author JAH
 \date   2009-03-15
 */
const UINT32 BAL_GetMajorVersion() {
    return BAL_MAJOR_VERSION;
} // end BAL_GetMajorVersion



/*!
 \brief  Returns the minor version of the BA Logix DLL
 \return the minor version
 \author JAH
 \date   2009-03-15
 */
const UINT32 BAL_GetMinorVersion() {
    return BAL_MINOR_VERSION;
} // end BAL_GetMinorVersion



/*!
 \brief  Returns the release serial version of the BA Logix DLL
 \return the release serial
 \author JAH
 \date   2009-03-15
 */
const UINT32 BAL_GetReleaseSerial() {
    return BAL_RELEASE_SERIAL;
} // end BAL_GetReleaseSerial





/*!
 \brief  Allocates a new image structure
 \param  height  - the pixel height of the image
 \param  width   - the pixel width of the image
 \param  width   - the length in bytes of the scanline (padding is often used to achieve faster bit alignment)
 \param  bpp     - the number of bytes per pixel
 \return A pointer to the allocate image structure, NULL if unsuccessful
 \author JAH
 \date   2009-02-19
 */
BAL_sImage* BAL_NewImage(UINT32 height, UINT32 width, UINT32 bpp) {
BAL_sImage *pImg = NULL;

// Validate image pixel dimension
if ((height == 0) || (width == 0)) {
fprintf(stderr, "Image pixel dimension must be greater than 0!\n");
return NULL;
}

// Validate the number of bytes per pixel
if (bpp == 0) {
fprintf(stderr, "A pixel must be greater than equal to 1 byte\n");
return NULL;
}

// Validate the width is greater than equal to the image width in bytes

pImg = (BAL_sImage *) malloc(sizeof(BAL_sImage));
// Allocate the image structure
if (pImg == NULL) {
fprintf(stderr, "Failed to allocate image structure\n");
return NULL;
}

pImg->height = height;
pImg->width  = width;
pImg->bpp    = bpp;
pImg->scan0 = (BYTE*)calloc(pImg->height * pImg->width, bpp);
memset(pImg->scan0, 0, pImg->height * pImg->width * pImg->bpp);
// Allocate the scanline buffer
//pImg->scan0 = (BYTE *) calloc(pImg->height*pImg->width, sizeof(BYTE));

if (pImg->scan0 == NULL) {
fprintf(stderr, "Failed to allocate scanline buffer\n");
free(pImg);
return NULL;
}

return pImg;
} // end BAL_NewImage



/*!
 \brief  Deallocates an image structure
 \param  ptr - a pointer to the image structure to be deallocated
 \return None
 \author JAH
 \date   2009-02-19
 */
void BAL_DeleteImage(BAL_sImage *ptr) {
    // Validate user's parameter
    if (ptr == NULL)
        return;

    // Destroy the pixel buffer
    if (ptr->scan0 != NULL) {
        free(ptr->scan0);
        ptr->scan0 = NULL;
    }

    // Destroy the image structure
    free(ptr);
} // end BAL_DeleteImage




/*!
 \brief  Performs the fast Partial Derivative of binary data for a 2x2 window
 \par    Detailed Description:
  		 This function performs a fast implementation of the BOOLean partial derivative.
 		 This function assumes the caller precomputed all possible results for a given window size.
		 This precomputed data will serve as a fast lookup instead of calculating the partial derivative
		 over an entire BOOLean image.  For a 2x2 window, there are 2 variables and 16 possible partial values
		 within the 2x2 window for each variable.  The lookup table for a 2x2 window will utilize:
		 2 variables * 4 bytes per window * 16 solutions per variable = 128 bytes of storage.
 \param  pPlane - 2D array containing the binary image plane (pPlane[height][width])
 \param  pEdge  - The partial derivative results (pEdge[height][width])
 \return true if successful, false otherwise
 \remark This function expects the input and output image structures to have 1 byte per pixel.
         This function expects the input and output image structures to have the same width.
		 This function expects the input and output image structures to have the same dimensions.
 \author JAH
 \date   2009-02-10
*/
BOOL BAL_PartialDerivativesBinary2x2(BAL_sImage *pPlane, BAL_sImage * pLikelihood) {
    UINT32 row, col;
    UINT32 row2, col2;
    UINT32 index, maxIndex;
    UINT32 winSize = 2;
    UINT32 offset;
    UINT32 rowEnd, colEnd;
    BAL_sImage* block;
#pragma region TD
// Validate the input image structure
    if (pPlane == NULL) {
        fprintf(stderr, "The plane image structure is NULL\n");
        return false;
    }

    // Validate the output image structure
    if (pLikelihood == NULL) {
        fprintf(stderr, "Edge image structure is NULL\n");
        return false;
    }

    // Validate the input image buffer
    if (pPlane->scan0 == NULL) {
        fprintf(stderr, "The plane image buffer is NULL\n");
        return false;
    }

    // Validate the output image buffer
    if (pLikelihood->scan0 == NULL) {
        fprintf(stderr, "The edge image buffer is NULL\n");
        return false;
    }

    // Validate the input image structure has valid dimensions
    if ((pPlane->height == 0) || (pPlane->width == 0)) {
        fprintf(stderr, "The plane data must have dimensions greater than 0\n");
        return false;
    }

    // Validate the input and output image structures are of the same dimensions
    if ((pPlane->height != pLikelihood->height) || (pPlane->width != pLikelihood->width)) {
        fprintf(stderr, "The plane and edge image structures must have the same dimensions!\n");
        return false;
    }
#pragma endregion


    // Initialize
    rowEnd = pPlane->height - winSize + 1;
    colEnd = pPlane->width - winSize + 1;
    maxIndex = (1 << (winSize*winSize)) - 1;
    block = BAL_NewImage(winSize, winSize, pPlane->bpp);

    // Initialize the edge map to all zero's
    memset(pLikelihood->scan0, 0, pLikelihood->height* pLikelihood->width);

    //row = 100;
    for (row = 0; row < pPlane->height-1; row++) {
        for (col = 0; col < pPlane->width-1; col++) {
            //pLikelihood->scan0[col * pLikelihood->height + row] = pPlane->scan0[(col + 1) * pPlane->height + row + 1];
            memset(block->scan0, 0, block->height * block->width);
            index = 0;
            // Determine the lookup index
            for (row2 = 0; row2 < winSize; row2++)
                for (col2 = 0; col2 < winSize; col2++)
                    index = (index << 1) + pPlane->scan0[(col + col2) * pPlane->height + row + row2];
            // Skip the all zeros case
            if ((index == 0) || (index == maxIndex))
                continue;

            for (col2 = 0; col2 < winSize; col2++)
                for (row2 = 0; row2 < winSize; row2++)
                    block->scan0[col2 * winSize + row2] = pLikelihood->scan0[(col + col2) * pLikelihood->height + row + row2];

            for (row2 = 0; row2 < winSize; row2++)
                for (col2 = 0; col2 < winSize; col2++)
                    block->scan0[col2 * winSize + row2] = (block->scan0[col2 * winSize + row2] || lookupTable2x2[index][row2][col2]);
            for (row2 = 0; row2 < winSize; row2++)
                for (col2 = 0; col2 < winSize; col2++)
                    pLikelihood->scan0[(col + col2) * pLikelihood->height + row + row2] = block->scan0[col2 * winSize + row2];
        }
    }
    BAL_DeleteImage(block);

    return true;
} // end BAL_PartialDerivativesBinary

BOOL BAL_PartialDerivativesBinary(BAL_sImage* pPlane, BAL_sImage* pLikelihood, UINT32 winsize) {
    UINT32 row, col;
    UINT32 row2, col2;
    UINT32 index, maxIndex;
    UINT32 winSize = 2;
    UINT32 offset;
#pragma region DT
    // Validate the input image structure
    if (pPlane == NULL) {
        fprintf(stderr, "The plane image structure is NULL\n");
        return false;
    }

    // Validate the output image structure
    if (pLikelihood == NULL) {
        fprintf(stderr, "Edge image structure is NULL\n");
        return false;
    }

    // Validate the input image buffer
    if (pPlane->scan0 == NULL) {
        fprintf(stderr, "The plane image buffer is NULL\n");
        return false;
    }

    // Validate the output image buffer
    if (pLikelihood->scan0 == NULL) {
        fprintf(stderr, "The edge image buffer is NULL\n");
        return false;
    }

    // Validate the input image structure has valid dimensions
    if ((pPlane->height == 0) || (pPlane->width == 0)) {
        fprintf(stderr, "The plane data must have dimensions greater than 0\n");
        return false;
    }

    // Validate the input and output image structures are of the same dimensions
    if ((pPlane->height != pLikelihood->height) || (pPlane->width != pLikelihood->width)) {
        fprintf(stderr, "The plane and edge image structures must have the same dimensions!\n");
        return false;
    }

#pragma endregion
    // Initialize
    maxIndex = (1 << (winSize * winSize)) - 1;

    // Initialize the edge map to all zero's
    memset(pLikelihood->scan0, 0, pLikelihood->height * pLikelihood->width);

    // Calculate the partial derivatives via a lookup table
    for (col = 0; col < pPlane->width - 1; col++) {
        for (row = 0; row < pPlane->height - 1; row++) {
            // Determine the lookup index
            index = 0;
            for (row2 = 0; row2 < winSize; row2++) {
                offset = (row + row2) * pPlane->width;
                for (col2 = 0; col2 < winSize; col2++) {
                    index = (index << 1) + pPlane->scan0[offset + col + col2];
                }
            }

            // Skip the all zeros and all ones cases
            if ((index == 0) || (index == maxIndex))
                continue;
            else
                index--;

            // Logically fuse the lookup results to the full partial results
            for (row2 = 0; row2 < winSize; row2++) {
                offset = (row + row2) * pLikelihood->width;
                for (col2 = 0; col2 < winSize; col2++) {
                    pLikelihood->scan0[offset + col + col2] = pLikelihood->scan0[offset + col + col2] || lookupTable2x2[index][row2][col2];
                } // end for col2
            } // end for row2
        } // end for col
    } // end for row

    return true;
} // end BAL_PartialDerivativesBinary



/*!
 \brief  Decomposes the image data into a binary plane
 \param  pData  - the image data (8 bit) (pData[height][width])
 \param  nPlane - the plane number (15 MSB - 0 LSB)
 \param  pPlane - The decomposed image plane (pPlane[height][width])
 \return true if successful, false otherwise
 \author JAH
 \date   2009-02-10
 */
BOOL BAL_BitplaneDecomposition(BAL_sImage *pData, UINT32 nPlane, BAL_sImage *pPlane) {
    BYTE   mask    = 1;
    BYTE   maxBits = 8;
    UINT32 row,col;

    // Validate the data image structure
    if (pData == NULL) {
        fprintf(stderr, "Data image structure is NULL\n");
        return false;
    }

    // Validate teh plane image structure
    if (pPlane == NULL) {
        fprintf(stderr, "Plane image structure is NULL\n");
        return false;
    }

    // Validate the data image buffer
    if (pData->scan0 == NULL) {
        fprintf(stderr, "The data image buffer is NULL\n");
        return false;
    }

    // Validate the plane image buffer
    if (pPlane->scan0 == NULL) {
        fprintf(stderr, "The plane image buffer is NULL\n");
        return false;
    }

    // Validate the data image structure has valid dimensions
    if ((pData->height == 0) || (pData->width == 0)) {
        fprintf(stderr, "Data image structure must have dimensions greater than 0\n");
        return false;
    }

    // Validate the input and output image structures are of the same dimensions
    if ((pData->height != pPlane->height) || (pData->width != pPlane->width)) {
        fprintf(stderr, "The data and plane image structures must have the same dimensions\n");
        return false;
    }

    // Validate parameters
    if ((nPlane < 0) && (nPlane >= maxBits)) {
        fprintf(stderr, "The decomposed plane number is out of range\n");
        return false;
    }

    // Initialize mask used in splitting the image into bit planes
    mask = 1 << nPlane;

    // Decompose the image into the selected bit-plane
    for (row=0; row < pData->height; row++)
        for (col=0; col < pData->width; col++)
            pPlane->scan0[col*pPlane->height + row] = ((pData->scan0[col*pData->height + row] & mask) != 0);

    return true;
} // end BAL_BitplaneDecomposition*



BOOL BAL_FibonacciDecomposition(BAL_sImage* pData,UINT32 nPlane, BAL_sImage* pPlane, UINT32 lutIndex) {
    UINT32	mask = 1;
    UINT32  maxBits = nPlane;
    UINT32  value = 0;
    const UINT32* lut = NULL;
    UINT32 row, col;
    mask = 1 << nPlane;

    lut = fibonacciLUTs[lutIndex];
    // Decompose the image into the selected bit-plane
    for (row = 0; row < pData->height; row++)
        for (col = 0; col < pData->width; col++)
        {
            pPlane->scan0[col * pPlane->height + row] = ((lut[pData->scan0[col * pData->height + row]] & mask) != 0);
        }
    return true;
}





/*!
 \brief  Fuses the likelihood map for a specified decomposed plane based upon the specifified decomposition method
 \par    Detailed Description:
         This function fused the likelihood map for a specified decomposed plane into the fused results.  The fusion
		 techniqued used is based upon the decomposition method.
		 BIT_PLANE_DECOMPOSITION:
		     fusion is weighted by the decomposed bit-plane bit-position
		 OTHER:
		     fusion is an equal weight additive fusion
 \param  pData   - the resultant decomposed plane data to be fused
 \param  pFusion - the fusion output result (8 bit) (pFusion[height][width])
 \param  nPlane  - the decomposed plane number
 \return true if successful, false otherwise
 \author JAH
 \date   2009-02-10
 */
BOOL BAL_BitplaneFusion(BAL_sImage* pLikelihood, BAL_sImage* pFusion, UINT32 nPlane) {
    UINT32 nBits = 0;
    UINT32 row, col;

    if ((pLikelihood->height != pFusion->height) || (pLikelihood->width != pFusion->width) || (pLikelihood->height == 0) || (pLikelihood->width == 0)) {
        fprintf(stderr, "The data and fusion image structures must have the same dimensions\n");
        return false;
    }

    nBits = 8 * pFusion->bpp;
    UINT32* pfscan0 = (UINT32*)pFusion->scan0;

    if (nPlane >= nBits) {
        fprintf(stderr, "The bit plane number exceeds the number of available bit planes\n");
        return false;
    }
    for (row=0; row < pLikelihood->height; row++) {
        for (col=0; col < pLikelihood->width; col++) {
            pfscan0[col * pFusion->height + row] |= (pLikelihood->scan0[col * pLikelihood->height + row] << nPlane);
            //pfscan0[row*pFusion->width + col] |= (pLikelihood->scan0[row* pLikelihood->width + col] << (nBits-nPlane-1));
        } // end for col
    } // end for row
    return true;
}

//
//BAL_API BAL_sImage* BAL_Reconstruction(BAL_sImage* img, BAL_sImage* edgeMap) {
//	BAL_sImage* pRecon =  Reconstruction(img, edgeMap);
//	return pRecon;
//}



/*!
 \brief  Performs the BALogix's edge detection algorithm
 \par    Detailed Description:
         This function performs the BALogix's edge detection algorithm.
		 This function uses bitplane decomposition and beta thresholding.
 \param  pImg         - the input grayscale image data
 \param  planeMask    - A 1 in a bit-position indicates that decomposed plane is to be processed, 0 indicates skip the decomposed plane
 \param  enableMorph  - true enables morphological thinning, false disables it.
 \return The resultant edge map results
 \author JAH
 \date   2009-02-17
 */
BAL_sImage* BAL_EdgeDetection(BAL_sImage* pImg, UINT32 nbitplanes, BAL_DECOMPOSITION_METHOD decompositionMethod, BOOL enableMorph,
        UINT32 p_code, UINT32 n_code, double beta, UINT32 winsize, double sigma, UINT32 kernelsize, BOOL usegausian, char** msg) {
    UINT32 nPlane = 0;
    BOOL status = true;
    BAL_sImage * pPlane = NULL;
    BAL_sImage * pLikelihood = NULL;
    BAL_sImage * pFusedLikelihood = NULL;
    BAL_sImage * pEdge = NULL;
    BAL_sImage * pImgG = NULL;
    BAL_sImage * pInput = NULL;
    UINT32 row, col;
    BOOL errorFlag = false;
    UINT32 bitplaneMask = 0;
    UINT32 lutIndex = 0;

    clock_t begin = clock();


    if (pImg == NULL) {
        fprintf(stderr,
                "Input image structure is NULL\n");
        return NULL;
    }

    if (pImg->scan0 == NULL) {
        fprintf(stderr,
                "Input image buffer is NULL\n");
        return NULL;
    }

    if ((pImg->height == 0) || (pImg->width == 0)) {
        fprintf(stderr,
                "Data buffer must have dimension greater than 0\n");
        return NULL;
    }

    pImgG = BAL_NewImage(pImg->height, pImg->width, sizeof(BYTE));
    if (usegausian == 1) {
        guassian(kernelsize, sigma, pImg, pImgG
        );
        pInput = pImgG;
    } else {
        pInput = pImg;
    }


    pEdge = BAL_NewImage(pImg->height, pImg->width, sizeof(BYTE));
    if (pEdge == NULL) {
        fprintf(stderr,
                "Failed to allocate edge image structure\n");
        return NULL;
    }

    pPlane = BAL_NewImage(pImg->height, pImg->width, sizeof(BYTE));
    if (pPlane == NULL) {
        fprintf(stderr,
                "Failed to allocate plane image structure\n");

        BAL_DeleteImage(pEdge);

        return NULL;
    }

    pLikelihood = BAL_NewImage(pImg->height, pImg->width, sizeof(BYTE));

    if (pLikelihood == NULL) {
        fprintf(stderr,
                "Failed to allocate likelihood image structure\n");

        BAL_DeleteImage(pEdge);
        BAL_DeleteImage(pPlane);

        return NULL;
    }

    pFusedLikelihood = BAL_NewImage(pImg->height, pImg->width, sizeof(UINT32));
    if (pFusedLikelihood == NULL) {
        fprintf(stderr,
                "Failed to allocate fused likelihood image structure\n");

        BAL_DeleteImage(pEdge);
        BAL_DeleteImage(pPlane);
        BAL_DeleteImage(pLikelihood);

        return NULL;
    }

    if (decompositionMethod == BAL_BITPLANE_DECOMPOSITION) {
        bitplaneMask = 256 - 1;
        UINT32 i;
        for (
                i = 0;
                i < 8 -
                    nbitplanes;
                i++) {
            bitplaneMask -= (UINT32) pow(2, i);
        }
    }

    if (decompositionMethod == BAL_FIBONACCI_DECOMPOSITION) {
        char lutName[10];
        sprintf(lutName,
                "%u,%u", n_code, p_code);
        while ((
                       strcmp(lutName, fibonacciLUTNames[lutIndex]
                       ) != 0) && (lutIndex < 136))
            lutIndex++;
        bitplaneMask = (UINT32) pow(2, n_code) - 1;
        UINT32 i;
        for (
                i = 0;
                i < n_code -
                    nbitplanes;
                i++)
            bitplaneMask -= (UINT32) pow(2, i);
    }

    while ((bitplaneMask != 0) && (!errorFlag)) {
        if ((bitplaneMask & 0x0001) == 1) {

            switch (decompositionMethod) {
                case BAL_BITPLANE_DECOMPOSITION:
                    status = BAL_BitplaneDecomposition(pInput, nPlane, pPlane);
                    if (status == false) {
                        fprintf(stderr,
                                "Failed to decompose image plane\n");
                        errorFlag = true;
                    }
                    break;
                case BAL_FIBONACCI_DECOMPOSITION:
                    status = BAL_FibonacciDecomposition(pInput, nPlane, pPlane, lutIndex);
                    if (status == false) {
                        fprintf(stderr,
                                "Failed to decompose image plane\n");
                        errorFlag = true;
                    }
                    break;
                default:
                    fprintf(stderr,
                            "Unknown decomposition method\n");
                    errorFlag = true;
                    break;
            }

            if (errorFlag)
                break;


            status = BAL_PartialDerivativesBinary2x2(pPlane, pLikelihood);
            for (
                    row = 0;
                    row < pEdge->
                            height;
                    row++)

                if (status == false) {
                    fprintf(stderr,
                            "Failed to obtain a likelihood map\n");
                    errorFlag = true;
                }

            switch (decompositionMethod) {
                case BAL_BITPLANE_DECOMPOSITION:
                    status = BAL_BitplaneFusion(pLikelihood, pFusedLikelihood, nPlane);
                    if (status == false) {

                        errorFlag = true;
                    }
                    break;
                case BAL_FIBONACCI_DECOMPOSITION:
                    status = BAL_BitplaneFusion(pLikelihood, pFusedLikelihood, nPlane);
                    if (status == false) {

                        errorFlag = true;
                    }
                    break;
                default:
                    fprintf(stderr,
                            "Unknown decomposition method\n");
                    errorFlag = true;
                    break;
            }

            if (errorFlag)
                break;
        }

        nPlane++;
        bitplaneMask = bitplaneMask >> 1;
    } // end while

    memset(pPlane
                   ->scan0, 0, pPlane->
            height * pPlane
                                       ->width);
    if (status) {

        status = WindowThresholdBeta(pFusedLikelihood, beta, pPlane);

        if (status == true) {
            if (enableMorph) {

                status = BW_Thin(pPlane, pEdge);
            } else {

                for (
                        row = 0;
                        row < pEdge->
                                height;
                        row++)
                    for (
                            col = 0;
                            col < pEdge->
                                    width;
                            col++)
                        pEdge->scan0[
                                col * pEdge
                                        ->height + row] = pPlane->scan0[
                                col * pPlane
                                        ->height + row];
            }
        }
    }

    clock_t end = clock();

    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    * msg = (char *) malloc(100);

    sprintf(*msg,
            "%f", time_spent);

    BAL_DeleteImage(pPlane);
    BAL_DeleteImage(pLikelihood);
    BAL_DeleteImage(pFusedLikelihood);
    BAL_DeleteImage(pImgG);
    return pEdge;
}// end BAL_EdgeDetection