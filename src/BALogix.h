//
// Created by golf on 20-5-16.
//

#ifndef PDBF_BALOGIX_H
#define PDBF_BALOGIX_H


typedef unsigned char          BYTE;
typedef unsigned char          UINT8;
typedef short int              INT16;
typedef unsigned short         UINT16;
typedef int                    INT32;
typedef unsigned int           UINT32;
typedef long long int          INT64;
typedef unsigned long long int UINT64;
typedef long                   BOOL;


#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


#ifndef MAXDOUBLE
#define MAXDOUBLE 1.7976931348623157e+308
#endif

#ifndef MINDOUBLE
#define MINDOUBLE -1.7976931348623157e+308
#endif


typedef struct _BAL_sImage {
    UINT32  height;  //!< The height of the image
    UINT32  width;   //!< The width of the image
    UINT32  bpp;     //!< bytes per pixel
    BYTE*   scan0;   //!< A pointer to the start of the first scanline
} BAL_sImage;

typedef struct _BAL_PixelCoordinate {
    UINT32  row;
    UINT32  col;
} BAL_PixelCoordinate;

/*!
 \brief  This enumerated type defines the image decomposition method
 \author JAH
 \date   2010-10-19
 */
typedef enum _BAL_DECOMPOSITION_METHOD {
    BAL_BITPLANE_DECOMPOSITION  = 0,
    BAL_FIBONACCI_DECOMPOSITION = 1
} BAL_DECOMPOSITION_METHOD;


const char*  BAL_GetVersion();
const UINT32 BAL_GetMajorVersion();
const UINT32 BAL_GetMinorVersion();
const UINT32 BAL_GetReleaseSerial();
BAL_sImage*  BAL_NewImage(UINT32 height, UINT32 width, UINT32 bpp);
void         BAL_DeleteImage(BAL_sImage *ptr);
BAL_sImage*  BAL_EdgeDetection(BAL_sImage *pImg, UINT32 nbitplanes, BAL_DECOMPOSITION_METHOD decompositionMethod, BOOL enableMorph, UINT32 p_code,
UINT32 n_code, double beta, UINT32 winsize, double sigma, UINT32 kernelsize, BOOL usegausian, char **msg);

#endif //PDBF_BALOGIX_H
