#include "gaussian.h"
#include <string.h>

BOOL guassian(UINT32 kernel_size, double sigma, BAL_sImage* in, BAL_sImage* out) {
	if (kernel_size % 2 == 0){
        fprintf(stderr, "kernel_size must be odd number!\n");
        return false;
    }
    memset(out->scan0, 0, out->height * out->width);
	// create gaussian filter
	int length = kernel_size;
	double * filter = (double*)calloc(kernel_size*kernel_size, sizeof(double));
    int i, j, row, col;
    double sum = 0;
    for (i = 0; i < kernel_size; i++) {
        for (j = 0; j < kernel_size; j++) {
            double x = i - (kernel_size - 1) / 2.0;
            double y = j - (kernel_size - 1) / 2.0;
            filter[i*kernel_size+j] = exp(((pow(x, 2) + pow(y, 2)) / ((2 * pow(sigma, 2)))) * (-1));
            sum += filter[i * kernel_size + j];
        }
    }

    for (i = 0; i < kernel_size; i++) {
        for (j = 0; j < kernel_size; j++) {
            filter[i * kernel_size + j] /= sum;
        }
    }
    
	//smoothing
    for (col = 0; col < in->width; col++) {
        for (row = 0; row < in->height; row++) {
            double intensity = 0;
            for (i = -length /2; i < length /2+1; i++) {
                for (j = -length / 2; j < length /2+1; j++) {
                    intensity += (((col + i) < 0) || ((row + j) < 0) || ((row + j) >= in->height) || ((col + i) > in->width)) ?
                        0 : (double)(in->scan0[(col + i) * in->height + row + j]) * filter[(i + length / 2) * length + j + length / 2];
                }
            }
            out->scan0[col * out->height + row] = (BYTE)intensity;
        }
    }
    free(filter);
}