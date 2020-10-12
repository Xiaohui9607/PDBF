#include "reconstruction.h"

BOOL Walk2Edge(BAL_sImage* recon, BAL_PixelCoordinate* p_start, BAL_PixelCoordinate* p_end, UINT32 rowIncr, UINT32 colIncr) {
	BOOL found = false;
	p_end->col = p_start->col;
	p_end->row = p_start->row;
	INT16* p_scan0 = (INT16*)recon->scan0;

	while (!found) {
		if ((p_end->row < 0) || (p_end->row >= (INT32)recon->height) ||
			(p_end->col < 0) || (p_end->col >= (INT32)recon->width)) 
			break;
		if (p_scan0[p_end->col * recon->height + p_end->row] > -1) {
			found = true;
		}
		else {
			p_end->row += rowIncr;
			p_end->col += colIncr;
		}
	}
	return found;
}



BAL_sImage* Reconstruction(BAL_sImage* img, BAL_sImage* edgeMap) {
	// assert img and edgemap has the same size
	INT16* rscan0;
	UINT32 row, col, nRoi;
	UINT32 length = 8;
	UINT32 n = 0;
	BOOL found, isBoolFlag;
	BAL_PixelCoordinate* rois = (BAL_PixelCoordinate*)calloc(length, sizeof(BAL_PixelCoordinate));
	double distance;
	double value;
	double sumValues, sumWeightedValues, sumDistances;
	UINT32 nDistances;
	INT32 rowIncr, colIncr;

	// Check if the unsigned char (Byte) array is black/white
	UINT16 max = 0;
	for (row = 0; row < img->height * img->width; row++) {
		max = (img->scan0[row] < max) ? max : img->scan0[row];
	}
	isBoolFlag = max == 1;

	// Dilate the edge map
	BAL_sImage* temp = Bw_Dilate(edgeMap);
	BAL_sImage* recon = BAL_NewImage(img->height, img->width, sizeof(INT16));
	BAL_sImage* result = BAL_NewImage(img->height, img->width, sizeof(BYTE));
	rscan0 = (INT16*)(recon->scan0);

	// recon = img * edgemap  type is INT16
	for (row = 0; row < img->height; row++) {
		for (col = 0; col < img->width; col++) {
			rscan0[col * img->height + row] = (INT16)img->scan0[col * img->height + row] * (INT16)temp->scan0[col * img->height + row];
		}
	}
	int i = 0;
	// Get roi list of all non-edges
	for (row = 0; row < recon->height; row++) {
		for (col = 0; col < recon->width; col++) {
			if (temp->scan0[col * temp->height + row] == 0) {
				i++;
				rscan0[col * recon->height + row] = -1;
				BAL_PixelCoordinate roi;
				roi.col = col;
				roi.row = row;
				if (n == length) {
					length *= 2;
					rois = realloc(rois, length * sizeof(BAL_PixelCoordinate));
				}
				rois[n].col = col;
				rois[n].row = row;
				n++;
			}
		}
	}
	//return temp;
	//if (n >256) {
	//	result->scan0[0] = n;
	//	//return result;
	//}

	// Matrix2Byte(recon, result);
	for (row = 0; row < recon->height; row++) {
		for (col = 0; col < recon->width; col++) {
			//result->scan0[col * result->height + row] = 111;
			result->scan0[col * result->height + row] = (BYTE)rscan0[col * recon->height + row];
		}
	}
	//return result;
	
	for (nRoi = 0; nRoi < n; nRoi++) {
		BAL_PixelCoordinate roi = rois[nRoi];
		BAL_PixelCoordinate end;
		//Initialize
		result->scan0[roi.col * result->height + roi.row] = 0;
		sumValues = 0.0;
		sumWeightedValues = 0.0;
		sumDistances = 0.0;
		nDistances = 0;

		for (rowIncr = -1; rowIncr < 2; rowIncr++) {
			for (colIncr = -1; colIncr < 2; colIncr++) {
				if ((rowIncr == 0) && (colIncr == 0))
					continue;
				found = Walk2Edge(recon, &roi, &end, rowIncr, colIncr);
				if (found) {
					nDistances++;
					distance = sqrt((double)((end.row - roi.row) * (end.row - roi.row) + 
						(end.col - roi.col) * (end.col - roi.col)));
					sumDistances += distance;
					sumValues += rscan0[end.col * recon->height + end.row];
					sumWeightedValues += rscan0[end.col * recon->height + end.row] * distance;
				}
			}
		}

		if (nDistances == 0) {
			value = 0.0;
		}
		else if (nDistances == 1) {
			value = sumValues;
		}
		else if (nDistances > 1) {
			value = (sumValues - sumWeightedValues / sumDistances) / (double)(nDistances - 1);
		}

		if (isBoolFlag) {
			result->scan0[roi.col*result->height+roi.row] = (value > 0.5) ? 1 : 0;
		}
		else {
			result->scan0[roi.col * result->height + roi.row] = (BYTE)value;
		}
	}
	
	free(rois);
	free(temp);
	free(recon);

	return result;
}