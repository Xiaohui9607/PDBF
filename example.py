import numpy as np
from ctypes import *
import cv2

lib = cdll.LoadLibrary('./build/libPDBF.so')

def detect(input, nbitplanes, beta, winsize, sigma, kernelsize, use_gaussian):
    use_gaussian = 1 if use_gaussian else 0
    res = np.zeros_like(input)

    lib.edgedetect(
        c_void_p(input.ctypes.data),  # input
        c_uint32(input.shape[0]),  # width
        c_uint32(input.shape[1]),  # height
        c_uint32(nbitplanes),  # nbitplanes
        c_float(beta),  # beta
        c_uint32(winsize),  # winsize
        c_float(sigma),  # sigma (smoothing)
        c_uint32(kernelsize),  # kernel size (smoothing)
        c_uint32(use_gaussian),  # use/not use gaussian
        c_void_p(res.ctypes.data)
    )
    return res

if __name__ == '__main__':
    input = cv2.cvtColor(cv2.imread("/home/golf/code/data/ED/Thermal/groundtruth/img_2.jpg"), cv2.COLOR_BGR2GRAY)
    edgemap = detect(input, nbitplanes=2, beta=0.0, winsize=2, sigma=1.0, kernelsize=5, use_gaussian=False)

    cv2.imshow("input", input)
    cv2.imshow("edgemap", edgemap*255)
    cv2.waitKey(0)
