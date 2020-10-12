cmake_minimum_required(VERSION 3.14)
project(PDBF C)

set(CMAKE_C_STANDARD 99)

add_library(PDBF SHARED BALogix.h BALogix.c bwmorph.h bwmorph.c graphics.h
        windowThresholding.h windowThresholding.c gaussian.h gaussian.c
        fibonaccilut.h reconstruction.c reconstruction.h edgedetect.c)