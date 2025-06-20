#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <x86intrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif
/// mex CFLAGS='-march=native -o4 -ffast-math -fopenmp -mavx2 -Wall' -lgomp backProjLInv_BoxChunk.c

#define ERR_HEAD "*** forwardProjLInv_BoxChunk[mex]: "
#define ERR_ID   "CRossS:forwardProjLInv_BoxChunk:"

#define NUM_THREADS 8
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Check number and type of inputs and outputs: -----------------------
    if (nrhs != 7) {
        mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                          ERR_HEAD "7 inputs required.");
    }
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]) || !mxIsSingle(prhs[2]) ||
        !mxIsSingle(prhs[3]) || !mxIsSingle(prhs[4]) || !mxIsInt32(prhs[5])  || !mxIsInt32(prhs[6])) {
        mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput",
                          ERR_HEAD "Inputs must be a full real float arrays +2x int32.");
    }

    const float *dataIn  = (float *)mxGetPr(prhs[0]);
    const float *img     = (float *)mxGetPr(prhs[1]);
    const float *vXY     = (float *)mxGetPr(prhs[2]);
    const float *vYX     = (float *)mxGetPr(prhs[3]);
    const float *mm2vox  = (float *)mxGetPr(prhs[4]);
    const int   *boxSize = (int   *)mxGetPr(prhs[5]);
    const int   *nCore   = (int   *)mxGetPr(prhs[6]);

    const mwSize *sizeD      = mxGetDimensions(prhs[0]);
    mwSize sizeW[2];
    sizeW[0]  = 1;
    sizeW[1]  = sizeD[1];
    const mwSize *sizeImgX   = mxGetDimensions(prhs[2]);
    const mwSize *sizeImgY   = mxGetDimensions(prhs[3]);
    mwSize sizeImg[2];
    sizeImg[0] = sizeImgX[1];
    sizeImg[1] = sizeImgY[1];

    const int halfBox = *boxSize/2;
    const int halfImg = *sizeImg/2;

    plhs[0] = mxCreateNumericMatrix(sizeW[0], sizeW[1], mxSINGLE_CLASS, mxREAL);
    float *projData = (float *)mxGetPr(plhs[0]);
    for(int iD = 0; iD < sizeD[1] ; iD++)
        projData[iD] = 0;

    // Perform back-projection
    float xTmp;
    float yTmp;
    float nTmp;
    float iV11Tmp;
    float iV12Tmp;
    float iV21Tmp;
    float iV22Tmp;
    float dXY, dYX, arg;
    float expArg;

    int xStart = 0;
    int xEnd   = sizeImg[0]-1;
    int yStart = 0;
    int yEnd   = sizeImg[1]-1;

    // Divide the data into chunks for each thread
    int chunkSize = (int)ceil((float)sizeD[1] / (float)*nCore);
    int threadNum = *nCore;
    int remainder = sizeD[1] - threadNum*(chunkSize-1);
    // printf("\nThreads: %d(Chunk Size = %d, Remainder = %d)\n",threadNum,chunkSize,remainder);

    // Allocate thread arrays
    float** projVecs = (float**)malloc(threadNum * sizeof(float*));
    int arraySize = chunkSize;
    for (int iT = 0; iT < threadNum; iT++)
    {
        projVecs[iT] = (float*)malloc(arraySize * sizeof(float));
        for (int iP = 0; iP < arraySize; iP++)
            projVecs[iT][iP] = 0;
    }

    // Create and start threads
    #pragma omp parallel for num_threads(threadNum)
    for (int iT = 0; iT < threadNum; iT++)
    {
        int chunkStart = iT*chunkSize;
        int chunkEnd   = chunkStart + chunkSize;
        if (iT>=remainder)
        {
            chunkStart = chunkStart - (iT-remainder);
            chunkEnd   = chunkStart + chunkSize-1;
        }

        for( long iD=chunkStart ; iD<chunkEnd ; iD++){
            expArg = 0;
            xTmp = dataIn[sizeD[0]*iD+1];
            yTmp = dataIn[sizeD[0]*iD+2];
            nTmp = dataIn[sizeD[0]*iD+3];
            iV11Tmp = -0.5*dataIn[sizeD[0]*iD+4];
            iV12Tmp = -0.5*dataIn[sizeD[0]*iD+5];
            iV21Tmp = -0.5*dataIn[sizeD[0]*iD+6];
            iV22Tmp = -0.5*dataIn[sizeD[0]*iD+7];

            xStart = (int)round(yTmp*(*mm2vox)+(float)halfImg+0.5)-halfBox;
            yStart = (int)round(xTmp*(*mm2vox)+(float)halfImg+0.5)-halfBox;
            xEnd = xStart + *boxSize;
            yEnd = yStart + *boxSize;
            if ( xStart < 0 )
                xStart = 0;
            if ( yStart < 0 )
                yStart = 0;
            if ( xEnd >= sizeImg[0] )
                xEnd  = sizeImg[0]-1;
            if ( yEnd >= sizeImg[1] )
                yEnd  = sizeImg[1]-1;

            for(int iX=xStart; iX<xEnd ; ++iX){
                dYX = vYX[iX] - yTmp;
                for(int iY=yStart; iY<yEnd ; ++iY){
                    dXY = vXY[iY] - xTmp;
                    arg = iV11Tmp*dXY*dXY + iV22Tmp*dYX*dYX + (iV12Tmp + iV21Tmp)*dXY*dYX;
                    expArg += nTmp*img[sizeImg[0]*iX+iY]*exp(arg);
                } // end iY
            }     // end iX
            projVecs[iT][iD-chunkStart] = expArg;
        }         // end iD
    }             // end iT

    for (int iT = 0; iT < threadNum; iT++)
    {
        int chunkStart = iT*chunkSize;
        int chunkEnd   = chunkStart + chunkSize;
        if (iT>=remainder)
        {
            chunkStart = chunkStart - (iT-remainder);
            chunkEnd   = chunkStart + chunkSize-1;
        }        
        for( long iD=chunkStart ; iD<chunkEnd ; iD++){
            projData[iD] = projVecs[iT][iD-chunkStart];
        }
    }

    // Free allocated memory
    for (int iT = 0; iT < threadNum; iT++)
        free(projVecs[iT]);
    free(projVecs);

}


