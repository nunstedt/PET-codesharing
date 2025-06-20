/* getAttenuationMEX.c */

/* 
 * Takes the two matrices X1,X2 making up a list of detector pairs and computes
 * attenuation coefficients for each corresponding LOR over the image Im (mu map).
 * 
 * INPUT: X1,X2,FOVx,FOVy,slicesz,Im,nx,ny,nLines
 *
 *      X1      (single)    - n-by-2 vector of (x,y)-coordinates for detector 1
 *      X2      (single)    - n-by-2 vector of (x,y)-coordinates for detector 2
 *      FOVx    (single)    - Field-of-view in x-direction
 *      FOVy    (single)    - Field-of-view in y-direction
 *      slizesz (single)    - Thickness of slice in z-direction
 *      Im      (single)    - Image with positive y-values going north in image (normally use flipud(Im) as input from Matlab)
 *      nx      (int32)     - n.o. pixels in x-direction
 *      ny      (int32)     - n.o. pixels in y-direction
 *      nLines  (int32)     - n.o. LORs (n.o. columns in X1,X2)
 *
 * OUTPUT: outPut
 * 
 *      outPut  (single)    - n-by-1 vector of attenuation coefficients (exp[-int mu(s) ds])
 *                            corresponding to the order of X1,X2.
 *
 */


# include <mex.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <immintrin.h>
# include "matrix.h"
// # include <dispatch/dispatch.h>
// #ifdef _OPENMP
//     #include<omp.h>
// #endif


/* gateway function, MEX */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
    /* check input types */
    if (!mxIsSingle(prhs[0])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input X1 must be of type single.");
    if (!mxIsSingle(prhs[1])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input X2 must be of type single.");
    if (!mxIsSingle(prhs[2])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input FOVx must be of type single.");
    if (!mxIsSingle(prhs[3])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input FOVy must be of type single.");
    if (!mxIsSingle(prhs[4])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input slizesz must be of type single.");
    if (!mxIsSingle(prhs[5])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input Im must be of type single.");
    if (!mxIsInt32(prhs[6])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input nx must be of type int32.");
    if (!mxIsInt32(prhs[7])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input ny must be of type int32.");    
    if (!mxIsInt32(prhs[8])) 
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputType", "Input nLines must be of type int32.");
    
    
    /* declare variables */
    float*  X1;      // detection points X1,X2
    float*  X2;
    float*  FOVx;
    float*  FOVy;
    float*  slizesz; // slice thickness, z direction
    int*    nx;        // # pixels x,y,z dir
    int*    ny;
    int*    nz;    
    float*  Im;      // Image, flipped up-down for each slice
    int*    nLines;    // # "LORs"    
    float   x1;
    float   x2;
    float   y1;
    float   y2;
    
    /* get inputs */
    X1      = mxGetSingles(prhs[0]);
    X2      = mxGetSingles(prhs[1]);
    FOVx    = mxGetSingles(prhs[2]);
    FOVy    = mxGetSingles(prhs[3]);
    slizesz = mxGetSingles(prhs[4]);
    Im      = mxGetSingles(prhs[5]);
    nx      = mxGetInt32s(prhs[6]);
    ny      = mxGetInt32s(prhs[7]);
    nLines  = mxGetInt32s(prhs[8]);    
    
    /* create outputs */
    plhs[0] = mxCreateNumericMatrix(nLines[0],1,mxSINGLE_CLASS,mxREAL);
    float *outPut;
    outPut = mxGetSingles(plhs[0]);    
    
    /* center image and set normalized step sizes */
    float nx2 = (float) *nx/2;
    float ny2 = (float) *ny/2;
    float dx  = (float) *FOVx / *nx; // [mm/pixel]
    float dy  = (float) *FOVy / *ny; // [mm/pixel]
        
    
    /* do the stuff */
    
    
    // for each line in lines
    for (int lx=0; lx<*nLines; ++lx){
        
        // normalize inputs to unit [pixel]
        x1 = X1[lx]/dx;                 // [mm/(mm/pixel) = pixel]
        x2 = X2[lx]/dx;
        y1 = X1[lx + *nLines]/dy;
        y2 = X2[lx + *nLines]/dy;
        

        /* TRAVERSE COLUMNS CASE */
        
        if (dy/dx > fabs(X1[lx + *nLines]-X2[lx + *nLines]) / fabs(X1[lx]-X2[lx]) ) {
                        
            float ky = (y2-y1)/(x2-x1);            
            
            /* DO THIS IN PARALLEL! */
            for (int ix=0; ix<*nx; ++ix){
                
                float yy1;
                float yy2;
                
                float xx1 = ix-nx2;
                float xx2 = xx1+1;

                if (ky >= 0){
                    yy1 = y1+ky*(xx1-x1)+ny2;
                    yy2 = y1+ky*(xx2-x1)+ny2;
                    }
                else {
                    yy1 = y1+ky*(xx2-x1)+ny2;
                    yy2 = y1+ky*(xx1-x1)+ny2;
                    }
                
                float c1 = floorf(yy1);
                float c2 = floorf(yy2);

                if (c2 == c1){
                    if (c1 >= 0 && c1 <= *ny-1)
                        {
                        int iy  = (int) c1;
                        outPut[lx] += sqrt(dx*dx+(dy*ky)*(dy*ky)) * Im[ix*(*ny)+iy];
                        }
                    }
                else
                    {
                    if (c1 >= -1 && c1 <= *ny-1)
                        {
                        float len = sqrtf(dx*dx +(dy*ky)*(dy*ky));
                        if (c1 >= 0)
                            {
                            int iy  = (int) c1;
                            outPut[lx] += (c2-yy1)/(yy2-yy1)*len * Im[ix*(*ny)+iy];
                            }
                        if (c2 <= *ny-1)
                            {
                            int iy  = (int) c2;
                            outPut[lx] += (yy2-c2)/(yy2-yy1)*len * Im[ix*(*ny)+iy];
                            }
                        }
                     }
                
                } // end parallel
        }
        
        
        
        
        
        /* TRAVERSE ROWS CASE */
        
        
        if (dy/dx <= fabs(X1[lx + *nLines]-X2[lx + *nLines]) / fabs(X1[lx]-X2[lx]) ) {
            
            float kx = (x2-x1)/(y2-y1);
            
            /* DO THIS IN PARALLEL! */
            for (int iy=0; iy<*ny; ++iy){
                
                float xx1;
                float xx2;
                
                float yy1 = iy-ny2;
                float yy2 = yy1+1;

                if (kx >= 0){
                    xx1 = x1+kx*(yy1-y1)+nx2;
                    xx2 = x1+kx*(yy2-y1)+nx2;
                    }
                else {
                    xx1 = x1+kx*(yy2-y1)+nx2;
                    xx2 = x1+kx*(yy1-y1)+nx2;
                    }
                
                float c1 = floorf(xx1);
                float c2 = floorf(xx2);
                
                if (c2 == c1){
                    if (c1 >= 0 && c1 <= *nx-1)
                        {
                        int ix  = (int) c1;
                        outPut[lx] += sqrt(dy*dy+(dx*kx)*(dx*kx)) * Im[ix*(*ny)+iy];
                        }
                    }
                else
                    {
                    if (c1 >= -1 && c1 <= *nx-1)
                        {
                        float len = sqrtf(dy*dy +(dx*kx)*(dx*kx));
                        if (c1 >= 0)
                            {
                            int ix  = (int) c1;
                            outPut[lx] += (c2-xx1)/(xx2-xx1)*len * Im[ix*(*ny)+iy];
                            }
                        if (c2 <= *ny-1)
                            {
                            int ix  = (int) c2;
                            outPut[lx] += (xx2-c2)/(xx2-xx1)*len * Im[ix*(*ny)+iy];
                            }
                        }
                     }
            
        
            }
            
        } // end columns case, parallel
            
    } // end for line in lines
    
    
    for (int ix=0; ix<*nLines; ++ix){
        outPut[ix] = expf(-outPut[ix]);
        }
   

} // end mexFunction





