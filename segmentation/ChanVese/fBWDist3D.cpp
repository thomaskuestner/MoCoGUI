// ========================================================================
// ***
// *** fBWDist3D.cpp
// ***
// *** A fast 2D/3D region growing algorithm for Matlab. See the Matlab
// *** wrapper file 'RegionGrowing.m' for details on its usage.
// ***
// *** Compile this file by making the directiory containing this file
// *** your current Matlab working directory and typing 
// ***
// *** >> mex fBWDist3D.cpp
// ***
// *** in the Matlab console.
// ***
// *** Copyright 2013 Christian Wuerslin, University of Tuebingen and
// *** University of Stuttgart, Germany.
// *** Contact: christian.wuerslin@med.uni-tuebingen.de
// ***
// ========================================================================

#include "mex.h"
#include <queue>
#include <cmath>

#define NQUEUES             500	// number of FIFO queues

using namespace std;

long             lNZ, lNX, lNY;     // The image dimensions;
bool            *pbMask;             // pointer to the image
queue<long>     *aqQueue;			// vector of FIFO queues
long             lQueueWriteInd = 0;

// ========================================================================
// Inline function to determin minimum of two numbers
template<class cType>
inline cType ifMin(cType a, cType b)
{
    return a < b ? a : b;
}
// ========================================================================



// ========================================================================
// Inline function to determin maximum of two numbers
template<class cType>
inline cType ifMax(cType a, cType b)
{
    return a > b ? a : b;
}
// ========================================================================



// ========================================================================
// ***
// *** FUNCTION fPop
// ***
// *** Function that pops a voxel location from the highest priority
// *** non-empty queue which fulfills the region growing criterion.
// ***
// ========================================================================
long fPop() {
	long lInd; // Index of the voxel
    
    // --------------------------------------------------------------------
    // Loop over the queues, start with highest priority (0)
    for (int iQueueInd = 0; iQueueInd < NQUEUES; iQueueInd++) {
        
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // While there are still entries in the queue, pop and determine
        // whether it fullfills the region growing criterion.
        while (!aqQueue[iQueueInd].empty()) {
            lInd = aqQueue[iQueueInd].front(); aqQueue[iQueueInd].pop();
            return lInd;// Return if valid entry found
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    }
    // --------------------------------------------------------------------
    
	return -1; // if all queues are empty
}
// ========================================================================
// *** END OF FUNCTION fPop
// ========================================================================



// ========================================================================
// ***
// *** FUNCTION fGetNHood
// ***
// *** Get the 6-neighbourhood of voxel lLinInd
// ***
// ========================================================================
void fGetNHood(const long lLinInd, long& lNHoodSize, long* lNHood) {
	long lX, lY, lZ, lTemp;
    
    lNHoodSize = 0;
    
    lY = lLinInd % lNY; // get y coordinate, add +/-1 if in image range
	if (lY >       0) lNHood[lNHoodSize++] = lLinInd - 1;
	if (lY < lNY - 1) lNHood[lNHoodSize++] = lLinInd + 1;
    
    
	lTemp = lLinInd/lNY; // That is a floor() operation in c
	lX = lTemp % lNX; // get x coordinate, add +/-1 if in image range. X increment is lNY
	if (lX >       0) lNHood[lNHoodSize++] = lLinInd - lNY;
	if (lX < lNX - 1) lNHood[lNHoodSize++] = lLinInd + lNY;
    
    if (lNZ > 1) { // 3D case
        lZ = lTemp/lNX; // z coordinate, add +/-1 if in image range. Z increment is lNX*lNY
        if (lZ >       0) lNHood[lNHoodSize++] = lLinInd - lNX*lNY;
        if (lZ < lNZ - 1) lNHood[lNHoodSize++] = lLinInd + lNX*lNY;
    }
}
// ========================================================================
// *** END OF FUNCTION fGetNHood
// ========================================================================



// ========================================================================
// ***
// *** MAIN MEX FUNCTION RegionGrowing_mex
// ***
// *** See m-file for description
// ***
// ========================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    // --------------------------------------------------------------------
    // Check the number of the input and output arguments.
    if(nrhs != 1) mexErrMsgTxt("Exactly one input argument required.");
    if(nlhs != 1) mexErrMsgTxt("Exactly one output argument required.");
    // --------------------------------------------------------------------
    
    // --------------------------------------------------------------------
    // Get pointer/values to/of the input and outputs objects
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // 1st input: Image (get dimensions as well)
    if (!mxIsLogical(prhs[0])) mexErrMsgTxt("First input argument must be of type logical.");
    pbMask = (bool*) mxGetData(prhs[0]);
    
    const int* pSize = mxGetDimensions(prhs[0]);
    long lNDims = mxGetNumberOfDimensions(prhs[0]);
    lNY = long(pSize[0]);
	lNX = long(pSize[1]);
    if (lNDims == 3) lNZ = long(pSize[2]); else lNZ = 1;
    long lImSize = lNX*lNY*lNZ;
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    // Get pointer to output arguments and allocate memory for the corresponding objects
    plhs[0] = mxCreateNumericArray(lNDims, pSize, mxDOUBLE_CLASS, mxREAL);	// create output array
    bool *pdDist = (bool*) mxGetData(plhs[0]);						// get data pointer to distance tranform
    bool *pbCandidate = (bool*) mxGetData(mxCreateNumericArray(lNDims, pSize, mxLOGICAL_CLASS, mxREAL));

    // --------------------------------------------------------------------
    // Start of the real functionality
    
    aqQueue = new queue<long>[NQUEUES];
    
    long    lNHoodSize;
    long    alNHood[26];
    long    lInd;
    
    // --------------------------------------------------------------------
    // Initialize
    for (long lI = 0; lI < lImSize; lI++) {
        if (pbMask[lI]) {
            pbCandidate[lI] = true;
            for (long lY = 
                
            }
        }
    }
    // --------------------------------------------------------------------
    
    // --------------------------------------------------------------------
   /* while (lRegSize < lImSize) {

        fGetNHood(lLinInd, lNHoodSize, alNHood);
        for (int iI = 0; iI < lNHoodSize; iI++) {
            lLinInd = alNHood[iI];
            if (pbCandidate[lLinInd]) continue;
            
            pbCandidate[lLinInd] = true;
            if (lQueueInd > NQUEUES - 1) lQueueInd = NQUEUES - 1;
            aqQueue[lQueueInd].push(lLinInd);
        }
        
        lLinInd = fPop();
        if (lLinInd < 0) return;
        
        pbMask[lLinInd] = true;
        dRegMean = (dRegMean*double(lRegSize) + pdImg[lLinInd]);
        lRegSize++;
        dRegMean = dRegMean/double(lRegSize);
    }*/
    // End of while loop
    // --------------------------------------------------------------------
}
// ========================================================================
// *** END OF MAIN MEX FUNCTION RegionGrowing_mex
// ========================================================================