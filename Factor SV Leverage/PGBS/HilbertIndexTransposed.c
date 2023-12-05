#include "mex.h"

#include <stdint.h> 
#include <string.h>


void HilbertIndexTransposed(uint32_t y[], int nbits, int ndimensions, uint32_t X[]) {
  /* in matlab the function is called like this
   * initialize the coordinates using uint32 - very important 
   * x5 = [uint32(5),uint32(10),uint32(20)]
   * HilbertIndexTransposed(5,3,x5)
   * Number of inputs 3 number of dimensions 3 number of bits 5
   * The three coordinates of the input are 5 10 20
   * ans =
   *
   *       10          14          27
   *
   * To compile the function  within the matlab commandline
   * mex -v HilbertIndexTransposed.c
   */
  
	uint32_t M,P, Q, t;
	int i;

  M=1<<(nbits-1);

  /* Inverse undo */
  for( Q = M; Q > 1; Q >>= 1 ) {
   
    P = Q - 1;
    for( i = 0; i < ndimensions; i++ ) {

      if( X[i] & Q ) 
        X[0] ^= P;
      else {
        t = (X[0]^X[i]) & P; 
        X[0] ^= t; 
        X[i] ^= t; 
      } 
    } 
  } /* exchange */

  /* Gray encode */
  for( i = 1; i < ndimensions; i++ ) 
        X[i] ^= X[i-1]; 
      
  t = 0;
  for( Q = M; Q > 1; Q >>= 1 ) {

      if( X[ndimensions-1] & Q ) 
        t ^= Q-1; 
  }
  for( i = 0; i < ndimensions; i++ ) 
    X[i] ^= t;

  for( i = 0; i < ndimensions; i++ ) 
    y[i] = X[i];

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  uint32_t *X;
  int nbits,ndimensions;
  uint32_t *hilbertAxes, *output;

  
	if (!mxIsClass(prhs[2], "uint32")) {
		mexErrMsgTxt("Error in HilbertIndexTransposed: expected first argument to be type uint32.");
	}
  
  /* first input.. the number of bits */
  nbits = mxGetScalar(prhs[0]);
  /* second input is the number of dimensions */
  ndimensions = mxGetScalar(prhs[1]);
  /* the actual coordinates are the third argument */
	hilbertAxes = (uint32_t *)mxGetPr(prhs[2]);

  /*mexPrintf("Number of inputs %d number of dimensions %d number of bits %d\n",nrhs,ndimensions,nbits);*/
  /*mexPrintf("The three coordinates of the input are %d %d %d\n",(int)hilbertAxes[0],(int)hilbertAxes[1],(int)hilbertAxes[2]);*/
  
	X = (uint32_t *)malloc(ndimensions * sizeof(uint32_t));
	memcpy(X, hilbertAxes, ndimensions * sizeof(uint32_t));

  /* allocate memory for the answer */
 	plhs[0] = mxCreateNumericMatrix(1,ndimensions,mxUINT32_CLASS, mxREAL);

  output = (uint32_t *)mxGetPr(plhs[0]);
  HilbertIndexTransposed(output,nbits,ndimensions,X);

	free(X);
}
