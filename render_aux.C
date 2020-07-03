////////////////////////////////////////////////////////
//
//   Han-Wei Shen 
//   01/25/1996
//

#ifndef RENDER_AUX
#define RENDER_AUX

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>

#include "render.h"
#include "Map.h"
#include "Trans_Stack.h"

static Matrix identity = {1.0, 0.0, 0.0, 0.0,
		   0.0, 1.0, 0.0, 0.0,
		   0.0, 0.0, 1.0, 0.0,
		   0.0, 0.0, 0.0, 1.0};

//////////////////////////////////////////
//
//  Some auxilary functions 
//
//

// Vector Dot product
//
REAL Dot(uvw *v1, uvw *v2)
{
  return (v1->u*v2->u + v1->v*v2->v + v1->w*v2->w);
}

//////////////////////////////////
//
// Normalize a vector, avoiding work if possible 
// (e.g. already normalized, or zero vector)
//
void Normalize(uvw *v)
{
REAL magn = Dot(v,v);
if (magn > EPS && fabs(magn-1.0) > EPS)
  {
      magn = sqrt(magn);
      v->u /=magn;
      v->v /=magn;
      v->w /=magn;
   }
}
///////////////////////////////////
//
// Clamp a value  to ensure it lies between min and max
//
REAL clamp(REAL value, REAL min, REAL max)
{
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

///////////////////////////////////
//
// Interpolate t of the way between y0 and y1.
//
REAL lerp(REAL t, REAL y0, REAL y1)
{
  return y0 + t*(y1-y0);
}

///////////////////////////////////
//
// Raise a REAL to an integer power (for phong highlight) 
//
REAL ipow( REAL x, int n)
{
int i;
REAL result = 1, factor = x;

// loop through bits of n to give O(log(n)) time bound 
for( i=1; i<= n; i<<= 1)
 {
   if( i & n) result *= factor;   // i is part of n's base 2 representation
   factor *= factor;              // factor is now X^( log(i)+1)
 }
  return result;
}

////////////////////////////////////////
//
//  Matrix-Vector multiplication
//
void matrix_mult(Matrix m, REAL opoint[], REAL npoint[]) 
{
  
  npoint[0] = 
    (opoint[0] * m[0][0]) + 
      (opoint[1] * m[1][0]) +
	(opoint[2] * m[2][0]) +
	  (opoint[3] * m[3][0]);

  npoint[1] = 
    (opoint[0] * m[0][1]) +
      (opoint[1] * m[1][1]) +
	(opoint[2] * m[2][1]) +
	  (opoint[3] * m[3][1]);

  npoint[2] = 
    (opoint[0] * m[0][2]) +
      (opoint[1] * m[1][2]) +
	(opoint[2] * m[2][2]) +
	  (opoint[3] * m[3][2]);

  npoint[3] = 
    (opoint[0] * m[0][3]) +
      (opoint[1] * m[1][3]) +
	(opoint[2] * m[2][3]) +
	  (opoint[3] * m[3][3]);
}

//////////////////////////////////////////
//
//  Inverse a Matrix 
//
int vrlib_invert_matrix(Matrix cm, Matrix inv)
{
  Matrix  src, result;
  int c, r, j, k;	
  int R[4];
  int i;
  void pivot(Matrix m, int c, int *R,int min, int max);   
//
//  We initialize the result matrix to the identity, then
//  do Gaussian elimination on the source matrix, applying
//  the row operations to the result matrix at the same time.
//  When we're done, the source contains the identity and the
//  result is the inverse.
//

// Initialize result 
memcpy(result,identity, sizeof(Matrix) );  
  
// Copy source, since this is destructive 
  memcpy(src, cm, sizeof(Matrix));
  
for(i=0; i<4; i++) R[i] = i;	/* Initialize row table */
  
  
  /* This is the forward elimination part */
  for(r=0, c=0; c<4 && r<4; r++, c++)   /* For each row ... */
   {
     /* Find a pivot row and swap rows */
     pivot(src,c,R,r,4);
     
     /* If the diagonal entry is non-zero then use it to eliminate entries in
       this column in lower rows.
      */
#define epsilon 1E-06
     if (fabs(src[R[r]][c]) > epsilon)
      {
	int    row = R[r];
	double div = src[row][c];
	
	for(j=c+1; j<4; j++)
	  src[row][j] /= div;
	for(j=0; j<4; j++)
	  result[row][j] /= div;
	src[row][c] = 1.;
	for(j=r+1; j<4; j++)
	 {
	   double mult = src[R[j]][c];
	   
	   if (fabs(mult) > 0) {
	    for(k=c+1; k<4; k++)
	      src[R[j]][k] -= mult*src[row][k];
	    src[R[j]][c] = 0.0;
	    for(k=0; k<4; k++)
	      result[R[j]][k] -= mult*result[row][k];
	  }
	 }
      }
     else
      {
	fprintf(stderr,"vrlib_invert_matrix: Attempt to invert a badly conditioned transform\n");
	return 0;
      }
   }
  
  /* Now we do back substitution */
  for(r=3,c=3; r>=0; r--,c--)
   {
     for(j=r-1; j>=0; j--)
      {
	double mult = src[R[j]][c];
	if (fabs(mult) > 0.)
	  for(k=0; k<4; k++)
	    result[R[j]][k] -= mult*result[R[r]][k];
      }
   }
  /* Reorder the rows */
  for(r=0; r<4; r++)
   {
     inv[r][0] = result[R[r]][0];
     inv[r][1] = result[R[r]][1];
     inv[r][2] = result[R[r]][2];
     inv[r][3] = result[R[r]][3];
   }
  
  return 1;
}

/////////////////////////////////////////////////////
//
// Locate the row in A with the largest value in the 
//  c'th column.  Only columns between min and max-1 are
//  considered.
//
void 
  pivot(Matrix m, int c, int* R,int min, int max) 
{
  int r, piv, tmp;
  
  piv = min;
  for(r=min+1; r<max; r++)
    if (fabs(m[R[r]][c]) > fabs(m[R[piv]][c]))
      piv=r;
  if (piv != min) 
   {
     tmp = R[piv]; R[piv] = R[min]; R[min] = tmp;
   }
  return;
}

#endif
