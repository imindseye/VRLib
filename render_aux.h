////////////////////////////////////////////////////////
//
//   Han-Wei Shen 
//   01/25/1996
//

#ifndef RENDER_AUX_H
#define RENDER_AUX_H

REAL Dot(uvw *v1, uvw *v2); 

void Normalize(uvw *v); 

REAL clamp(REAL value, REAL min, REAL max); 

REAL lerp(REAL t, REAL y0, REAL y1); 

REAL ipow( REAL x, int n); 

void matrix_mult(Matrix m, REAL opoint[], REAL npoint[]) ; 

int vrlib_invert_matrix(Matrix cm, Matrix inv); 

void 
  pivot(Matrix m, int c, int* R,int min, int max) ; 

#endif
