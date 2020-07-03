/*
 * Trans_Stack.h - This is the Trans_Stack class file.  It is used by
 * other classes to provide transformations.
 *
 * Written by Mike Krogh, NCSA, 3/9/92
 *
 */

#include <stdio.h>
#include <math.h>
#include "Trans_Stack.h"


Trans_Stack::Trans_Stack() {
/* This subroutine allocates memory for the transformation stack and fills in
 * the top with identity (for safety's sake).  The user should call loadmatrix
 * after this subroutine to place something on the stack.
 */

  register int i,j;

  /* define the current location in the stack; */
  TransStackCurrent = 0;

  
  /* get the memory for the stack */
  /*
  if ((TransStack = new Matrix[TransStackMax]) == NULL) {
    printf("Trans_Stack: error, not enough memory for TransStack\n");
    return;
  }
  */

  /* identity for the current matrix */
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      if (i==j)
	(TransStack[TransStackCurrent])[i][j] = 1.0;
      else
	(TransStack[TransStackCurrent])[i][j] = 0.0;

}




Trans_Stack::~Trans_Stack(void) {

  // delete TransStack;

}



/* returns the top matrix */
void Trans_Stack::getmatrix(Matrix m) {

  register int i,j;

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      m[i][j] = (TransStack[TransStackCurrent])[i][j];

}




/* replaces the top matrix with m*top */
void Trans_Stack::multmatrix(Matrix m) {

  register int i,j;
  Matrix tmp;

  /* make a copy of the top of the stack */
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      tmp[i][j] = (TransStack[TransStackCurrent])[i][j];

/*
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      (*(TransStack+TransStackCurrent))[i][j]
        = m[i][0] * tmp[0][j] +
	  m[i][1] * tmp[1][j] +
	  m[i][2] * tmp[2][j] +
	  m[i][3] * tmp[3][j];
*/

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      (TransStack[TransStackCurrent])[i][j]
        = tmp[i][0] * m[0][j] +
	  tmp[i][1] * m[1][j] +
	  tmp[i][2] * m[2][j] +
	  tmp[i][3] * m[3][j];

}




/* multiplies a 3D point by the Matrix */
void Trans_Stack::multpoint3d(REAL opoint[3], REAL npoint[3]) {

  int i;
  REAL tmp[4];

  for (i=0;i<4;i++) 
    tmp[i] = opoint[0] * (TransStack[TransStackCurrent])[0][i]
      + opoint[1] * (TransStack[TransStackCurrent])[1][i]
	+ opoint[2] * (TransStack[TransStackCurrent])[2][i]
	  + (TransStack[TransStackCurrent])[3][i];

  for (i=0;i<3;i++)
    npoint[i] = tmp[i] / tmp[3];

}




/* multiplies a 4D point by the Matrix */
void Trans_Stack::multpoint4d(REAL opoint[4], REAL npoint[4]) {

  int i;
  REAL tmp[4];

  for (i=0;i<4;i++) 
    tmp[i] = opoint[0] * (TransStack[TransStackCurrent])[0][i]
      + opoint[1] * (TransStack[TransStackCurrent])[1][i]
	+ opoint[2] * (TransStack[TransStackCurrent])[2][i]
	  + opoint[3] * (TransStack[TransStackCurrent])[3][i];

  for (i=0;i<4;i++)
    npoint[i] = tmp[i];

}




/* replaces the top matrix */
void Trans_Stack::loadmatrix(Matrix m) {

  register int i,j;

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      (TransStack[TransStackCurrent])[i][j] = m[i][j];

}




/* pushs the top matrix */
void Trans_Stack::pushmatrix(void) {

  register int i,j;

  if (TransStackCurrent >= TransStackMax-1) {
    printf("Trans_Stack: error, can't pushmatrix, stack is full.\n");
    return;
  }

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      (TransStack[TransStackCurrent+1])[i][j]
	= (TransStack[TransStackCurrent])[i][j];
  
  TransStackCurrent++;

}




/* pops the top matrix */
void Trans_Stack::popmatrix(void) {

  if (TransStackCurrent <= 0) {
    printf("Trans_Stack: error, can't popmatrix, stack is empty.\n");
    return;
  }
  
  TransStackCurrent--;

}




/* replaces the top matrix with a perspective matrix */
void Trans_Stack::perspective(short fovy, REAL aspect, REAL near,
			      REAL far) {
  Matrix m;
  REAL *m1;
  register int i;
  register double angle;

  angle = (double)DEGTORAD(((REAL)fovy / 10.0));

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  m[1][1] = (REAL) (cos(angle/2.0) / sin(angle/2.0));
  m[0][0] = m[1][1] / aspect;
  m[2][2] = -(far+near) / (far-near);
  m[2][3] = -1.0;
  m[3][2] = (REAL) (-(double)(2.0*far*near) / (double)(far-near));
  loadmatrix(m);

}




/* replaces the top matrix with a perspective matrix */
void Trans_Stack::window(REAL left, REAL right, REAL bottom, 
			 REAL top, REAL near, REAL far) {

  Matrix m;
  REAL *m1;
  register int i;

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  m[0][0] = (2.0*near) / (right-left);
  m[1][1] = (2.0*near) / (top-bottom);
  m[2][0] = (right+left) / (right-left);
  m[2][1] = (top+bottom) / (top-bottom);
  m[2][2] = -(far+near) / (far-near);
  m[2][3] = -1.0;
  m[3][2] = -(2.0*far*near) / (far-near);
  loadmatrix(m);

}




/* replaces the top matrix with a 3D orthographic matrix */
void Trans_Stack::ortho(REAL left, REAL right, REAL bottom,
			REAL top, REAL near, REAL far) {

  Matrix m;
  REAL *m1;
  register int i;

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  m[0][0] = 2.0 / (right-left);
  m[1][1] = 2.0 / (top-bottom);
  m[2][2] = -2.0 / (far-near);
  m[3][0] = -(right+left) / (right-left);
  m[3][1] = -(top+bottom) / (top-bottom);
  m[3][2] = -(far+near) / (far-near);
  m[3][3] = 1.0;
  loadmatrix(m);

}




/* replaces the top matrix with a 2D orthographic matrix */
void Trans_Stack::ortho2(REAL left, REAL right, REAL bottom, REAL top) {

  Matrix m;
  REAL *m1;
  register int i;

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  m[0][0] = 2.0 / (right-left);
  m[1][1] = 2.0 / (top-bottom);
  m[2][2] = -1.0;
  m[3][0] = -(right+left) / (right-left);
  m[3][1] = -(top+bottom) / (top-bottom);
  m[3][3] = 1.0;
  loadmatrix(m);

}




/* This subroutine defines a viewing transformation with the eye located in
 * polar coordinates looking at the world origin with twist rotation about
 * the line of site.  Precisely, polarview does:
 * polarview = rot(-azim,z)*rot(-inc,x)*rot(-twist,z)*trans(0.0,0.0,-dist)
 */
void Trans_Stack::polarview(REAL dist, short azim, short inc, short twist) {

  /* premultiply the stack by trans(0.0,0.0,-dist) */
  translate(0.0,0.0,-dist);

  /* premultiply the stack by rotate(-twist,z) */
  rotate(-twist,'z');

  /* premultiply the stack by rotate(-inc,x) */
  rotate(-inc,'x');

  /* premultiply the stack by rotate(-axim,z) */
  rotate(-azim,'z');

}




/* This subroutine defines a viewing transformation with the eye at the point
 * (vx,vy,vz) looking at the point (px,py,pz).  Twist is the right-hand
 * rotation about this line.  The resultant matrix is multiplied with
 * the top of the transformation stack and then replaces it.  Precisely,
 * lookat does:
 * lookat = trans(-vx,-vy,-vz)*rotate(theta,y)*rotate(phi,x)*rotate(-twist,z)
 */
void Trans_Stack::lookat(REAL vx, REAL vy, REAL vz, REAL px, REAL py,
			 REAL pz, short twist) {

  Matrix m;
  REAL *m1;
  register int i;
  register double tmp;

  /* pre multiply stack by rotate(-twist,z) */
  rotate(-twist,'z');
  
  /* pre multiply by rotate (phi,x) */
  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;
  m[0][0] = 1.0;
  tmp = sqrt((double) ((px-vx)*(px-vx) + (py-vy)*(py-vy) + (pz-vz)*(pz-vz)));
  m[1][1] = (REAL)(sqrt((double) ((px-vx)*(px-vx) + (pz-vz)*(pz-vz))) / tmp );
  m[1][2] = (REAL)( (double)(vy-py) / tmp );
  m[2][1] = -m[1][2];
  m[2][2] = m[1][1];
  m[3][3] = 1.0;
  multmatrix(m);

  /* premultiply by rotate(theta,y) */
  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;
  tmp = sqrt((double) ((px-vx)*(px-vx) + (pz-vz)*(pz-vz)));
  m[0][0] = (REAL) ((double)(vz-pz) / tmp );
  m[0][2] = (REAL) -((double)(px-vx) / tmp );
  m[1][1] = 1.0;
  m[2][0] = -m[0][2];
  m[2][2] = m[0][0];
  m[3][3] = 1.0;
  multmatrix(m);

  /* premultiply by trans(-vx,-vy,-vz) */
  translate(-vx,-vy,-vz);

}




/* performs a rotation around an axis */
void Trans_Stack::rotate(short a, char axis) {

  Matrix m;
  REAL *m1;
  register int i;
  register double angle;

  angle = (double)DEGTORAD(((REAL)a / 10.0));

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  if (axis == 'x') {
    m[0][0] = 1.0;
    m[1][1] = (REAL)cos(angle);
    m[1][2] = (REAL)sin(angle);
    m[2][1] = -m[1][2];
    m[2][2] = m[1][1];
    m[3][3] = 1.0;
  } else if (axis == 'y') {
    m[0][0] = (REAL)cos(angle);
    m[0][2] = (REAL) -sin(angle);
    m[1][1] = 1.0;
    m[2][0] = -m[0][2];
    m[2][2] = m[0][0];
    m[3][3] = 1.0;
  } else if (axis == 'z') {
    m[0][0] = (REAL)cos(angle);
    m[0][1] = (REAL)sin(angle);
    m[1][0] = -m[0][1];
    m[1][1] = m[0][0];
    m[2][2] = 1.0;
    m[3][3] = 1.0;
  } else {
    printf("Trans_Stack: error, illegal axis in rotate.\n");
    return;
  }
  
  multmatrix(m);

}




/* performs a rotation around an axis */
void Trans_Stack::rot(REAL a, char axis) {

  Matrix m;
  REAL *m1;
  register int i;
  register double angle;

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  angle = (double)DEGTORAD(a);

  if (axis == 'x') {
    m[0][0] = 1.0;
    m[1][1] = (REAL)cos(angle);
    m[1][2] = (REAL)sin(angle);
    m[2][1] = -m[1][2];
    m[2][2] = m[1][1];
    m[3][3] = 1.0;
  } else if (axis == 'y') {
    m[0][0] = (REAL)cos(angle);
    m[0][2] = (REAL) -sin(angle);
    m[1][1] = 1.0;
    m[2][0] = -m[0][2];
    m[2][2] = m[0][0];
    m[3][3] = 1.0;
  } else if (axis == 'z') {
    m[0][0] = (REAL)cos(angle);
    m[0][1] = (REAL)sin(angle);
    m[1][0] = -m[0][1];
    m[1][1] = m[0][0];
    m[2][2] = 1.0;
    m[3][3] = 1.0;
  } else {
    printf("Trans_Stack: error, illegal axis in rot.\n");
    return;
  }

  multmatrix(m);

}




/* performs a translation */
void Trans_Stack::translate(REAL x, REAL y, REAL z) {

  Matrix m;
  REAL *m1;
  register int i;

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  m[0][0] = 1.0;
  m[1][1] = 1.0;
  m[2][2] = 1.0;
  m[3][0] = x;
  m[3][1] = y;
  m[3][2] = z;
  m[3][3] = 1.0;

  multmatrix(m);

}




/* performs scaling */
void Trans_Stack::scale(REAL x, REAL y, REAL z) {

  Matrix m;
  REAL *m1;
  register int i;

  m1 = &m[0][0];
  for (i=0;i<16;i++)
    *(m1+i) = 0.0;

  m[0][0] = x;
  m[1][1] = y;
  m[2][2] = z;
  m[3][3] = 1.0;

  multmatrix(m);

}




/* prints the matrix */
void Trans_Stack::printit() {

  register int i,j;

  printf("\nTrans_Stack matrix = \n");
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      printf("%f ",(TransStack[TransStackCurrent])[i][j]);
    printf("\n");
  }
  printf("\n");

}
