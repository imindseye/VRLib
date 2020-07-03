/*
 * Trans_Stack.h - This is the Trans_Stack header file.  It is used by
 * other classes to provide transformations.
 *
 * Written by Mike Krogh, NCSA, 3/9/92
 *
 */

#ifndef TRANSSTACK_H
#define TRANSSTACK_H

#ifndef PI
#define PI 3.141592654
#endif

// Degree-to-Radians and Radians-to-Degrees Conversion macros 
#define DEGTORAD(a)  (a*PI/180.0)
#define RADTODEG(a)  (a*180.0/PI)

/* the 'thing' stored in the stack */
#ifndef __GL_GL_H__
typedef REAL Matrix[4][4];
#define XMAXSCREEN 1279
#define YMAXSCREEN 1023
#endif

#define TransStackMax 32

class Trans_Stack {

private:
  Matrix TransStack[TransStackMax];
  int TransStackCurrent;

public:

  Trans_Stack(void);

  ~Trans_Stack(void);

  /* returns the top matrix */
  void getmatrix(Matrix);

  /* replaces the top matrix with m*top */
  void multmatrix(Matrix);

  /* multiplies a 3D point by the Matrix */
  void multpoint3d(REAL[3], REAL[3]);

  /* multiplies a 4D point by the Matrix */
  void multpoint4d(REAL[4], REAL[4]);

  /* replaces the top matrix */
  void loadmatrix(Matrix);

  /* pushs the top matrix */
  void pushmatrix(void);

  /* pops the top matrix */
  void popmatrix(void);

  /* replaces the top matrix with a perspective matrix */
  void perspective(short, REAL, REAL, REAL);

  /* replaces the top matrix with a perspective matrix */
  void window(REAL, REAL, REAL, REAL, REAL, REAL);

  /* replaces the top matrix with a 3D orthographic matrix */
  void ortho(REAL, REAL, REAL, REAL, REAL, REAL);

  /* replaces the top matrix with a 2D orthographic matrix */
  void ortho2(REAL, REAL, REAL, REAL);

  /* This subroutine defines a viewing transformation with the eye located in
   * polar coordinates looking at the world origin with twist rotation about
   * the line of site.  Precisely, polarview does:
   * polarview = rot(-azim,z)*rot(-inc,x)*rot(-twist,z)*trans(0.0,0.0,-dist)
   */
  void polarview(REAL, short, short, short);

  /* This subroutine defines a viewing transformation with the eye at the point
   * (vx,vy,vz) looking at the point (px,py,pz).  Twist is the right-hand
   * rotation about this line.  The resultant matrix is multiplied with
   * the top of the transformation stack and then replaces it.  Precisely,
   * lookat does:
   * lookat = trans(-vx,-vy,-vz)*rotate(theta,y)*rotate(phi,x)*rotate(-twist,z)
   */
  void lookat(REAL, REAL, REAL, REAL, REAL, REAL, short);

  /* performs a rotation around an axis */
  void rotate(short, char);

  /* performs a rotation around an axis */
  void rot(REAL, char);

  /* performs a translation */
  void translate(REAL, REAL, REAL);

  /* performs scaling */
  void scale(REAL, REAL, REAL);

  /* prints the matrix */
  void printit();

};

#endif

