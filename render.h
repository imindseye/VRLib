/////////////////////////////////////////////////////////////////////
//
//                       Volume Renderer Class
//                            Han-Wei Shen 
//                        hwshen@nas.nasa.gov 
//
//         MRJ Technology Solutions/NASA Ames Research Center 
//
#ifndef RENDER_H 
#define RENDER_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <iostream>

#include "image.h"
#include "Map.h"
#include "Trans_Stack.h"
#include "minmax.h"

#define EPS 1.0E-6

#define TRUE 1
#define FALSE 0

//#define Ka 0.5
//#define Kd 0.5

#define Ka 0.6
#define Kd 0.4

#define Ks 0.1
#define Alpha 2.0

//#define ambient_light 0.3
//#define light_strength 0.9
#define ambient_light  1
#define light_strength 5

static REAL eye_W[4] =  {0,0,1,0};	/* eye vector */
static REAL light_W[4] = {0.40824829,0.40824829,0.816496,0};  

void Image_to_File(image_type*, char*); 

///////////////////////////////////////////////////////
//
//  Used by those vectors such as gradient, light, etc. 
//
typedef struct _uvw
{
  REAL u,v,w; 
} uvw; 

///////////////////////////////////////////////////////
//
// This structure saves some state shared 
// between get_value and get_normal. This 
// makes get_normal faster since it can reuse 
// some of the work done by get_value.
//
struct interpolation_state
{
  unsigned long offsets[8]; 
  REAL tx,ty,tz; 
}; 

union VolumePtr {
  REAL* fVolume; 
}; 

//////////////////////////////////////////////////////
class volumeRender {
public: 

  enum VolumeType {
    RAW    = 1
  }; 

protected:

  VolumePtr vptr; 

  void get_bounds(); 
  void update_bounds(REAL q[4]); 

  uvw  *gradient;           // volume gradient
  int has_gradient; 
				 
  Map  *map;                // color map

  float* lookup; 
  int lookupSize; 

  void set_color_map(char* mapname);  // old color map file
  void set_color_map(Map*);           // old color map file

  // viewing bounding box  (data space)
  int vxmin, vxmax, vymin, vymax, vzmin, vzmax; 				 
  int vxdim, vydim, vzdim; 

  // actual in core data bounding box (data space)
  int lxmin, lxmax, lymin, lymax, lzmin, lzmax; 				 
  int lxdim, lydim, lzdim, lxdimlydim; 

  // rendering range, i.e.  bounding box (data space)
  int rxmin, rxmax, rymin, rymax, rzmin, rzmax; 				 
  int rxdim, rydim, rzdim; 

  int udim, vdim;           // image dimensions

  // image bound 
  int umin, umax, vmin, vmax, wmin, wmax; 

  REAL xangle, yangle, zangle; 

  uvw eye, light, h, eye_nonN;   // eye, light, and half vector
                            // in the data coordinate system

  Matrix screen_to_data;    // mapping between different 
  Matrix data_to_world;     // coordinate systems
  Matrix world_to_data; 

  int UNIFORM_FLAG; 
  REAL UNIFORM_VAL; 

  float curMin, curMax; 

  int get_value(REAL p[4], REAL*,
                interpolation_state*); 

  int get_normal(uvw*, interpolation_state*); 

  void local_lighting(REAL*, interpolation_state*, 
                      REAL obj_color[4], REAL result[3]); 

  void local_lighting(REAL*, 
                      REAL obj_color[4], REAL result[3]); 

  void depth_lighting(REAL*, interpolation_state*, 
                      REAL obj_color[4], REAL result[3]); 

  void compute_gradient(REAL*, int); // compute gradients for 
                                     // a slice of the data



  int mapLookup(float, float*); 

  int check_inbound(REAL[4]); 				     

  void update_transform(); 
  void update_viewing(); 

  void render();  // regular volume rendering 





  // set the range of data that will be actually accessible 
  // (i.e. in core) 
  //
  void set_data_and_bbx(int imin, int imax, int jmin, 
		    int jmax, int kmin, int kmax, 
		    void* data, 
                    int computeGrd = 1, uvw* grad = NULL); 

  // a simple setup routine that will call set_vewing_bbx(), 
  // set_data_and_bbx(), and set_clipping_bbx() using 
  // default parameters 
  //
  void set_volume_simple(int imin, int imax, int jmin, 
		  int jmax, int kmin, int kmax, 
		  void* volume); 



  void get_opacity(float, float*, 
		   interpolation_state*); 

public:

  volumeRender(int xdim, int ydim, int zdim, 
	       int udim, int vdim, 
	       void *volume); 

  ~volumeRender(); 

  void execute(int is_uniform = 0, 
	       REAL mean = 0.0 ); 

  int  readCmapFile(char *filename);

  Matrix data_to_screen;    // transformation matrixes for 

  void setColorMap(int size, float* ctable); 

  // rotation angle: alpha, beta, gamma 
  void set_view(float xA, float yA, float zA);

  void get_viewing_bbx(int& imin,int& imax,int& jmin,
		       int& jmax,int& kmin,int& kmax) {
       imin = vxmin; imax = vxmax;  jmin = vymin; jmax= vymax; 
       kmin = vzmin; kmax = vzmax; 				 
  }


  void set_clipping_bbx(int imin, int imax, int jmin, 
		      int jmax, int kmin, int kmax); 


  void out_to_image(char* fname);  // to be implemented 			
  //-----------------------------------------
  // following are just some query functions for the 
  // rendering setup -- not essential for performing 
  // the actual rendering 
  // 
  void get_view_vector(REAL&, REAL&, REAL&);  // in data space  
  void get_eye_vector(REAL&, REAL&, REAL&);  // in data space  

  float distance_to_viewplane(float, float, float); 

  void set_min_max(float min, float max) {
    curMin = min; curMax = max; }

  void get_min_max(float& min, float& max) {
    min = curMin; max = curMax; }

  void set_image_size(int usize, int vsize); 

  void get_image_dims(int& usize, int& vsize) {
                      usize = udim; vsize = vdim; }

  // set the viewing bounding box in data space 
  //
  void set_viewing_bbx(int imin, int imax, int jmin, 
		       int jmax, int kmin, int kmax); 

  uvw* get_gradient() {return gradient;}

  void computeVolMinMax(float*, int, float &vol_min, float &vol_max);

  image_type *image;               // output image

}; 

#endif
