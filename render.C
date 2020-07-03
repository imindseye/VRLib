////////////////////////////////////////////////////////////////
//
//             Definitions of the volume renderer class
//                           Han-Wei Shen 
//                        hwshen@nas.nasa.gov 
//
//         MRJ Technology Solutions/NASA Ames Research Center 
//
//     ** Note: This code evolves from the C programs that I used 
//              in 1993 at Los Alamos National Laboratory. The very 
//              original author I believe is Jamie Painter. Since then, 
//              various changes have been made, but nothing dramatic. 
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "render.h"
#include "image.h"
#include "Map.h"
#include "Trans_Stack.h"
#include "minmax.h"
#include "render_aux.h"

#define  STEPSIZE  0.9

/////////////////////////////////////////////////////////////
//
//                      Constructor
//            The user provides image dimensions 
//
volumeRender::volumeRender(int xsize, int ysize, int zsize, 
			   int usize, int vsize, 
			   void* volume):
  udim(usize),  vdim(vsize),  xangle(0), 
  yangle(0), zangle(0), gradient(NULL), image(NULL)
{
  set_volume_simple(0,xsize-1,0,ysize-1,0,zsize-1, 
		    volume); 
  //  image = image_new(0,udim-1,0,vdim-1);
}

/////////////////////////////////////////////////////////////
//
//                       Destructor
//
volumeRender::~volumeRender()
{
  //  ** note: Rather not to free the data. Pass the ball back 
  //  to the user. 
  //  if (data!=NULL) free(data); 

  if (gradient!=NULL) delete[]gradient; 
  if (image!=NULL) free(image); 
}

/////////////////////////////////////////////////////////////
//
//             Specify a Color Map : either pass a 
//             color map that has been set up or 
//             provide the color map file name. 
//
void volumeRender::set_color_map(Map* newmap)
{
  map = newmap; 
}
void volumeRender::set_color_map(char* mapname)
{
  map = new Map; 
  map->read_file(mapname); 
}

/////////////////////////////////////////////////////////////
//
//  Compute gradients, only needs to be done if 
//  lighting is calculated 
//
void volumeRender::compute_gradient(REAL* data, int z)
{
  int x, y;
  int idx; 
  uvw *normal;

  if (z < lzmin || z > lzmax)
    return;  // wrong z value

    for (y=lymin; y<=lymax; y++)
      for (x=lxmin; x<=lxmax; x++) {

      idx = (x-lxmin) + (y-lymin)*lxdim + (z-lzmin)*lxdimlydim;
      normal = gradient + idx; 

      if (x==lxmin)
	normal->u= (data[idx+1] - data[idx]);
      else if (x== lxmax)
	normal->u= (data[idx] - data[idx-1]);
      else
	normal->u= (data[idx+1] - data[idx-1]) / 2.0;
	  
      if (y==lymin)
	normal->v= (data[idx+lxdim] - data[idx]);
      else if (y==lymax)
	normal->v= (data[idx] - data[idx-lxdim]);
      else
	normal->v= (data[idx+lxdim] - data[idx-lxdim]) / 2.0;

      if (z==lzmin)
	normal->w= (data[idx+lxdimlydim] - data[idx]);
      else if (z== lzmax)
	normal->w= (data[idx] - data[idx-lxdimlydim]);
      else
	normal->w= (data[idx+lxdimlydim] - data[idx-lxdimlydim]) / 2.0;

      // Normalize the gradient, if necessary 
      Normalize(normal);
      }
}

////////////////////////////////////////////////////
//
//  Check if the given point is in the volume 
//
int volumeRender::check_inbound(REAL p[4])
{
  int x1, y1, z1;
  REAL x,y,z;

  // Make the coordinates relative to our subvolume 
  x = p[0];
  y = p[1];
  z = p[2];

  // Integerize the coordinates 
  x1 = (int)floor((double)x);
  y1 = (int)floor((double)y);
  z1 = (int)floor((double)z);

  if (x1 < rxmin || x1 >=rxmax || 
      y1 < rymin || y1 >=rymax ||
      z1 < rzmin || z1 >=rzmax) {

    return FALSE;
  }
  else if (x1 < lxmin || x1 > lxmax || 
      y1 < lymin || y1 > lymax ||
      z1 < lzmin || z1 > lzmax) {

    return FALSE;
  }
  return TRUE; 
}


//////////////////////////////////////////////////////////////////
//
//  get opacity of a value: implement Levoy's opacity ramp 
//  (code written by xinyue li) 
//
void volumeRender::get_opacity(float val, float *alpha, interpolation_state *is)
{
  float iso_opacity,gradient_magnitude,temp;
  uvw normal;
  float iso_value = 1.0, thickness = 2.0;

  iso_opacity = 1.0; 

  get_normal(&normal,is);
  gradient_magnitude = sqrt(normal.u*normal.u+
			    normal.v*normal.v+normal.w*normal.w);
  if((gradient_magnitude==0)&&(val==iso_value)) 
    *alpha = iso_opacity;
  else if((gradient_magnitude>0)&&
	  ((val-thickness*gradient_magnitude)<=iso_value)&&
	  ((val+thickness*gradient_magnitude)>=iso_value)){
    temp = fabs((iso_value-val)/(gradient_magnitude*thickness));
    temp = 1 - temp;
    *alpha = iso_opacity * temp;
  }
  else *alpha = 0.0;
}


////////////////////////////////////////////////////
//
//  Given a point in the volume, 
//  interpolate the scalar value 
//  on that point
//
int volumeRender::get_value(REAL p[4] ,REAL *val, 
              interpolation_state *is) 
{
  int x1, y1, z1;
  int z1offset, y1offset;
  REAL x,y,z;

  // temps for trilinear interpolation 
  REAL top, bot, front, back;

  *val = 0.0;

  // Make the coordinates relative to our subvolume 
  x = p[0];
  y = p[1];
  z = p[2];

  // Integerize the coordinates 
  x1 = (int)floor((double)x);
  y1 = (int)floor((double)y);
  z1 = (int)floor((double)z);

  if (x1 < rxmin || x1 >= rxmax ||   
      y1 < rymin || y1 >= rymax ||
      z1 < rzmin || z1 >= rzmax) {

    return FALSE;   // no your business
  }
  else if (x1 < lxmin || x1 >= lxmax || 
      y1 < lymin || y1 >= lymax ||
      z1 < lzmin || z1 >= lzmax) {

    return FALSE;   // no data 
  }

  //
  // These are the 't' values used for the x, y 
  // and z interpolation steps ,i.e. the fraction 
  // in x, y, z we are into the voxel
  //

  is->tx = x - x1;
  is->ty = y - y1;
  is->tz = z - z1;

  // convert coordinates to be relative to this volume 
  x1 -= lxmin;
  y1 -= lymin;
  z1 -= lzmin;

  // Compute offsets to the eight sournding voxels

  z1offset = lxdimlydim*z1;
  y1offset = y1*lxdim;

  is->offsets[0] = z1offset + y1offset + x1;   /* [x1,y1,z1] */
  is->offsets[1] = is->offsets[0] + 1;         /* [x2,y1,z1] */
  is->offsets[2] = is->offsets[1] + lxdim;     /* [x2,y2,z1] */
  is->offsets[3] = is->offsets[0] + lxdim;     /* [x1,y2,z1] */

  is->offsets[4] = is->offsets[0] + lxdimlydim; /* [x1,y1,z2] */
  is->offsets[5] = is->offsets[1] + lxdimlydim; /* [x2,y1,z2] */
  is->offsets[6] = is->offsets[2] + lxdimlydim; /* [x2,y2,z2] */
  is->offsets[7] = is->offsets[3] + lxdimlydim; /* [x1,y2,z2] */

  REAL* data = vptr.fVolume; 
  // Interpolate in the z=z1 plane first 
  bot   = lerp(is->tx,data[is->offsets[0]],data[is->offsets[1]]);
  top   = lerp(is->tx,data[is->offsets[3]],data[is->offsets[2]]);
  front = lerp(is->ty,bot,top);

  // now in the z=z2 plane 
  bot  = lerp(is->tx,data[is->offsets[4]],data[is->offsets[5]]);
  top  = lerp(is->tx,data[is->offsets[7]],data[is->offsets[6]]);
  back = lerp(is->ty,bot,top);

  // finally, interpolate between the two z planes
  *val = lerp(is->tz,front,back);

  return TRUE;
}

//////////////////////////////////////////////////////////////
//
// Trilinearly interpolate the normal vector
// Assumes the interpolation state is already setup in 'is'.
//
int volumeRender::get_normal(uvw *val, interpolation_state *is) 
{

  assert(gradient!=NULL); 

// temps for trilinear interpolation 
  REAL top, bot, front, back;

// Interpolate the u,v and w fields of the gradient and put the result in val
#define LERP_1_FIELD(FIELD) \
  bot = lerp(is->tx,gradient[is->offsets[0]].FIELD,  \
	     gradient[is->offsets[1]].FIELD);        \
  top = lerp(is->tx,gradient[is->offsets[3]].FIELD,  \
	     gradient[is->offsets[2]].FIELD);        \
  front = lerp(is->ty,bot,top);                      \
  bot = lerp(is->tx,gradient[is->offsets[4]].FIELD,  \
	     gradient[is->offsets[5]].FIELD);        \
  top = lerp(is->tx,gradient[is->offsets[7]].FIELD,  \
	     gradient[is->offsets[6]].FIELD);        \
  back = lerp(is->ty,bot,top);                       \
  val->FIELD = lerp(is->tz,front,back); 

  LERP_1_FIELD(u);
  LERP_1_FIELD(v);
  LERP_1_FIELD(w);

// Normalize the normal vector, if necessary
  Normalize(val);
  return TRUE;
}

//////////////////////////////////////////////////////
//
// Evaluate a simple lighting model at a point.
// This is, roughly, equation 16.13 from Foley, vanDam, Feiner & 
// Hughes, 2nd edition.
//
// Their is no distance attenuation since the light is considered to be
// at infinity.
//
void volumeRender::local_lighting( REAL *point, interpolation_state *is,
		 REAL obj_color[4], REAL result[3] )
{
  uvw normal;
  REAL sign = 1.0;
  REAL NdotV, NdotL, NdotH;
  REAL ambient, diffuse, specular;

  get_normal(&normal,is);

  // Logically flip the normal if it is backfacing relative to the eye vector.
  NdotV = Dot(&normal,&eye);
  if (NdotV < 0) sign = -1;	  // flip normal

  // Ambient term
  ambient = Ka * ambient_light;

  // Diffuse term
  NdotL = sign*Dot(&normal,&light);
  if (NdotL < 0) 
    diffuse = 0;
  else
    diffuse = light_strength * Kd * NdotL;
  
  // Specular highlight
  NdotH = sign*Dot(&normal,&h);
  if (NdotH < 0) 
    specular = 0;
  else
    specular = light_strength * Ks * ipow( NdotH, 30) ;

  result[0] = clamp((ambient + diffuse)*obj_color[0] + specular,0.0,1.0);
  result[1] = clamp((ambient + diffuse)*obj_color[1] + specular,0.0,1.0);
  result[2] = clamp((ambient + diffuse)*obj_color[2] + specular,0.0,1.0);
}

///////////////////////////////////////////////////////////////
//
//               User pass in the normal 
//
void volumeRender::local_lighting( REAL *n, 
		 REAL obj_color[4], REAL result[3] )
{
  uvw normal;
  REAL sign = 1.0;
  REAL NdotV, NdotL, NdotH;
  REAL ambient, diffuse, specular;

  normal.u = n[0];   normal.v = n[1];   normal.u = n[2]; 

  // Logically flip the normal if it is backfacing relative to the eye vector.
  NdotV = Dot(&normal,&eye);
  if (NdotV < 0) sign = -1;	  // flip normal

  // Ambient term
  ambient = Ka * ambient_light;

  // Diffuse term
  NdotL = sign*Dot(&normal,&light);
  if (NdotL < 0) 
    diffuse = 0;
  else
    diffuse = light_strength * Kd * NdotL;
  
  // Specular highlight
  NdotH = sign*Dot(&normal,&h);
  if (NdotH < 0) 
    specular = 0;
  else
    specular = light_strength * Ks * ipow( NdotH, 30) ;

  result[0] = clamp((ambient + diffuse)*obj_color[0] + specular,0.0,1.0);
  result[1] = clamp((ambient + diffuse)*obj_color[1] + specular,0.0,1.0);
  result[2] = clamp((ambient + diffuse)*obj_color[2] + specular,0.0,1.0);
}


///////////////////////////////////////////////////////////////////
//
// This is some sort of depth attenuation trick that 
// Mike came up with as a substitute for a real lighting 
// model.
//
void volumeRender::depth_lighting( REAL *point, interpolation_state *is,
		 REAL obj_color[4], REAL result[3])
		  
{
  REAL zn = (point[2]-vzmin)/(vzdim-1);
  REAL weight = exp(-Alpha*zn);

  result[0] = weight*obj_color[0];
  result[1] = weight*obj_color[1];
  result[2] = weight*obj_color[2];
}

/////////////////////////////////////////////////////////////////
//
// transform each corner from world coordinates to 
// image coordinates and keep the min and max values.
//
void volumeRender::get_bounds() 
{
  REAL p[4],q[4];
  extern  void matrix_mult(Matrix,REAL*,REAL*);

  p[0] = (REAL)rxmin; p[1] = (REAL)rymin; p[2] = (REAL)rzmin; p[3] = 1.0;

  matrix_mult(data_to_screen,p,q);
  umin = (int)rint((double)q[0]);
  umax = umin-1;
  vmin = (int)rint((double)q[1]);
  vmax = vmin-1;
  wmin = (int)rint((double)q[2]);
  wmax = wmin-1;

  p[0] = (REAL)rxmax;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  p[1] = (REAL)rymax;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  p[0] = (REAL)rxmin;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  p[0] = (REAL)rxmin;
  p[1] = (REAL)rymin;
  p[2] = (REAL)rzmax;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  p[0] = (REAL)rxmax;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  p[1] = (REAL)rymax;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  p[0] = (REAL)rxmin;
  matrix_mult(data_to_screen,p,q);
  update_bounds(q);

  if (umin < 0)     umin = 0;
  if (umin >= udim) umin = udim-1;
  if (umax < 0)     umax = 0;
  if (umax >= udim) umax = udim-1;

  if (vmin < 0)     vmin = 0;
  if (vmin >= vdim) vmin = vdim-1;
  if (vmax < 0)     vmax = 0;
  if (vmax >= vdim) vmax = vdim-1;
}  

/////////////////////////////////////////////////////
// Update the bounds of the volume in screen space given one of the corner
// verticies.
//
// Note that the integer pixel lattice passes through the the screen coordinate
// system at the 0.5,0.5 marks in continous screen coordinates.  That is a 
// pixel spans between integer coordinates with the centers at the 0.5,0.5
// marks.
//
// Some care is taken here to ensure that ray samples are never counted twice.
//
void volumeRender::update_bounds(REAL q[4]) 
{
  int u = (int)rint((double)q[0]);
  int v = (int)rint((double)q[1]);
  int w = (int)rint((double)q[2]);

  if (u < umin) umin = u;
  if (u > umax) umax = u-1;
  if (v < vmin) vmin = v;
  if (v > vmax) vmax = v-1;
  if (w < wmin) wmin = w;
  if (w > wmax) wmax = w-1;
}

/////////////////////////////////////////////////////
//
//  Output to an ASCII ppm file 
//
void volumeRender::out_to_image(char* filename)
{
  unsigned char* red   = new unsigned char[udim*vdim]; 
  unsigned char* green = new unsigned char[udim*vdim]; 
  unsigned char* blue  = new unsigned char[udim*vdim]; 

  for (int i=0; i<udim; i++)
    for (int j=0; j<vdim; j++) {
      int idx = j*udim+i; 
      pixel*p = image_index(image, i, j); 
      red[idx]   = (unsigned char) (p->bp.r); 
      green[idx] = (unsigned char) (p->bp.g); 
      blue[idx]  = (unsigned char) (p->bp.b); 
    }

  FILE* out = fopen(filename, "w"); 
  fprintf(out, "P3\n"); 
  fprintf(out, "%d  %d\n", udim, vdim); 
  fprintf(out, "%d\n", 256); 

  for (int i=0; i<udim; i++)
    for (int j=vdim-1; j>=0; j--) {
      int idx = j*udim+i; 
      fprintf(out, "%d %d %d\n", red[idx], 
	      green[idx], blue[idx]); 
    }
  fclose(out); 

  delete []red;   delete []green;   delete []blue; 
}

///////////////////////////////////////////////////////////////////
//
// Here is where the rendering work is actually done.
//
void volumeRender::render() 
{
  extern void matrix_mult(Matrix,REAL in[],REAL out[]); 
  REAL sumred, sumgreen, sumblue, sumalpha;
  pixel *p;
  REAL p1[4],p2[4],inc[4];
  REAL rgba[4];
  REAL val1;
  REAL outcolor[3];
  REAL alpha;
  interpolation_state is;
  int use_uniform = 0; 
  int step_count = 0; 
  float* alphalut; 

  // compute the bounding volume 
  get_bounds();

  // reset the image 

  if (image!=NULL) free(image); 
  image = image_new(umin, umax, vmin, vmax); 

  zero_rect(image, umin, umax, vmin, vmax); 

  p1[3] = 1.0;        

  if (UNIFORM_FLAG) {
    //       use_uniform = map->lookup(UNIFORM_VAL, rgba); 
    use_uniform = mapLookup(UNIFORM_VAL, rgba); 
    alphalut = new float[wmax-wmin+1]; 
    sumalpha = 0.0; 
    for (int i=0; i<wmax-wmin+1; i++) {
      alpha = rgba[3]*(1.0 - sumalpha);         // add in to result
      sumalpha += alpha;
      alphalut[i] = sumalpha; 
    }
  }

  for (int u=umin; u<=umax; u++) {     // loop over each pixel 
    p1[0] = (REAL)u + 0.5;             //pixels centers are at the 0.5 mark 
    for (int v=vmin; v<=vmax; v++) {
      p1[1] = (REAL)v + 0.5;
      sumred = sumgreen = sumblue = sumalpha = 0.0;

      // Determine basepoint of ray 
      p1[2] = (REAL) wmin + 0.5;
      matrix_mult(screen_to_data,p1,p2);
      
      // Determine increment between sample steps
      p1[2] += 1.0;
      matrix_mult(screen_to_data, p1,inc);
      inc[0] -= p2[0];  inc[1] -= p2[1];   inc[2] -= p2[2];
      step_count = 0; 

      for (int z=wmin; z<=wmax; z+=1, p2[0] += inc[0], 
	                        p2[1] +=inc[1], p2[2] += inc[2]) {
	if (use_uniform) {
	  //	  if (check_inbound(p2) == FALSE) 
	    if (p2[0] < rxmin || p2[0] >=rxmax || 
		p2[1] < rymin || p2[1] >=rymax ||
		p2[2] < rzmin || p2[2] >=rzmax) 
	      continue; 
	  else 
	    step_count++; 
	  if (alphalut[step_count-1]>=0.99) break; 
	}
	else if (get_value(p2,&val1,&is)) {// get the data value
	  //get_opacity(val1,&alpha,&is);	  
	  // if (map->lookup(val1,rgba)) {// lookup corresponding RGBA
	  if (mapLookup(val1,rgba)) {// lookup corresponding RGBA
	    if (rgba[3] > EPS) {                        // partly opaque?
	      if (has_gradient) {
		local_lighting(p2,&is,rgba,outcolor);     // compute lighting
	      }
	      else {
		//depth_lighting(p2,&is,rgba,outcolor);  // compute lighting
		outcolor[0] = rgba[0]; 
		outcolor[1] = rgba[1]; 
		outcolor[2] = rgba[2]; 
	      }
	      alpha = rgba[3]*(1.0 - sumalpha);	  // add in to result
	      //alpha = alpha*(1.0 - sumalpha);
	      sumred   += (outcolor[0]*alpha);
	      sumgreen += (outcolor[1]*alpha);
	      sumblue  += (outcolor[2]*alpha);
	      sumalpha += alpha;
	    }
	  }
	if (sumalpha >= 0.99)  break;
	}
      }
      if (use_uniform&& step_count!=0) {
	sumred   = (rgba[0]*alphalut[step_count-1]);
	sumgreen   = (rgba[1]*alphalut[step_count-1]);
	sumblue   = (rgba[2]*alphalut[step_count-1]);
	sumalpha = alphalut[step_count-1]; 
      }
      p = image_index(image,u,v);
      p->bp.r = (unsigned char)clamp(rint((double)(sumred*255.0)),0,255);
      p->bp.g = (unsigned char)clamp(rint((double)(sumgreen*255.0)),0,255);
      p->bp.b = (unsigned char)clamp(rint((double)(sumblue*255.0)),0,255);
      p->bp.a = (unsigned char)clamp(rint((double)(sumalpha*255.0)),0,255);
    }
  }
  if (use_uniform) delete[]alphalut; 
}

/////////////////////////////////////////////////////////////////
//
//   The volume rendering main routine
//   After the execution, image can be obtained 
//   from the 'image' member field.
//
void volumeRender::execute(int is_uniform, REAL val)
{
  UNIFORM_FLAG = is_uniform; 

  if (is_uniform) 
    UNIFORM_VAL = val; 

  render(); 
}
///////////////////////////////////////////////////////////////////
//
//   Specify the in-core data and its bounding box
//
void volumeRender::set_clipping_bbx(int imin, int imax, int jmin, 
				    int jmax, int kmin, int kmax)
{
  if (imin < lxmin) rxmin = lxmin; // data is not available.....
  else rxmin = imin; 
  if (imax > lxmax) rxmax = lxmax; 
  else rxmax = imax; 

  if (jmin < lymin) rymin = lymin; 
  else rymin = jmin; 
  if (jmax > lymax) rymax = lymax; 
  else rymax = jmax; 

  if (kmin <lzmin) rzmin = lzmin; 
  else rzmin = kmin; 
  if (kmax >lzmax) rzmax = lzmax; 
  else rzmax = kmax; 

  rxdim = rxmax-rxmin+1; 
  rydim = rymax-rymin+1; 
  rzdim = rzmax-rzmin+1; 
}
////////////////////////////////////////////////////////////////////
//
//   Specify the in-core data and its bounding box
//
void volumeRender::set_data_and_bbx(int imin, int imax, int jmin, 
				    int jmax, int kmin, int kmax,
				    void* data, 
				    int computeGradient, uvw* grad)
{
  lxmin = imin; lxmax = imax; 
  lymin = jmin; lymax = jmax; 
  lzmin = kmin; lzmax = kmax; 
  lxdim = lxmax-lxmin+1; 
  lydim = lymax-lymin+1; 
  lzdim = lzmax-lzmin+1; 
  lxdimlydim = lxdim*lydim; 

  has_gradient = 0; 

  vptr.fVolume = (REAL*)data; 

  if (grad !=NULL) {
    gradient = grad; 
    has_gradient = 1; 
  }
  else if (computeGradient) { 
    int size;
    int z;

    // allocate space and compute gradient 
    //    size = lxdim * lydim * lzdim * sizeof(uvw);
    size = lxdim * lydim * lzdim; 
    printf(" allocating %d uvws for gradient field\n", size); 

    if (gradient !=NULL) delete[]gradient; 

    gradient = new uvw[size]; 
    if (gradient == NULL) 
    //    if ((gradient = (uvw *) malloc(size)) == NULL)
      printf("error allocating gradient memory in node\n");

    for (z=lzmin; z<=lzmax; z++) {
      //           printf(" compute gradient for z = %d\n", z); 
      //            fflush(stdout); 
      compute_gradient(vptr.fVolume, z); 
    }
    has_gradient = 1; 
  }
}
////////////////////////////////////////////////////////////////////
//
//  Specify a viewing bounding volume. This can be used to 
//  zoom in the volume data. 
//
void volumeRender::set_viewing_bbx(int imin, int imax, int jmin, 
				   int jmax, int kmin, int kmax)
{
  assert(imax >=imin && jmax>=jmin && kmax>=kmin); 
  vxmin = imin; vxmax = imax; 
  vymin = jmin; vymax = jmax; 
  vzmin = kmin; vzmax = kmax; 
  vxdim = vxmax-vxmin+1; 
  vydim = vymax-vymin+1; 
  vzdim = vzmax-vzmin+1; 
}
///////////////////////////////////////////////////////////////////
//
//   A simple function to set all the bounding boxes.
//   This function is used when the user just want to render
//   the whole volume using standard transformation.
//
void volumeRender::set_volume_simple(int imin, int imax, int jmin, 
				     int jmax, int kmin, int kmax, 
				     void* data)
{
  assert(imax >=imin && jmax>=jmin && kmax>=kmin); 
  set_viewing_bbx(imin, imax, jmin, jmax, kmin, kmax); 

  set_data_and_bbx(imin, imax, jmin, jmax, kmin, kmax, data); 

  set_clipping_bbx(imin, imax, jmin, jmax, kmin, kmax); 

}
////////////////////////////////////////////////////////////////////
//
//             Form the transformation matrices 
//

////////////////////////////////////////////////////////////////////
//
//  Updating the tranformation matrixes according to the 
//  current viewing angle and bounding box
//
//  Calculate the matrix that goes between data space to 
//  world space. This matrix accounts for the rotation 
//  of the data and scaling of the data.
//
void volumeRender::update_transform()
{
  Trans_Stack s;
  extern  int vrlib_invert_matrix(Matrix,Matrix);
  int dmax;
  REAL  sc, zscale;

  REAL vxcenter = (vxmin+vxmax)/2.0; 
  REAL vycenter = (vymin+vymax)/2.0; 
  REAL vzcenter = (vzmin+vzmax)/2.0; 

  s.translate(-vxcenter, -vycenter, -vzcenter); 
  dmax = MAX(vxdim,MAX(vydim,vzdim));
  sc = 0.8 * 2.0/(REAL)dmax;
  s.scale( sc, sc, sc );

  while(xangle < 0)   xangle += 360.0;
  while(xangle > 360) xangle -= 360.0;
  while(yangle < 0)   yangle += 360.0;
  while(yangle > 360) yangle -= 360.0;
  while(zangle < 0)   zangle += 360.0;
  while(zangle > 360) zangle -= 360.0;

  s.rotate( (short) (xangle*10 + 0.5), 'x');
  s.rotate( (short) (yangle*10 + 0.5), 'y');
  s.rotate( (short) (zangle*10 + 0.5), 'z');
  s.getmatrix(data_to_world);
  vrlib_invert_matrix(data_to_world,world_to_data);

// Now add in the additional parts of the transform that
// thansform to screen space.  The Z axis is scaled so that
// the sample spacing is +1.0 in screen space.  
//
  s.scale(1.0,-1.0,-1.0);  // flip over the image, and reverse Z
  s.translate(1.0,1.0,1.0);
  zscale = 1.0/(sc * STEPSIZE);
  s.scale((REAL)udim/2.0,(REAL)vdim/2.0,zscale);
  
  s.getmatrix(data_to_screen);
  vrlib_invert_matrix(data_to_screen,screen_to_data);
}
////////////////////////////////////////////////////////////////////
//
// Recompute viewing information for the next frame
//
void volumeRender::update_viewing() 
{
  Trans_Stack s; 
  REAL dlight[4]; 
  REAL deye[4]; 
  
  // Transform the light vector from world space to data space 
  s.loadmatrix(world_to_data);
  s.multpoint4d(light_W,dlight);

  // Transform the eye vector from world space to data space 
  s.multpoint4d(eye_W,deye);

  light.u = dlight[0];
  light.v = dlight[1];
  light.w = dlight[2];
  Normalize(&light);

  eye_nonN.u = eye.u = deye[0];
  eye_nonN.v = eye.v = deye[1];
  eye_nonN.w = eye.w = deye[2];
  Normalize(&eye);

  // Compute halfway vector 
  h.u = eye.u + light.u;
  h.v = eye.v + light.v;
  h.w = eye.w + light.w;
  Normalize(&h);
}

//////////////////////////////////////////////////////////////////////
//
//                  set the rendering view 
//
void volumeRender::set_view(float xA, float yA, float zA)
{
  xangle = xA; 
  yangle = yA; 
  zangle = zA; 
  update_transform();  // create new transformation matrixes
  update_viewing();    // create new eye, light, and half vectors
}
//////////////////////////////////////////////////////////////////////
//
//            get the viewing vector in data space 
//
void volumeRender::get_view_vector(float& vx, float& vy, float& vz)
{
  REAL vector[4]; 
  REAL p1[4], p2[4], q1[4], q2[4]; 
  REAL resB[4]; 

  vector[0] = 0; 
  vector[1] = 0; 
  vector[2] = 1.0; 
  vector[3] = 0.0; 

  Trans_Stack s; 
  s.loadmatrix(screen_to_data);

  s.multpoint4d(vector,resB); 

  p1[0] = 0; 
  p1[1] = 0; 
  p1[2] = 0.0; 
  p1[3] = 1.0; 

  s.multpoint4d(p1,q1); 

  p2[0] = 0; 
  p2[1] = 0; 
  p2[2] = 1.0; 
  p2[3] = 1.0; 

  s.multpoint4d(p2,q2); 

  vx = q2[0]-q1[0];  
  vy = q2[1]-q1[1];  
  vz = q2[2]-q1[2];  

  printf(" **** compare %f %f %f to %f %f %f ******\n", 
	 resB[0], resB[1], resB[2], vx, vy, vz); 


}
float volumeRender::distance_to_viewplane(float x, float y, float z)
{
  REAL p[4],q[4];
  extern  void matrix_mult(Matrix,REAL*,REAL*);

  p[0] = (REAL)x; p[1] = (REAL)y; p[2] = (REAL)z; p[3] = 1.0;

  matrix_mult(data_to_screen,p,q);
  printf(" %f %f %f %f \n",q[0], q[1], q[2], q[3]);  
  return(q[2]); 
}
void volumeRender::get_eye_vector(float &ex, float &ey, float&ez)
{
  ex = eye_nonN.u; 
  ey = eye_nonN.v; 
  ez = eye_nonN.w; 
}

////////////////////////////////////////////////////////////////////
int volumeRender::mapLookup(float val, float rgba[4])
{
  int id =  (int)(((val - curMin) * lookupSize)/(curMax-curMin)); 
  if (id <0) id = 0; 
  else if (id >=lookupSize) id = lookupSize-1; 
  int cid = id*4; 
  rgba[0] = lookup[cid];   rgba[1] = lookup[cid+1]; 
  rgba[2] = lookup[cid+2]; rgba[3] = lookup[cid+3]; 
  return(1); 
}
///////////////////////////////////////////////////////////////////

void volumeRender::setColorMap(int table_size, float *table)
{
  lookup = table;
  lookupSize = table_size; 
}
///////////////////////////////////////////////////////////////////

int volumeRender::readCmapFile(char *filename)
{
    FILE *in;

    in = fopen(filename, "r");
    if (!in)
        return (0);

    char buf[81], cname[81];
    int type, numRegions, numCP;
    int i, j, total;
    float norm_min, norm_max, vol_min, vol_max;
    float red, green, blue, opac;

    if (fgets(buf, 81, in))
        sscanf(buf, "%s%d%d%d%f%f%f%f", cname, &type, &numRegions,
                    &lookupSize, &norm_min, &norm_max, &vol_min, &vol_max);


    float range = vol_max - vol_min;
    curMin = vol_min + norm_min * range;
    curMax = vol_min + norm_max * range;

    float *ctable = new float[4*lookupSize];

    // read numCP for index transfer function
    fgets(buf, 81, in);
    sscanf(buf, "%s%d%d", cname, &type, &numCP);

    // skip control points for index transfer function
    for ( i=0; i < numCP; i++)
         fgets(buf, 81, in);

    // skip control points for color transfer functions
    for ( i=0; i < 4; i++) {
         fgets(buf, 81, in);
         sscanf(buf, "%d%d", &type, &numCP);
         for (j = 0; j < numCP; j++)
              fgets(buf, 81, in);
    }

    i = 0;
    total = 4 * lookupSize;
    while (fgets(buf, 81, in) && i < total) {
         sscanf(buf, "%f%f%f%f", &red, &green, &blue, &opac);
         ctable[i++] = red;
         ctable[i++] = green;
         ctable[i++] = blue;
         ctable[i++] = opac;
    }
    fclose(in);

    setColorMap(lookupSize, ctable);

    return (1);
}

void volumeRender::set_image_size(int usize, int vsize)
{
  udim = usize; vdim = vsize; 
  if (image!=NULL) free(image); 
  image = image_new(0,udim-1,0,vdim-1);
  set_view(xangle, yangle, zangle); 
}

//////////////////////////////////////////////////////////////////////
//
//         Compute minmax values from volume data
//
void volumeRender::computeVolMinMax(float* volume, int size,
                                   float &vol_min, float &vol_max)
{
    assert(volume!=NULL);

    vol_min = vol_max = volume[0];

    for (int i=1; i < size; i++) {
         if (volume[i] < vol_min)
             vol_min = volume[i];
         else if (volume[i] > vol_max)
             vol_max = volume[i];
    }
}

