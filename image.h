/*
 *  cm5 volume renderer: image.h
 *
 * Purpose:
 *
 * Definitions for image data types and operations
 *
 * James S. Painter (painter@cs.utah.edu)     August 15 1992
 * Univesity of Utah, Computer Science Department 
 *
 * Copyright (c) 1992 University of Utah
 * Copying, use and development for non-commercial purposes permitted.
 *                  All rights for commercial use reserved.
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <stdlib.h>

/*
 * Data Types 
 */
typedef unsigned char byte;	/* a convenient data type */

typedef union pixel_u	/* A pixel can be treated as a long or R,G,B,A */
{
  unsigned long lp;	
  struct pixel_b { byte a,b,g,r; } bp;   /* order of bytes must matches X11 */
} pixel;

				/* boundaries of an image */
typedef struct image_bounds_tag 
{
  int umin, umax, vmin, vmax;
} image_bounds_type;

				/* an image  */
typedef struct image_tag
{
  image_bounds_type b;		/* bounds for image */
  int rowsize;		     /* rowsize in pixels (for indexing) */ 
  int totalsize;	     /* total size in bytes (for msg passing) */
  pixel image[1]; 	     /* the image itself. allocated dynamically */
                             /*  immediately after the image header so */  
                             /*  images can be passed in one contiguous msg */ 
} image_type;



/* Macros */

  /* Free memory for an image */
#define image_free(img) (free(img))

  /* bounds check whether a u,v coordinate is inside an image */
#define in_image_p(img,u,v)     \
  ((u) >= (img)->b.umin && u <= (img)->b.umax && \
   (v) >= (img)->b.vmin && v <= (img)->b.vmax)


/* 
 *  Return a pointer to the pixel at coordinates (u,v) within an image 
 *   Asumptions:  Assumes u and v are in the range of the image.
 *                If in doubt, use in_image_p! 
 */
#define image_index(img,u,v) \
  ((img)->image + (((u)-(img)->b.umin) + ((v)-(img)->b.vmin)*(img)->rowsize))


/*
 * See if an image is empty
 */
#define image_empty_p(img)  \
   ((img)->b.umin > (img)->b.umax || (img)->b.vmin > (img)->b.vmax)


/* Prototypes */

  /* Allocate memory for a new image */
extern image_type *image_new( int umin, int umax, int vmin, int vmax );

  /* Same as above but use a static buffer given by the address in buffer */
extern image_type *image_initialize( int umin, int umax, int vmin, int vmax, void  *buffer );
 

/* Composite two images into a third */
enum op_code {OVER,UNDER};
extern  void image_composite( image_type *imgA, image_bounds_type *Abnds,
			      image_type *imgB, enum op_code op);

/* Extract a subimage */
extern image_type *
image_extract
  (image_type *in,		/* Input image */
   image_bounds_type *b		/* image bounds for output */
  );

/* Duplicate an image (quickly) */			     
extern image_type * image_dup( image_type *in );

void zero_rect (image_type *img, int umin, int umax,  int vmin, int vmax );

void copy_rect( image_type *imgA, image_type *imgB, 
	       int umin, int umax, int vmin, int vmax );

#endif /* IMAGE_H */
