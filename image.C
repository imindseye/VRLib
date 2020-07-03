/*
 *  cm5 volume renderer: image.c
 *
 * Purpose:
 *
 * Basic operations on images
 *
 * James S. Painter (painter@cs.utah.edu)     August 15 1992
 * Univesity of Utah, Computer Science Department 
 *
 * Copyright (c) 1992 University of Utah
 * Copying, use and development for non-commercial purposes permitted.
 *                  All rights for commercial use reserved.
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "image.h"
#include "minmax.h"

image_type *image_initialize( int umin, int umax, int vmin, int vmax, void *buf)
{
  int rowsize, npixels;
  image_type *result;
  size_t totalsize;

  if (umin > umax || vmin > vmax)
    {
      /* null images are allowed */
      npixels = 0;
      rowsize = 0;
    }
  else
    {
      rowsize = (umax-umin+1);
      npixels = rowsize*(vmax-vmin+1);
    }

  result = (image_type *) buf;

  result->b.umin = umin;
  result->b.umax = umax;
  result->b.vmin = vmin;
  result->b.vmax = vmax;
  result->rowsize = rowsize;
  result->totalsize = sizeof(image_type) + sizeof(pixel)*(npixels);

  /* The image data follows directly after the image */

  return result;
}

image_type *image_new( int umin, int umax, int vmin, int vmax)
{
  int rowsize, npixels;
  void *result;
  size_t totalsize;

  if (umin > umax || vmin > vmax)
    {
      /* null images are allowed */
      npixels = 0;
      rowsize = 0;
    }
  else
    {
      rowsize = (umax-umin+1);
      npixels = rowsize*(vmax-vmin+1);
    }

  totalsize = sizeof(image_type) + sizeof(pixel)*(npixels);
  result =  malloc( totalsize );
  if (result == NULL)
    printf( "no memory for image_new" );
  return image_initialize( umin, umax, vmin, vmax, result );
}


  
/* Extract a subimage */
extern image_type *
image_extract
  (image_type *in,		/* Input image */
   image_bounds_type *b		/* bounds for output image */
   )
{
  int umin, umax, vmin, vmax, v;
  image_type *out;

  /* Determine output image bounds */
  umin = MAX(b->umin, in->b.umin);
  umax = MIN(b->umax, in->b.umax);
  vmin = MAX(b->vmin, in->b.vmin);
  vmax = MIN(b->vmax, in->b.vmax);

  /* Create the image */
  out = image_new( umin, umax, vmin, vmax );

  /* Fill it in row-by-row */
  for(v=vmin; v<=vmax; v++)
    {
      memcpy( image_index(out,umin,v), image_index(in,umin,v), 
	     out->rowsize*sizeof(pixel) );
    }
  return out;
}


/*
 * Duplicate an image (quickly)
 */
image_type *image_dup(image_type *img)
{
  image_type *result = (image_type *) malloc(img->totalsize);

  memcpy( result, img, img->totalsize );
  return result;
}


void copy_rect( image_type *imgA, image_type *imgB, 
	       int umin, int umax, int vmin, int vmax )
{
  int rowsize = (umax - umin + 1)*sizeof(pixel);
  int columns = (vmax - vmin + 1);
  int instride  = imgB->rowsize;
  int outstride = imgA->rowsize;
  register pixel *pin, *pout, *pend;

  if (rowsize < 0 || columns < 0)
    return;

  pin  = image_index( imgB, umin, vmin );
  pout = image_index( imgA, umin, vmin );
  pend = pin + columns*instride;
  for( ;pin < pend; pin += instride, pout += outstride)
    memcpy( pout, pin, rowsize );
}

void zero_rect( image_type *img, 
		       int umin, int umax, int vmin, int vmax )
{
  register int i,v;
  register pixel *p, *pend;
  register int rowsize = (umax - umin + 1)*sizeof(pixel);

  if (rowsize < 0) return;

  p = image_index( img, umin, vmin );
  pend = p + (vmax-vmin+1)*img->rowsize;
  while (p <pend)
    {
      memset( p, '\0', rowsize );
      p += img->rowsize;
    }
}

