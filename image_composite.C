/*
 *  cm5 volume renderer: image_composite
 *
 * Purpose:
 *
 * Composite two local images using the Cook/Porter/Duff OVER
 * operator.         
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
#include "image.h"
#include "minmax.h"

static image_bounds_type union_bounds
  ( image_bounds_type *b1, image_bounds_type *b2)
{
  image_bounds_type result;


  result.umin = MIN(b1->umin,b2->umin);
  result.umax = MAX(b1->umax,b2->umax);

  result.vmin = MIN(b1->vmin,b2->vmin);
  result.vmax = MAX(b1->vmax,b2->vmax);
  return result;
}

/* Zero out any part of img within the bounding box of the union of 
   b1 and b2 but not actually in either bound.
 */
static image_bounds_type 
  zero_empty_part(image_type *img, 
		  image_bounds_type *b1,
		  image_bounds_type *b2 )
{
  image_bounds_type b, *btmp;
  int vmin, vmax;

  b = union_bounds(b1,b2);

  /* Arrange it so that the top of b1 is the same or higher than the top 
   * of b2 
   */
  if (b1->vmax < b2->vmax)
    { btmp = b1;      b1 = b2;      b2 = btmp;  }
			    
  /* Zero out the part between the top of b1 and the top of b2 */
  vmin = MAX(b2->vmax+1,b1->vmin);
  zero_rect( img, b.umin, b1->umin-1, vmin, b1->vmax );
  zero_rect( img, b1->umax+1, b.umax, vmin, b1->vmax );

  /* Zero out the part between the bottom of b2 and the bottom of b1 */
  if (b2->vmin < b1->vmin)
    {
      vmax = MIN(b2->vmax,b1->vmin-1);
      zero_rect( img, b.umin, b2->umin-1, b2->vmin, vmax );
      zero_rect( img, b2->umax+1, b.umax, b2->vmin, vmax );
    }
  else
    {
      vmax = MIN(b1->vmax,b2->vmin-1);
      zero_rect( img, b.umin, b1->umin-1, b1->vmin, vmax );
      zero_rect( img, b1->umax+1, b.umax, b1->vmin, vmax );
    }
  
  /* Handle the overlap scanlines */
  if (b2->umax < b1->umin)
    zero_rect( img, b2->umax+1, b1->umin-1, b1->vmin, b2->vmax );
  else
    zero_rect( img, b1->umax+1, b2->umin-1, b1->vmin, b2->vmax );

  /* Handle the scan lines between the two regions */
  if (b2->vmax < b1->vmin)
    zero_rect( img, b.umin, b.umax, b2->vmax+1, b1->vmin-1 );

  return b;
}

void
  image_composite( image_type * imgA, image_bounds_type *bndsA,
 		   image_type * imgB, enum op_code op )
{
  image_bounds_type unionB;
  register pixel *pxlA, *pxlB;
  register REAL one_minus_alpha;
  int u,v; 
  int umin, umax, vmin, vmax;
  int leftsize, rightsize; 

/*
  printf("%d %d %d %d %d %d %d %d\n",
	 bndsA->umin, bndsA->umax, bndsA->vmin, bndsA->vmax,
	 imgB->b.umin, imgB->b.umax, imgB->b.vmin, imgB->b.vmax );
  fflush(stdout);
*/
  /* Check for quick kill case: one image is empty */
  if (bndsA->umin > bndsA->umax || bndsA->vmin > bndsA->vmax)
    { 
      /* A empty: Copy B into A */
      copy_rect( imgA, imgB, imgB->b.umin, imgB->b.umax,
		 imgB->b.vmin, imgB->b.vmax );
      *bndsA = imgB->b;  
      return; 
    }

  if (imgB->b.umin > imgB->b.umax || imgB->b.vmin > imgB->b.vmax)
    return;			/* B empty: nothing to do */



  /* Find the bounds of the overlap region and 
     Zero the empty part of the resulting area */

  unionB =  zero_empty_part( imgA, bndsA, &imgB->b );

  /* Copy scanlines in img B which are outside the bounds of the "live"
   *  part of A.
   */

  copy_rect( imgA, imgB, imgB->b.umin, imgB->b.umax, 
	    imgB->b.vmin, bndsA->vmin-1 );

  copy_rect( imgA, imgB, imgB->b.umin, imgB->b.umax, 
	    bndsA->vmax+1, imgB->b.vmax );

  /* Figure out the overlap region between the two image bounds*/
  umin = MAX(bndsA->umin,imgB->b.umin);
  umax = MIN(bndsA->umax,imgB->b.umax);
  vmin = MAX(bndsA->vmin,imgB->b.vmin);
  vmax = MIN(bndsA->vmax,imgB->b.vmax);

  leftsize = (umin - imgB->b.umin)*sizeof(pixel);
  rightsize = (imgB->b.umax - umax)*sizeof(pixel);

  /* Step through all the pixels in image A compositing them over B */
  for(v=vmin; v<=vmax; v++)
    {
      if (leftsize > 0)
	memcpy( image_index( imgA, imgB->b.umin, v),
	        image_index( imgB, imgB->b.umin, v),
	        leftsize );
      if (rightsize > 0)
	memcpy( image_index( imgA, umax+1, v),
	        image_index( imgB, umax+1, v),
	        rightsize );
      pxlA = image_index( imgA, umin, v );
      pxlB = image_index( imgB, umin, v );

      /* Composite A over B with the result back to A */
      if (op == OVER)
	for(u=umin; u<=umax; u++, pxlA++, pxlB++)
	  {

	    if (pxlA->bp.a == 0)
	      *pxlA  = *pxlB;
/*
	    else if (pxlA->bp.a == 255)
	      continue;
*/
	    else
	      {
		one_minus_alpha = (REAL) (255 - pxlA->bp.a) / 255.0;

		pxlA->bp.r = (byte) ((REAL) pxlA->bp.r + 
				     (REAL) pxlB->bp.r*one_minus_alpha);
		pxlA->bp.g = (byte) ((REAL) pxlA->bp.g + 
				     (REAL) pxlB->bp.g*one_minus_alpha);
		pxlA->bp.b = (byte) ((REAL) pxlA->bp.b + 
				     (REAL) pxlB->bp.b*one_minus_alpha);
		pxlA->bp.a = (byte) ((REAL) pxlA->bp.a + 
				     (REAL) pxlB->bp.a*one_minus_alpha);
	      }
	  }
      else
	for(u=umin; u<=umax; u++, pxlA++, pxlB++)
	  {
/*
	    if (pxlB->bp.a == 255)
	      *pxlA  = *pxlB;
	    else
*/
	    if (pxlB->bp.a == 0)
	      continue;
	    else
	      {

		one_minus_alpha = (REAL) (255 - pxlB->bp.a)  / 255.0;
		
		pxlA->bp.r = (byte) ((REAL) pxlB->bp.r + 
				     (REAL) pxlA->bp.r*one_minus_alpha);
		pxlA->bp.g = (byte) ((REAL) pxlB->bp.g + 
				     (REAL) pxlA->bp.g*one_minus_alpha);
		pxlA->bp.b = (byte) ((REAL) pxlB->bp.b + 
				     (REAL) pxlA->bp.b*one_minus_alpha);
		pxlA->bp.a = (byte) ((REAL) pxlB->bp.a + 
				     (REAL) pxlA->bp.a*one_minus_alpha);
	      }
	  }
    }
  /* Return the new bounds */
  *bndsA = unionB;
}









