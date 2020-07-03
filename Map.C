// standard include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Map.h"

#define DEBUG_write(msg)     { fputs(msg,stderr) ;  exit(-1); }

// Local data type for a Map message
//
typedef struct {
  int entries;
  REAL data[ITEMSPERENTRY*MAP_TABLE_SIZE];
} Map_msg;

/////////////////////////////////////////////
//
//         Default Constructor 
//
Map::Map() {

  data[0] = 0.0;
  func[0] = MAP_LINEAR;
  rvalue[0] = 0.0;
  gvalue[0] = 0.0;
  bvalue[0] = 0.0;
  opacity[0] = 0.0;
  data[1] = 1.0;
  func[1] = MAP_LINEAR;
  rvalue[1] = 1.0;
  gvalue[1] = 1.0;
  bvalue[1] = 1.0;
  opacity[1] = 1.0;
  entries = 2;
		
}
/////////////////////////////////////////////
//
//   Constructor with two initialized
//   entries
//
Map::Map(REAL min, REAL fnc1, REAL rgbo1[4],
	 REAL max, REAL fnc2, REAL rgbo2[4]) {

  data[0] = min;
  func[0] = fnc1;
  rvalue[0] = rgbo1[0];
  gvalue[0] = rgbo1[1];
  bvalue[0] = rgbo1[2];
  opacity[0] = rgbo1[3];
  data[1] = max;
  func[1] = fnc2;
  rvalue[1] = rgbo2[0];
  gvalue[1] = rgbo2[1];
  bvalue[1] = rgbo2[2];
  opacity[1] = rgbo2[3];
  entries = 2;
		
}
///////////////////////////////////////////
//
//    Destructor, nothing special
//    Let the compiler to worry about it
//
Map::~Map() {
	
}
//////////////////////////////////////////
//
// reads in a map table from an HDF file, 
// which is the map format used by Han-Wei 
//
void Map::read_file(char *name) {

  FILE *fp;
  int i;
  
  if ((fp=fopen(name,"r")) == NULL) {
    DEBUG_write("error, can't open map table for reading\n");
    return;
  }
  fscanf(fp,"%d",&i);
  if ((i < 2) || (i > MAP_TABLE_SIZE)) {
    DEBUG_write("error, too many entries in map table file\n");
    return;
  }
  entries = i;
  for (i=0;i<entries;i++) {
    double d,f,r,g,b,a;		/* must be doublees */
    fscanf(fp,"%lf %lf %lf %lf %lf %lf", &d, &f, &r, &g, &b, &a);
    data[i] = d;         /* Assignments here may be converting to float, */
    func[i] = f;         /* depends on how REAL is defined */
    rvalue[i] = r;
    gvalue[i] = g;
    bvalue[i] = b;
    opacity[i] = a;
  }
  fclose(fp);
}
////////////////////////////////////////////
//
// write the current table to an HDF file
//
void Map::write_file(char *name) {

  FILE *fp;
  int i;

  if ((fp=fopen(name,"w")) == NULL) {
    DEBUG_write("error, can't open map table for writing\n");
    return;
  }
  fprintf(fp,"%d\n",entries);
  for (i=0;i<entries;i++) {
    fprintf(fp,"%f %f %f %f %f %f\n",data[i],func[i],
	    rvalue[i],gvalue[i],bvalue[i],opacity[i]);
  }

  fclose(fp);
}
/////////////////////////////////////////////
//
// looks up the rgb values for a data value; 
// return false if out of bounds
//
// This is really slow and a little brain-dead 
// for volume rendering, will implement a 
// quick version -- Han-Wei 
//
//
int Map::lookup(REAL value, REAL rgbo[4]) {

  int base, start, end, i, lim;

  // find which two values the value of interest lies between
  // check for out of bounds

  if (value < data[0])
    return FALSE;
  if (value > data[entries-1])
    return FALSE;

  // Do a binary search looking for the two entries that bracket the
  // value we want.
  base=0;
  for(lim=entries; lim != 0; lim >>=1)
    {
      i = base + (lim >> 1);
      if (value == data[i])
	break;
      if (value > data[i])
	{			// move right
	  base = i + 1;		// ...new base is above our last query
	  lim--;		// ...necessary if lim was odd.
	}
      // else move left
    }
  // i holds the last value we checked.  It is one of the two bracketting
  // values
  if (value >= data[i])	// i and i+1 bracket it
    {
      start = i;
      end = i+1;
    }
  else
    {				// i-1 and i bracket it
      start = i-1;
      end = i;
    }
  if (func[start] == MAP_CONST ||      // constant entry
      start == (entries-1)     ||      // ... or last entry
      value == data[start])            // ... or perfect match
    {
      // No need to interpolate
      rgbo[0] = rvalue[start];
      rgbo[1] = gvalue[start];
      rgbo[2] = bvalue[start];
      rgbo[3] = opacity[start];
      return TRUE;
    }
  // handle linear interp case
  if (func[start] == MAP_LINEAR) 
    {
      rgbo[0] = ((value - data[start]) * (rvalue[end] - rvalue[start])
		 / (data[end] - data[start])) + rvalue[start];
      rgbo[1] = ((value - data[start]) * (gvalue[end] - gvalue[start])
		 / (data[end] - data[start])) + gvalue[start];
      rgbo[2] = ((value - data[start]) * (bvalue[end] - bvalue[start])
		 / (data[end] - data[start])) + bvalue[start];
      rgbo[3] = ((value - data[start]) * (opacity[end] - opacity[start])
		 / (data[end] - data[start])) + opacity[start];
      return TRUE;
    }

  return FALSE;
}

/////////////////////////////////////////////////
//
// adds a new range to the map map using a constant
//
int Map::add(REAL value, REAL rgbo[4], int *pos) {

  int i,j;

  // make sure there is space in the table
  if (entries >= MAP_TABLE_SIZE-1) {
    DEBUG_write("error, no more space in map table\n");
    return FALSE;
  }
	
  // check for out of bounds
  if (value <= data[0]) {
    DEBUG_write("error, below minimum in map table\n");
    return FALSE;
  }
  if (value >= data[entries-1]) {
    DEBUG_write("error, above maximum in map table\n");
    return FALSE;
  }

  // find which two values the value of interest lies between
  i=0;
  while (i<entries-2) {
    if (value < data[i+1])
      break;
    i++;
  }

  // make sure this is not a duplicate entry
  if (data[i] == value) {
    DEBUG_write("error, duplicate map entry\n");
    return FALSE;
  }

  // insert at i+1
  i++;
  for (j=entries;j>=i;j--) {
    data[j] = data[j-1];
    func[j] = func[j-1];
    rvalue[j] = rvalue[j-1];
    gvalue[j] = gvalue[j-1];
    bvalue[j] = bvalue[j-1];
    opacity[j] = opacity[j-1];
  }
  entries++;

  data[i] = value;
  rvalue[i] = rgbo[0];
  gvalue[i] = rgbo[1];
  bvalue[i] = rgbo[2];
  opacity[i] = rgbo[3];

  *pos = i;
  return TRUE;

}
//////////////////////////////////////////////////
//
// adds a new range to the map map using a constant
//
void Map::add_constant(REAL value, REAL rgbo[4]) {

  int i;

  if (add(value,rgbo,&i) == FALSE)
    return;

  func[i] = MAP_CONST;

}
///////////////////////////////////////////////////
//
// adds a new range to the map map using a linear interp
//
void Map::add_linear(REAL value, REAL rgbo[4]) {

  int i;

  if (add(value,rgbo,&i) == FALSE)
    return;

  func[i] = MAP_LINEAR;

}
///////////////////////////////////////////////////
//
// delete an entry
//
int Map::delete_entry(REAL value) {

  int i,j;

  // can't delete end points
  if (value == data[0]) {
    DEBUG_write("error, end point can't be deleted\n");
    return FALSE;
  }
  if (value == data[entries-1]) {
    DEBUG_write("error, end point can't be deleted\n");
    return FALSE;
  }

  for (i=1;i<=entries-2;i++) {
    if (value < data[i]) {
      DEBUG_write("error, map entry doesn't exist\n");
      return FALSE;
    }

    if (value == data[i]) {
      for (j=i;j<entries-1;j++) {
	data[j] = data[j+1];
	func[j] = func[j+1];
	rvalue[j] = rvalue[j+1];
	gvalue[j] = gvalue[j+1];
	bvalue[j] = bvalue[j+1];
	opacity[j] = opacity[j+1];
      }
      entries--;
      return TRUE;
    }
  }

  DEBUG_write("error, map entry doesn't exist\n");
  return FALSE;
      
}
////////////////////////////////////////////////////
//
// change an entry
//
int Map::change_entry(REAL value, REAL fnc, REAL rgbo[4]) {

  int i;

  for (i=0;i<=entries-1;i++) {
    if (value < data[i]) {
      DEBUG_write("error, map entry doesn't exist\n");
      return FALSE;
    }

    if (value == data[i]) {
      func[i] = fnc;
      rvalue[i] = rgbo[0];
      gvalue[i] = rgbo[1];
      bvalue[i] = rgbo[2];
      opacity[i] = rgbo[3];
      return TRUE;
    }
  }

  DEBUG_write("error, map entry doesn't exist\n");
  return FALSE;
      
}
