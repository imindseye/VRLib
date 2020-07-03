#ifndef MAP_H
#define MAP_H

// standard include files
#include <stdio.h>

#define TRUE 1
#define FALSE 0

// project include files

// local defines
#define MAP_TABLE_SIZE 256
#define MAP_CONST        1.0
#define MAP_LINEAR       2.0

// number of REALs per entry (data,func,rgb values,opacity)
#define ITEMSPERENTRY       6

////////////////////////////////////////////////
//
//    A color map object class
//
class Map {

public:
  // the data values
  REAL data[MAP_TABLE_SIZE];
  // mode of transformation
  REAL func[MAP_TABLE_SIZE];
  // the map values
  REAL rvalue[MAP_TABLE_SIZE];
  REAL gvalue[MAP_TABLE_SIZE];
  REAL bvalue[MAP_TABLE_SIZE];
  REAL opacity[MAP_TABLE_SIZE];
  // the number of entries in the table
  int entries;
  // debugging message buffer
  char message[256];

  // adds a new range to the map
  int add(REAL value, REAL rgbo[4], int *pos);

  Map();

  Map(REAL min, REAL fnc1, REAL rgbo1[4],
      REAL max, REAL fnc2, REAL rgbo2[4]);

  ~Map(void);
  
  void read_file(char*);       // reads the map function from a file
  void write_file(char*);      // writes the map function to a file

  // looks up the rgb values for a data value; returns false if out of range
  int lookup(REAL value, REAL rgbo[4]);

  // adds a new range to the map map using a constant
  void add_constant(REAL value, REAL rgbo[4]);
  
  // adds a new range to the map map using a linear interp
  void add_linear(REAL value, REAL rgbo[4]);

  // delete an entry
  int delete_entry(REAL value);

  // change an entry
  int change_entry(REAL value, REAL fnc, REAL rgbo[4]);

};

#endif
