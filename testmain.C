#include<stdio.h>
#include<render.h>

void usage(char* prgm) {
  printf(" usage: %s udim vdim volume colormap alpha beta gamma out \n", 
	 prgm); 
  exit(0); 
}

int main(int argc, char* argv[]) {

  if (argc!= 9) usage(argv[0]); 

  int udim = atoi(argv[1]); 
  int vdim = atoi(argv[2]); 

  float alpha = (float) atoi(argv[5]); 
  float beta = (float) atoi(argv[5]); 
  float gamma = (float) atoi(argv[5]); 

  FILE* in;
  in = fopen(argv[3],"r"); 
  if (in == NULL) {
    printf(" can't open file %s\n", argv[3]); 
    exit(0);
  }

  printf(" read file %s ....\n", argv[3]); 
  int xdim, ydim, zdim; 
  fread(&xdim, sizeof(int), 1, in);
  fread(&ydim, sizeof(int), 1, in);
  fread(&zdim, sizeof(int), 1, in);

  int junk =1; 
  //  fread(&junk, sizeof(int), 1, in);
  printf(" %d %d %d %d\n", xdim, ydim, zdim, junk); 

  int size = xdim*ydim*zdim;
  float *volume = new float[size];
  fread(volume, sizeof(float), size, in);

  volumeRender vr(xdim,ydim,zdim,udim,vdim,volume); 
  vr.readCmapFile(argv[4]); 
  vr.set_view(alpha, beta, gamma); 
  vr.execute(); 
  vr.out_to_image(argv[8]); 
  

}
