#include "box_collision.hpp"
#include <stdio.h>
bool box_collision(float *min_A, float *max_A, float *min_B, float *max_B)
{
  float shift = 1.0;
  bool bc_A = !( (min_A[0]>max_B[0])||(max_A[0]<min_B[0])||(min_A[1]>max_B[1])||(max_A[1]<min_B[1])||(min_A[2]>max_B[2])||(max_A[2]<min_B[2]));
  bool bc_B = !( (min_A[0]-shift>max_B[0])||(max_A[0]-shift<min_B[0])||(min_A[1]-shift>max_B[1])||(max_A[1]-shift<min_B[1])||(min_A[2]-shift>max_B[2])||(max_A[2]-shift<min_B[2]));
  bool bc_C = !( (min_A[0]+shift>max_B[0])||(max_A[0]+shift<min_B[0])||(min_A[1]+shift>max_B[1])||(max_A[1]+shift<min_B[1])||(min_A[2]+shift>max_B[2])||(max_A[2]+shift<min_B[2]));
  bool bc_D = !( (min_A[0]>max_B[0]-shift)||(max_A[0]<min_B[0]-shift)||(min_A[1]>max_B[1]-shift)||(max_A[1]<min_B[1]-shift)||(min_A[2]>max_B[2]-shift)||(max_A[2]<min_B[2]-shift));
  bool bc_E = !( (min_A[0]>max_B[0]+shift)||(max_A[0]<min_B[0]+shift)||(min_A[1]>max_B[1]+shift)||(max_A[1]<min_B[1]+shift)||(min_A[2]>max_B[2]+shift)||(max_A[2]<min_B[2]+shift));
  return (bc_A||bc_B||bc_C||bc_D||bc_E);
  //return bc_A;
}

bool region_collision(float *min_A, float *max_A, float *min_B, float *max_B, float b)
{
  float shift[3] = {0.0, 0.0, 0.0};

  for(int k=0;k<3;k++)
  {

    if(min_A[k] - min_B[k] <= -0.5 )
      shift[k] = -1.0;

    if(min_A[k] - min_B[k] >=   0.5 )
      shift[k] =  1.0;

    shift[k] = get_shift(min_A[k], min_B[k]);
  }


  bool bc_A = !( (min_A[0]-b>max_B[0]+shift[0])||(max_A[0]+b<min_B[0]+shift[0])||(min_A[1]-b>max_B[1]+shift[1])||(max_A[1]+b<min_B[1]+shift[1])||(min_A[2]-b>max_B[2]+shift[2])||(max_A[2]+b<min_B[2]+shift[2]));

  return bc_A;
}

bool adjacent_subvolume(float *min_A, float *max_A, float *min_B, float *max_B)
{
  double dx[3]; //half width of a subvolume in each dimension

  for(int k=0;k<3;k++)
    dx[k] = (max_A[k] - min_A[k]);

  bool bc_A = false;

  //consider cells within min_A[k]-1.2*dx[k], max_A[k]+1.2*dx[k]
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;

  //shift by - 1 in x direction
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;
  //shift by + 1 in x direction
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;

  //shift by - 1 in y direction
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;
  //shift by + 1 in y direction
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;

  //shift by - 1 in z direction
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;
  //shift by + 1 in z direction
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;

  //shift by - 1 in x,y direction
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;
  //shift by + 1 in x,y direction
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;

  //shift by - 1 in x,z direction
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;
  //shift by + 1 in x,z direction
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;  

  //shift by - 1 in y,z direction
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;
  //shift by + 1 in y,z direction
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;  

  //shift by - 1 in x,y,z direction
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;

  //shift by + 1 in x,y,z direction
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;  

  //shift by -1 in x, +1 in y
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;
  //shift by +1 in x, -1 in y
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]-1.2*dx[2] && max_B[2] < max_A[2]+1.2*dx[2])
    bc_A = true;

  //shift by -1 in x, +1 in z
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;
  //shift by +1 in x, -1 in z
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]-1.2*dx[1] && max_B[1] < max_A[1]+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;

  //shift by -1 in y, +1 in z
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;

  //shift by +1 in y, -1 in z
  if(min_B[0]>min_A[0]-1.2*dx[0] && max_B[0]<max_A[0]+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;

  //shift by -1 in x, -1 in y, +1 in z
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;

  //shift by +1 in x, +1 in y, -1 in z
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;

  //shift by +1 in x, -1 in y, +1 in z
  if(min_B[0]>min_A[0]+1-1.2*dx[0] && max_B[0]<max_A[0]+1+1.2*dx[0] && min_B[1] > min_A[1]-1-1.2*dx[1] && max_B[1] < max_A[1]-1+1.2*dx[1] && min_B[2] > min_A[2]+1-1.2*dx[2] && max_B[2] < max_A[2]+1+1.2*dx[2])
    bc_A = true;

  //shift by -1 in x, +1 in y, -1 in z
  if(min_B[0]>min_A[0]-1-1.2*dx[0] && max_B[0]<max_A[0]-1+1.2*dx[0] && min_B[1] > min_A[1]+1-1.2*dx[1] && max_B[1] < max_A[1]+1+1.2*dx[1] && min_B[2] > min_A[2]-1-1.2*dx[2] && max_B[2] < max_A[2]-1+1.2*dx[2])
    bc_A = true;

  return bc_A;
}


bool adjacent_region(float *min_A, float *max_A, float *min_B, float *max_B, float b)
{
  double ff = 1.2; //safety
  double dx[3]; //half width of a subvolume in each dimension

  for(int k=0;k<3;k++)
    dx[k] = b;
    //dx[k] = (max_A[k] - min_A[k]);

  bool bc_A = false;

  //consider regions that overlap with no shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true;

  //consider regions that overlap with -1 x-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true;  

  //consider regions that overlap with +1 x-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true; 

  //consider regions that overlap with -1 y-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true;  

  //consider regions that overlap with +1 y-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true; 

  //consider regions that overlap with -1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;  

  //consider regions that overlap with +1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true; 

  //consider regions that overlap with -1 x,y-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true;  

  //consider regions that overlap with +1 x,y-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true;  

  //consider regions that overlap with -1 x,z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;  

  //consider regions that overlap with +1 x,z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true;  

  //consider regions that overlap with -1 y,z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;  

  //consider regions that overlap with +1 y,z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true;  

  //consider regions that overlap with -1 x,y,z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;  

  //consider regions that overlap with +1 x,y,z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true; 

  //consider regions that overlap with -1 x, +1 y-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true;  

  //consider regions that overlap with +1 x, -1 y-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2] || max_B[2]<min_A[2]-ff*dx[2] ))
    bc_A = true; 

  //consider regions that overlap with -1 x +1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true; 

  //consider regions that overlap with +1 x -1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1] || max_B[1]<min_A[1]-ff*dx[1] ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true; 

  //consider regions that overlap with -1 y, +1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true;  

  //consider regions that overlap with +1 y, -1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0] || max_B[0]<min_A[0]-ff*dx[0] ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;  

  //consider regions that overlap with -1 x, -1 y, +1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true;

  //consider regions that overlap with +1 x, +1 y, -1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;

  //consider regions that overlap with +1 x, -1 y, +1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]+1 || max_B[0]<min_A[0]-ff*dx[0]+1 ||
        min_B[1]>max_A[1]+ff*dx[1]-1 || max_B[1]<min_A[1]-ff*dx[1]-1 ||
        min_B[2]>max_A[2]+ff*dx[2]+1 || max_B[2]<min_A[2]-ff*dx[2]+1 ))
    bc_A = true;

  //consider regions that overlap with -1 x, +1 y, -1 z-direction shift
  if(!( min_B[0]>max_A[0]+ff*dx[0]-1 || max_B[0]<min_A[0]-ff*dx[0]-1 ||
        min_B[1]>max_A[1]+ff*dx[1]+1 || max_B[1]<min_A[1]-ff*dx[1]+1 ||
        min_B[2]>max_A[2]+ff*dx[2]-1 || max_B[2]<min_A[2]-ff*dx[2]-1 ))
    bc_A = true;


  return bc_A;
}

float get_shift(float min_A, float min_B)
{
  //get the periodic box wrap shift for B
  //if the min of region B >= min_A + 0.5, then shift=-1
  //if the min of region B <= min_A - 0.5, then shift=+1
  //otherwise, shift = 0.0

  //min_A 3.749165e-01 7.499136e-01 2.499516e-01 max_A 5.001888e-01 1.000055e+00 5.002480e-01 
  //min_B 4.997529e-01 7.499428e-01 2.499662e-01 max_B 6.251784e-01 1.000160e+00 5.001441e-01
  //shift 1.000000e+00 1.000000e+00 1.000000e+00 bc_A 0
  float shift = 0.0;

  if(min_B >= min_A + 0.5)
    shift = -1.0;

  if(min_B <= min_A - 0.5)
    shift =  1.0;
 
  return shift; 
}
