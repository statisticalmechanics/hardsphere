/***********************************************************
* determine which cell does particle i with
* coordinate x,y,z belong to
* be careful when particles are on the boundary of cells
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include <math.h>

int CellDetermine(int i)
{
 int ncell;
 int nx,ny,nz;

 nx = (int)(position[i].x/rcell);
 ny = (int)(position[i].y/rcell);
 nz = (int)(position[i].z/rcell);
 
 //avoid n = LCell
 if(nx == -1) nx = 0;
 if(ny == -1) ny = 0;
 if(nz == -1) nz = 0;

 if(nx == LCell) nx = LCell-1;
 if(ny == LCell) ny = LCell-1;
 if(nz == LCell) nz = LCell-1;


 if(fabs(position[i].x/rcell - (double)(nx)) < tolerance)
 {
  if(velocity[i].x < 0.)
    nx = (nx-1+LCell)%LCell;
 }
 else if(fabs(position[i].x/rcell - (double)(nx+1)) < tolerance)
 {
  if(velocity[i].x > 0.)
    nx = (nx+1)%LCell;
 }
 
 if(fabs(position[i].y/rcell - (double)(ny)) < tolerance)
 {
  if(velocity[i].y < 0.)
    ny = (ny-1+LCell)%LCell;
 }
 else if(fabs(position[i].y/rcell - (double)(ny+1)) < tolerance)
 {
  if(velocity[i].y > 0.)
    ny = (ny+1)%LCell;
 }
 
 if(fabs(position[i].z/rcell - (double)(nz)) < tolerance)
 {
  if(velocity[i].z < 0.)
    nz = (nz-1+LCell)%LCell;
 }
 else if(fabs(position[i].z/rcell - (double)(nz+1)) < tolerance)
 {
  if(velocity[i].z > 0.)
    nz = (nz+1)%LCell;
 }
 
 ncell = (nz*LCell+ny)*LCell+nx;

 return ncell;
}
