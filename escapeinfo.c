/***************************************************
* find escape time of each  particle i
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"


void EscapeInfo(int i) //
{
 double xold,yold,zold;

 int cello,celln;

 int nx,ny,nz;
 double xa,xb;
 double ya,yb;
 double za,zb;
 double tx,ty,tz;
 int onboundx,onboundy,onboundz;

if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)
{
 //indicator of particles on boundary of cells
 onboundx = 0;
 onboundy = 0;
 onboundz = 0;

 nx = (int)(position[i].x/rcell);
 ny = (int)(position[i].y/rcell);
 nz = (int)(position[i].z/rcell);

//avoid nx = LCell
 if(nx == -1) nx = 0;
 if(ny == -1) ny = 0;
 if(nz == -1) nz = 0;

 if(nx == LCell) nx = LCell-1;
 if(ny == LCell) ny = LCell-1;
 if(nz == LCell) nz = LCell-1;
 
 if(fabs(position[i].x/rcell - (double)(nx)) < tolerance)
 {
  if(velocity[i].x < 0.)
  {
    nx = (nx-1+LCell)%LCell;
    onboundx = 1;
  }
 }
 else if(fabs(position[i].x/rcell - (double)(nx+1)) < tolerance)
 {
  if(velocity[i].x > 0.)
  {
    nx = (nx+1)%LCell;
    onboundx = 1;
  }
 }
 
 if(fabs(position[i].y/rcell - (double)(ny)) < tolerance)
 {
  if(velocity[i].y < 0.)
  {
    ny = (ny-1+LCell)%LCell;
    onboundy = 1;
  }
 }
 else if(fabs(position[i].y/rcell - (double)(ny+1)) < tolerance)
 {
  if(velocity[i].y > 0.)
  {
    ny = (ny+1)%LCell;
    onboundy = 1;
  }
 }
 
 if(fabs(position[i].z/rcell - (double)(nz)) < tolerance)
 {
  if(velocity[i].z < 0.)
  {
    nz = (nz-1+LCell)%LCell;
    onboundz = 1;
  }
 }
 else if(fabs(position[i].z/rcell - (double)(nz+1)) < tolerance)
 {
  if(velocity[i].z > 0.)
  {
    nz = (nz+1)%LCell;
    onboundz = 1;
  }
 }
 
 xa = nx*rcell;
 xb = (nx+1)*rcell;
 ya = ny*rcell;
 yb = (ny+1)*rcell;
 za = nz*rcell;
 zb = (nz+1)*rcell;
 
 if(velocity[i].x > 0.) //->
 {
   if(onboundx == 1)
    tx = rcell/velocity[i].x;
   else
    tx = (xb - position[i].x)/velocity[i].x;
 }
 else // <-
 {
   if(onboundx == 1)
    tx = -rcell/velocity[i].x;
   else
    tx = (xa - position[i].x)/velocity[i].x;
 }
 
 if(velocity[i].y > 0.) //->
 {
   if(onboundy == 1)
    ty = rcell/velocity[i].y;
   else
    ty = (yb - position[i].y)/velocity[i].y;
 }
 else // <-
 {
   if(onboundy == 1)
    ty = -rcell/velocity[i].y;
   else
    ty = (ya - position[i].y)/velocity[i].y;
 }
 
 if(velocity[i].z > 0.) //->
 {
   if(onboundz == 1)
    tz = rcell/velocity[i].z;
   else
    tz = (zb - position[i].z)/velocity[i].z;
 }
 else // <-
 {
   if(onboundz == 1)
    tz = -rcell/velocity[i].z;
   else
    tz = (za - position[i].z)/velocity[i].z;
 }

 EscapeTime[i] = MIN(tx,ty);
 EscapeTime[i] = MIN(EscapeTime[i],tz);

} // end if cell list is used

 return;
}


