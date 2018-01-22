/***************************************************
* return distance r, u(r), du/dr, u_tail, w_tail
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

double Distance(int i,int j) // return rij^2
{
  double r2,dx,dy,dz;
  
  // miminum image convention under periodic boundary condition
  dx = position[i].x-position[j].x;
  dx = dx - L*round(dx/L);

  dy = position[i].y-position[j].y;
  dy = dy - L*round(dy/L);

  dz = position[i].z-position[j].z;
  dz = dz - L*round(dz/L);
 
  r2 = SQR(dx)+SQR(dy)+SQR(dz); // reduced unit
  //r = sqrt(SQR(dx)+SQR(dy)+SQR(dz));

  return r2;
}

double SigmaIJ(int iID,int jID) // input identity[i]!
{
 double sigmaij;

 if(iID == jID)
 {
  if(iID == 1) // A-A
  sigmaij = sigmaA;
  else if (iID == 2) // B-B
  sigmaij = sigmaB;
  else // C-C
  sigmaij = sigmaC;
 }
 else
 {
  if(iID+jID == 3) // A-B (1,2)
  sigmaij = sigmaAB;
  else if(iID+jID == 5) // B-C (2,3)
  sigmaij = sigmaBC;
  else // C-A (3,1)
  sigmaij = sigmaCA;
 }

 return sigmaij;
}

/******
// A: ID = 1 B:ID = 2
double Potential(int iID,int jID,double r) // u(rij)
{
  double u;
 
  return u;
}

double Potential_dr(int iID,int jID,double r) // du(rij)/dr
{
  double dudr;
 
  return dudr;
}
*********/
