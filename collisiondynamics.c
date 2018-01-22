/***************************************************
* collision dyanmics involving particle i and j
* i and j are in contact at the moment
* virial is calculated
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"


void Collision(int i, int j) //
{
 VECTOR rij,vij,pij;
 double sigmaij2,bij;
 double mu;

 sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
 mu = mass[i]*mass[j]/(mass[i]+mass[j]);

 rij.x = position[i].x - position[j].x;
 rij.y = position[i].y - position[j].y;
 rij.z = position[i].z - position[j].z;
  //pbc
 MinimumImage(&(rij.x));
 MinimumImage(&(rij.y));
 MinimumImage(&(rij.z));
 
 vij.x = velocity[i].x - velocity[j].x;
 vij.y = velocity[i].y - velocity[j].y;
 vij.z = velocity[i].z - velocity[j].z;

 bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z; 

 pij.x = -rij.x*bij/sigmaij2*2.*mu;
 pij.y = -rij.y*bij/sigmaij2*2.*mu;
 pij.z = -rij.z*bij/sigmaij2*2.*mu;

 velocity[i].x = velocity[i].x + pij.x/mass[i];
 velocity[i].y = velocity[i].y + pij.y/mass[i];
 velocity[i].z = velocity[i].z + pij.z/mass[i];
 
 velocity[j].x = velocity[j].x - pij.x/mass[j];
 velocity[j].y = velocity[j].y - pij.y/mass[j];
 velocity[j].z = velocity[j].z - pij.z/mass[j];

 Virial = pij.x*rij.x + pij.y*rij.y + pij.z*rij.z;

 return;
}


