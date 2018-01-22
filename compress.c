/***************************************************
* compression
* 
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Compress(void) // 
{
 int i,j;
 double ratio;
 double scaling;

 //ratio = compressfactor*compressratio;
 
 ratio = 1.0 + compressfactor*(compressratio-1.0);
// printf("compress ratio = %lf\t",compressratio);
// printf("ratio = %.15lf\n",ratio);

 if(ratio - 1.0 < tolerance) 
 {
  printf("compression cause particles contact\n");
 //return;
  exit(1);
 }
 else
 {

 Vo = V;
 Lo = L;
// packingfraction = M_PI/6.*rho*(fA*CUBIC(sigmaA*ratio)+fB*CUBIC(sigmaB*ratio));
 packingfraction = M_PI/6.*rho*(fA*CUBIC(sigmaA*ratio)+fB*CUBIC(sigmaB*ratio)+fC*CUBIC(sigmaC*ratio));
// rho = packingfraction/M_PI*6./(fA*CUBIC(sigmaA)+fB*CUBIC(sigmaB));
 rho = packingfraction/M_PI*6./(fA*CUBIC(sigmaA)+fB*CUBIC(sigmaB)+fC*CUBIC(sigmaC));
 V = NumberOfParticles/rho*CUBIC(sigmaA);

 rhoA = NA/V*CUBIC(sigmaA);
 rhoB = NB/V*CUBIC(sigmaA);
 rhoC = NC/V*CUBIC(sigmaA);
 L = pow(V,1.0/3.0); 
 scaling = L/Lo;


 for(i=0;i<NumberOfParticles;i++)
 {
  position[i].x *= scaling;
  position[i].y *= scaling;
  position[i].z *= scaling;
 }
 
 if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)
 {
 rcell *= scaling;
 if(rcell < rc) MakeCell();
// problem? on boundary when compressed???
 }

 return;
 }
}
