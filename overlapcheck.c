/***************************************************
* check overlaps between every pair
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"


void OverlapCheck(void) //
{
 int i,j;
 double rij2;
 double sigmaij2;
 double tor;

 tor = 10.E-10;

 for(i=0;i<NumberOfParticles-1;i++)
 for(j=i+1;j<NumberOfParticles;j++)
 {
  rij2 = Distance(i,j);
  sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));

  if(rij2 - sigmaij2 < -tor)
  {
   printf("overlap occurs between %d and %d, rij = %.30lf\tsigmaij = %.30lf\n",i,j,sqrt(rij2),sqrt(sigmaij2));
//   exit(1);
  }
 }

 return;
}
