/***************************************************
* find collision time and partner of particle i
* and particle j
* can use cell list
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"


void CollisionInfo(int i) //
{
 int j;
 VECTOR rij,vij;
 double bij; // rij*vij
 double rij2,vij2,sigmaij2; // square
 double discr; // discriminate b^2 - 4ac
 double tij;
 int iCell,cell,jCell;
 int l,m,n; // pbc box index
 double xtemp,ytemp,ztemp;

 CollisionTime[i] = TimeBig;

if(NumberOfCells<=NumberOfNeighborCells || CellSwitch == 0) // not use cell list
{
 for(j=0;j<NumberOfParticles;j++) 
 {
 if(j != i)
 {
  vij.x = velocity[i].x - velocity[j].x;
  vij.y = velocity[i].y - velocity[j].y;
  vij.z = velocity[i].z - velocity[j].z;
  for(l=-1;l<=1;l++)
  for(m=-1;m<=1;m++)
  for(n=-1;n<=1;n++)
  {
   xtemp = position[j].x + l*L;
   ytemp = position[j].y + m*L;
   ztemp = position[j].z + n*L;

   rij.x = position[i].x - xtemp;
   rij.y = position[i].y - ytemp;
   rij.z = position[i].z - ztemp;
  
   bij = rij.x*vij.x + rij.y*vij.y + rij.z*vij.z; 
   
 //  rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
 //  if(rij2 < sigmaij2)
 //  printf("r2(%d, %d) = %.7lf\n",i,j,rij2);
   
   if(bij < 0.)
   {
   rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
   vij2 = SQR(vij.x) + SQR(vij.y) + SQR(vij.z); 
   sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
   discr = SQR(bij)- vij2*(rij2-sigmaij2);
  
   if(discr > 0.)
   {
    tij = (-bij - sqrt(discr))/vij2;
 //   tij = (-bij - sqrt(discr));
 //   tij = MIN(tij/vij2,(rij2-sigmaij2)/tij);

    if(tij < CollisionTime[i])
    {
     CollisionTime[i] = tij;
     CollisionPartner[i] = j;
    }
    if(tij < CollisionTime[j]) // update j simutaneously
    {
     CollisionTime[j] = tij;
     CollisionPartner[j] = i;
    }
   } // endif discriminant > 0
  } //endif bij < 0
  }// end loop periodic images

 } //end if j is not i
 }//end loop j
}//end if cell list not used
else // use cell list when CellSwitch == 1 && NumberOfCells>27
{
 iCell = CellTrack[i].WhichCell;
 for(cell=0;cell<NumberOfNeighborCells;cell++)
 {
  jCell = NeighborCellList[iCell][cell];
  j=HeadOfChain[jCell];
  while(j != -1)
  {
 if(j != i)
 {
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

  if(bij < 0.)
  {
   rij2 = SQR(rij.x) + SQR(rij.y) + SQR(rij.z); 
   vij2 = SQR(vij.x) + SQR(vij.y) + SQR(vij.z); 
   sigmaij2 = SQR(SigmaIJ(identity[i],identity[j]));
   discr = SQR(bij)- vij2*(rij2-sigmaij2);

   if(discr > 0.)
   {
    tij = (-bij - sqrt(discr))/vij2;

    if(tij < CollisionTime[i])
    {
     CollisionTime[i] = tij;
     CollisionPartner[i] = j;
    }
    if(tij < CollisionTime[j])
    {
     CollisionTime[j] = tij;
     CollisionPartner[j] = i;
    }
   } // endif discriminant > 0
  } //endif bij < 0
 } //end if j is not i

    j = CellTrack[j].Next;
  }//end while j != -1
 }//end loop cell
}//end cell list is used

 return;
}


