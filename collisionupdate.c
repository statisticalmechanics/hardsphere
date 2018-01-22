/**********************************************************
* update collision info involving the transferred particle
* do not need to loop all the 27 cells
* the new 9 cells are enough, ideally
* need to update both i and j collision time
***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void CollisionUpdate(int i, int iCell) //
{
 int j;
 VECTOR rij,vij;
 double bij; // rij*vij
 double rij2,vij2,sigmaij2; // square
 double discr; // discriminate b^2 - 4ac
 double tij;
 int cell,jCell;

 //CollisionTime[i] = TimeBig; do not set to TimeBig

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
  } //end if j != i

    j = CellTrack[j].Next;
  }//end while j != -1
 }//end loop cell

 return;
}
