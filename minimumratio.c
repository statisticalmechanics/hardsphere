/***************************************************
* find particle i and j whose surface distance
* is minimum dmin 
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

void MinRatio(void) // update i,j,minimum compress ratio
{
 double r,ratio;
 double ri,rj;
 int i,j;
 double sigmaij;

 int iCell,cell,jCell;

 compressratio = 10.E10;

if(CellSwitch == 0 || NumberOfCells<=NumberOfNeighborCells) // not use cell
{
 for(i=0;i<NumberOfParticles-1;i++)
 for(j=i+1;j<NumberOfParticles;j++)
 {
  r = Distance(i,j);
  r = sqrt(r); 
  sigmaij = SigmaIJ(identity[i],identity[j]);
  ratio = r/sigmaij;

  if(ratio<compressratio)
  {
   compressratio = ratio;
   icompress = i;
   jcompress = j;
  }
 }//end loop i,j
}//end if not use cell
else // use cell list
{
 for(i=0;i<NumberOfParticles;i++)
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
    r = Distance(i,j);
    r = sqrt(r); 
  sigmaij = SigmaIJ(identity[i],identity[j]);
  ratio = r/sigmaij;
    if(ratio<compressratio)
    {
     compressratio = ratio;
     icompress = i;
     jcompress = j;
    }
    }//end if j != i
    j = CellTrack[j].Next;
   }//end linked list
  }//end loop 27 cell
 }//end loop particle i
}

 return;
}
